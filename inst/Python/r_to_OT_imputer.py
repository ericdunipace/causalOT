#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# from https://github.com/BorisMuzellec/MissingDataOT/blob/master/imputers.py
# adding in non-uniform masses

import numpy as np
import torch
from geomloss import SamplesLoss

# from OTutils import nanmean, MAE, RMSE

import logging


class OTimputer():
    """
    'One parameter equals one imputed value' model (Algorithm 1. in the paper)

    Parameters
    ----------

    eps: float, default=0.01
        Sinkhorn regularization parameter.
        
    lr : float, default = 0.01
        Learning rate.

    opt: torch.nn.optim.Optimizer, default=torch.optim.Adam
        Optimizer class to use for fitting.
        
    max_iter : int, default=10
        Maximum number of round-robin cycles for imputation.

    niter : int, default=15
        Number of gradient updates for each model within a cycle.

    batchsize : int, defatul=128
        Size of the batches on which the sinkhorn divergence is evaluated.

    n_pairs : int, default=10
        Number of batch pairs used per gradient update.

    tol : float, default = 0.001
        Tolerance threshold for the stopping criterion.

    weight_decay : float, default = 1e-5
        L2 regularization magnitude.

    order : str, default="random"
        Order in which the variables are imputed.
        Valid values: {"random" or "increasing"}.

    unsymmetrize: bool, default=True
        If True, sample one batch with no missing 
        data in each pair during training.

    scaling: float, default=0.9
        Scaling parameter in Sinkhorn iterations
        c.f. geomloss' doc: "Allows you to specify the trade-off between
        speed (scaling < .4) and accuracy (scaling > .9)"


    """
    def __init__(self, 
                 power=2,
                 eps=0.01, 
                 lr=1e-2,
                 opt=torch.optim.RMSprop, 
                 niter=2000,
                 batchsize=128,
                 n_pairs=1,
                 noise=0.1,
                 scaling=.9):
        self.eps = eps
        self.lr = lr
        self.opt = opt
        self.niter = niter
        self.batchsize = batchsize
        self.n_pairs = n_pairs
        self.noise = noise
        self.sk = SamplesLoss("sinkhorn", p=power, blur=eps, scaling=scaling, backend="tensorized")

    def fit_transform(self, X, mass, verbose=True, report_interval=500, X_true=None):
        """
        Imputes missing values using a batched OT loss

        Parameters
        ----------
        X : torch.DoubleTensor or torch.cuda.DoubleTensor
            Contains non-missing and missing data at the indices given by the
            "mask" argument. Missing values can be arbitrarily assigned
            (e.g. with NaNs).

        mask : torch.DoubleTensor or torch.cuda.DoubleTensor
            mask[i,j] == 1 if X[i,j] is missing, else mask[i,j] == 0.

        verbose: bool, default=True
            If True, output loss to log during iterations.

        X_true: torch.DoubleTensor or None, default=None
            Ground truth for the missing values. If provided, will output a
            validation score during training, and return score arrays.
            For validation/debugging only.

        Returns
        -------
        X_filled: torch.DoubleTensor or torch.cuda.DoubleTensor
            Imputed missing data (plus unchanged non-missing data).


        """

        X = X.clone()
        n, d = X.shape
        mass = mass.clone()
        mass = mass / mass.sum()
        
        if self.batchsize > n // 2:
            e = int(np.log2(n // 2))
            self.batchsize = 2**e
            if verbose:
                logging.info(f"Batchsize larger that half size = {len(X) // 2}. Setting batchsize to {self.batchsize}.")

        mask = torch.isnan(X).double()
        imps = (self.noise * torch.randn(mask.shape).double() + nanmean(X, 0))[mask.bool()]
        imps.requires_grad = True

        optimizer = self.opt([imps], lr=self.lr)

        if verbose:
            logging.info(f"batchsize = {self.batchsize}, epsilon = {self.eps:.4f}")

        if X_true is not None:
            maes = np.zeros(self.niter)
            rmses = np.zeros(self.niter)

        for i in range(self.niter):

            X_filled = X.detach().clone()
            X_filled[mask.bool()] = imps
            loss = 0
            
            for _ in range(self.n_pairs):

                idx1 = np.random.choice(n, self.batchsize, replace=False)
                idx2 = np.random.choice(n, self.batchsize, replace=False)
    
                X1 = X_filled[idx1]
                X2 = X_filled[idx2]
                
                mass1 = mass[idx1]
                mass1 = mass1 / mass1.sum()
                
                mass2 = mass[idx2]
                mass2 = mass2 / mass2.sum()
    
                loss = loss + self.sk(mass1, X1, mass2, X2)

            if torch.isnan(loss).any() or torch.isinf(loss).any():
                ### Catch numerical errors/overflows (should not happen)
                logging.info("Nan or inf loss")
                break

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            if X_true is not None:
                maes[i] = MAE(X_filled, X_true, mask).item()
                rmses[i] = RMSE(X_filled, X_true, mask).item()

            if verbose and (i % report_interval == 0):
                if X_true is not None:
                    logging.info(f'Iteration {i}:\t Loss: {loss.item() / self.n_pairs:.4f}\t '
                                 f'Validation MAE: {maes[i]:.4f}\t'
                                 f'RMSE: {rmses[i]:.4f}')
                else:
                    logging.info(f'Iteration {i}:\t Loss: {loss.item() / self.n_pairs:.4f}')

        X_filled = X.detach().clone()
        X_filled[mask.bool()] = imps

        if X_true is not None:
            return X_filled, maes, rmses
        else:
            return X_filled


class RRimputer():
    """
    Round-Robin imputer with a batch sinkhorn loss

    Parameters
    ----------
    models: iterable
        iterable of torch.nn.Module. The j-th model is used to predict the j-th
        variable using all others.

    eps: float, default=0.01
        Sinkhorn regularization parameter.
        
    lr : float, default = 0.01
        Learning rate.

    opt: torch.nn.optim.Optimizer, default=torch.optim.Adam
        Optimizer class to use for fitting.
        
    max_iter : int, default=10
        Maximum number of round-robin cycles for imputation.

    niter : int, default=15
        Number of gradient updates for each model within a cycle.

    batchsize : int, defatul=128
        Size of the batches on which the sinkhorn divergence is evaluated.

    n_pairs : int, default=10
        Number of batch pairs used per gradient update.

    tol : float, default = 0.001
        Tolerance threshold for the stopping criterion.

    weight_decay : float, default = 1e-5
        L2 regularization magnitude.

    order : str, default="random"
        Order in which the variables are imputed.
        Valid values: {"random" or "increasing"}.

    unsymmetrize: bool, default=True
        If True, sample one batch with no missing 
        data in each pair during training.

    scaling: float, default=0.9
        Scaling parameter in Sinkhorn iterations
        c.f. geomloss' doc: "Allows you to specify the trade-off between
        speed (scaling < .4) and accuracy (scaling > .9)"

    """
    def __init__(self,
                 models, 
                 power = 2,
                 eps= 0.01, 
                 lr=1e-2, 
                 opt=torch.optim.Adam, 
                 max_iter=10,
                 niter=15, 
                 batchsize=128,
                 n_pairs=10, 
                 tol=1e-3,
                 noise=0.1,
                 weight_decay=1e-5, 
                 order='random',
                 unsymmetrize=True, 
                 scaling=.9):

        self.models = models
        self.sk = SamplesLoss("sinkhorn", p=power, blur=eps,
                              scaling=scaling, backend="auto")
        self.lr = lr
        self.opt = opt
        self.max_iter = max_iter
        self.niter = niter
        self.batchsize = batchsize
        self.n_pairs = n_pairs
        self.tol = tol
        self.noise = noise
        self.weight_decay=weight_decay
        self.order=order
        self.unsymmetrize = unsymmetrize

        self.is_fitted = False

    def fit_transform(self, X, mass, verbose=True,
                      report_interval=1, X_true=None):
        """
        Fits the imputer on a dataset with missing data, and returns the
        imputations.

        Parameters
        ----------
        X : torch.DoubleTensor or torch.cuda.DoubleTensor, shape (n, d)
            Contains non-missing and missing data at the indices given by the
            "mask" argument. Missing values can be arbitrarily assigned 
            (e.g. with NaNs).

        mask : torch.DoubleTensor or torch.cuda.DoubleTensor, shape (n, d)
            mask[i,j] == 1 if X[i,j] is missing, else mask[i,j] == 0.

        verbose : bool, default=True
            If True, output loss to log during iterations.
            
        report_interval : int, default=1
            Interval between loss reports (if verbose).

        X_true: torch.DoubleTensor or None, default=None
            Ground truth for the missing values. If provided, will output a 
            validation score during training. For debugging only.

        Returns
        -------
        X_filled: torch.DoubleTensor or torch.cuda.DoubleTensor
            Imputed missing data (plus unchanged non-missing data).

        """

        X = X.clone()
        n, d = X.shape
        mask = torch.isnan(X).double()
        normalized_tol = self.tol * torch.max(torch.abs(X[~mask.bool()]))
        mass = mass.clone()
        mass = mass / mass.sum()

        if self.batchsize > n // 2:
            e = int(np.log2(n // 2))
            self.batchsize = 2**e
            if verbose:
                logging.info(f"Batchsize larger that half size = {len(X) // 2}."
                             f" Setting batchsize to {self.batchsize}.")

        order_ = torch.argsort(mask.sum(0))

        optimizers = [self.opt(self.models[i].parameters(),
                               lr=self.lr, weight_decay=self.weight_decay) for i in range(d)]

        imps = (self.noise * torch.randn(mask.shape).double() + nanmean(X, 0))[mask.bool()]
        X[mask.bool()] = imps
        X_filled = X.clone()

        if X_true is not None:
            maes = np.zeros(self.max_iter)
            rmses = np.zeros(self.max_iter)

        for i in range(self.max_iter):

            if self.order == 'random':
                order_ = np.random.choice(d, d, replace=False)
            X_old = X_filled.clone().detach()

            loss = 0

            for l in range(d):
                j = order_[l].item()
                n_not_miss = (~mask[:, j].bool()).sum().item()

                if n - n_not_miss == 0:
                    continue  # no missing value on that coordinate

                for k in range(self.niter):

                    loss = 0

                    X_filled = X_filled.detach()
                    X_filled[mask[:, j].bool(), j] = self.models[j](X_filled[mask[:, j].bool(), :][:, np.r_[0:j, j+1: d]]).squeeze()

                    for i in range(self.n_pairs):
                        
                        idx1 = np.random.choice(n, self.batchsize, replace=False)
                        X1 = X_filled[idx1]
                        mass1 = mass[idx1]
                        mass1 = mass1/mass1.sum()

                        if self.unsymmetrize:
                            n_miss = (~mask[:, j].bool()).sum().item()
                            idx2 = np.random.choice(n_miss, self.batchsize, replace= self.batchsize > n_miss)
                            X2 = X_filled[~mask[:, j].bool(), :][idx2]
                            mass2 = mass[idx2]
                            mass2 = mass2/mass2.sum()

                        else:
                            idx2 = np.random.choice(n, self.batchsize, replace=False)
                            X2 = X_filled[idx2]

                        loss = loss + self.sk(mass1, X1, 
                                              mass2, X2)

                    optimizers[j].zero_grad()
                    loss.backward()
                    optimizers[j].step()

                # Impute with last parameters
                with torch.no_grad():
                    X_filled[mask[:, j].bool(), j] = self.models[j](X_filled[mask[:, j].bool(), :][:, np.r_[0:j, j+1: d]]).squeeze()

            if X_true is not None:
                maes[i] = MAE(X_filled, X_true, mask).item()
                rmses[i] = RMSE(X_filled, X_true, mask).item()

            if verbose and (i % report_interval == 0):
                if X_true is not None:
                    logging.info(f'Iteration {i}:\t Loss: {loss.item() / self.n_pairs:.4f}\t'
                                 f'Validation MAE: {maes[i]:.4f}\t'
                                 f'RMSE: {rmses[i]: .4f}')
                else:
                    logging.info(f'Iteration {i}:\t Loss: {loss.item() / self.n_pairs:.4f}')

            if torch.norm(X_filled - X_old, p=np.inf) < normalized_tol:
                break

        if i == (self.max_iter - 1) and verbose:
            logging.info('Early stopping criterion not reached')

        self.is_fitted = True

        if X_true is not None:
            return X_filled, maes, rmses
        else:
            return X_filled

    def transform(self, X, mask, verbose=True, report_interval=1, X_true=None):
        """
        Impute missing values on new data. Assumes models have been previously 
        fitted on other data.
        
        Parameters
        ----------
        X : torch.DoubleTensor or torch.cuda.DoubleTensor, shape (n, d)
            Contains non-missing and missing data at the indices given by the
            "mask" argument. Missing values can be arbitrarily assigned 
            (e.g. with NaNs).

        mask : torch.DoubleTensor or torch.cuda.DoubleTensor, shape (n, d)
            mask[i,j] == 1 if X[i,j] is missing, else mask[i,j] == 0.

        verbose: bool, default=True
            If True, output loss to log during iterations.
            
        report_interval : int, default=1
            Interval between loss reports (if verbose).

        X_true: torch.DoubleTensor or None, default=None
            Ground truth for the missing values. If provided, will output a 
            validation score during training. For debugging only.

        Returns
        -------
        X_filled: torch.DoubleTensor or torch.cuda.DoubleTensor
            Imputed missing data (plus unchanged non-missing data).

        """

        assert self.is_fitted, "The model has not been fitted yet."

        n, d = X.shape
        normalized_tol = self.tol * torch.max(torch.abs(X[~mask.bool()]))

        order_ = torch.argsort(mask.sum(0))

        X[mask] = nanmean(X)
        X_filled = X.clone()

        for i in range(self.max_iter):

            if self.order == 'random':
                order_ = np.random.choice(d, d, replace=False)
            X_old = X_filled.clone().detach()

            for l in range(d):

                j = order_[l].item()

                with torch.no_grad():
                    X_filled[mask[:, j].bool(), j] = self.models[j](X_filled[mask[:, j].bool(), :][:, np.r_[0:j, j+1: d]]).squeeze()

            if verbose and (i % report_interval == 0):
                if X_true is not None:
                    logging.info(f'Iteration {i}:\t '
                                 f'Validation MAE: {MAE(X_filled, X_true, mask).item():.4f}\t'
                                 f'RMSE: {RMSE(X_filled, X_true, mask).item():.4f}')

            if torch.norm(X_filled - X_old, p=np.inf) < normalized_tol:
                break

        if i == (self.max_iter - 1) and verbose:
            logging.info('Early stopping criterion not reached')

        return X_filled

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# https://github.com/BorisMuzellec/MissingDataOT/blob/master/utils.py

import torch
import numpy as np

from scipy import optimize

def nanmean(v, *args, **kwargs):
    """
    A Pytorch version on Numpy's nanmean
    """
    v = v.clone()
    is_nan = torch.isnan(v)
    v[is_nan] = 0
    return v.sum(*args, **kwargs) / (~is_nan).float().sum(*args, **kwargs)


#### Quantile ######
def quantile(X, q, dim=None):
    """
    Returns the q-th quantile.

    Parameters
    ----------
    X : torch.DoubleTensor or torch.cuda.DoubleTensor, shape (n, d)
        Input data.

    q : float
        Quantile level (starting from lower values).

    dim : int or None, default = None
        Dimension allong which to compute quantiles. If None, the tensor is flattened and one value is returned.


    Returns
    -------
        quantiles : torch.DoubleTensor

    """
    return X.kthvalue(int(q * len(X)), dim=dim)[0]


#### Automatic selection of the regularization parameter ####
def pick_epsilon(X, quant=0.5, mult=0.05, max_points=2000):
    """
        Returns a quantile (times a multiplier) of the halved pairwise squared distances in X.
        Used to select a regularization parameter for Sinkhorn distances.

    Parameters
    ----------
    X : torch.DoubleTensor or torch.cuda.DoubleTensor, shape (n, d)
        Input data on which distances will be computed.

    quant : float, default = 0.5
        Quantile to return (default is median).

    mult : float, default = 0.05
        Mutiplier to apply to the quantiles.

    max_points : int, default = 2000
        If the length of X is larger than max_points, estimate the quantile on a random subset of size max_points to
        avoid memory overloads.

    Returns
    -------
        epsilon: float

    """
    means = nanmean(X, 0)
    X_ = X.clone()
    mask = torch.isnan(X_)
    X_[mask] = (mask * means)[mask]

    idx = np.random.choice(len(X_), min(max_points, len(X_)), replace=False)
    X = X_[idx]
    dists = ((X[:, None] - X) ** 2).sum(2).flatten() / 2.
    dists = dists[dists > 0]

    return quantile(dists, quant, 0).item() * mult


#### Accuracy Metrics ####
def MAE(X, X_true, mask):
    """
    Mean Absolute Error (MAE) between imputed variables and ground truth. Pytorch/Numpy agnostic
    
    Parameters
    ----------
    X : torch.DoubleTensor or np.ndarray, shape (n, d)
        Data with imputed variables.

    X_true : torch.DoubleTensor or np.ndarray, shape (n, d)
        Ground truth.

    mask : torch.BoolTensor or np.ndarray of booleans, shape (n, d)
        Missing value mask (missing if True)

    Returns
    -------
        MAE : float

    """
    if torch.is_tensor(mask):
        mask_ = mask.bool()
        return torch.abs(X[mask_] - X_true[mask_]).sum() / mask_.sum()
    else: # should be an ndarray
        mask_ = mask.astype(bool)
        return np.absolute(X[mask_] - X_true[mask_]).sum() / mask_.sum()



def RMSE(X, X_true, mask):
    """
    Root Mean Squared Error (MAE) between imputed variables and ground truth. Pytorch/Numpy agnostic

    Parameters
    ----------
    X : torch.DoubleTensor or np.ndarray, shape (n, d)
        Data with imputed variables.

    X_true : torch.DoubleTensor or np.ndarray, shape (n, d)
        Ground truth.

    mask : torch.BoolTensor or np.ndarray of booleans, shape (n, d)
        Missing value mask (missing if True)

    Returns
    -------
        RMSE : float

    """
    if torch.is_tensor(mask):
        mask_ = mask.bool()
        return (((X[mask_] - X_true[mask_]) ** 2).sum() / mask_.sum()).sqrt()
    else: # should be an ndarray
        mask_ = mask.astype(bool)
        return np.sqrt(((X[mask_] - X_true[mask_])**2).sum() / mask_.sum())

##################### MISSING DATA MECHANISMS #############################

##### Missing At Random ######

def MAR_mask(X, p, p_obs):
    """
    Missing at random mechanism with a logistic masking model. First, a subset of variables with *no* missing values is
    randomly selected. The remaining variables have missing values according to a logistic model with random weights,
    re-scaled so as to attain the desired proportion of missing values on those variables.

    Parameters
    ----------
    X : torch.DoubleTensor or np.ndarray, shape (n, d)
        Data for which missing values will be simulated. If a numpy array is provided,
        it will be converted to a pytorch tensor.

    p : float
        Proportion of missing values to generate for variables which will have missing values.

    p_obs : float
        Proportion of variables with *no* missing values that will be used for the logistic masking model.

    Returns
    -------
    mask : torch.BoolTensor or np.ndarray (depending on type of X)
        Mask of generated missing values (True if the value is missing).

    """

    n, d = X.shape

    to_torch = torch.is_tensor(X) ## output a pytorch tensor, or a numpy array
    if not to_torch:
        X = torch.from_numpy(X)

    mask = torch.zeros(n, d).bool() if to_torch else np.zeros((n, d)).astype(bool)

    d_obs = max(int(p_obs * d), 1) ## number of variables that will have no missing values (at least one variable)
    d_na = d - d_obs ## number of variables that will have missing values

    ### Sample variables that will all be observed, and those with missing values:
    idxs_obs = np.random.choice(d, d_obs, replace=False)
    idxs_nas = np.array([i for i in range(d) if i not in idxs_obs])

    ### Other variables will have NA proportions that depend on those observed variables, through a logistic model
    ### The parameters of this logistic model are random.

    ### Pick coefficients so that W^Tx has unit variance (avoids shrinking)
    coeffs = pick_coeffs(X, idxs_obs, idxs_nas)
    ### Pick the intercepts to have a desired amount of missing values
    intercepts = fit_intercepts(X[:, idxs_obs], coeffs, p)

    ps = torch.sigmoid(X[:, idxs_obs].mm(coeffs) + intercepts)

    ber = torch.rand(n, d_na)
    mask[:, idxs_nas] = ber < ps

    return mask

##### Missing not at random ######

def MNAR_mask_logistic(X, p, p_params =.3, exclude_inputs=True):
    """
    Missing not at random mechanism with a logistic masking model. It implements two mechanisms:
    (i) Missing probabilities are selected with a logistic model, taking all variables as inputs. Hence, values that are
    inputs can also be missing.
    (ii) Variables are split into a set of intputs for a logistic model, and a set whose missing probabilities are
    determined by the logistic model. Then inputs are then masked MCAR (hence, missing values from the second set will
    depend on masked values.
    In either case, weights are random and the intercept is selected to attain the desired proportion of missing values.

    Parameters
    ----------
    X : torch.DoubleTensor or np.ndarray, shape (n, d)
        Data for which missing values will be simulated.
        If a numpy array is provided, it will be converted to a pytorch tensor.

    p : float
        Proportion of missing values to generate for variables which will have missing values.

    p_params : float
        Proportion of variables that will be used for the logistic masking model (only if exclude_inputs).

    exclude_inputs : boolean, default=True
        True: mechanism (ii) is used, False: (i)

    Returns
    -------
    mask : torch.BoolTensor or np.ndarray (depending on type of X)
        Mask of generated missing values (True if the value is missing).

    """

    n, d = X.shape

    to_torch = torch.is_tensor(X) ## output a pytorch tensor, or a numpy array
    if not to_torch:
        X = torch.from_numpy(X)

    mask = torch.zeros(n, d).bool() if to_torch else np.zeros((n, d)).astype(bool)

    d_params = max(int(p_params * d), 1) if exclude_inputs else d ## number of variables used as inputs (at least 1)
    d_na = d - d_params if exclude_inputs else d ## number of variables masked with the logistic model

    ### Sample variables that will be parameters for the logistic regression:
    idxs_params = np.random.choice(d, d_params, replace=False) if exclude_inputs else np.arange(d)
    idxs_nas = np.array([i for i in range(d) if i not in idxs_params]) if exclude_inputs else np.arange(d)

    ### Other variables will have NA proportions selected by a logistic model
    ### The parameters of this logistic model are random.

    ### Pick coefficients so that W^Tx has unit variance (avoids shrinking)
    coeffs = pick_coeffs(X, idxs_params, idxs_nas)
    ### Pick the intercepts to have a desired amount of missing values
    intercepts = fit_intercepts(X[:, idxs_params], coeffs, p)

    ps = torch.sigmoid(X[:, idxs_params].mm(coeffs) + intercepts)

    ber = torch.rand(n, d_na)
    mask[:, idxs_nas] = ber < ps

    ## If the inputs of the logistic model are excluded from MNAR missingness,
    ## mask some values used in the logistic model at random.
    ## This makes the missingness of other variables potentially dependent on masked values

    if exclude_inputs:
        mask[:, idxs_params] = torch.rand(n, d_params) < p

    return mask

def MNAR_self_mask_logistic(X, p):
    """
    Missing not at random mechanism with a logistic self-masking model. Variables have missing values probabilities
    given by a logistic model, taking the same variable as input (hence, missingness is independent from one variable
    to another). The intercepts are selected to attain the desired missing rate.

    Parameters
    ----------
    X : torch.DoubleTensor or np.ndarray, shape (n, d)
        Data for which missing values will be simulated.
        If a numpy array is provided, it will be converted to a pytorch tensor.

    p : float
        Proportion of missing values to generate for variables which will have missing values.

    Returns
    -------
    mask : torch.BoolTensor or np.ndarray (depending on type of X)
        Mask of generated missing values (True if the value is missing).

    """

    n, d = X.shape

    to_torch = torch.is_tensor(X) ## output a pytorch tensor, or a numpy array
    if not to_torch:
        X = torch.from_numpy(X)

    ### Variables will have NA proportions that depend on those observed variables, through a logistic model
    ### The parameters of this logistic model are random.

    ### Pick coefficients so that W^Tx has unit variance (avoids shrinking)
    coeffs = pick_coeffs(X, self_mask=True)
    ### Pick the intercepts to have a desired amount of missing values
    intercepts = fit_intercepts(X, coeffs, p, self_mask=True)

    ps = torch.sigmoid(X * coeffs + intercepts)

    ber = torch.rand(n, d) if to_torch else np.random.rand(n, d)
    mask = ber < ps if to_torch else ber < ps.numpy()

    return mask


def MNAR_mask_quantiles(X, p, q, p_params, cut='both', MCAR=False):
    """
    Missing not at random mechanism with quantile censorship. First, a subset of variables which will have missing
    variables is randomly selected. Then, missing values are generated on the q-quantiles at random. Since
    missingness depends on quantile information, it depends on masked values, hence this is a MNAR mechanism.

    Parameters
    ----------
    X : torch.DoubleTensor or np.ndarray, shape (n, d)
        Data for which missing values will be simulated.
        If a numpy array is provided, it will be converted to a pytorch tensor.

    p : float
        Proportion of missing values to generate for variables which will have missing values.

    q : float
        Quantile level at which the cuts should occur

    p_params : float
        Proportion of variables that will have missing values

    cut : 'both', 'upper' or 'lower', default = 'both'
        Where the cut should be applied. For instance, if q=0.25 and cut='upper', then missing values will be generated
        in the upper quartiles of selected variables.
        
    MCAR : bool, default = True
        If true, masks variables that were not selected for quantile censorship with a MCAR mechanism.
        
    Returns
    -------
    mask : torch.BoolTensor or np.ndarray (depending on type of X)
        Mask of generated missing values (True if the value is missing).

    """
    n, d = X.shape

    to_torch = torch.is_tensor(X) ## output a pytorch tensor, or a numpy array
    if not to_torch:
        X = torch.from_numpy(X)

    mask = torch.zeros(n, d).bool() if to_torch else np.zeros((n, d)).astype(bool)

    d_na = max(int(p_params * d), 1) ## number of variables that will have NMAR values

    ### Sample variables that will have imps at the extremes
    idxs_na = np.random.choice(d, d_na, replace=False) ### select at least one variable with missing values

    ### check if values are greater/smaller that corresponding quantiles
    if cut == 'upper':
        quants = quantile(X[:, idxs_na], 1-q, dim=0)
        m = X[:, idxs_na] >= quants
    elif cut == 'lower':
        quants = quantile(X[:, idxs_na], q, dim=0)
        m = X[:, idxs_na] <= quants
    elif cut == 'both':
        u_quants = quantile(X[:, idxs_na], 1-q, dim=0)
        l_quants = quantile(X[:, idxs_na], q, dim=0)
        m = (X[:, idxs_na] <= l_quants) | (X[:, idxs_na] >= u_quants)

    ### Hide some values exceeding quantiles
    ber = torch.rand(n, d_na)
    mask[:, idxs_na] = (ber < p) & m

    if MCAR:
    ## Add a mcar mecanism on top
        mask = mask | (torch.rand(n, d) < p)

    return mask


def pick_coeffs(X, idxs_obs=None, idxs_nas=None, self_mask=False):
    n, d = X.shape
    if self_mask:
        coeffs = torch.randn(d)
        Wx = X * coeffs
        coeffs /= torch.std(Wx, 0)
    else:
        d_obs = len(idxs_obs)
        d_na = len(idxs_nas)
        coeffs = torch.randn(d_obs, d_na)
        Wx = X[:, idxs_obs].mm(coeffs)
        coeffs /= torch.std(Wx, 0, keepdim=True)
    return coeffs


def fit_intercepts(X, coeffs, p, self_mask=False):
    if self_mask:
        d = len(coeffs)
        intercepts = torch.zeros(d)
        for j in range(d):
            def f(x):
                return torch.sigmoid(X * coeffs[j] + x).mean().item() - p
            intercepts[j] = optimize.bisect(f, -50, 50)
    else:
        d_obs, d_na = coeffs.shape
        intercepts = torch.zeros(d_na)
        for j in range(d_na):
            def f(x):
                return torch.sigmoid(X.mv(coeffs[:, j]) + x).mean().item() - p
            intercepts[j] = optimize.bisect(f, -50, 50)
    return intercepts


