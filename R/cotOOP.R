# COT general form using object oriented code and R6 classes

#' @name Measure_
#' @title An R6 object for measures
#' @rdname Measure_-class
#' @description Internal R6 class object for Measure objects
#' @keywords internal
Measure_ <- R6::R6Class("Measure", # change name later for Roxygen purposes
 public = {list(
   # data
   
   #' @field balance_functions the functions of the data that 
   #' we want to adjust towards the targets
   balance_functions = "torch_tensor",
   
   #' @field balance_target the values the balance_functions are targeting
   balance_target = "vector",
   
   # values
   #' @field adapt What aspect of the data will be adapted. One of "none","weights", or "x".
   adapt  = "character",
   
   #' @field device the [torch::torch_device()] of the data.
   device = "torch_device",
   
   #' @field dtype the [torch::torch_dtype] of the data.
   dtype  = "torch_dtype",
   
   #' @field n the rows of the covariates, x.
   n      = "integer",
   
   #' @field d the columns of the covariates, x.
   d      = "integer",
   
   #' @field probability_measure is the measure a probability measure?
   probability_measure = "logical",
   
   # functions
   #' @description 
   #' generates a deep clone of the object without gradients.
   detach = function() { #removes gradient calculation in copy
     orig_adapt <- self$adapt
     
     if (orig_adapt != "none") {
       if (orig_adapt == "weights") {
         private$mass_$requires_grad <- FALSE
       } else if (orig_adapt == "x") {
         private$data_$requires_grad <- FALSE
       }
     }
     
     temp_obj <- self$clone(deep = TRUE)
     
     if (orig_adapt != "none") {
       temp_obj$adapt <- "none"
       if (orig_adapt == "weights") {
         private$mass_$requires_grad <- TRUE
       } else if (orig_adapt == "x") {
         private$data_$requires_grad <- TRUE
       }
     }
     return(temp_obj)
   },
   
   #' @description 
   #' Makes a copy of the weights parameters.
   get_weight_parameters = function() {
     private$mass_$clone()
   },
   
   #' @description prints the measure object
   #' @param ... Not used
   print = function(...) {
     cat("Measure: ",  rlang::obj_address(self), "\n", sep = "")
     cat("  x      : a ", self$n, "x", self$d, " matrix \n", sep = "")
     if(self$d > 5) {
     cat("           " , paste(round(as_matrix(private$data_[1,1:5]), digits = 2), collapse = ", "), ", \u2026", "\n", sep = "")
       if(self$n > 1) cat("           " , paste(round(as_matrix(private$data_[2,1:5]), digits = 2), collapse = ", "), ", \u2026", "\n", sep = "")
       if(self$n > 2) cat("           " , paste(round(as_matrix(private$data_[3,1:5]), digits = 2), collapse = ", "), ", \u2026", "\n", sep = "") 
       if(self$n > 3) cat("           " , paste(round(as_matrix(private$data_[4,1:5]), digits = 2), collapse = ", "), ", \u2026", "\n", sep = "")   
       if(self$n > 4) cat("           " , paste(round(as_matrix(private$data_[5,1:5]), digits = 2), collapse = ", "), ", \u2026", "\n", sep = "")   
     } else {
     cat("           " , paste(round(as_matrix(private$data_[1,1:5]), digits = 2), collapse = ", "), "\n", sep = "")
       if(self$n > 1) cat("           " , paste(round(as_matrix(private$data_[2,1:5]), digits = 2), collapse = ", "),  "\n", sep = "")
       if(self$n > 2) cat("           " , paste(round(as_matrix(private$data_[3,1:5]), digits = 2), collapse = ", "),  "\n", sep = "") 
       if(self$n > 3) cat("           " , paste(round(as_matrix(private$data_[4,1:5]), digits = 2), collapse = ", "),  "\n", sep = "")   
       if(self$n > 4) cat("           " , paste(round(as_matrix(private$data_[5,1:5]), digits = 2), collapse = ", "),  "\n", sep = "")  
     }
     if(self$n > 5) {
     cat("           \u22EE \n")
     # cat("           .\n")
     # cat("           .\n")
     }
     if(self$n > 5) { 
     cat("  weights: ", paste(round(t(as_numeric(self$weights[1:5])), digits = 2), collapse = ", "), ", \u2026", "\n", sep = "")
     } else {
     cat("  weights: ", paste(round(t(as_numeric(self$weights[1:5])), digits = 2), collapse = ", "), "\n", sep = "")
     }
     if(all(as_logical(!is.na(self$balance_target)))) {
     cat("  balance: ", "\n", sep = "")
     cat("   funct.: ",  paste(round(as_matrix(self$balance_functions[1,1:5]), digits = 2), collapse = ", "), "\n", sep = "")
     cat("   target: ",  paste(round(t(as_numeric(self$balance_target[1:5])), digits = 2), collapse = ", "), " \u2026", "\n", sep = "")
     }
     cat("  adapt  : ", self$adapt, "\n", sep = "")
     cat("  dtype  : ", capture.output(self$dtype), "\n", sep = "")
     cat("  device : ", capture.output(self$device), "\n", sep = "")
   },
   
   #' @description Constructor function
   #' @param x The data points
   #' @param weights The empirical measure. If NULL, assigns equal weight to each observation
   #' @param probability.measure Is the empirical measure a probability measure? Default is TRUE.
   #' @param adapt Should we try to adapt the data ("x"), the weights ("weights"), or neither ("none"). Default is "none".
   #' @param balance.functions A matrix of functions of the covariates to target for mean balance. If NULL and `target.values` are provided, will use the data in `x`.
   #' @param target.values The targets for the balance functions. Should be the same length as columns in `balance.functions.`
   #' @param dtype The [torch::torch_dtype] or NULL.
   #' @param device The device to have the data on. Should be result of [torch::torch_device()] or NULL.
   initialize = function(x, weights = NULL, 
                         probability.measure = TRUE, 
                         adapt = c("none","weights", "x"), 
                         balance.functions = NA_real_,
                         target.values = NA_real_,
                         dtype = NULL, device = NULL) {
     # browser()
     if(!is.matrix(x)) x <- as.matrix(x)
     self$d <- ncol(x)
     self$n <- nrow(x)
     
     self$adapt <- match.arg(adapt)
     self$probability_measure <- isTRUE(probability.measure)
     
     
     self$device <- cuda_device_check(device)
     self$dtype  <- cuda_dtype_check(dtype, self$device)
     # device
     # if (is.null(device)) {
     #   cuda_opt <- torch::cuda_is_available() && torch::cuda_device_count() >= 1
     #   if (cuda_opt) {
     #     self$device <-  torch::torch_device("cuda")
     #   } else {
     #     self$device <-  torch::torch_device("cpu")
     #   }
     # } else {
     #   stopifnot("device argument must be NULL or an object of class 'torch_device'" = torch::is_torch_device(device))
     #   self$device <- device
     # }
     # 
     # #dtype
     # if ( is.null(dtype) ) {
     #   if (grepl("cuda", capture.output(print(self$device)) ) ) {
     #     dtype <- torch::torch_float()
     #   } else {
     #     dtype <- torch::torch_double()
     #   }
     # }
     # stopifnot("Argument 'dtype' must be of class 'torch_dtype'. Please see '?torch_dtype' for more info." = torch::is_torch_dtype(dtype))
     # 
     # set data
     private$data_ <- torch::torch_tensor(x, dtype = self$dtype, device = self$device)$contiguous()
     
     if (self$adapt == "x") {
       if(!private$data_$requires_grad) private$data_$requires_grad <- TRUE
     }
     
     # browser()
     self$balance_target <- as.numeric(target.values)
     if(self$adapt == "none") self$balance_target <- NA_real_
     if(missing(balance.functions)) balance.functions <- NA_real_
     if(all(is.na(self$balance_target)) || self$adapt == "none") balance.functions <- NA_real_
     if(all(is.na(balance.functions)) && !all(is.na(self$balance_target))) balance.functions <- private$data_
     # browser()
     if (!inherits(balance.functions, "torch_tensor") && !all(is.na(balance.functions)) ) {
       self$balance_functions <- torch::torch_tensor(as.matrix(balance.functions), 
                                                     dtype = self$dtype, device = self$device)$contiguous()
     } else {
       self$balance_functions <- balance.functions
     }
     if (!inherits( self$balance_target, "torch_tensor") && !all(is.na(self$balance_target)) )  {
       self$balance_target <- torch::torch_tensor(self$balance_target, 
                                                     dtype = self$dtype, device = self$device)$contiguous()
     } else if (!all(is.na(self$balance_target))) {
       self$balance_target <- self$balance_target$to(device = self$device, dtype = self$dtype)$contiguous()
     }
     
     if (self$adapt == "x") {
       if (!all(is.na(self$balance_functions)) && !self$balance_functions$requires_grad) self$balance_functions$requires_grad <- TRUE
       if (!all(is.na(self$balance_functions))) {
         if(rlang::obj_address(private$data_) != rlang::obj_address(self$balance_functions)) warning("x and balance.functions may not come from same data. Adapting the two together may not make sense. If you would like x and balance.functions to have the same underlying data, you can either feed the same torch tensor into both or leave balance.functions blank and it will use the values in data by default if target.values are supplied.")
       }
     }
     
     stopifnot("Argument 'target.values' must be NA or the same length as the number of columns in 'balance.functions'." = all(is.na(self$balance_target)) || length(self$balance_target) == ncol(self$balance_functions))
 
     # adjust be Std Dev
     if (inherits( self$balance_target, "torch_tensor") && !as.logical(all(is.na(self$balance_target$to(device = "cpu")))) ) {
       if (!self$adapt == "x") {
         sds <- self$balance_functions$std(1)
         self$balance_functions <- self$balance_functions/sds
         self$balance_target <- self$balance_target /sds
         non_zero_sd <-  as.logical(0 < sds$to(device = "cpu"))
         if (sum(non_zero_sd) == 0) {
           warning("All columns of balance.functions have zero variance. Balance functions not used.")
           self$balance_functions <- NA_real_
           self$balance_target <- NA_real_
         } else {
           sel_nz <- which(non_zero_sd)
           self$balance_functions <- self$balance_functions[,sel_nz, drop = FALSE]
           self$balance_target <- self$balance_target[sel_nz]
         }
         
       }
     }
     
     weights <- check_weights(weights, private$data_, self$device)
     
     if (self$adapt == "weights"){
       # browser()
       private$mass_ <- torch::torch_zeros(length(weights)-1L,
                                            dtype = self$dtype,
                                           device = self$device)
       private$get_mass_ <- function() {
         private$transform_(private$mass_)
       }
       private$assign_mass_ <- function(value) {
         # torch::autograd_set_grad_mode(enabled = FALSE)
         torch::with_no_grad({
         private$mass_$copy_(private$inv_transform_(value)$to(device = self$device))
         if(! torch::is_undefined_tensor(private$mass_$grad) ) private$mass_$grad$copy_(0.0)
         })
         # torch::autograd_set_grad_mode(enabled = TRUE)
       }
       
       private$mass_$requires_grad <- TRUE
     } else {
       private$get_mass_ <- function() {
         return(private$mass_)
       }
       private$assign_mass_ <- function(value) {
         private$mass_ <- value$to(device = self$device)
       }
     }
     
     if (self$probability_measure) {
       private$inv_transform_ <- function(value) {
         min_neg <- round(log(.Machine$double.xmin)) - 50.0
         logs <- torch::torch_log(value/sum(value))
         logs[logs< min_neg] <- min_neg
         logs <- logs[2:length(logs)] - logs[1]
         return(logs)
       }
       
       private$transform_ <- function(value) {
         full_param <- torch::torch_cat(
           list( torch::torch_zeros(1L, device = self$device, dtype = value$dtype),
                value)
         )
         return(full_param$log_softmax(1)$exp())
       }
       
     } else {
       private$inv_transform_ <- torch::torch_log
       private$transform_ <- torch::torch_exp
     }
     private$assign_mass_(weights)
     
     # assign initial values
     private$init_weights_ <- self$weights$detach()
     if (self$adapt == "none" || self$adapt == "weights") {
       private$init_data_ <- private$data_
     } else if (self$adapt == "x") {
       private$init_data_ <- private$data_$detach()$clone()
     } else {
       stop("adapt hasn't been properly set!")
     }
     if (torch::cuda_is_available()) torch::cuda_empty_cache()
     return(invisible(self))
   } 
 )},
 active = {list(
   #' @field grad gets or sets gradient
   grad = function(value) {
     
     # return grad values as appropriate
     if (missing(value)) {
       if(self$adapt == "none") {
         return(NULL)
       } else if (self$adapt == "x") {
         return(private$data_$grad)
       } else if (self$adapt == "weights") {
         return(private$mass_$grad)
       }
     }
     
     
     # save grad values
     torch::with_no_grad({
       if (!is_torch_tensor(value)) {
         value <- torch::torch_tensor(value, dtype = self$dtype,
                                      device = self$device)
       }
       if(self$adapt == "none") {
         stop("No elements of this Measure object have gradients")
       } else if (self$adapt == "x") {
         
         l_v <- length(value)
         dim_v <- dim(value)
         
         if(any(dim_v != c(self$n,self$d)) && l_v != self$d && l_v != self$d * self$n) {
           stop(sprintf("Length of input for the `x` gradients must be of length %s or %s. Alternatively, better to supply a matrix of dimension %s by %s directly.", self$d, self$d*self$n, self$n, self$d))
         }
         
         private$data_$grad <- value$to(device = private$data_$device)
         
       } else if (self$adapt == "weights") {
         
         l_v <- length(value)
         
         if (l_v != (self$n - 1)) {
           stop(sprintf("Input value must be of length %s for the weight gradients. The first value is fixed to make the vector identifiable and thus does not have a gradient.", self$n - 1))
         }
         
         private$mass_$grad <- value$to(device = private$mass_$device)
       }
     })
   },
   
   #' @field init_weights returns the initial value of the weights
   init_weights = function(value) {
     if(!missing(value)) stop("Can't change the initial weights. Try setting the weights by using the `$weights` operator.")
     return(private$init_weights_$clone())
   },
   
   #' @field init_data returns the initial value of the data
   init_data = function(value) {
     if(!missing(value)) stop("Can't change the initial values. Try setting the data by using the `$x` operator.")
     return(private$init_data_$clone())
   },
   
   #' @field requires_grad checks or turns on/off gradient
   requires_grad = function(value) {
     if (missing(value)) {
       rg <- switch(self$adapt,
              "none" = FALSE,
              TRUE)
       return(rg)
     }
     value <- match.arg(value, c("none","weights","x"))
     if (value == "none") {
       if(self$adapt == "weights") {
         private$mass_$requires_grad <- FALSE
       } else if (self$adapt == "x") {
         private$data_$requires_grad <- FALSE
       }
       self$adapt <- "none"
     } else if (value == "x") {
       if (self$adapt == "weights") {
         warning("Turning off gradients for weights and turning on for x values.")
         private$mass_$requires_grad <- FALSE
       } else if (self$adapt == "x") {
         warning("Gradients already on for the x values.")
       }
       private$data_$requires_grad <- TRUE
       self$adapt <- "x"
     } else if (value == "weights") {
       if (self$adapt == "x") {
         warning("Turning off gradients for x values and turning on for weights")
         private$data_$requires_grad <- FALSE
         
       } else if (self$adapt == "weights") {
         warning("Gradients already on for the weights.")
       }
       if(self$adapt != "weights") {
         private$get_mass_ <- function() {
           private$transform_(private$mass_)
         }
         private$assign_mass_ <- function(value) {
           # torch::autograd_set_grad_mode(enabled = FALSE)
           torch::with_no_grad({
             private$mass_$copy_(private$inv_transform_(value)$to(device = self$device))
           })
           # torch::autograd_set_grad_mode(enabled = TRUE)
         }
         
         if (self$probability_measure) {
           private$inv_transform_ <- function(value) {
             min_neg <- round(log(.Machine$double.xmin)) - 50.0
             logs <- torch::torch_log(value/sum(value))
             logs[logs< min_neg] <- min_neg
             logs <- logs[2:length(logs)] - logs[1]
             return(logs)
           }
           
           private$transform_ <- function(value) {
             full_param <- torch::torch_cat(
               list( torch::torch_zeros(1L, dtype = value$dtype),
                     value)
             )
             return(full_param$log_softmax(1)$exp())
           }
           
         } else {
           private$inv_transform_ <- torch::torch_log
           private$transform_ <- torch::torch_exp
         }
         
         private$assign_mass_(self$weights)
       }
       private$mass_$requires_grad <- TRUE
       self$adapt <- "weights"
     }
   },
   
   #' @field weights gets or sets weights
   weights = function(value) {
     if(missing(value)) {
       return(private$get_mass_())
     }
     stopifnot("Input is NA" = !isTRUE(all(is.na(value))))
     stopifnot("Input is NULL" = !isTRUE(is.null(value)))
     stopifnot("Input value is not same length as nrows of data" = (length(value) == self$n) )
     
     # browser()
     if(!inherits(value, "torch_tensor")) {
       value <- torch::torch_tensor(value, dtype = private$mass_$dtype, device = self$device)$contiguous()
     } else {
       stopifnot("Input tensor and original weights have different dtypes! " = isTRUE(value$dtype == private$mass_$dtype))
       if (isFALSE(value$device == private$mass_$device) ) {
         value <- value$to(device = private$mass_$device)
       }
     }
     
     if (self$probability_measure) {
       stopifnot("supplied weights must be >=0" = all(as.logical((value >=0)$to(device = "cpu"))))
       if(as.logical((sum(value) != 1)$to(device = "cpu"))) value <- (value/sum(value))$detach()
     }
     
     private$assign_mass_(value)
     
   },
   
   #' @field x Gets or sets the data.
   x = function(value) {
     
     # return data tensor if no value provided
     if (missing(value)) {
       return(private$data_)
     }
     
     # check input data
     stopifnot("Input is NA" = !all(is.na(value)))
     stopifnot("Input is NULL" = !is.null(value))
     
     # check if input is tensor
     if (!inherits(value, "torch_tensor")) {
       value <- torch::torch_tensor(value, device = self$device, dtype = private$data_$dtype)
     } else {
       stopifnot("Input tensor and original data have different dtypes! " = isTRUE(value$dtype == private$data_$dtype))
     }
     
     # make sure dimensions are correct
     l_value <- length(value)
     
     if(l_value != (self$n * self$d) ) stop(sprintf("Input must either be a matrix of dimension %s by %s or a vector with total length %s", self$n, self$d, self$n * self$d))
     
     dims_v <- dim(value)
     
     if (length(dims_v) != 2) {
       if(length(dims_y) >  2) warning("Tensor being reshaped to a two dimensional tensor")
       value <- value$view(c(self$n, self$d))
     }
     
     # set data
     # torch::autograd_set_grad_mode(enabled = FALSE)
     torch::with_no_grad({
     private$data_$copy_(value$to(device = self$device))
     })
     # torch::autograd_set_grad_mode(enabled = TRUE)
     
     # check if balance.function data is equal to data
     bf_not_equal_data <- isTRUE(rlang::obj_address(self$balance_functions) != rlang::obj_address(private$data_))
     
     if (bf_not_equal_data && !all(is.na(self$balance_functions))) {
       warning("Measure data reset but not the balance_functions. You may need to manually reset this as well.")
     }
     
   }
 )},
 private = {list(
   # values
   data_ = "torch_tensor",
   init_data_ = "torch_tensor",
   init_weights_ = "torch_tensor",
   mass_ = "torch_tensor",
   
   # functions
   assign_mass_ = "function", #transforms if needed
   deep_clone = function(name, value) {
     if (inherits(value, "torch_tensor")) {
       value$clone()
     } else {
       value
     }
   },
   get_mass_ = "function", #transforms if needed
   inv_transform_ = "function", # needs inv log_softmax for prob measure
   transform_ = "function" # needs log_softmax
   
 )}
)

#' @name Measure
#' @title An R6 Class for setting up measures
#'
#' @param x The data points
#' @param weights The empirical measure. If NULL, assigns equal weight to each observation
#' @param probability.measure Is the empirical measure a probability measure? Default is TRUE.
#' @param adapt Should we try to adapt the data ("x"), the weights ("weights"), or neither ("none"). Default is "none".
#' @param balance.functions A matrix of functions of the covariates to target for mean balance. If NULL and `target.values` are provided, will use the data in `x`.
#' @param target.values The targets for the balance functions. Should be the same length as columns in `balance.functions.`
#' @param dtype The torch_tensor dtype or NULL.
#' @param device The device to have the data on. Should be result of [torch::torch_device()] or NULL.
#' @return Returns a Measure object
#' 
#' @details # Public fields
#'   \if{html}{\out{<div class="r6-fields">}}
#'   \describe{
#'     \item{\code{balance_functions}}{the functions of the data that
#'       we want to adjust towards the targets}
#'     \item{\code{balance_target}}{the values the balance_functions are targeting}
#'     \item{\code{adapt}}{What aspect of the data will be adapted. One of "none","weights", or "x".}
#'     \item{\code{device}}{the \code{\link[torch:torch_device]{torch::torch_device}} of the data.}
#'     \item{\code{dtype}}{the \link[torch:torch_dtype]{torch::torch_dtype} of the data.}
#'     \item{\code{n}}{the rows of the covariates, x.}
#'     \item{\code{d}}{the columns of the covariates, x.}
#'     \item{\code{probability_measure}}{is the measure a probability measure?}
#'   }
#'   \if{html}{\out{</div>}}
#' @details # Active bindings
#'   \if{html}{\out{<div class="r6-active-bindings">}}
#'   \describe{
#'     \item{\code{grad}}{gets or sets gradient}
#'     \item{\code{init_weights}}{returns the initial value of the weights}
#'     \item{\code{init_data}}{returns the initial value of the data}
#'     \item{\code{requires_grad}}{checks or turns on/off gradient}
#'     \item{\code{weights}}{gets or sets weights}
#'     \item{\code{x}}{Gets or sets the data}
#'   }
#'   \if{html}{\out{</div>}}
#' @details # Methods
#' \subsection{Public methods}{
#' \itemize{
#' \item \href{#method-Measure-detach}{\code{Measure$detach()}}
#' \item \href{#method-Measure-get_weight_parameters}{\code{Measure$get_weight_parameters()}}
#' \item \href{#method-Measure-clone}{\code{Measure$clone()}}
#' }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-Measure-detach"></a>}}
#' \if{latex}{\out{\hypertarget{method-Measure-detach}{}}}
#' \subsection{Method \code{detach()}}{
#' generates a deep clone of the object without gradients.
#' \subsection{Usage}{
#' \if{html}{\out{<div class="r">}}\preformatted{Measure$detach()}\if{html}{\out{</div>}}
#' }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-Measure-get_weight_parameters"></a>}}
#' \if{latex}{\out{\hypertarget{method-Measure-get_weight_parameters}{}}}
#' \subsection{Method \code{get_weight_parameters()}}{
#' Makes a copy of the weights parameters.
#' \subsection{Usage}{
#' \if{html}{\out{<div class="r">}}\preformatted{Measure$get_weight_parameters()}\if{html}{\out{</div>}}
#' }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-Measure-clone"></a>}}
#' \if{latex}{\out{\hypertarget{method-Measure-clone}{}}}
#' \subsection{Method \code{clone()}}{
#' The objects of this class are cloneable with this method.
#' \subsection{Usage}{
#' \if{html}{\out{<div class="r">}}\preformatted{Measure$clone(deep = FALSE)}\if{html}{\out{</div>}}
#' }
#' \subsection{Arguments}{
#' \if{html}{\out{<div class="arguments">}}
#' \describe{
#'   \item{\code{deep}}{Whether to make a deep clone.}
#' }
#' \if{html}{\out{</div>}}
#' }
#' }
#' @examples 
#' if(torch::torch_is_installed()) {
#' m <- Measure(x = matrix(0, 10, 2), adapt = "none")
#' print(m)
#' m$x
#' m$x <- matrix(1,10,2) # must have same dimensions
#' m$x
#' m$weights
#' m$weights <- 1:10/sum(1:10)
#' m$weights
#' 
#' # with gradients
#' m <- Measure(x = matrix(0, 10, 2), adapt = "weights")
#' m$requires_grad # TRUE
#' m$requires_grad <- "none" # turns off
#' m$requires_grad # FALSE
#' m$requires_grad <- "x"
#' m$requires_grad # TRUE
#' m <- Measure(matrix(0, 10, 2), adapt = "none")
#' m$grad # NULL
#' m <- Measure(matrix(0, 10, 2), adapt = "weights")
#' loss <- sum(m$weights * 1:10)
#' loss$backward()
#' m$grad
#' # note the weights gradient is on the log softmax scale
#' #and the first parameter is fixed for identifiability
#' m$grad <- rep(1,9)  
#' m$grad
#' }
#' @export
Measure <- function(x, 
                    weights = NULL, 
                    probability.measure = TRUE, 
                    adapt = c("none","weights", "x"), 
                    balance.functions = NA_real_,
                    target.values = NA_real_,
                    dtype = NULL, 
                    device = NULL) {
  
  return(Measure_$new(x = x, 
                      weights = weights,
                      probability.measure = probability.measure,
                      adapt = adapt,
                      balance.functions = balance.functions,
                      target.values = target.values,
                      dtype = dtype, 
                      device = device))
  
}

#' Internal function to select appropriate loss function
#'
#' @description Selects sinkhorn or energy distance losses depending on value
#' of penalty parameter
#' 
#' @param ot an OT object
#' @keywords internal
oop_loss_select <- function(ot) {
  lambda <- ot$penalty
  if (is.finite(lambda)) {
    return(sinkhorn_dist(ot))
  } else if ( is.infinite(lambda) ) {
    return(inf_sinkhorn_dist(ot))
  }
}

#' @name OTProblem_-class
#' @title An R6 class to construct OTProblems
#' @rdname OTProblem_-class
#' @description OTProblem R6 class
#' @keywords internal
OTProblem_ <- R6::R6Class("OTProblem",
 public = {list(
   # objects
   
   #' @field device the [torch::torch_device()] of the data.
   device = "torch_device",
   
   #' @field dtype the [torch::torch_dtype] of the data.
   dtype  = "torch_dtype",
   
   #' @field selected_delta the delta value selected after `choose_hyperparameters`
   selected_delta = "numeric", # final delta
   
   #' @field selected_lambda the lambda value selected after `choose_hyperparameters`
   selected_lambda = "numeric", # final lambda

   # functions
   
   #' @param o2 A number or object of class OTProblem
   #' @description adds `o2` to the OTProblem
   add = function(o2) {
     private$unaryop(o2, "+")
   },
   
   #' @param o2 A number or object of class OTProblem
   #' @description subtracts `o2` from OTProblem
   subtract = function(o2) {
     private$unaryop(o2, "-")
   },
   
   #' @param o2 A number or object of class OTProblem
   #' @description multiplies OTProblem by `o2`
   multiply = function(o2) {
     private$unaryop(o2, "*")
   },
   
#' @param o2 A number or object of class OTProblem
#' @description divides OTProblem by `o2`
   divide = function(o2) {
     private$unaryop(o2, "/")
   },

#' @description prints the OT problem object
#' @param ... Not used
print = function(...) {
  obj <- rlang::expr_text(private$objective)
  obj <- gsub("oop_loss_select", "OT", obj)
  obj <- gsub('private$ot_objects[[\"', "", obj, fixed = TRUE)
  obj <- gsub('\"]]', "", obj, fixed = TRUE)
  obj <- gsub("\n    ", "", obj)
  obj <- gsub("+ ", "+\n  ", obj, fixed = TRUE)
  obj <- gsub("- ", "-\n  ", obj, fixed = TRUE)
  cat("OT Problem: \n")
  cat("  ", obj, "\n", sep = "")
},

#' @description Constructor method
#' @param measure_1 An object of class [Measure]
#' @param measure_2 An object of class [Measure]
#' @param ... Not used at this time 
#'
#' @return An R6 object of class "OTProblem"
initialize = function(measure_1, measure_2) {
  # browser()
  add_1 <- rlang::obj_address(measure_1)
  add_2 <- rlang::obj_address(measure_2)
  addresses <- c(add_1, add_2)
  address_names <- paste0(addresses, collapse = ", ")
  
  stopifnot("argument 'measure_1' must be of class 'Measure'" = inherits(measure_1, "Measure"))
  stopifnot("argument 'measure_2' must be of class 'Measure'" = inherits(measure_2, "Measure"))
  
  dtype <- measure_1$dtype
  device <- measure_1$device
  
  if(isFALSE(measure_2$dtype == dtype) ) {
    stop(sprintf("Measures must have same data type! measure_1 is of type %s, while measure_2 is of type %s.", dtype, measure_2$dtype))
  }
  if (!(measure_2$device == device) ) { # can't use != with torch device
    stop(sprintf("Measures should be on same device! measure_1 is is on device %s, while measure_2 is on device %s.", device, measure_2$device) )
  }
  if (measure_1$d != measure_2$d) {
    stop(sprintf("Measures should have the same number of columns! measure_1 has %s columns, while measure_2 has %s columns.", measure_1$d, measure_2$d))
  }
  
  self$dtype <- dtype
  self$device <- device
  
  #environment with names as obj_add, and measures as elements of environment
  private$measures <- rlang::env(!!add_1 := measure_1, !!add_2 := measure_2)
  
  # env with names as obj_add1, obj_add2 (sorted),
  #contains vector with c(obj_add1, obj_add2)
  private$problems <- rlang::env(!!address_names := addresses)
  
  # envionrment with names as obj_add1, obj_add2 (sorted), then a list with f, g duals
  # self$duals <- rlang::env(!!address_names := list(!!addresses[1] := torch::torch_zeros(private$measures[[addresses[1] ]]$n, device = device, dtype = dtype)),
  #                          !!addresses[2] := torch::torch_zeros(private$measures[[addresses[2] ]]$n, device = device, dtype = dtype) )
  
  # ot_objects
  # envionrment with names as obj_add1, obj_add2 (sorted), with OT class objects
  private$ot_objects <- rlang::env()
  
  # target_objects
  # envionrment with names as obj_add1, obj_add2 (sorted), with list of balance.functions, means and delta values
  private$target_objects <- rlang::env()
  
  #penalty list
  private$penalty_list <- list(lambda = NA_real_, delta = NA_real_)
  
  # parameter list initialize
  private$parameters <- list()
  
  
  private$objective <- rlang::expr(
    oop_loss_select(private$ot_objects[[!!address_names]])
  )
  
  private$args_set <- FALSE
  private$opt <- private$sched <- NULL
  
  return(invisible(self))
},

#' @param lambda The penalty parameters to try for the OT problems. If not provided, function will select some
#' @param delta The constraint paramters to try for the balance function problems, if any
#' @param grid.length The number of hyperparameters to try if not provided
#' @param cost.function The cost function for the data. Can be any function that takes arguments `x1`, `x2`, `p`. Defaults to the Euclidean distance
#' @param p The power to raise the cost matrix by. Default is 2
#' @param cost.online Should online costs be used? Default is "auto" but "tensorized" stores the cost matrix in memory while "online" will calculate it on the fly.
#' @param debias Should debiased OT problems be used? Defaults to TRUE
#' @param diameter Diameter of the cost function.
#' @param ot_niter Number of iterations to run the OT problems
#' @param ot_tol The tolerance for convergence of the OT problems
#'
#' @return NULL
setup_arguments = function(lambda, delta, 
                           grid.length = 7L,
                           cost.function = NULL, 
                           p = 2,
                           cost.online = "auto",
                           debias = TRUE,
                           diameter = NULL, ot_niter = 1000L,
                           ot_tol = 1e-3) {
  
  prob_names <- ls(private$problems)
  if(private$args_set) warning("OT problems already set up. This function will erase previous objects")
  problem_1 <- problem_2 <- measure_1 <- measure_2 <- NULL
  device_vector <- NULL
  not_warned <- FALSE
  
  for(v in prob_names) {
    problem_1 <- private$problems[[v]][[1]]
    problem_2 <- private$problems[[v]][[2]]
    measure_1 <- private$measures[[problem_1]]
    measure_2 <- private$measures[[problem_2]]
    private$ot_objects[[v]] <- OT$new(x = measure_1$x$detach(),
                                      y = measure_2$x$detach(),
                                      a = measure_1$weights$detach(),
                                      b = measure_2$weights$detach(),
                                      penalty = 10.0, 
                                      cost_function = cost.function, 
                                      p = p, 
                                      debias = debias, 
                                      tensorized = cost.online,
                                      diameter = diameter,
                                      device = self$device,
                                      dtype = self$dtype)
    if(not_warned && isTRUE(!(device_vector == measure_1$device)) ){
      warning("All measures not on same device. This could slow things down.")
      not_warned <- FALSE
    } else {
      device_vector <- measure_1$device
    }
    
    
  }
  
  
  # ot opt param
  private$ot_niter <- as.integer(ot_niter)
  private$ot_tol <- as.numeric(ot_tol)
  
  # targets
  measure_addresses <- ls(private$measures)
  delta_set <- NA_real_
  not_na_bt <- not_na_bf <- FALSE
  meas <- NULL
  for (v in measure_addresses) {
    meas <- private$measures[[v]]
    not_na_bt <- !all(is.na(meas$balance_target) )
    not_na_bf <- !all(is.na(meas$balance_functions) )
    if (not_na_bt && not_na_bf) {
      if(meas$adapt != "none") {
        if (meas$balance_functions$requires_grad) {
          delta_set <- 0.0
        } else {
          delta_set <- NA_real_
        }
        private$target_objects[[v]] <- list(
          bf = meas$balance_functions,
          bt = meas$balance_target,
          delta = delta_set)
      }
      
    }
  }
  
  # parameters
  measure_addresses <- ls(private$measures)
  adapt <- NULL
  meas <- NULL
  for (v in measure_addresses) {
    meas <- private$measures[[v]]
    adapt <- meas$adapt
    if(adapt != "none") {
      private$parameters[[v]] <-
        switch(adapt,
               "x" = meas$.__enclos_env__$private$data_,
               "weights" = meas$.__enclos_env__$private$mass_)
    }
  }
  
  # make sure ot_args ok
  stopifnot("'ot_niter' must be > 0" = private$ot_niter > 0)
  stopifnot("'ot_tol' must be > 0" = private$ot_tol > 0)
  
  # ot penalty
  
  diameters <- length(private$ot_objects)
  names_ot <- ls(private$ot_objects)
  v <- NULL
  for(i in seq_along(private$ot_objects)) {
    v <- names_ot[i]
    diameters[i] <- private$ot_objects[[v]]$diameter
  }
  max_diameter <- max(diameters)
  stopifnot("The maximum diameter of the OT problem is not finite!" = is.finite(max_diameter))
  l_md <-log(max_diameter)
  stopifnot("The log of the maximum diameter of the OT problem is not finite!" = is.finite(l_md))
  if ( missing(lambda) || is.null(lambda) || all(is.na(lambda)) ) {
    lambda <- c( #log(max_diam * 1e-6) to log(max_diam * 1e4), about
      exp(seq(l_md - 13.82, l_md + 9.21, length.out = grid.length)),
      Inf)
  } 
  private$penalty_list$lambda <- sort(lambda, decreasing = TRUE)
  private$penalty_list$lambda[lambda == 0] <- max_diameter / 1e9
  
  if ( length(private$target_objects) != 0) {
    
    if ( missing(delta) || is.null(delta) || all(is.na(delta)) ) {
      
      diffs    <- numeric(length(private$target_objects))
      names_TO <- ls(private$target_objects)
      measure <- NULL
      for(i in seq_along(private$target_objects)) {
        v <- names_TO[i]
        measure <- private$measures[[v]]
        if(measure$adapt == "weights") {
          diffs[i] <- ((private$target_objects[[v]]$bf * measure$weights$detach()$view(c(measure$n,1)))$sum(1) - private$target_objects[[v]]$bt)$abs()$max()$item()
        }
        
      }
      max_diffs <- max(diffs)
      
      delta <- c(seq(1e-4, max_diffs, length.out = grid.length))
      
    }
    
    
    private$penalty_list$delta <- sort(as.numeric(delta), decreasing = TRUE)
    
    stopifnot("'delta' values must be >=0."=all(private$penalty_list$delta >= 0))
    
  }
  
  # flag to warn about overwrite next time function called
  private$args_set <- TRUE
  
  return(invisible(self))
},

#' @description Solve the OTProblem at each parameter value. Must run setup_arguments first.
#' @param niter The nubmer of iterations to run solver at each combination of hyperparameter values 
#' @param tol The tolerance for convergence
#' @param optimizer The optimizer to use. One of "torch" or "frank-wolfe"
#' @param torch_optim The `torch_optimizer` to use. Default is [torch::optim_lbfgs]
#' @param torch_scheduler The [torch::lr_scheduler] to use. Default is [torch::lr_reduce_on_plateau]
#' @param torch_args Arguments passed to the torch optimizer and scheduler
#' @param osqp_args Arguments passed to [osqp::osqpSettings()] if appropriate
#' @param quick.balance.function Should [osqp::osqp()] be used to select balance function constraints (delta) or not. Default true.
  solve = function(niter = 1000L, tol = 1e-5, optimizer = c("torch", "frank-wolfe"),
                   torch_optim = torch::optim_lbfgs,
                   torch_scheduler = torch::lr_reduce_on_plateau,
                   torch_args = NULL,
                   osqp_args = NULL,
                   quick.balance.function = TRUE) {
    
    # check that everything setup already
    stopifnot("arguments not set! Run '$setup_arguments' first." = private$args_set)
    # check that niter and tol arguments provided
    # stopifnot("`niter` argument must be provided" = !missing(niter))
    stopifnot("Argument `niter` must be > 0" = (niter > 0))
    # stopifnot("Argument `tol` must be provided" = !missing(tol))
    stopifnot("Argument `tol` must be >=0" = (tol >= 0))
    
    # collect osqp args
    private$osqp_args <- osqp_args[names(osqp_args) %in% methods::formalArgs(osqp::osqpSettings)] 
    
    # check feasibility of deltas
    # can also do a quick, approximate selection of deltas
    private$delta_values_setup(run.quick = quick.balance.function, osqp_args = private$osqp_args)
    
    # setup optimizer
    optimizer <- match.arg(optimizer)
    stopifnot("Optimizer must be one of 'torch' or 'frank-wolfe'." = (optimizer %in% c("torch", "frank-wolfe") ))
    opt_call <- opt <- scheduler_call <- opt_sched <- NULL
    
    # assign optimizer to `private$optimization_step` holder
    if (optimizer == "torch") {
      
      private$torch_optim_setup(torch_optim,
                                torch_scheduler,
                                torch_args)
      
    } else if (optimizer == "frank-wolfe") {
      private$frankwolfe_setup()
      private$optimization_step <- private$frankwolfe_step
    } else {
      stop("Optimizer must be one of 'torch' or 'frank-wolfe'.")
    }
    
    # strategy for this function
    # outer loop iterates over lambda
    # inner loop iterates over delta
    # for mirror descent ONLY
    # add BF violations
    # calculate loss
    # torch_optim step
    # for frank-wolfe
    # run ot opt
    # get results from LP
    # step (armijo line search)
    
    # setup holder for weights
    private$weights_list <- vector("list", length(private$penalty_list$lambda))
    names(private$weights_list) <- as.character(private$penalty_list$lambda)
    
    # setup diagnostic collection of variables
    private$iterations_run <- vector("list", length(private$penalty_list$lambda)) |>
      setNames(as.character(private$penalty_list$lambda))
    
    private$final_loss <- vector("list", length(private$penalty_list$lambda)) |>
      setNames(as.character(private$penalty_list$lambda))
    
    
    # optimize over lambda values
    private$iterate_over_lambda(niter, tol)
    torch_cubic_reassign()
    return(invisible(self))
  },

  #' @param n_boot_lambda The number of bootstrap iterations to run when selecting lambda
  #' @param n_boot_delta The number of bootstrap iterations to run when selecting delta
  #' @param lambda_bootstrap The penalty parameter to use when selecting lambda. Higher numbers run faster.
  #'
  #' @description Selects the hyperparameter values through a bootstrap algorithm
  choose_hyperparameters =  function(n_boot_lambda = 100L, n_boot_delta = 1000L, lambda_bootstrap = Inf) {
      
      # check arguments
      stopifnot("n_boot_lambda must be >= 0"= n_boot_lambda>=0)
      stopifnot("n_boot_delta must be >= 0"= n_boot_delta>=0)
      stopifnot("lambda_bootstrap must be > 0"=lambda_bootstrap>0)
      
      # alter lists if needed for inherited classes
      res <- private$setup_choose_hyperparameters()
      
      # pull out current delta and lambda values
      delta_values <- res$delta
      lambda_values <- res$lambda
      
      # vector parameters for temporary holders
      n_delta <- length(delta_values)
      n_lambda<- length(lambda_values)
      
      # current weights list
      weights_list  <- res$weights_list
      
      # setup final metrics list
      private$weights_metrics <- list(delta = vector("list", n_lambda),
                                      lambda= vector("list", n_lambda))
      
      # check if function already ran before
      if(is.numeric(self$selected_lambda)) {
        warning(sprintf("Lambda value of %s previously selected. This function will erase previous bootstrap selection", self$selected_lambda))
      }
      
      if (n_delta > 1) { # begin delta select
        # setup temporary vector to hold delta evaluations
        delta_temp_metric <- vector("numeric", n_delta) 
        class(delta_temp_metric) <- c("weightEnv", class(delta_temp_metric))
        
        lambda_temp_metric <- vector("list", n_lambda)
        
        for (l in seq_along(lambda_values)) {
          lambda_temp_metric[[l]] <- vector("numeric", n_delta)
        }
        
        self$selected_delta <- vector("numeric", n_lambda) 
        
        # boot measure holder
        boot_measure <- NULL
        
        # running delta evaluations
        for ( i in 1:n_boot_delta ) {
          boot_measure <- private$draw_boot_measure()
          for ( l in seq_along(lambda_values) ) {
            delta_temp_metric <- vapply(weights_list[[l]],
                                        FUN = private$eval_delta,
                                        FUN.VALUE = 0.0,
                                        boot = boot_measure
            )
            lambda_temp_metric[[l]] <- lambda_temp_metric[[l]] + delta_temp_metric/n_boot_delta
          }
        }
        
        
        # assign metrics back to final location
        private$weights_metrics$delta <- lambda_temp_metric
        
        # select final deltas
        delta_eval_holder <- vector("numeric", n_delta)
        d_idx <- NULL
        targ_names <- ls(private$target_objects)
        param_names<- ls(private$parameters)
        
        # run through the parameter list and assign all to temp_wt
        # for (address in param_names ) {
        #   for(l in l_names) {
        #     # only look at delta eval if has target
        #     if (address %in% targ_names) {
        #       for (d in d_names) {
        #         delta_eval_holder[[d]] <- private$weights_metrics$delta[[l]][[d]][[address]]
        #       }
        #       d_idx <- which.min(delta_eval_holder)
        #     } else {
        #       d_idx <- 1L
        #     }
        #     
        #     # assign final selected weights
        #     temp_wt[[l]][[t_address]] <- private$weights_list[[l]][[d_idx]][[address]]
        #   }
        # }
        min_d_idx <- NULL
        for(l in seq_along(lambda_values)) {
          min_d_idx <- which.min(lambda_temp_metric[[l]])[1]
          weights_list[[l]] <-  weights_list[[l]][[min_d_idx]]
          self$selected_delta[[l]] <- delta_values[[min_d_idx]]
        }
        
      } else {
        new_weights_list <- vector("list", length(weights_list))
        for (l in seq_along(lambda_values)) {
          new_weights_list[[l]] <- weights_list[[l]][[1L]]
        }
        weights_list <- new_weights_list
        
        if(length(private$target_objects) > 0) {
          self$selected_delta <- list(rep(NA_real_, length(private$target_objects)))
          names_target <- ls(private$target_objects)
          for (v in seq_along(private$target_objects) ) {
            self$selected_delta[[1L]][[v]] <- private$target_objects[[names_target[v] ]]$delta
          }
          names(self$selected_delta[[1L]]) <- names_target
        }
      }# end of delta selection
      
      if (n_lambda > 1) {
        
        # use Energy Dist
        private$set_lambda(lambda_bootstrap)
        
        # setup holder
        lambda_metrics <- rep(0.0, n_lambda)
        boot_measure <- NULL
        
        for (i in 1:n_boot_lambda ) {
          boot_measure <- private$draw_boot_measure()
          lambda_metrics <- vapply(X = weights_list,
                                   FUN = private$eval_lambda,
                                   FUN.VALUE = 0.0,
                                   boot = boot_measure)/n_boot_lambda +
            lambda_metrics
        }
        # choose final lambda value
        idx_lambda <- which.min(lambda_metrics)
        selected_lambda <- as.numeric(lambda_values)[idx_lambda]
        
        # set metrics
        private$weights_metrics$lambda <- lambda_metrics
        
        # pull out final wts
        weights_list <- weights_list[[idx_lambda]]
        
        # save selected lambda
        self$selected_lambda <- selected_lambda
        private$set_lambda(selected_lambda)
        
        if(length(self$selected_delta) > 1) self$selected_delta <- self$selected_delta[[idx_lambda]]
      } else {
        weights_list <- weights_list[[1L]]
        self$selected_lambda <- as.numeric(lambda_values)[1L]
        private$set_lambda(self$selected_lambda)
      }
      
      # set weights back to the measures
      private$parameters_get_set(value = weights_list)
      
      
      # private$ot_update(only_params = FALSE, get_weights = TRUE, use_grad = FALSE)
      
    },

#' @description Provides diagnostics after solve and choose_hyperparameter methods have been run.
#'
#' @return a list with slots
#' \itemize{
#' \item `loss` the final loss values
#' \item `iterations` The number of iterations run for each combination of parameters
#' \item `balance.function.differences` The final differences in the balance functions
#' \item `hyperparam.metrics` A list of the bootstrap evalustion for delta and lambda values}
   info = function(){
     losses <- if (is.list(private$final_loss)) {
       do.call("rbind", private$final_loss)
     } else {
       NULL
     }
     
     metrics <- private$weights_metrics
     if(!is.character(metrics)) {
       delta_df <- do.call("cbind", metrics$delta)
       metrics["delta"] <- list(delta_df)
     } else {
       metrics <- "Hyperparameters not selected yet"
     }
     
     
     bal <- as.list(private$balance_check())
     
     iter <- if(is.list(private$iterations_run)) {
       do.call("rbind", private$iterations_run)
     } else {
       NULL
     }
     
     return(list(loss = losses,
                 iterations = iter,
                 balance.function.differences = bal,
                 hyperparam.metrics = metrics))
   }
   
 )},
 active = {list(
#' @field loss prints the current value of the objective. Only availble after the solve method has been run
   loss = function() {
     private$ot_update()
     return(eval(private$objective)$to(device = self$device))
   },
   
   #' @field penalty Returns a list of the lambda and delta penalities that will be iterated through. To set these values, use the `setup_arguments` function.
   penalty = function() {
     return(private$penalty_list)
   }
 )},
 private = {list(
   # objects
   args_set = "logical",
   final_loss = "list",
   iterations_run = "list",
   lbfgs_reset = 0L,
   lbfgs_count = 0L,
   # @field measures An an environment of measure objects named by address
   measures = "env", # environment of measure objects, named by address
   # @field objective An [rlang::expr()] giving the objective function
   objective = "rlang",
   opt_calls = "list", #saves arguments for the optimizers to reset
   opt = "R6", # store optimizer so easier to reset for LBFGS
   osqp_args = "list",
   ot_niter = "integer",
   ot_objects = "env", # environment with ot R6 classes
   ot_tol = "numeric",
   parameters = "list", # list of the parameters of model
   # @field problems An environment giving the object addresses for each measure in the OT problems
   problems = "env", # character vector of object addresses crossed
   penalty_list = "list", # list of penalties to try
   sched = "R6", # store scheduler
   target_objects = "env", # balance target objects
   weights_metrics = "list", # list of bootstrap metric for each penalty
   weights = "weightEnv", # holder of transformed parameters that are weights
   weights_list = "list", # list of estimated weights for each penalty
   
   # functions
   # adds new element to specified environment
   append = function(o, env_name) {
     listE1 = ls(private[[env_name]])
     listE2 = ls(o$.__enclos_env__$private[[env_name]])
     for(v in listE2) {
       if(v %in% listE1) {
         next
       } else {
         private[[env_name]][[v]] <- o$.__enclos_env__$private[[env_name]][[v]]
       }
     }
   },
   
   # update balance constraint parameters
   bal_param_update = function(osqp_args = NULL, tol = 1e-7) {
     # update balance constraint parameters and then returns
     # the langrangian terms to add to the loss
     
     # save weights incase not already done
     private$weights <- private$get_weights()
     
     l_to <- length(private$target_objects)
     
     loss <- torch::torch_tensor(0.0, dtype = self$dtype, device = self$device)
     # calc_deriv <- FALSE
     
     machine_tol <- .Machine$double.xmin
     coef <- torch::torch_tensor(10000.0, dtype = self$dtype, device = self$device)
     
     if (length(l_to) > 0) {
       
       # variables in for loop
       delta <- bal <- ng_bal <- n <- d <- v <- w <- NULL
       problem <- q <- res <- l <- u <- A <- gamma <- NULL
       abs_bal <- const.viol <- which.max <- which.sign <- NULL
       
       
       # addresses of objects
       target_addresses <- ls(private$target_objects)
      
       # don't print everything by default
       # if(missing(osqp_args) || is.null(osqp_args)) osqp_args <- list(verbose = FALSE)
       
       # run through target means
       for (adds in target_addresses) {
         # get objects
         v   <- private$target_objects[[adds]]
         w   <- private$weights[[adds]]
         
         # setup values
         delta <- v$delta
         # n   <- nrow(v$bf)
         # d   <- ncol(v$bf)
         bal <-  v$bf$transpose(2,1)$matmul(w) - v$bt
         if (v$bf$requires_grad) { # center v$bf
           d   <- ncol(v$bf)
           torch::with_no_grad(
             v$bf$subtract_(bal$view(c(1L,d)))
             )
           next
         }
         
         # #Linear program variables
         # q  <- c(-as.numeric(bal), delta * rep(1, 2 * d)) #linear terms
         # A  <- rbind(cbind(matrix(0,d * 2,d), Matrix::Diagonal(d*2, 1)), cbind(Matrix::Diagonal(d,1),Matrix::Diagonal(d,-1),Matrix::Diagonal(d,1)) ) #constraints
         # l  <- rep(0, 3 * d) #constraint lower bounds
         # u  <- c(rep(Inf, 2 * d), rep(0,d)) #constraint upper bounds
         # 
         # # Linear program
         # problem <- osqp::osqp(q = q, A = A, l = l, u = u, pars = osqp_args)
         # res <- problem$Solve()
         # 
         # # dual variables for targets
         # gamma <- res$x[1:d]
         
         # approx answer if infeasible
         # if (res$info$status_val == -3 || res$info$status_val == -4) {
           ng_bal <- bal$detach()
           abs_bal   <- ng_bal$abs()
           const.viol<- abs_bal > delta
           which.max <- abs_bal$argmax()
           which.sign<- ng_bal[which.max]$sign()
           gamma <- torch::torch_zeros_like(ng_bal)
           #   if ( length(which.max) == 0 ) {
           #   browser()
           #   torch::torch_zeros_like(bal[1L], dtype = bal$dtype)
           # } else {
           #   bal$detach()[which.max]$sign()
           # }
           gamma[which.max] <- which.sign * const.viol[which.max]
         # }
         
         # update loss
           # cur_loss <- bal$dot(gamma) -  sum(abs(gamma)) * delta 
           # loss <- loss + cur_loss * 10000.0
           cur_loss <- bal$dot(gamma) - sum(abs(gamma)) * delta 
           if ( cur_loss$item() / (machine_tol + delta) > tol ) {
             loss <- loss + cur_loss * coef
           }
           
       }

     } 
     return(loss)
     # return(list(loss = loss, calc_grad = calc_deriv))
   },
   
   # check balance constraints
   balance_check = function() {
     l_to <- length(private$target_objects)
     
     ret <- NULL
     
     if (length(l_to) > 0) {
       delta <- bal <- n <- d <- v <- w <- NULL
       
       
       # addresses of objects
       target_addresses <- ls(private$target_objects)
       ret <- rlang::env()
       
       # don't print everything by default
       # if(missing(osqp_args) || is.null(osqp_args)) osqp_args <- list(verbose = FALSE)
       
       # run through target means
       for (adds in target_addresses) {
         # get objects
         v   <- private$target_objects[[adds]]
         w   <- private$measures[[adds]]$weights
         
         # setup values
         delta <- v$delta
         # n   <- nrow(v$bf)
         # d   <- ncol(v$bf)
         bal <-  v$bf$transpose(2,1)$matmul(w) - v$bt
         if(v$bf$requires_grad) bal <- bal/v$bf$std(1)
         
         ret[[adds]] <- list(balance = bal, delta = delta)
         
       }
     }
     return(ret)
   },
   
   # deep clone function
   deep_clone = function(name, value) {
     if (name %in% c("problems", "target_objects", "measures", "ot_objects")){
       list2env(as.list.environment(value, all.names = TRUE),
                parent = emptyenv())
     } else {
       value
     }
   },
   
   delta_values_setup = function(run.quick = TRUE, osqp_args = NULL) {
     
     # check if any target_objects
     n_tf    <- length(private$target_objects)
     
     # check if ndelta>1
     n_delta <- length(private$penalty_list$delta)
     
     # default to not verbose
     if (is.null(osqp_args)) {
       osqp_args <- list(verbose = FALSE)
     }
     
     # run if are target_objects
     if (n_tf > 0 && n_delta > 1) {
       
       w       <- vector("list", n_delta)
       
       # if check all deltas for all lambdas
       if (isFALSE(run.quick)) {
         names_d  <- as.character(private$penalty_list$delta)
         d <- bf <- w <- v <- m <- NULL
         
         target_addresses <- ls(private$target_objects)
         successful.deltas <- NULL
         
         for (addy in  target_addresses) {
           
           v  <- private$target_objects[[addy]]
           m  <- private$measures[[addy]]
           
           if(v$bf$requires_grad) next
           
           bf <- SBW_4_oop$new(source = v$bf,
                               target = v$bt,
                               # a      = m$weights,
                               prob.measure = m$probability_measure,
                               osqp_opts = osqp_args)
           
           for (nd in names_d) {
             d <- eval(rlang::parse_expr(nd))
             check <- bf$solve(d)
             if (!is.null(check)) {
               successful.deltas <- c(successful.deltas, d)
             } else {
               next
             }
           }
           successful.deltas <- sort(unique(successful.deltas), decreasing = TRUE)
           names_d  <- as.character(successful.deltas)
         }
         
         # set feasible deltas
         private$penalty_list$delta <- successful.deltas
         
       } else { # quick run
         
         # selects best delta for each target
         names_d  <- as.character(private$penalty_list$delta)
         d <- bf <- w <- v <- m <- NULL
         means <- NULL
         w_star <- a_star <- means <- NULL
         
         target_addresses <- ls(private$target_objects)
         
         nboot <- 1000L
         
         # check for each target
         for (addy in  target_addresses) {
           w        <- vector("list", n_delta)
           names(w) <- names_d
           
           v  <- private$target_objects[[addy]]
           m  <- private$measures[[addy]]
           
           if(v$bf$requires_grad) next
           
           bf <- SBW_4_oop$new(source = v$bf$to(device = "cpu"),
                               target = v$bt$to(device = "cpu"),
                               # a      = m$weights,
                               prob.measure = m$probability_measure,
                               osqp_opts = osqp_args)
           
           # w  <- lapply(names_d, function(nd) {
           #   d <- eval(rlang::parse_expr(nd))
           #   bf$solve(d)
           #   }) |> 
           #   setNames(names_d)
           for (nd in names_d) {
             d <- eval(rlang::parse_expr(nd))
             w[[nd]] <- bf$solve(d)
           }
           
           means <- sbw_oop_bs_(w, 
                                nboot, 
                                as.matrix(bf$source_scale$to(device = "cpu")), 
                                as.numeric(bf$target_scale$to(device = "cpu")), 
                                as.numeric(m$init_weights$to(device = "cpu")))
           
           # means <- rep(0.0, length(w))
           # for (i in 1L:1000L) {
           #   a_star <- rmultinom(1L, bf$n, as.numeric(m$init_weights))
           #   w_star <- lapply(w, function(ww) ww * a_star)
           #   means  <- means + vapply(w_star, bf$evalBoot, FUN.VALUE = 0.0)/1000.0
           # }
           
           private$target_objects[[addy]]$delta <- as.numeric(names(w)[which.min(means)[1L]])
           
         }
         
         
         # set global delta to NA, which will keep the individualized ones
         private$penalty_list$delta <- NA_real_ 
       }
     }
   },
   
   # make sure no barycenters!!
   frankwolfe_setup = function() {
     
     addresses <- ls(private$parameters)
     meas <- NULL
     
     for (a in addresses) {
       meas <- private$measures[[a]]
       if (meas$adapt == "x") stop("Our Frank-Wolfe algorithm can't handle barycenters. Depending on your problem, you can optimize the barycenter and then optimize the weights with Frank-Wolfe. Alternatively, you can use the torch optimizers directly.")
     }
     
   },
   
   # frankwolfe optimization algorithm
   frankwolfe_step = function(opt, osqp_args, tol) {
     
     # get weights and retain grad
     # can probably be deleted
     weight_setup <- function() {
       private$weights <- private$get_weights()
       pw <- NULL
       for(w in ls(private$weights)) {
         pw <- private$weights[[w]]
         if(pw$requires_grad) pw$retain_grad()
       }
     }
     
     # retain grad
     weight_retain_grad <- function() {
       pw <- NULL
       for(w in ls(private$weights)) {
         pw <- private$weights[[w]]
         if(pw$requires_grad) pw$retain_grad()
       }
     }
     
     # get grad of weights
     weight_grad <- function() {
       w_names <- ls(private$weights)
       out_grad <- rlang::env()
       class(out_grad) <- c("weightEnv", class(out_grad))
       pw <- NULL
       for(w in w_names) {
         pw <- private$weights[[w]]
         if(pw$requires_grad) out_grad[[w]] <- pw$grad
       }
       
       return(out_grad)
     }
     
     # get addresses
     bf_address    <- ls(private$target_objects) # measures with balance functions
     param_address <- ls(private$parameters) # all parameters/ meas with grads
     
     # setup holders
     old_weights <- private$parameters_get_set(clone = TRUE)
     # class(old_weights) <- c("weightEnv", class(old_weights))
     
     new_weights <- rlang::env()
     class(new_weights) <- c("weightEnv", class(new_weights))
     
     ## Run functions ##
     
     # zero out gradients
     private$zero_grad()
     
     # eval loss and get gradients
     # weight_setup()
     # private$ot_update() # now in loss fun
     init_loss          <- self$loss
     weight_retain_grad()
     init_loss$backward()
     
     # setup osqp args
     osqp_arg_call <- rlang::call2(osqp::osqpSettings, !!!osqp_args)
     osqp_args <- eval(osqp_arg_call)
     
     res <- cur_env <- meas <- osqp_opt <- n <- prob.measure <- sums <- sum_bounds <- bt_cpu <- NULL
     
     
     # get solutions most correlated with the gradients
     for (addy in param_address) {
       meas <- private$measures[[addy]]
       n <- meas$n
       prob.measure <- meas$probability_measure
       sums <- switch(1L + prob.measure,
                      Matrix::Matrix(data = 0, nrow = 1, ncol = n),
                      Matrix::Matrix(data = 1, nrow = 1, ncol = n))
       sum_bounds <- switch(1L + prob.measure,
                            c(0, Inf),
                            c(1, 1))
       if (addy %in% bf_address) {
         cur_env <- private$target_objects[[addy ]]
         bt_cpu <- as.numeric(cur_env$bt$to(device = "cpu"))
         osqp_opt <- osqp::osqp(q = as.numeric(private$weights[[addy]]$grad$to(device = "cpu")), 
                                A = rbind(
                                  t(as.matrix(cur_env$bf$to(device = "cpu"))),
                                  Matrix::Diagonal(n, x = 1),
                                  sums
                                ),
                                l = c(-cur_env$delta + bt_cpu, rep(0,n), sum_bounds[1]),
                                u = c(cur_env$delta + bt_cpu, rep(Inf,n), sum_bounds[2]),
                                pars = osqp_args)
       } else {
         osqp_opt <- osqp::osqp(q = as.numeric(private$weights[[addy]]$grad$to(device = "cpu")), 
                                A = rbind(
                                  Matrix::Diagonal(n, x = 1),
                                  sums
                                ),
                                l = c(rep(0,n), sum_bounds[1]),
                                u = c(rep(Inf,n), sum_bounds[2]),
                                pars = osqp_args)
       }
       
       res <- osqp_opt$Solve()
       res$x[res$x < 0] <- 0
       
       #save into env of weights
       new_weights[[addy]] <- switch(1L + prob.measure,
         res$x,
         renormalize(res$x))
     }
     
     deriv <- weight_grad()
     # class(deriv) <- c("weightEnv", class(deriv))
     
     deltaG <- new_weights - old_weights
     derphi0 <- sum(deltaG * deriv)
     
     # need function to update weights and run update on OT
     armijo_loss_fun <- function(x, dx, alpha, ...) {
       
       if(inherits(alpha, "torch_tensor")) {
         alpha <- as.numeric(alpha$item())
       }
       # assign linearly shifted weights
       private$weights <- x + dx * alpha
       
       # update OT problems
       # private$ot_update()
       
       # evaluate loss (which updates ot problems)
       loss <- self$loss$detach()
       
       # return value
       return(loss)
     }
     
     if (-derphi0$item() >  tol) {
       search_res <-  scalar_search_armijo(phi = armijo_loss_fun,
                                           phi0 = init_loss$detach()$item(),
                                           derphi0 = derphi0$item(),
                                           x = old_weights,
                                           dx = deltaG,
                                           c1 = 1e-4, alpha0 = 1.0, 
                                           amin = 0)
       if ( !is.null(search_res$alpha)) {
         if(inherits(search_res$alpha, "torch_tensor")) search_res$alpha <- as.numeric(search_res$alpha$item())
         search_res$alpha = min(max(search_res$alpha,0), 1)
         new_weights <- old_weights +  deltaG * search_res$alpha
         private$parameters_get_set(new_weights)
         loss <- search_res$phi1$detach()
       } else {
         private$parameters_get_set(old_weights)
         loss <- init_loss$detach()
       }
       
     } else {
       private$parameters_get_set(old_weights)
       loss <- init_loss$detach()
     }
     
     return(loss)
   },
   
   get_weights = function(clone = FALSE) {
     
     out <- rlang::env()
     class(out) <- c("weightEnv", class(out))
     
     addresses <- ls(private$measures)
     for (v in addresses) {
       out[[v]] <- 
         private$measures[[v]]$weights
     }
     return(out)
   },
   
   # returns curent gradient value for parameters
   grad = function() {
     param_address <- ls(private$parameters)
     out <- rlang::env()
     for ( p in param_address ) {
       out[[p]] <- private$parameters[[p]]$grad
     }
     return(out)
   },
   
   iterate_over_delta = function(niter, tol) {
     
     # counter for delta
     n_delta <- length(private$penalty_list$delta)
     names_delta <- as.character(private$penalty_list$delta)
     
     # set weights_delta holder
     weights_delta <- vector("list", n_delta)
     names(weights_delta) <- names_delta
     
     iter_delta <- vector("numeric", n_delta) |> setNames(names_delta)
     losses_delta <- vector("numeric", n_delta) |> setNames(names_delta)
     
     # run loop over delta values
     for (d in private$penalty_list$delta) {
       # set bf constraint if needed
       private$set_delta(d)
       
       # run inner loop that actually optimizes
       res  <- private$optimization_loop(niter, tol)
       
       # copy estimated weights to weights_delta holder
       weights_delta[[paste(d)]] <- private$parameters_get_set(clone = TRUE)
       
       # save iterations and final losses
       losses_delta[[paste(d)]] <- res$loss
       iter_delta[[paste(d)]] <- res$iter
       
     }
     
     return(list(loss = losses_delta, 
                 iter = iter_delta,
                 weights = weights_delta))
     
   },
   
   iterate_over_lambda = function(niter, tol) {
     for (l in private$penalty_list$lambda) {  
       # set penalty for OT problems
       private$set_lambda(l)
       
       # initialize ot dual
       private$ot_update(only_params = FALSE, get_weights = TRUE, use_grad = FALSE)
       
       # optimize over delta values
       res <- private$iterate_over_delta(niter, tol) # calls main workhorse functions
       
       # assign weights_delta to weights_list permanent object
       private$weights_list[[as.character(l)]] <- res$weights
       
       # save diagnostic values
         # save iterations run
         private$iterations_run[[as.character(l)]] <- res$iter
         
         # save final loss
         private$final_loss[[as.character(l)]] <- res$loss
         
       
       
     }
   },
   
   lr_reduce = function(loss) { #complicatd lr reduce to fix bugs in torch
     sched <- private$sched
     if (is.null(sched)) return(TRUE)
     
     check <- TRUE
     
     if (inherits(sched, "lr_reduce_on_plateau")) {
       get_lr <- function() {
         tryCatch(
           sched$get_lr(),
           error = function(e) sapply(sched$optimizer$state_dict()$param_groups, function(p) p[["lr"]])
         )
       }
       check <- FALSE
       
       old_mode <- sched$threshold_mode
       if (as.logical(loss == sched$best) ) sched$threshold_mode <- "abs"
       
       init_lr <- sched$optimizer$defaults$lr
       lr <- get_lr()
       min_lr <- sched$min_lrs[[1]]
       
       improved <- as.logical(sched$.is_better(loss, sched$best))
       sched$step(loss)
       
       if (inherits(private$opt, "optim_lbfgs") && isTRUE(private$opt$defaults$line_search_fn == "strong_wolfe") ) {
         if (private$lbfgs_count == sched$patience) {
           if (private$lbfgs_reset > 0) check <- TRUE
           best <- sched$best
           private$torch_optim_reset()
           private$lbfgs_count <- 0L
           private$lbfgs_reset <- private$lbfgs_reset + 1L
           private$sched$best <- best
         } else if (!improved) {
           private$lbfgs_count <- private$lbfgs_count + 1L
         } else if (improved) {
           private$lbfgs_count <- 0L
           private$lbfgs_reset <- 0L
         }
         
       } else {
         # if ( sched$num_bad_epochs == sched$patience && init_lr != lr) {
         if (abs(lr - min_lr)/(lr + .Machine$double.eps) < 1e-3) {
           check <- TRUE 
         } else {
           check <- FALSE
         }
       }
       
       sched$threshold_mode <- old_mode
     } else {
       sched$step()
     }
     
     
     return(check)
   },
   # holder variable for selected optimization method
   optimization_step = "function",
   
   # optimization inner loop
   optimization_loop = function(niter, tol) {
     # set initial loss
     loss_old <- self$loss$detach()$item()
     
     # check convergence variable initialization  
     check <- TRUE
     
     # run optimization steps for given penalties
     for ( i in 1:niter ) {
       # take opt step and return loss
       loss <- private$optimization_step(private$opt, private$osqp_args, tol = tol)$detach()$to(device = "cpu")$item()
       
       # reduce lr
       check <- private$lr_reduce(loss)
       # check <- lr_reduce(opt_sched, loss$detach())
       
       # see if converged
       if ( check && (i >1) && converged(loss, loss_old, tol) ) break
       
       # if not converged, save old loss
       loss_old <- loss
     }
     
     #reset torch optimizer if present
     private$torch_optim_reset()
     private$lbfgs_count <- private$lbfgs_reset <- 0L
     
     return(list(loss = loss,
                 iter = i))
     
   },
   
   #update ot problems
   ot_update = function(only_params = TRUE, get_weights = TRUE, use_grad = TRUE) {
     ot_adds    <- ls(private$ot_objects)
     param_adds <- ls(private$parameters)
     
     # set weights
     if (isTRUE(get_weights)) {
       private$weights <- private$get_weights()
     }
     
     # variables for the loop
     has_grad <- weight_grad <- cost_grad <- FALSE
     problem_addy <- NULL
     measure_1 <- measure_2 <- NULL
     
     cost_forward <- function(add_1, add_2, ot) {
       measure_1 <- private$measures[[add_1]]
       measure_2 <- private$measures[[add_2]]
       a_1 <- measure_1$adapt
       a_2 <- measure_2$adapt
       
       if(a_1 == "x" || a_2 == "x") {
         x <- measure_1$.__enclos_env__$private$data_
         y <- measure_2$.__enclos_env__$private$data_
         
         if(x$requires_grad) {
           update_cost(ot$C_xy, x, y$detach())
           if(ot$debias) update_cost(ot$C_xx, x, x$detach())
         }
         if(y$requires_grad) {
           update_cost(ot$C_yx, y, x$detach())
           if(ot$debias) update_cost(ot$C_yy, y, y$detach())
         }
         has_grad <- TRUE
       } else{
         has_grad <- FALSE
       }
       
       return(has_grad)
     }
     
     weights_forward <- function(add_1, add_2, ot) {
       
       has_grad <- FALSE
       if(private$measures[[add_1]]$adapt == "weights") {
         ot$a <- private$weights[[add_1]]
         has_grad <- TRUE
       }
       if(private$measures[[add_2]]$adapt == "weights") {
         ot$b <-  private$weights[[add_2]]
         has_grad <- TRUE
       }
       return(has_grad)
     }
     
     # loop over OT problems
     ot_prob_loop <- function() {
       cur_ot <- NULL
       for (addy in ot_adds) {
         cur_ot <- private$ot_objects[[addy]]
         
         # need to update weights used
         problem_addy <- private$problems[[addy]]
         
         # update cost if needed
         cost_grad <- cost_forward(problem_addy[[1L]],
                                   problem_addy[[2L]],
                                   cur_ot)
         
         # update weights if needed
         wt_grad <- weights_forward(problem_addy[[1L]],
                                    problem_addy[[2L]],
                                    cur_ot)
         
         has_grad <- (cost_grad || wt_grad)
         
         # run sinkhorn if needed
         if ( is.finite(cur_ot$penalty) && (has_grad || !only_params) ) {
           cur_ot$sinkhorn_opt(niter = private$ot_niter, tol = private$ot_tol)
         }
       }
     }
     # differential run with/without grad
     if (use_grad) {
       ot_prob_loop() # grad, if needed
     } else {
       torch::with_no_grad(ot_prob_loop) # no grad
     }
     
   },
   # return or set parameter weights
   parameters_get_set = function(value, clone = FALSE) {
     
     # return a clone of the weights
     if (missing(value)) {
       out <- rlang::env()
       class(out) <- c("weightEnv", class(out))
       
       param_addresses <- ls(private$parameters)
       for (v in param_addresses) {
         out[[v]] <- if (isTRUE(clone)) {
           if(private$measures[[v]]$adapt == "weights") {
             private$measures[[v]]$weights$detach()$clone()
           } else if (private$measures[[v]]$adapt == "x") {
             private$measures[[v]]$x$detach()$clone()
           }
         } else {
           if(private$measures[[v]]$adapt == "weights") {
             private$measures[[v]]$weights
           } else if (private$measures[[v]]$adapt == "x") {
             private$measures[[v]]$x
           }
         }
       }
       return(out)
     }
     
     # set weights
     param_addresses <- ls(private$parameters)
     
     if(!rlang::is_environment(value)){
       names(value) <- value_addresses <- 1:length(value)
     } else {
       value_addresses <- ls(value)
     }
     
     if(length(value_addresses) == length(param_addresses)) {
       v <- NULL
       u <- NULL
       torch::with_no_grad(
         for (i in seq_along(param_addresses)) {
           v <- param_addresses[[i]]
           u <- value_addresses[[i]]
           if(private$measures[[v]]$adapt == "weights") {
             private$measures[[v]]$weights <- value[[u]]
           } else if (private$measures[[v]]$adapt == "x") {
             private$measures[[v]]$x <- value[[u]]
           } else {
             stop("Error in assignment. Tried to assign to a measure without any gradients.")
           }
         }
       )
     } else {
       stop("Input must have same number of groups as do the parameters.")
     }
   },
   
   torch_optim_step = function(opt, osqp_args = NULL, tol) {
     
     closure <- function() {
       opt$zero_grad()
       # self$forward()
       loss <- private$bal_param_update(osqp_args = osqp_args, tol) +
         self$loss  
       
       # only run ot if no bal constraint violations
       # if (loss$item() == 0) {
         # private$ot_update()
         # loss <- self$loss + loss
       # }
       
       loss$backward()
       return(loss$to(device = "cpu"))
     }
     
     if( inherits(opt, "optim_lbfgs") ) {
       loss <- opt$step(closure)
     } else {
       loss <- closure()
       opt$step()
     }
     
     return(loss)
     
   },
   torch_optim_setup = function(torch_optim, 
                                torch_scheduler,
                                torch_args) {
     
     names_args <- names(torch_args)
     optim_args_names <- names(formals(torch_optim))
     opt_args <- torch_args[match(optim_args_names,
                                  names_args, nomatch = 0L)]
     opt_call <- rlang::call2(torch_optim,
                              params = private$parameters,
                              !!!opt_args)
     private$opt <- eval(opt_call)
     
     torch_lbfgs_check(private$opt)
     
     private$optimization_step <- private$torch_optim_step
     
     if (!is.null(torch_scheduler)) {
       scheduler_args_names <- names(formals(torch_scheduler))
       sched_args <- torch_args[match(scheduler_args_names, 
                                      names_args, 
                                      nomatch = 0L)]
       if(inherits(torch_scheduler, "lr_reduce_on_plateau") &&
          is.null( sched_args$patience ) ) {
         sched_args$patience <- 1L
       }
       if(inherits(torch_scheduler, "lr_reduce_on_plateau") && is.null( sched_args$min_lr ) && inherits(private$opt, "optim_lbfgs") ) {
         sched_args$min_lr <- private$opt$defaults$lr * 1e-3
       }
       scheduler_call <- rlang::call2(torch_scheduler,
                                      optimizer = private$opt,
                                      !!!sched_args)
       opt_sched <- eval(scheduler_call)
     } else {
       opt_sched <- scheduler_call <- NULL
     }
     
     private$opt_calls <- list(opt = opt_call,
                               sched = scheduler_call)
     
     private$sched <- opt_sched
     
     # return(list(opt = opt, opt_call = opt_call,
     #             sched = opt_sched, sched_call = scheduler_call))
     
   },
   torch_optim_reset = function(lr = NULL) {
     if(!is.null(private$opt)) {
       opt_call <- private$opt_calls$opt
       if(is.null(lr)) {
         private$opt <- eval(rlang::call_modify(opt_call, 
                                                params = private$parameters))
       } else {
         # browser()
         # def <- rlang::call_match(opt_call[[1]], opt_call, defaults = TRUE)
         private$opt <- eval(rlang::call_modify(opt_call, 
                                                params = private$parameters,
                                                lr = lr))
       }
       
     }
     
     if(!is.null(private$sched)) {
       private$sched <- eval(rlang::call_modify(
         private$opt_calls$sched, 
                               optimizer = private$opt))
     }
     
   },
   set_lambda = function(lambda) {
     stopifnot("lambda value must be >= 0" = (lambda >= 0))
     ot_prob_names <- ls(private$ot_objects)
     for (v in ot_prob_names) {
       private$ot_objects[[v]]$penalty <- lambda
     }
   },
   set_delta = function(delta) {
     if(length(private$target_objects) > 0 && !is.na(delta)) {
       stopifnot("delta value must be >= 0" = (delta >= 0))
       target_names <- ls(private$target_objects)
       for (v in target_names) {
         if(!private$target_objects[[v]]$bf$requires_grad) private$target_objects[[v]]$delta <- delta
       }
     }
   },
   set_penalties = function(lambda, delta) {
     if(!missing(lambda)) {
       private$set_lambda(lambda)
     }
     
     if(!missing(delta)) {
       private$set_delta(delta)
     }
   },
   setup_choose_hyperparameters = function() {
     return(
       list(delta = private$penalty_list$delta,
            lambda = private$penalty_list$lambda,
            weights_list = private$weights_list)
            )
   },
   
   unaryop = function(o, fun) {
     if(! is_ot_problem(o)) {
       private$objective <- rlang::parse_expr(
         paste(rlang::expr_text(private$objective), fun, o)
         )
     } else {
       stopifnot(self$device == o$device)
       stopifnot(self$dtype == o$dtype)
       private$append(o, "measures")
       private$append(o, "problems")
       private$objective <- rlang::parse_expr(paste(
        rlang::expr_text(private$objective),
        fun,
        rlang::expr_text(o$.__enclos_env__$private$objective)
       ))
     }
   },
  
   zero_grad = function() {
     for ( p in private$parameters ) {
       if(torch::is_undefined_tensor(p$grad)) next
       if(!p$requires_grad) {
         warning("One of the objects in parameters doesn't require a grad. Report this bug!")
         next
       }
       torch::with_no_grad( p$grad$copy_(0.0) )
     }
   }
 )}
)

#' Object Oriented OT Problem
#'
#' @param measure_1 An object of class [Measure]
#' @param measure_2 An object of class [Measure]
#' @param ... Not used at this time 
#'
#' @return An R6 object of class "OTProblem"
#' @details # Public fields
#'   \if{html}{\out{<div class="r6-fields">}}
#'   \describe{
#'     \item{\code{device}}{the \code{\link[torch:torch_device]{torch::torch_device()}} of the data.}
#'     \item{\code{dtype}}{the \link[torch:torch_dtype]{torch::torch_dtype} of the data.}
#'     \item{\code{selected_delta}}{the delta value selected after \code{choose_hyperparameters}}
#'     \item{\code{selected_lambda}}{the lambda value selected after \code{choose_hyperparameters}}
#'   }
#'   \if{html}{\out{</div>}}
#' @details # Active bindings
#'   \if{html}{\out{<div class="r6-active-bindings">}}
#'   \describe{
#'     \item{\code{loss}}{prints the current value of the objective. Only availble after the \href{#method-OTProblem-solve}{\code{OTProblem$solve()}} method has been run}
#'     \item{\code{penalty}}{Returns a list of the lambda and delta penalities that will be iterated through. To set these values, use the \href{#method-OTProblem-setup_arguments}{\code{OTProblem$setup_arguments()}} function.}
#'   }
#'   \if{html}{\out{</div>}}
#' @details # Methods
#'   \subsection{Public methods}{
#'     \itemize{
#'     \item \href{#method-OTProblem-add}{\code{OTProblem$add()}}
#'     \item \href{#method-OTProblem-subtract}{\code{OTProblem$subtract()}}
#'     \item \href{#method-OTProblem-multiply}{\code{OTProblem$multiply()}}
#'     \item \href{#method-OTProblem-divide}{\code{OTProblem$divide()}}
#'     \item \href{#method-OTProblem-setup_arguments}{\code{OTProblem$setup_arguments()}}
#'     \item \href{#method-OTProblem-solve}{\code{OTProblem$solve()}}
#'     \item \href{#method-OTProblem-choose_hyperparameters}{\code{OTProblem$choose_hyperparameters()}}
#'     \item \href{#method-OTProblem-info}{\code{OTProblem$info()}}
#'     \item \href{#method-OTProblem-clone}{\code{OTProblem$clone()}}
#'     }
#'     }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-OTProblem-add"></a>}}
#' \if{latex}{\out{\hypertarget{method-OTProblem-add}{}}}
#' \subsection{Method \code{add()}}{
#'   adds \code{o2} to the OTProblem
#'   \subsection{Usage}{
#'     \if{html}{\out{<div class="r">}}\preformatted{OTProblem$add(o2)}\if{html}{\out{</div>}}
#'   }
#'   \subsection{Arguments}{
#'     \if{html}{\out{<div class="arguments">}}
#'     \describe{
#'       \item{\code{o2}}{A number or object of class OTProblem}
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-OTProblem-subtract"></a>}}
#' \if{latex}{\out{\hypertarget{method-OTProblem-subtract}{}}}
#' \subsection{Method \code{subtract()}}{
#'   subtracts \code{o2} from OTProblem
#'   \subsection{Usage}{
#'     \if{html}{\out{<div class="r">}}\preformatted{OTProblem$subtract(o2)}\if{html}{\out{</div>}}
#'   }
#'   \subsection{Arguments}{
#'     \if{html}{\out{<div class="arguments">}}
#'     \describe{
#'       \item{\code{o2}}{A number or object of class OTProblem}
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-OTProblem-multiply"></a>}}
#' \if{latex}{\out{\hypertarget{method-OTProblem-multiply}{}}}
#' \subsection{Method \code{multiply()}}{
#'   multiplies OTProblem by \code{o2}
#'   \subsection{Usage}{
#'     \if{html}{\out{<div class="r">}}\preformatted{OTProblem$multiply(o2)}\if{html}{\out{</div>}}
#'   }
#'   \subsection{Arguments}{
#'     \if{html}{\out{<div class="arguments">}}
#'     \describe{
#'       \item{\code{o2}}{A number or an object of class OTProblem}
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-OTProblem-divide"></a>}}
#' \if{latex}{\out{\hypertarget{method-OTProblem-divide}{}}}
#' \subsection{Method \code{divide()}}{
#'   divides OTProblem by \code{o2}
#'   \subsection{Usage}{
#'     \if{html}{\out{<div class="r">}}\preformatted{OTProblem$divide(o2)}\if{html}{\out{</div>}}
#'   }
#'   \subsection{Arguments}{
#'     \if{html}{\out{<div class="arguments">}}
#'     \describe{
#'       \item{\code{o2}}{A number or object of class OTProblem}
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-OTProblem-setup_arguments"></a>}}
#' \if{latex}{\out{\hypertarget{method-OTProblem-setup_arguments}{}}}
#' \subsection{Method \code{setup_arguments()}}{
#'   \subsection{Usage}{
#'     \if{html}{\out{<div class="r">}}\preformatted{OTProblem$setup_arguments(
#'       lambda,
#'       delta,
#'       grid.length = 7L,
#'       cost.function = NULL,
#'       p = 2,
#'       cost.online = "auto",
#'       debias = TRUE,
#'       diameter = NULL,
#'       ot_niter = 1000L,
#'       ot_tol = 0.001
#'     )}\if{html}{\out{</div>}}
#'   }
#'   \subsection{Arguments}{
#'     \if{html}{\out{<div class="arguments">}}
#'     \describe{
#'       \item{\code{lambda}}{The penalty parameters to try for the OT problems. If not provided, function will select some}
#'       \item{\code{delta}}{The constraint paramters to try for the balance function problems, if any}
#'       \item{\code{grid.length}}{The number of hyperparameters to try if not provided}
#'       \item{\code{cost.function}}{The cost function for the data. Can be any function that takes arguments \code{x1}, \code{x2}, \code{p}. Defaults to the Euclidean distance}
#'       \item{\code{p}}{The power to raise the cost matrix by. Default is 2}
#'       \item{\code{cost.online}}{Should online costs be used? Default is "auto" but "tensorized" stores the cost matrix in memory while "online" will calculate it on the fly.}
#'       \item{\code{debias}}{Should debiased OT problems be used? Defaults to TRUE}
#'       \item{\code{diameter}}{Diameter of the cost function.}
#'       \item{\code{ot_niter}}{Number of iterations to run the OT problems}
#'       \item{\code{ot_tol}}{The tolerance for convergence of the OT problems}
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#'   \subsection{Returns}{
#'     NULL
#'   }
#'   \subsection{Examples}{
#'     \if{html}{\out{<div class="r example copy">}}
#'     \preformatted{ ot$setup_arguments(lambda = c(1000,10))
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-OTProblem-solve"></a>}}
#' \if{latex}{\out{\hypertarget{method-OTProblem-solve}{}}}
#' \subsection{Method \code{solve()}}{
#'   Solve the OTProblem at each parameter value. Must run setup_arguments first.
#'   \subsection{Usage}{
#'     \if{html}{\out{<div class="r">}}\preformatted{OTProblem$solve(
#'       niter = 1000L,
#'       tol = 1e-05,
#'       optimizer = c("torch", "frank-wolfe"),
#'       torch_optim = torch::optim_lbfgs,
#'       torch_scheduler = torch::lr_reduce_on_plateau,
#'       torch_args = NULL,
#'       osqp_args = NULL,
#'       quick.balance.function = TRUE
#'     )}\if{html}{\out{</div>}}
#'   }
#'   \subsection{Arguments}{
#'     \if{html}{\out{<div class="arguments">}}
#'     \describe{
#'       \item{\code{niter}}{The nubmer of iterations to run solver at each combination of hyperparameter values}
#'       \item{\code{tol}}{The tolerance for convergence}
#'       \item{\code{optimizer}}{The optimizer to use. One of "torch" or "frank-wolfe"}
#'       \item{\code{torch_optim}}{The \code{torch_optimizer} to use. Default is \link[torch:optim_lbfgs]{torch::optim_lbfgs}}
#'       \item{\code{torch_scheduler}}{The \link[torch:lr_scheduler]{torch::lr_scheduler} to use. Default is \link[torch:lr_reduce_on_plateau]{torch::lr_reduce_on_plateau}}
#'       \item{\code{torch_args}}{Arguments passed to the torch optimizer and scheduler}
#'       \item{\code{osqp_args}}{Arguments passed to \code{\link[osqp:osqpSettings]{osqp::osqpSettings()}} if appropriate}
#'       \item{\code{quick.balance.function}}{Should \code{\link[osqp:osqp]{osqp::osqp()}} be used to select balance function constraints (delta) or not. Default true.}
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#'   \subsection{Examples}{
#'     \if{html}{\out{<div class="r example copy">}}
#'     \preformatted{ ot$solve(niter = 1, torch_optim = torch::optim_rmsprop)
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-OTProblem-choose_hyperparameters"></a>}}
#' \if{latex}{\out{\hypertarget{method-OTProblem-choose_hyperparameters}{}}}
#' \subsection{Method \code{choose_hyperparameters()}}{
#'   Selects the hyperparameter values through a bootstrap algorithm
#'   \subsection{Usage}{
#'     \if{html}{\out{<div class="r">}}\preformatted{OTProblem$choose_hyperparameters(
#'       n_boot_lambda = 100L,
#'       n_boot_delta = 1000L,
#'       lambda_bootstrap = Inf
#'     )}\if{html}{\out{</div>}}
#'   }
#'   \subsection{Arguments}{
#'     \if{html}{\out{<div class="arguments">}}
#'     \describe{
#'       \item{\code{n_boot_lambda}}{The number of bootstrap iterations to run when selecting lambda}
#'       \item{\code{n_boot_delta}}{The number of bootstrap iterations to run when selecting delta}
#'       \item{\code{lambda_bootstrap}}{The penalty parameter to use when selecting lambda. Higher numbers run faster.}
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#'   \subsection{Examples}{
#'     \if{html}{\out{<div class="r example copy">}}
#'     \preformatted{ ot$choose_hyperparameters(n_boot_lambda = 10, 
#'                                              n_boot_delta = 10, 
#'                                              lambda_bootstrap = Inf)
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-OTProblem-info"></a>}}
#' \if{latex}{\out{\hypertarget{method-OTProblem-info}{}}}
#' \subsection{Method \code{info()}}{
#'   Provides diagnostics after solve and choose_hyperparameter methods have been run.
#'   \subsection{Usage}{
#'     \if{html}{\out{<div class="r">}}\preformatted{OTProblem$info()}\if{html}{\out{</div>}}
#'   }
#'   \subsection{Returns}{
#'     a list with slots
#'     \itemize{
#'       \item \code{loss} the final loss values
#'       \item \code{iterations} The number of iterations run for each combination of parameters
#'       \item \code{balance.function.differences} The final differences in the balance functions
#'       \item \code{hyperparam.metrics} A list of the bootstrap evalustion for delta and lambda values}
#'   }
#'   \subsection{Examples}{
#'     \if{html}{\out{<div class="r example copy">}}
#'     \preformatted{ ot$info()
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#' }
#' \if{html}{\out{<hr>}}
#' \if{html}{\out{<a id="method-OTProblem-clone"></a>}}
#' \if{latex}{\out{\hypertarget{method-OTProblem-clone}{}}}
#' \subsection{Method \code{clone()}}{
#'   The objects of this class are cloneable with this method.
#'   \subsection{Usage}{
#'     \if{html}{\out{<div class="r">}}\preformatted{OTProblem$clone(deep = FALSE)}\if{html}{\out{</div>}}
#'   }
#'   \subsection{Arguments}{
#'     \if{html}{\out{<div class="arguments">}}
#'     \describe{
#'       \item{\code{deep}}{Whether to make a deep clone.}
#'     }
#'     \if{html}{\out{</div>}}
#'   }
#' }
#' @examples
#' ## ------------------------------------------------
#' ## Method `OTProblem(measure_1, measure_2)`
#' ## ------------------------------------------------
#'
#' if (torch::torch_is_installed()) {
#'   # setup measures
#'   x <- matrix(1, 100, 10)
#'   m1 <- Measure(x = x)
#'   
#'   y <- matrix(2, 100, 10)
#'   m2 <- Measure(x = y, adapt = "weights")
#'   
#'   z <- matrix(3,102, 10)
#'   m3 <- Measure(x = z)
#'   
#'   # setup OT problems
#'   ot1 <- OTProblem(m1, m2)
#'   ot2 <- OTProblem(m3, m2)
#'   ot <- 0.5 * ot1 + 0.5 * ot2
#'   print(ot)
#'
#' ## ------------------------------------------------
#' ## Method `OTProblem$setup_arguments`
#' ## ------------------------------------------------
#'
#'   ot$setup_arguments(lambda = 1000)
#'
#' ## ------------------------------------------------
#' ## Method `OTProblem$solve`
#' ## ------------------------------------------------
#'
#'   ot$solve(niter = 1, torch_optim = torch::optim_rmsprop)
#'
#' ## ------------------------------------------------
#' ## Method `OTProblem$choose_hyperparameters`
#' ## ------------------------------------------------
#'
#'   ot$choose_hyperparameters(n_boot_lambda = 1,
#'                             n_boot_delta = 1, 
#'                             lambda_bootstrap = Inf)
#'
#' ## ------------------------------------------------
#' ## Method `OTProblem$info`
#' ## ------------------------------------------------
#'
#' ot$info()
#' }
#' @export
OTProblem <- function(measure_1, measure_2,...) {
  
  OTProblem_$new(measure_1 = measure_1, 
                  measure_2 = measure_2)
}

is_ot_problem <- function(obj) {
  isTRUE(inherits(obj, "OTProblem"))
}

binaryop <- function(e1, e2, fun) {
  UseMethod("binaryop")
}

#' @export
binaryop.OTProblem <- function(e1, e2, fun) {
  if ( !is_ot_problem(e1) ) {
    stopifnot("LHS must be numeric or OTProblem" = is.numeric(e1))
    o_new <- e2$clone(deep = TRUE)
    o_old <- e1
  } else if ( !is_ot_problem(e2)) {
    stopifnot("RHS must be numeric or OTProblem"=is.numeric(e2))
    o_new <- e1$clone(deep = TRUE)
    o_old <- e2
  } else if(is_ot_problem(e1) && is_ot_problem(e2)) {
    o_new <- e1$clone(deep = TRUE)
    o_old <- e2
  } else {
    stop("error in basic operation on 'OTProblem' objects")
  }
  
  o_new[[fun]](o_old)
  
  return(o_new)
  
}

#' @export
`+.OTProblem` <- function(e1, e2) {
  o_new <- binaryop.OTProblem(e1, e2, "add")
  return(invisible(o_new))
  
}

#' @export
`-.OTProblem` <- function(e1, e2) {
  o_new <- binaryop.OTProblem(e1, e2, "subtract")
  return(invisible(o_new))
  
}

#' @export
`*.OTProblem` <- function(e1, e2) {
  o_new <- binaryop.OTProblem(e1, e2, "multiply")
  return(invisible(o_new))
  
}

#' @export
`/.OTProblem` <- function(e1, e2) {
  o_new <- binaryop.OTProblem(e1, e2, "divide")
  return(invisible(o_new))
  
}


OTProblem_$set("private", 
               "draw_boot_measure",
function() {
  addresses <- ls(private$measures)
  
  out <- rlang::env()
  n <- NULL
  prob <- NULL
  meas <- NULL
  
  for (add in addresses) {
    meas <- private$measures[[add]]
    n <- meas$n
    prob <- meas$init_weights
    out[[add]] <- prob$multinomial(n,replacement = TRUE)$add(1L)$bincount(minlength=n)
  }
  
  return(out)
}               
)

OTProblem_$set("private", 
               "eval_delta",
function (wts, boot) {
  addresses <- ls(private$target_objects)
  means     <- rep(0.0, length(addresses)) |> setNames(addresses)
  ns        <- rep(0.0, length(addresses)) |> setNames(addresses)
  n <- NULL
  w <- NULL
  b <- NULL
  m <- NULL
  s <- NULL
  w_star <- NULL
  target_obj<- NULL
  
  for (add in addresses) {
    if (private$measures[[add]]$adapt == "weights") {
      w <- wts[[add]]
    } else {
      w <- self$measure[[add]]$init_weights 
    }
    b <- boot[[add]]
    
    target_obj <- private$target_objects[[add]]
    n <- nrow(target_obj$bf)
    
    w_star <- w * b
    if(!target_obj$bf$requires_grad) {
      m <- target_obj$bf$transpose(2,1)$matmul( w_star )
      targ <- target_obj$bt
    } else {
      s <- target_obj$bf$detach()$std(1)
      m <- target_obj$bf$detach()$mean(1)/s
      targ <- target_obj$bt/s
    }
    
    means[[add]] <- (m - targ)$abs()$mean()$item()
    ns[[add]] <- n
  }
  
  return(weighted.mean(means, ns))
}               
               
)

OTProblem_$set("private", "eval_lambda",
function (wts, boot) {
  addresses <- ls(private$measures)
  prob_adds <- ls(private$ot_objects)
  wt_adds   <- ls(wts)
  sel_probs <- NULL
  prob_hold <- NULL
  p1_add    <- p2_add <- NULL
  
  n         <- NULL
  w         <- NULL
  b         <- NULL
  m         <- NULL
  w_star    <- NULL
  ot_obj    <- NULL
  meas      <- NULL
  
  for (add in addresses) {
   b <- boot[[add]]
   meas <- private$measures[[add]]
   
   w <- if(add %in% wt_adds) {
     wts[[add]]$detach()
   } else {
     meas$init_weights$detach()
   }
   
   if (meas$adapt == "weights" || meas$adapt == "none") {
     w_star <- w * b
   } else if (meas$adapt == "x") {
     w_star <- meas$init_weights * b
   } else {
     stop("Bug in this if else statement!")
   }
   
   sel_probs <- prob_adds[grep(add, prob_adds)]
   for (i in sel_probs) {
     prob_hold <- private$problems[[i]]
     p1_add <- prob_hold[[1]]
     p2_add <- prob_hold[[2]]
     if (p1_add == add){
       private$ot_objects[[i]]$a <- w_star
       if (meas$adapt == "x") {
         update_cost(private$ot_objects[[i]]$C_xy, 
                     w$detach(), #w should be the data values in this case, not weights
                     private$measures[[p2_add]]$x$detach())
       }
     } else if (p2_add == add) {
       private$ot_objects[[i]]$b <- w_star
       if (meas$adapt == "x") {
         update_cost(private$ot_objects[[i]]$C_yx, 
                     w$detach(), #w should be the data values in this case, not weights
                     private$measures[[p1_add]]$x$detach())
       }
     } else {
       stop("Problem address not found. You found a bug!")
     }
   }
   
  }
  if(is.finite(private$ot_objects[[i]]$penalty ))  {
    for(o in ls(private$ot_objects)) {
      private$ot_objects[[o]]$sinkhorn_opt(20, 1e-3)
    }
  }
  dists <- self$loss$detach()$item()
  return(dists)
}               
               
)


# operators to add weight environments for OTProblems
setOldClass("weightEnv")
binaryop.weightEnv <- function(e1,e2, fun) {
  
  we1 <- inherits(e1, "weightEnv")
  we2 <- inherits(e2, "weightEnv")
  
  if(isFALSE(we1) || isFALSE(we2)) {
    if(we1 && !we2) {
      w1 <- e1
      w2 <- e2
    } else if (!we1 && we2) {
      w1 <- e2
      w2 <- e1
    } else {
      stop("You found a bug!")
    }
    return(unaryop.weightEnv(w1, w2, fun))
  }
  
  listE1 <- ls(e1)
  listE2 <- ls(e2)
  
  stopifnot("weightEnv objects must have the same length." = length(listE1) == length(listE2))
  
  out <- rlang::env()
  class(out) <- c("weightEnv",class(out))
  for(i in seq_along(e1)) {
    v <- listE1[[i]]
    u <- listE2[[i]]
    out[[v]] <- fun(e1[[u]], e2[[v]])
  }
  
  return(out)
}

unaryop.weightEnv <- function(e1,e2, fun) {
  listE1 <- ls(e1)
  
  out <- rlang::env()
  class(out) <- c("weightEnv",class(out))
  for(e in listE1) {
    out[[e]] <- fun(e1[[e]], e2)
  }
  
  return(out)
}

#' @export
`+.weightEnv` <- function(e1, e2) {
  binaryop(e1, e2, `+`)
}

#' @export
`-.weightEnv` <- function(e1, e2) {
  binaryop(e1, e2, `-`)
}

#' @export
`*.weightEnv` <- function(e1, e2) {
  binaryop(e1, e2, `*`)
}

#' @export
`/.weightEnv` <- function(e1, e2) {
  binaryop(e1, e2, `/`)
}


#' @export
sum.weightEnv <- function(..., na.rm = FALSE) {
  out <- 0.0
  l   <- list(...)
  if (length(l) == 1) {
    e1 <- l[[1]]
    listE1 <- ls(e1)
    for (e in listE1) {
      out <- sum(e1[[e]]) + out
    }
  } else {
    for (i in seq_along(l) ) {
      out <- sum(l[[i]]) + out
    }
  }
  
  return(out)
}


#### Optimizers relying on OTProblem class ####
#### NNM class ####
# maybe make a part of OTProblem...unclear how to unite several problems though...
NNM <- R6::R6Class(
  classname = "NNM",
  public = {list(
    setup_arguments = function(lambda = NULL, delta = NULL , 
                               grid.length = 7L,
                               cost.function = NULL, 
                               p = 2,
                               cost.online = "auto",
                               debias = TRUE,
                               diameter = NULL, ot_niter = 1000L,
                               ot_tol = 1e-5) {
      super$setup_arguments(lambda = 0, delta = NULL,
                            grid.length = 1,
                            cost.function = cost.function,
                            p = p,
                            cost.online = cost.online,
                            debias = FALSE,
                            diameter = diameter,
                            ot_niter = ot_niter,
                            ot_tol = ot_tol)
      ot <- private$ot_objects[[ls(private$ot_objects)[[1]] ]]
      # private$device <- ot$device
      private$C_xy <- ot$C_xy
      private$tensorized <- ot$tensorized
      private$b <- ot$b
      private$n <- as.integer(ot$n)
    },
    solve = function(...) {
      C_xy <- private$C_xy
      if (!private$tensorized) { 
        
        x = as_matrix(C_xy$data$x)
        y = as_matrix(C_xy$data$y)
        d = ncol(x)
        dim.red <- if(utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
          "0"
        } else {
          "1"
        }
        argmin_op <- rkeops::keops_kernel(
          formula = paste0("ArgMin_Reduction(", C_xy$fun, ",",dim.red,")"),
          args = c(
            paste0("X = Vi(",d,")"),
            paste0("Y = Vj(",d,")"))
        )
        mins = torch::torch_tensor(c(argmin_op(list(x,y))) + 1, 
                                   dtype = torch::torch_int64(),
                                   device = private$b$device)
        
      } else {
        mins = C_xy$data$argmin(1)
      }
      w_nnm = torch::torch_bincount(self = mins, weights = private$b, minlength = private$n)
      
      for (m in ls(private$measures) ) {
        if(private$measures[[m]]$adapt == "weights") private$measures[[m]]$weights <- w_nnm
      } 
      
      # return(w_nnm)
    },
    choose_hyperparameters = function(...) {
      self$selected_lambda <- 0.0
    }
  )},
  private = list(
    b = "tensor",
    C_xy = "cost",
    device = "torch_device",
    n = "integer",
    tensorized = "logical"
  ),
  inherit = OTProblem_
)

#### Primal optimizer ####
# this uses the cotOOP paradigm, OTProblem
# therefore, no individual functions

#### Dual function optimizer ####

# forward functions in torchscript
dual_forward_code_tensorized <- "

 def calc_w1(f: Tensor, C_xy: Tensor, a_log: Tensor, b_log: Tensor, lambda: float, n: int):
   f_minus_C_lambda = f.view([n,1]) + a_log.view([n,1]) - C_xy/lambda #  f/lambda - C/lambda + a_log
   g_lambda = b_log -(f_minus_C_lambda).logsumexp(0) # = g/lambda + b_log
   w1 = (f_minus_C_lambda + g_lambda).logsumexp(1).exp()
   return w1
   
 def calc_w2(f: Tensor, C_xx: Tensor, a_log: Tensor, lambda: float, n: int):
   f_lambda = f + a_log
   w2 = ( f_lambda.view([n,1]) + f_lambda  - C_xx / lambda).logsumexp(1).log_softmax(0).exp()
   return w2
   
 def cot_dual(gamma: Tensor, C_xy: Tensor, C_xx: Tensor, a_log: Tensor, b_log: Tensor, lambda: float, n: int):
    
   f_star = gamma.detach()
   
   w1 = calc_w1(f_star, C_xy, a_log, b_log, lambda, n)
   w2 = calc_w2(f_star, C_xx, a_log,        lambda, n)
   
   measure_diff = (w1-w2).detach()
   loss =  -1.0 * gamma.dot(measure_diff) 
   # mult by -1 because is a maximization and are turning into minimization
   # print(-loss.item())
   return {'loss' : loss, 'avg_diff' : measure_diff.detach().norm(), 'bf_diff' : torch.zeros(1,dtype=loss.dtype)}
   
 def cot_bf_dual(gamma: Tensor, C_xy: Tensor, C_xx: Tensor, a_log: Tensor, b_log: Tensor, lambda: float, n: int, beta: Tensor, bf: Tensor, bt: Tensor, delta: float):
   
   w1 = calc_w1(gamma.detach(), C_xy, a_log, b_log, lambda, n)
   w2 = calc_w2(gamma.detach(), C_xx, a_log,        lambda, n)
   
   measure_diff = (w1-w2).detach()
   loss_gamma = gamma.dot(measure_diff) 
   
   bf_diff = bf.transpose(0,1).matmul(w1) - bt
   
   beta_check = bf_diff * beta.detach() - delta * beta.detach().abs()
   
   loss_beta = bf_diff.dot(beta) - beta.abs().sum() * delta
   
   loss = (loss_gamma + loss_beta) * -1.0 # mult by neg 1 because is a maximization
                                   
   return {'loss' : loss, 'avg_diff' : measure_diff.detach().norm(), 'bf_diff' :  bf_diff.detach().abs().max(), 'beta_check' : beta_check}
   
"

# rkeops forward functions
dual_forwards_keops <- list(
  calc_w1 = function(f, C_xy, a_log, b_log, lambda, n) {
    xmat <- as.matrix(C_xy$data$x$to(device = "cpu"))
    ymat <- as.matrix(C_xy$data$y$to(device = "cpu"))
    f_lambda <- f + a_log
    # f_lambda <- f
    exp_sums_g <- C_xy$reduction( list(ymat, xmat,  
                                       as.numeric(f_lambda$to(device = "cpu")),
                                       1.0 / lambda) )
    
    g_lambda <- if (utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
       - c(exp_sums_g) + as_numeric(b_log)
    } else {
       - (log(exp_sums_g[,2]) + exp_sums_g[,1]) + as_numeric(b_log)
    }
    
    exp_sums_a1 <- C_xy$reduction( list(xmat, ymat, 
                                        g_lambda,
                                        1.0 / lambda) )
    a1_log <- if (utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
        torch::torch_tensor(c(exp_sums_a1), dtype = f$dtype, device = f$device) + f_lambda
      } else {
        torch::torch_tensor(log(exp_sums_a1[,2]) + exp_sums_a1[,1], dtype = f$dtype, device = f$device) + f_lambda
      }
    return(a1_log$log_softmax(1)$exp()$view(-1))
  },
  calc_w2 = function(f, C_xy, a_log, lambda, n) {
    f_lambda <- f + a_log
    # f_lambda <- f
    xmat <- as.matrix(C_xy$data$x$to(device = "cpu"))
    # ymat <- as.matrix(C_xy$data$x$to(device = "cpu"))
    
    if (utils::packageVersion("rkeops") >= pkg_vers_number("2.0")) {
      log_exp_sums_a2 <- C_xy$reduction( list(xmat, xmat,
                                          as.numeric(f_lambda$to(device = "cpu")),
                                          1.0 / lambda) )
      a2_log <-  torch::torch_tensor(c(log_exp_sums_a2), dtype = f$dtype, device = f$device) + f_lambda
    } else {
        exp_sums_a2 <- C_xy$reduction( list(xmat, xmat,
                                        as.numeric(f_lambda$to(device = "cpu")),
                                        1.0 / lambda) )
        a2_log <-  torch::torch_tensor(log(exp_sums_a2[,2]) + exp_sums_a2[,1], dtype = f$dtype, device = f$device) + f_lambda
      }
    
    return(a2_log$log_softmax(1)$exp()$view(-1))
  },
  cot_dual = function(gamma, C_xy, C_xx, a_log, b_log, lambda, n) {
    f_star = gamma$detach() #+ a_log
    
    w1 = dual_forwards_keops$calc_w1(f_star, C_xy, a_log, b_log, lambda, n)$detach()$to(device = gamma$device)
    w2 = dual_forwards_keops$calc_w2(f_star, C_xx, a_log, lambda, n)$detach()$to(device = gamma$device)
    
    measure_diff = w1-w2
    loss = -1.0 * gamma$dot(measure_diff) # mult by neg 1 because is a maximization
    
    return(list(
      loss = loss,
      avg_diff = measure_diff$norm(),
      bf_diff = torch::torch_zeros(1, dtype = loss$dtype)
    ))
  },
  cot_bf_dual = function(gamma, C_xy, C_xx, a_log, b_log, lambda, n,
                         beta, bf, bt, delta) {
    f_star = gamma$detach() #+ a_log
    beta_d = beta$detach()
    # f_star1 = f_star2 - bf$matmul(beta$detach()) 
    
    w1 = dual_forwards_keops$calc_w1(f_star, C_xy, a_log, b_log, lambda, n)$detach()$to(device = gamma$device)
    w2 = dual_forwards_keops$calc_w2(f_star, C_xx, a_log, lambda, n)$detach()$to(device = gamma$device)
    
    measure_diff = w1-w2
    loss_gamma = gamma$dot(measure_diff) # mult by neg 1 because is a maximization
    
    bf_diff = bf$transpose(1,2)$matmul(w1) - bt
    
    beta_check = bf_diff * beta_d - delta * beta_d$abs()
    
    loss_beta = bf_diff$dot(beta) - delta * beta$abs()$sum()
    
    loss = (loss_gamma + loss_beta) * -1.0
    
    return(list(
      loss = loss,
      avg_diff = measure_diff$abs()$mean(),
      bf_diff = bf_diff$abs()$max(),
      beta_check = beta_check
    ))
  }
  
)



# optimizer without bf
cotDualOpt <- torch::nn_module(
  classname = "cotDualOpt",
  initialize = function(n, d = NULL, device = NULL, dtype = NULL) {
    
    self$device <- cuda_device_check(device)
    self$dtype  <- cuda_dtype_check(dtype, self$device)
    
    self$n <- torch::jit_scalar(as.integer(n))
    self$gamma <- torch::nn_parameter(
      torch::torch_zeros(self$n,
                         dtype = self$dtype, 
                         device = self$device),
      requires_grad = TRUE)
    private$set_forward(bf = FALSE)
  },
  forward = function(C_xy, C_xx, a_log, b_log, lambda, bf=NULL, bt=NULL, delta = NULL) {
    private$ts_forward(self$gamma, C_xy$data, C_xx$data, a_log, b_log, torch::jit_scalar(lambda), self$n)
  },
  backward = function(res) {
    res$loss$backward()
  },
  clone_param = function(requires_grad = FALSE) {
    if(isTRUE(requires_grad) ) {
      return(self$parameters)
    } else {
      param <- self$parameters
      out_param <- vector("list", length(param))
      names(out_param) <- names(param)
      
      for(i in seq_along(param) ) {
        out_param[[i]] <- param[[i]]$detach()$clone()
      }
      return(out_param)
    }
  },
  converged = function(res, avg_diff_old, loss_old, old_param,
                        tol=1e-5, lambda, delta) {
    machine_tol <- .Machine$double.eps
    if(is.na(delta) || is.null(delta)) delta <- 0
    
    avg_diff <- res$avg_diff$item()
    loss <- res$loss$detach()$item()
    # param <- self$parameters
    bf_diff <- res$bf_diff$item()
    
    # rel_param_sq_sum <- torch::torch_zeros(1, dtype = torch::torch_double())
    # for (i in seq_along(param)) {
    #   rel_param_sq_sum$add_( ((param[[i]]$detach() - old_param[[i]])/(old_param[[i]] + machine_tol))$norm()$square() )
    # }
    
    sq_tol <- tol * tol
    
    abs_avg_diff  <- abs(avg_diff - avg_diff_old)
    abs_loss      <- abs(loss - loss_old)
    
    abs_loss_check <- abs_loss < sq_tol
    abs_avg_diff_check <- abs_avg_diff < sq_tol
    
    abs_bf_diff <- (bf_diff - delta)
    rel_bf_diff <- abs_bf_diff/(delta + machine_tol)
    
    rel_avg_diff_check <- abs_avg_diff/(avg_diff + machine_tol) < tol
    rel_loss_check <- abs_loss/(abs(loss) + machine_tol) < tol
    # rel_param_sq_check <- rel_param_sq_sum$item() < sq_tol
    
    
    must_pass <- avg_diff < 1/lambda &&  (rel_bf_diff <= tol || abs_bf_diff <= min(sq_tol, delta/1e4) )
    
    rel_checks <- rel_loss_check || rel_avg_diff_check #|| rel_param_sq_check
    abs_checks <- abs_loss_check || abs_avg_diff_check
    
    return ( must_pass && (rel_checks || abs_checks ) )
  },
  calc_w1 =  "torch_jit",
  calc_w2 =  "torch_jit",
  dtype = "torch_dtype",
  device = "torch_dtype",
  private = list(
    ts_forward = "function",
    set_forward = function(bf = FALSE) {
      dual_forwards <- private$dual_forwards()
      private$ts_forward = switch(bf +1L,
                                  dual_forwards$cot_dual,
                                  dual_forwards$cot_bf_dual)
      self$calc_w1 =  dual_forwards$calc_w1
      self$calc_w2 =  dual_forwards$calc_w2
    },
    dual_forwards = function() {
      torch::jit_compile(dual_forward_code_tensorized)
    }
  )
)

cotDualOpt_keops <- torch::nn_module(
  classname = "cotDualOpt_keops",
  inherit = cotDualOpt,
  forward = function(C_xy, C_xx, a_log, b_log, lambda, bf, bt, delta) {
    private$ts_forward(self$gamma, C_xy, C_xx, a_log, b_log, torch::jit_scalar(lambda), self$n)
  },
  private = list(
    set_forward = function(...) {
      private$ts_forward = dual_forwards_keops$cot_dual
      self$calc_w1 = dual_forwards_keops$calc_w1
      self$calc_w2 = dual_forwards_keops$calc_w2
    }
  )
)

# optimizer with bf
cotDualBfOpt <- torch::nn_module(
  classname="cotDualBfOpt",
  inherit = cotDualOpt,
  initialize = function(n, d, device = NULL, dtype = NULL) {
    self$device <- cuda_device_check(device)
    self$dtype  <- cuda_dtype_check(dtype, self$device)
    
    self$n <- torch::jit_scalar(as.integer(n))
    self$d <- torch::jit_scalar(as.integer(d))
    self$gamma <- torch::nn_parameter(
      torch::torch_zeros(self$n, 
                         dtype = self$dtype,
                         device = self$device), requires_grad = TRUE)
    self$beta <- torch::nn_parameter(
      torch::torch_zeros(self$d, 
                         dtype = self$dtype,
                         device = self$device), requires_grad = TRUE)
    private$set_forward(bf = TRUE)
  },
  forward = function(C_xy, C_xx, a_log, b_log, lambda, bf, bt, delta) {
    torch::with_no_grad(self$gamma$sub_(bf$matmul(self$beta)))
    res <- private$ts_forward(self$gamma, C_xy$data, C_xx$data, a_log, b_log, torch::jit_scalar(lambda), self$n,
                       self$beta, bf, bt, torch::jit_scalar(delta))
    return(res)
  },
  backward = function(res) {
    res$loss$backward()
    torch::with_no_grad(self$beta$mul_(res$beta_check > 0))
  }
)

cotDualBfOpt_keops <- torch::nn_module(
  classname = "cotDualBfOpt_keops",
  inherit = cotDualBfOpt,
  forward = function(C_xy, C_xx, a_log, b_log, lambda, bf, bt, delta) {
    torch::with_no_grad(self$gamma$sub_(bf$matmul(self$beta)))
    res <- private$ts_forward(self$gamma, C_xy, C_xx, a_log, b_log, torch::jit_scalar(lambda), self$n,
                       self$beta, bf, bt, torch::jit_scalar(delta))
    return(res)
  },
  private = list(
    set_forward = function(...) {
      private$ts_forward = dual_forwards_keops$cot_bf_dual
      self$calc_w1 = dual_forwards_keops$calc_w1
      self$calc_w2 = dual_forwards_keops$calc_w2
    }
  )
)

# main object to do training, inherit from OTproblem and makes changes where needed
cotDualTrain <- R6::R6Class(
  classname = "cotDualTrain",
  public = {list(
    ot = "OT",
    setup_arguments = function(lambda = NULL, delta = NULL , 
                               grid.length = 7L,
                               cost.function = NULL, 
                               p = 2,
                               cost.online = "auto",
                               debias = TRUE,
                               diameter = NULL, ot_niter = 1000L,
                               ot_tol = 1e-5) {
      super$setup_arguments(lambda, delta, 
                            grid.length,
                            cost.function, 
                            p,
                            cost.online ,
                            debias,
                            diameter , ot_niter,
                            ot_tol)
      
      m_add   <- ls(private$measures)
      p_add   <- ls(private$problems)
      adapt_m <- private$problems[[ p_add[[1]] ]][1]
      targ_m  <- private$problems[[ p_add[[1]] ]][2]
      
      self$ot <- private$ot_objects[[ p_add[1] ]]
      
      private$C_xy <- self$ot$C_xy
      private$C_xx <- self$ot$C_xx
      # private$a_log <- log_weights(self$ot$a$detach())
      # private$b_log <- log_weights(self$ot$b$detach())
      # makes sure the non- changing weights are used...
      stopifnot("Wrong measure is being fed for adapatation. Report this bug please!" = private$measures[[ adapt_m ]]$requires_grad)
      private$a_log <- log_weights(private$measures[[ adapt_m ]]$init_weights$detach())
      private$b_log <- log_weights(private$measures[[ targ_m  ]]$init_weights$detach())
      
      runbf <-  length(private$target_objects) >= 1
      if (runbf) {
        t_add <- ls(private$target_objects)
        targ  <- private$target_objects[[ t_add[1] ]]
        
        private$bf <- targ$bf
        private$bt <- targ$bt
        
      } else {
        private$bf <- private$bt <- NULL
      }
      
      tensorized <- self$ot$tensorized
      nn_fun <- switch(tensorized * 2 + runbf + 1L ,
                       cotDualOpt_keops,
                       cotDualBfOpt_keops,
                       cotDualOpt,
                       cotDualBfOpt
      )
      private$nn_holder  <- nn_fun$new(n = self$ot$n, 
                                       d = length(private$bt),
                                       device = private$measures[[m_add[1L]]]$device,
                                       dtype = private$measures[[m_add[1L]]]$dtype)
      private$parameters <- private$nn_holder$parameters
      
      private$penalty_list$lambda[ is.infinite(private$penalty_list$lambda) ] <- self$ot$diameter * 1e5
      private$penalty_list$lambda[private$penalty_list$lambda == 0] <- self$ot$diameter / 1e9
      # runs faster and more accurately when reversed (small -> large)
      private$penalty_list$lambda <- sort(private$penalty_list$lambda, decreasing = FALSE)
      
      private$lambda <- private$penalty_list$lambda[1L]
      
      # private$niter <- private$ot_niter
      # private$tol <- private$ot_tol
      private$prev_lambda <- private$lambda
      return(invisible(self))
      
    }
  )},
  active = {list(
    weights = function(value) {
      f1 <- private$nn_holder$gamma$detach()$clone() #+ private$a_log
      # f2 <- private$nn_holder$gamma$detach()
      if (!is.null(private$bf)) {
        f1$sub_(private$bf$matmul(private$nn_holder$beta$detach()))
      } 
      if (inherits(private$C_xy, "costTensor") ){
        C_xy <- private$C_xy$data
        C_xx <- private$C_xx$data
      } else {
        C_xy <- private$C_xy
        C_xx <- private$C_xx
      }
      w1 <- private$nn_holder$calc_w1(f1, C_xy, private$a_log, private$b_log, torch::jit_scalar(private$lambda), torch::jit_scalar(private$nn_holder$n))
      w2 <- private$nn_holder$calc_w2(f1, C_xx, private$a_log, torch::jit_scalar(private$lambda), torch::jit_scalar(private$nn_holder$n))
      w  <- (w1 + w2) * 0.5
      # return(list(a = a, a1 = a1, a2 = a2))
      return(w)
    }
  )},
  private = {list(
    eval_lambda = function (wts, boot) {
      addresses <- ls(private$measures)
      prob_adds <- ls(private$ot_objects)
      wt_adds   <- ls(wts)
      sel_probs <- NULL
      prob_hold <- NULL
      p1_add    <- p2_add <- NULL
      
      n         <- NULL
      w         <- NULL
      b         <- NULL
      m         <- NULL
      w_star    <- NULL
      ot_obj    <- NULL
      meas      <- NULL
      
      for (add in addresses) {
        b    <- boot[[add]]
        meas <- private$measures[[add]]
        
        w <- if(add %in% wt_adds) {
          wts[[add]]$detach()
        } else {
          meas$init_weights$detach()
        }
        
        if (meas$adapt == "weights" || meas$adapt == "none") {
          w_star <- w * b
        } else if (meas$adapt == "x") {
          w_star <- meas$init_weights * b
        } else {
          stop("Bug in this if else statement!")
        }
        
        sel_probs <- prob_adds[grep(add, prob_adds)]
        for (i in sel_probs) {
          prob_hold <- private$problems[[i]]
          p1_add <- prob_hold[[1]]
          p2_add <- prob_hold[[2]]
          if (p1_add == add){
            private$ot_objects[[i]]$a <- w_star
            if (meas$adapt == "x") {
              update_cost(private$ot_objects[[i]]$C_xy, w, private$measures[[p2_add]]$x$detach())
            }
          } else if (p2_add == add) {
            private$ot_objects[[i]]$b <- w_star
            if (meas$adapt == "x") {
              update_cost(private$ot_objects[[i]]$C_yx, w, private$measures[[p1_add]]$x$detach())
            }
          } else {
            stop("OT Problem address not found. You found a bug!")
          }
        }
        
      }
      # browser()
      # private$ot_objects[[i]]$penalty <- private$boot_lambda #0.05 #private$ot_objects[[i]]$diameter
      if(is.finite(private$ot_objects[[i]]$penalty ))  private$ot_objects[[i]]$sinkhorn_opt(20, 1e-3)
      dists <- self$loss$detach()$item()
      return(dists)
    },
    set_lambda = function (l) {
      stopifnot(l >= 0.0)
      stopifnot(l <= Inf)
      if(l == 0.0) l <- self$ot$diameter / 1e9
      super$set_lambda(l)
      if(is.infinite(l)) l <- self$ot$diameter * 1e4
     
      private$lambda <- torch::jit_scalar(l)
    },
    set_delta = function (d) {
      stopifnot(d >= 0.0 || is.na(d))
      private$delta <- if(!is.na(d)) {
        torch::jit_scalar(d)
      } else if(length(private$target_objects) >= 1L) {
         private$target_objects[[ls(private$target_objects)]]$delta
      } else {
        NA_real_
      }
      
    },
    set_penalties = function (value) {
      if (is.list(value)) {
        if (!all(names(value) %in% c("lambda", "delta")) ) {
          stop("If penalties are provided as a list, must be with names in c('lamba', 'delta')")
        }
        lambda <- value$lambda
        delta <- value$delta
      } else {
        if (any(!is.null(names(value))) ) {
          if (!all(names(value) %in% c("lambda", "delta")) ) {
            stop("If penalties are provided as a named vector, must be with names in c('lamba', 'delta')")
          }
          lambda <- value["lambda"]
          delta <- value["delta"]
          names(lambda) <- NULL
          names(delta) <- NULL
        } else {
          lambda <- value[1]
          delta <- NA_real_
          if(length(value) > 1) warning("For unnamed vectors, only the first number is used to set the OT penalty, lambda. Other values are ignored.")
        }
      }
      private$set_lambda(lambda)
      private$set_delta(delta)
    },
    optimization_loop = function (niter, tol) {
      
      #reset torch optimizer if present
      private$torch_optim_reset()
      
      # reset the parameters to 0, speeds up estimation for lower lambdas
      # torch::with_no_grad(private$nn_holder$gamma$mul_(private$lambda/private$prev_lambda))
      # torch::with_no_grad(private$nn_holder$gamma$copy_(0.0))
      if(!is.null( private$nn_holder$beta) ) torch::with_no_grad(private$nn_holder$beta$copy_(0.0))
      
      avg_diff_old <- loss_old <- 10.0
      old_param <- NULL #private$nn_holder$clone_param()
      
      for (i in 1:niter) {
        private$opt$zero_grad()
        res <- private$nn_holder$forward(private$C_xy, private$C_xx, 
                                         private$a_log, 
                                         private$b_log, 
                                         private$lambda, private$bf, private$bt, private$delta)
        # print(res$loss$detach()$item())
        if ((i > 2L) && private$nn_holder$converged(res, avg_diff_old, 
                                                    loss_old, old_param,
                                        tol,
                                        private$lambda, private$delta)) {
          break
        } else {
          avg_diff_old <- res$avg_diff$item()
          loss_old <- res$loss$detach()$item()
          # old_param <- private$nn_holder$clone_param()
        }
        
        private$nn_holder$backward(res)
        private$opt$step()
        check <- private$lr_reduce(avg_diff_old)
        # private$sched$step()
        
      }
      private$prev_lambda <- private$lambda
      return(list(loss = res$loss$detach()$item(),
                  iter = i))
      
    }, # overwrite super function
    parameters_get_set = function (value, clone = FALSE) {
      
      ifthenfun <- function(v) {
        if (is.list(private$parameters[[v]]) ) {
          private$parameters[[v]]$params
        } else {
          private$parameters[[v]]
        }
      }
      
      # return a clone of the weights
      if ( missing(value) ) {
        out <- rlang::env()
        class(out) <- c("weightEnv", class(out))
        
        param_addresses <- ls(private$measures)
        for (v in param_addresses) {
          # out[[v]] <- if (isTRUE(clone)) {
          #   ifthenfun(v)$detach()$clone()
          # } else {
          #   ifthenfun(v)
          # }
          if(private$measures[[v]]$adapt == "weights") {
            out[[v]] <- self$weights
          }
         
        }
        return(out)
      }
      
      # set weights
      param_addresses <- ls(private$measures)
      
      if(!rlang::is_environment(value)){
        stopifnot("value must be a list or environment" = is.list(value))
        names(value) <- value_addresses <- 1:length(value)
      } else {
        value_addresses <- ls(value)
      }
      
      if(length(value_addresses) == 1) {
        v <- NULL
        u <- NULL
        torch::with_no_grad(
          for ( i in seq_along(param_addresses) ) {
            v <- param_addresses[[i]]
            u <- value_addresses[[1L]]
            if (private$measures[[v]]$adapt == "weights")  {
              private$measures[[v]]$weights <- value[[u]]
            }
          }
        )
      } else {
        stop("Input must have same number of groups as do the parameters.")
      }
    }, #overwrite super
    torch_optim_setup = function (torch_optim, 
                                 torch_scheduler,
                                 torch_args) {
      
      names_args <- names(torch_args)
      optim_args_names <- names(formals(torch_optim))
      opt_args <- torch_args[match(optim_args_names,
                                   names_args, nomatch = 0L)]
      param <- private$parameters
      gamma_lr <- if(!is.null(opt_args$lr)) {
        opt_args$lr
        } else {
          1e-2 #private$lambda/100
        }
      param$gamma <- if(!is.list(param$gamma) ) {
        list(params = param$gamma, lr = gamma_lr)
      } else {
        list(params = param$gamma$params, lr = gamma_lr)
      }
      if(!is.null(param$beta)) {
        param$beta <- if(!is.list(param$beta)) {
          list(params = param$beta,
                           lr = min(private$delta, min(opt_args$lr, 1e-2)))
        } else {
          list(params = param$beta$params,
               lr = min(private$delta, min(opt_args$lr, 1e-2)))
        }
      }
      private$parameters <- param
      opt_call <- rlang::call2(torch_optim,
                               params = private$parameters,
                               !!!opt_args)
      private$opt <- eval(opt_call)
      
      torch_lbfgs_check(private$opt)
      if (inherits(private$opt, "optim_lbfgs")) {
        warning("Torch's LBFGS optimizer does not work well on the dual problem. Please use another optimizer.")
      }
      
      if (!is.null(torch_scheduler)) {
        scheduler_args_names <- names(formals(torch_scheduler))
        sched_args <- torch_args[match(scheduler_args_names, 
                                       names_args, 
                                       nomatch = 0L)]
        if(inherits(torch_scheduler, "lr_reduce_on_plateau") &&
           is.null( sched_args$patience )  && inherits(private$opt, "optim_lbfgs")) {
          sched_args$patience <- 1L
        }
        if(inherits(torch_scheduler, "lr_reduce_on_plateau") && is.null( sched_args$min_lr ) && inherits(private$opt, "optim_lbfgs") ) {
          sched_args$min_lr <- private$opt$defaults$lr * 1e-3
        }
        if(inherits(torch_scheduler, "lr_multiplicative") &&
           is.null( sched_args$lr_lambda ) ) {
          sched_args$lr_lambda <- function(epoch) {0.99}
        }
        scheduler_call <- rlang::call2(torch_scheduler,
                                       optimizer = private$opt,
                                       !!!sched_args)
        opt_sched <- eval(scheduler_call)
      } else {
        opt_sched <- scheduler_call <- NULL
      }
      
      private$opt_calls <- list(opt = opt_call,
                                sched = scheduler_call)
      
      private$sched <- opt_sched
      
      # return(list(opt = opt, opt_call = opt_call,
      #             sched = opt_sched, sched_call = scheduler_call))
      
    },
    torch_optim_reset = function (lr = NULL) {
      # browser()
      if(!is.null(private$opt)) {
        default_lr <- private$opt$defaults$lr
        opt_call <- private$opt_calls$opt
        
        if (!is.null(lr)) {
          private$parameters$gamma$lr <- lr
        # } else if(is.null(opt_call$lr)) {
        #   private$parameters$gamma$lr <- private$lambda/100
        } else if (!is.null(opt_call$lr)) {
          private$parameters$gamma$lr <- opt_call$lr
        } else {
          private$parameters$gamma$lr <- default_lr
        }
          
        if( is.null(lr)) lr <- default_lr
        
        if(!is.null(private$parameters$beta)) {
          if(!is.na(private$delta)) {
            lr_new <- min(private$delta, lr)
          } else {
            lr_new <- lr
          }
          private$parameters$beta$lr <- lr_new
        }
        
        if(is.null(lr)) {
          private$opt <- eval(rlang::call_modify(opt_call, params = private$parameters))
        } else {
          # browser()
          # def <- rlang::call_match(opt_call[[1]], opt_call, defaults = TRUE)
          private$opt <- eval(rlang::call_modify(opt_call, params = private$parameters,
                                                 lr = lr))
        }
        
      }
      
      if(!is.null(private$sched)) {
        private$sched <- eval(rlang::call_modify(private$opt_calls$sched, 
                                                 optimizer = private$opt))
      }
      
    },
    ot_update = function (...) {NULL},
    lambda = "numeric",
    delta = "numeric",
    # tol = "numeric",
    nn_holder = "dualCotOpt",
    # niter = "integer",
    optim = "optim",
    prev_lambda = "numeric",
    sched = "scheduler",
    C_xy = "torch_tensor",
    C_xx = "torch_tensor",
    a_log = "torch_tensor",
    b_log = "torch_tensor",
    bf = "torch_tensor",
    bt = "torch_tensor"
  )},
  inherit = OTProblem_
)
