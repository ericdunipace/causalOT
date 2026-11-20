sinkhorn_compiled <- local({
  compiled <- NULL
  function() {
    if (is.null(compiled)) {
      compiled <<- torch::jit_compile(
        "
def online_center(x: Tensor, y: Tensor, p: float) -> Tuple[Tensor, Tensor]:
    if p == 2.0:
        mu = (x.sum(0, keepdim=True) + y.sum(0, keepdim=True)) / float(x.size(0) + y.size(0))
        return x - mu, y - mu
    return x, y
      
def online_cost_block(x: Tensor, y: Tensor, p: float) -> Tensor:
    if p == 2.0:
        C = (x * x).sum(1).view(-1, 1) / 2.0 + (y * y).sum(1).view(1, -1) / 2.0 - x.matmul(y.t())
        return C.clamp_min(0.0)
    if p == 1.0:
        return torch.cdist(x, y, p=1.0)
    return torch.cdist(x, y, p=p).pow(p) / p

def online_ratio_ok(r: Tensor) -> bool:
    # need ratio neer 1; NaN should fail as well.
    return float(r.min()) > 1e-30 and float(r.max()) < 1e30

def ridge_scale(A: Tensor) -> float:
    # ridge at a few levels above machine precision
    if A.element_size() == 8:  # float64
        return 1e-10
    return 1e-5
      
def ridge(X: Tensor, y: Tensor, scale: float) -> Tensor:
    p = X.size(1)

    I = torch.eye(p, dtype=X.dtype, device=X.device)
    A = X.t().matmul(X) + (ridge_scale(X) * scale) * I
    return torch.linalg_solve(A, X.t().mv(y))


def anderson_step(Tf: Tensor, r: Tensor, F_hist: List[Tensor], R_hist: List[Tensor],
                  m: int) -> Tensor:
    F_hist.append(Tf) # f potentials
    R_hist.append(r) # diff between iteration n and n+1
    if len(R_hist) > m + 1: # remove history 
        F_hist.pop(0)
        R_hist.pop(0)
    n_hist = len(R_hist)
    if n_hist == 1:
        return Tf
    dR = torch.stack([R_hist[j] - R_hist[j - 1] for j in range(1, n_hist)], 1)
    dF = torch.stack([F_hist[j] - F_hist[j - 1] for j in range(1, n_hist)], 1)
    scale = float((dR * dR).sum(dim=0).max())
    if not (scale > 0.0):  # if residuals 0, return original tensor
        return Tf
    gamma = ridge(dR, r, scale)
    return Tf - dF.matmul(gamma)
      
def softmin_fixed_K(eps: float, K_xy: Tensor, y_potential: Tensor, b_log: Tensor):
            return -eps * (y_potential/eps + b_log - K_xy).logsumexp(1)

def softmin(eps: float, C_xy: Tensor, y_potential: Tensor, b_log: Tensor):
            return -eps * (y_potential/eps + b_log - C_xy/eps).logsumexp(1)

def online_softmin(x: Tensor, y: Tensor, G: Tensor, P: float, p: float, rows: int) -> Tensor:
    # out_i = logsumexp_j (G_j - P * C(x_i, y_j))
    # G = potential/epsilon + a_log
    # P = penalty
    xc = x
    yc = y
    # xc, yc = online_center(x, y, p) # only used for p = 2 accuracy; now moved into R
    n = xc.size(0)
    out = torch.empty(n, dtype=x.dtype, device=x.device)
    if p == 2.0: # p = 2 gets extra efficiency
        yt = (yc * P).t().contiguous()
        bias = (G - P * (yc * yc).sum(1) / 2.0).view(1, -1)
        xx = P * (xc * xc).sum(1) / 2.0
        for i in range(0, n, rows):
            E = torch.addmm(bias, xc[i:i + rows], yt)
            out[i:i + rows] = E.logsumexp(1) - xx[i:i + rows]
    else:
        Gv = G.view(1, -1)
        for i in range(0, n, rows):
            out[i:i + rows] = (Gv - P * online_cost_block(xc[i:i + rows], yc, p)).logsumexp(1)
    return  -out / P

def softmin_sweep(e: float, X: Tensor, Y: Tensor, f_xy: Tensor, a_log: Tensor,
                   b_log: Tensor, b: Tensor, p: float,
                   rows_x: int, rows_y: int, online: bool) -> Tuple[Tensor, Tensor]:
    # Dense: X, Y are K_xy, K_yx (cost / e). Online: X, Y are the points x, y.
    if online:
        g_yx = online_softmin(Y, X, f_xy / e + a_log, 1.0 / e, p, rows_y)
    else:
        g_yx = softmin_fixed_K(e, Y, f_xy, a_log)
    g_yx = g_yx - b.dot(g_yx)
    if online:
        Tf = online_softmin(X, Y, g_yx / e + b_log, 1.0 / e, p, rows_x)
    else:
        Tf = softmin_fixed_K(e, X, g_yx, b_log)
    return Tf, g_yx
      
def online_ratio_pass(x: Tensor, y: Tensor, f: Tensor, g: Tensor, a: Tensor, b: Tensor,
               e: float, p: float, rows: int, need_acc: bool) -> Tuple[Tensor, Tensor]:
    # K_ji = exp((f_i + g_j - C(x_i, y_j)) / e), one block of rows of y at a time.
    # r_j = sum_i K_ji a_i;  acc_i = sum_j K_ji b_j / r_j (only if need_acc).
    xc = x
    yc = y
    # xc, yc = online_center(x, y, p) # centering now down in R one time
    if p == 2.0:
        # (f_i + g_j - C_ij)/e = x_i.y_j/e + (f_i - |x_i|^2/2)/e + (g_j - |y_j|^2/2)/e
        xt = (xc / e).t().contiguous()
        fx = ((f - (xc * xc).sum(1) / 2.0) / e).view(1, -1)
        gy = ((g - (yc * yc).sum(1) / 2.0) / e).view(-1, 1)
    else:
        xt = xc
        fx = f.view(1, -1)
        gy = g.view(-1, 1)
    r = torch.empty(yc.size(0), dtype=x.dtype, device=x.device)
    acc = torch.zeros(xc.size(0), dtype=x.dtype, device=x.device)
    for j in range(0, yc.size(0), rows):
        if p == 2.0:
            K = torch.addmm(fx, yc[j:j + rows], xt).add_(gy[j:j + rows]).exp_()
        else:
            K = ((gy[j:j + rows] + fx - online_cost_block(yc[j:j + rows], xc, p)) / e).exp_()
        rj = K.mv(a)
        r[j:j + rows] = rj
        if need_acc:
            acc += K.t().mv(b[j:j + rows] / rj)
    return r, acc

def dense_ratio_pass(K_xy: Tensor, K_yx: Tensor, f: Tensor, g: Tensor, a: Tensor, b: Tensor,
               e: float, need_acc: bool) -> Tuple[Tensor, Tensor]:
    E = (f.view(-1, 1) / e + g.view(1, -1) / e - K_xy).exp_()
    r = E.t().mv(a)
    acc = torch.empty(0, dtype=K_xy.dtype, device=K_xy.device)
    if need_acc:
        acc = E.mv(b / r)        
    return r, acc

def ratio_sweep(e: float, 
      X: Tensor, Y: Tensor, 
      f: Tensor,  g: Tensor, 
      a: Tensor, b: Tensor, p: float,
      rows_y: int,
      online: bool, need_acc: bool) -> Tuple[Tensor, Tensor, bool]:
    # Dense: X, Y are K_xy, K_yx (cost / e). Online: X, Y are the points x, y.
    if online:
        r, acc = online_ratio_pass(X, Y, f, g, a, b, e, p, rows_y, need_acc)
    else:
        r, acc = dense_ratio_pass(X, Y, f, g, a, b, e, need_acc)
    ok = online_ratio_ok(r) and online_ratio_ok(acc)
    g = g - e * r.log()
    f = f - e * acc.log()
    norm = b.dot(g)
    return f + norm, g - norm, ok
  
def ratio_floor(f: Tensor, g: Tensor, e: float) -> float:
    # The ratio update can't resolve viol below its rounding noise, about u * |C| / e.
    u = 1.1e-16
    if f.element_size() == 4:
        u = 6e-8
    return 64.0 * u * (float(f.abs().max()) + float(g.abs().max())) / e

      
def sinkhorn_sweep(e: float, 
      X: Tensor, Y: Tensor, 
      f: Tensor,  g: Tensor, 
      a_log: Tensor, b_log: Tensor, a: Tensor, b: Tensor,  p: float,
      rows_x: int, rows_y: int,
      online: bool, use_ratio: bool) -> Tuple[Tensor, Tensor, bool]:
    if use_ratio:
        Tf, g_new, ok = ratio_sweep(e, X, Y, f, g, a, b, p, rows_y, online, True)
        if ok:
            return Tf, g_new, True
    Tf, g_new = softmin_sweep(e, X, Y, f, a_log, b_log, b, p, rows_x, rows_y, online)
      
    return Tf, g_new, False

def sinkhorn_xy(X: Tensor, Y: Tensor, 
                a_log: Tensor, b_log: Tensor, a: Tensor, b: Tensor,
                f_xy: Tensor, g_yx: Tensor, 
                eps_path: Tensor, eps: float, tol: float, m: int,
                p: float, rows_x: int, rows_y: int,
                online: bool, use_ratio: bool):
    e_prev = -1.0
    K_xy = X
    K_yx = Y
    iterations = 0
    converged = False
    # best = float('inf')
    # stall = 0
    F_hist: List[Tensor] = []
    R_hist: List[Tensor] = []
    viol_prev = -1.0
    Tf = f_xy
    for k in range(eps_path.numel()):
        e = float(eps_path[k])
        if not online and e != e_prev:
            K_xy = X / e
            K_yx = Y / e
            e_prev = e
        Tf, g_yx, use_ratio = sinkhorn_sweep(e, K_xy, K_yx, f_xy, g_yx, 
                                  a_log, b_log, a, b,
                                  p, rows_x, rows_y, online, use_ratio)
        iterations += 1
        r = Tf - f_xy
        viol = float(r.abs().max()) / e

        # stopping rules
        if e == eps and viol < tol:
            converged = True
            break
        if use_ratio and e == eps and viol < ratio_floor(f_xy, g_yx, e):
            use_ratio = False
            # stall = 0
        # if e == eps:
        #     if viol < 0.99 * best:
        #         best = viol
        #         stall = 0
        #     else:
        #         stall += 1
        #         if stall >= 10:
        #           break 
        if e != eps or m == 0:
            f_xy = Tf
            continue
        if viol_prev >= 0.0 and viol > 2.0 * viol_prev:
            F_hist.clear()
            R_hist.clear()
        viol_prev = viol
        f_xy = anderson_step(Tf, r, F_hist, R_hist, m)
    # return Tf because if code has taken an anderson accel step, then potential is not strictly correct
    return {'f': Tf, 'g': g_yx, 'iterations': torch.tensor(iterations), 'converged': torch.tensor(converged)}
      
def self_sweep(e: float, X: Tensor, f: Tensor, a_log: Tensor, a: Tensor,
               p: float, rows: int, online: bool, use_ratio: bool) -> Tuple[Tensor, bool]:
    # Dense: X is K_xx = C_xx / e. Online: X is the points x.
    if use_ratio:
        if online:
            r, _ = online_ratio_pass(X, X, f, f, a, a, e, p, rows, False)
        else:
            r, _ = dense_ratio_pass(X, X, f, f, a, a, e, False)
        if online_ratio_ok(r):
            return f - 0.5 * e * r.log(), True
    if online:
        S = online_softmin(X, X, f / e + a_log, 1.0 / e, p, rows)
    else:
        S = softmin_fixed_K(e, X, f, a_log)
    return 0.5 * (f + S), False


def sinkhorn_xx(X: Tensor, a_log: Tensor, a: Tensor, f_xx: Tensor,
                eps_path: Tensor, eps: float, tol: float, m: int,
                p: float, rows: int, online: bool, use_ratio: bool):
    # Symmetric (self) potential with the damped update f -> (f + T(f)) / 2.
    # Same structure as sinkhorn_xy; the returned f is the last damped sweep.
    e_prev = -1.0
    K_xx = X
    iterations = 0
    converged = False
    # best = float('inf')
    # stall = 0
    F_hist: List[Tensor] = []
    R_hist: List[Tensor] = []
    viol_prev = -1.0
    Tf = f_xx
    for k in range(eps_path.numel()):
        e = float(eps_path[k])
        if not online and e != e_prev:
            K_xx = X / e
            e_prev = e
        Tf, use_ratio = self_sweep(e, K_xx, f_xx, a_log, a, p, rows, online, use_ratio)
        iterations += 1
        r = Tf - f_xx
        viol = float(r.abs().max()) / e
        if e == eps and viol < tol:
            converged = True
            break
        if use_ratio and e == eps and viol < ratio_floor(f_xx, f_xx, e):
            use_ratio = False
            # stall = 0
        # if e == eps:
        #     if viol < 0.99 * best:
        #         best = viol
        #         stall = 0
        #     else:
        #         stall += 1
        #         if stall >= 10:
        #           break 
        if e != eps or m == 0:
            f_xx = Tf
            continue
        if viol_prev >= 0.0 and viol > 2.0 * viol_prev:
            F_hist.clear()
            R_hist.clear()
        viol_prev = viol
        f_xx = anderson_step(Tf, r, F_hist, R_hist, m)
    # return Tf because if code has taken an anderson accel step, then potential f_xx is not strictly correct
    return {'f': Tf, 'iterations': torch.tensor(iterations), 'converged': torch.tensor(converged)}
"
      )
    }
    compiled
  }
})
