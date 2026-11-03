ccp_forward_code_tensorized <- "

 def ccp_log_weights(w: Tensor):
   return torch.clamp(torch.log(w), min=-100000.0)

 def ccp_candidate_weights(log_w: Tensor, direction: Tensor, tau: float):
   return (log_w + direction / tau).log_softmax(0).exp()

 def ccp_model_terms(w_new: Tensor, w: Tensor, log_w: Tensor,
                     direction: Tensor):
   # linear term <grad, w_new - w> and Bregman term KL(w_new || w) of the
   # relative smoothness upper model; direction is the negative gradient
   log_new = torch.clamp(torch.log(w_new), min=-100000.0)
   lin = -direction.dot(w_new - w)
   kl = (w_new * (log_new - log_w)).sum()
   return {'lin' : lin, 'kl' : kl}

 def ccp_kkt(w: Tensor, f_xy: Tensor, f_xx: Tensor, support_tol: float):
   kkt = f_xy - f_xx
   active = w > support_tol
   active_w = w.masked_select(active)
   active_kkt = kkt.masked_select(active)
   if active_w.numel() > 0:
     kappa = active_w.dot(active_kkt) / active_w.sum()
     max_deviation = (active_kkt - kappa).abs().max()
   else:
     kappa = torch.tensor(float('nan'), dtype=w.dtype, device=w.device)
     max_deviation = torch.tensor(float('nan'), dtype=w.dtype, device=w.device)
   # Frank-Wolfe gap: covers both KKT conditions (equal gradient on the
   # support, no smaller gradient off it) and bounds f(w) - f(w*) by convexity
   fw_gap = w.dot(kkt) - kkt.min()
   return {'kkt_constant' : kappa, 'kkt_max_deviation' : max_deviation,
           'fw_gap' : fw_gap}

"

cotCCPTrain <- R6::R6Class(
  classname = "cotCCPTrain",
  inherit = OTProblem_,
  public = list(
    solve = function(
      niter = 1000L,
      tol = 1e-5,
      osqp_args = NULL,
      quick.balance.function = TRUE,
      ...
    ) {
      stopifnot(
        "arguments not set! Run '$setup_arguments' first." = private$args_set
      )
      stopifnot("Argument `niter` must be > 0" = (niter > 0))
      stopifnot("Argument `tol` must be >=0" = (tol >= 0))

      private$osqp_args <- osqp_args[
        names(osqp_args) %in% methods::formalArgs(osqp::osqpSettings)
      ]
      private$delta_values_setup(
        run.quick = quick.balance.function,
        osqp_args = private$osqp_args
      )
      private$optimizer_name <- "mirror"
      private$optimization_step <- private$ccp_step
      private$ccp_setup()
      private$fit_status_run <- list()
      private$iterate_over_lambda(niter, tol)
      private$check_fit_path()

      return(invisible(self))
    }
  ),
  private = list(
    ccp_ot_name = "character",
    ccp_adapt_address = "character",
    ccp_kkt_resid = "numeric",
    ccp_obj_delta = "numeric",
    ccp_tau = "numeric",
    ccp_tau_min = "numeric",
    ccp_tau_max = "numeric",
    ccp_tau_shrink = "numeric",
    ccp_tau_grow = "numeric",
    ccp_tau_multipliers = "numeric",
    ccp_backtracking = "logical",
    ccp_step_rule = "character",
    ccp_prev_log_x = NULL,
    ccp_prev_x = NULL,
    ccp_momentum_t = "numeric",
    ccp_tau_factor = "numeric",
    ccp_diagnostics_every = "integer",
    ccp_iter = "integer",
    ccp_state_ready = "logical",
    ccp_candidate_count = "integer",
    ccp_w = "torch_tensor",
    ccp_log_w = "torch_tensor",
    ccp_direction = "torch_tensor",
    ccp_loss = "numeric",
    ccp_loss_offset = "numeric",
    ccp_conv_scale = "numeric",
    ccp_potentials = "list",
    ccp_candidate_w = "torch_tensor",
    ccp_candidate_loss = "numeric",
    ccp_candidate_potentials = "list",
    ccp_candidate_direction = "torch_tensor",
    ccp_candidate_tau = "numeric",
    ccp_ts = "torch_jit",

    ccp_setup = function() {
      if (length(private$ot_objects) != 1L) {
        stop("solver = 'ccp' currently supports exactly one OT problem.")
      }
      if (length(private$target_objects) != 0L) {
        stop("solver = 'ccp' does not yet support balance constraints.")
      }
      if (
        length(private$penalty_list$lambda) != 1L ||
          !is.finite(private$penalty_list$lambda[[1L]]) ||
          private$penalty_list$lambda[[1L]] <= 0
      ) {
        stop("solver = 'ccp' requires a single finite positive lambda.")
      }
      ot_name <- ls(private$ot_objects)[[1L]]
      ot <- private$ot_objects[[ot_name]]
      if (!ot$debias) {
        stop("solver = 'ccp' requires debias = TRUE.")
      }
      if (!ot$tensorized) {
        stop("solver = 'ccp' currently requires tensorized costs.")
      }

      problem <- private$problems[[ot_name]]
      adapt <- vapply(
        problem,
        function(add) private$measures[[add]]$adapt,
        character(1L)
      )
      if (sum(adapt == "weights") != 1L) {
        stop(
          "solver = 'ccp' currently supports exactly one adapted weight measure."
        )
      }
      private$ccp_ot_name <- ot_name
      private$ccp_adapt_address <- problem[[which(adapt == "weights")]]
      lambda <- private$ot_objects[[ot_name]]$penalty
      if (!is.numeric(private$ccp_tau_factor)) {
        private$ccp_tau_factor <- 0.1
      }
      private$ccp_tau_min <- lambda * private$ccp_tau_factor
      private$ccp_tau_max <- lambda * 10.0
      private$ccp_tau <- private$ccp_tau_min
      private$ccp_tau_shrink <- 1.5
      private$ccp_tau_grow <- 3.0
      private$ccp_tau_multipliers <- c(1.0, 3.0, 10.0)
      if (!is.logical(private$ccp_backtracking)) {
        private$ccp_backtracking <- FALSE
      }
      if (
        !isTRUE(
          private$ccp_step_rule %in% c("fixed", "adaptive", "accelerated")
        )
      ) {
        private$ccp_step_rule <- "fixed"
      }
      private$ccp_prev_log_x <- NULL
      private$ccp_prev_x <- NULL
      private$ccp_momentum_t <- 1.0
      if (!is.numeric(private$ccp_diagnostics_every)) {
        private$ccp_diagnostics_every <- 10L
      }
      private$ccp_diagnostics_every <- as.integer(
        max(1L, private$ccp_diagnostics_every)
      )
      private$ccp_iter <- 0L
      private$ccp_state_ready <- FALSE
      private$ccp_candidate_count <- 0L
      private$ccp_kkt_resid <- Inf
      private$ccp_obj_delta <- Inf
      private$ccp_ts <- torch::jit_compile(ccp_forward_code_tensorized)
    },

    ccp_clone_potentials = function(ot) {
      lapply(ot$potentials, function(x) x$detach()$clone())
    },

    ccp_eval_sinkhorn = function(ot) {
      ot$sinkhorn_ccp_eval(niter = private$ot_niter, tol = private$ot_tol)
    },

    ccp_capture_state = function(ot, eval) {
      private$ccp_w <- ot$a$detach()
      private$ccp_log_w <- private$ccp_ts$ccp_log_weights(private$ccp_w)
      private$ccp_loss <- eval$loss$detach()$to(device = "cpu")$item()
      private$ccp_direction <- eval$direction$detach()
      private$ccp_state_ready <- TRUE
    },

    ccp_set_current_state = function(ot) {
      if (
        !isTRUE(private$ccp_state_ready) ||
          is.null(ot$potentials) ||
          length(ot$potentials) == 0
      ) {
        eval <- private$ccp_eval_sinkhorn(ot)
        private$ccp_capture_state(ot, eval)
      }
      private$ccp_potentials <- private$ccp_clone_potentials(ot)
    },

    ccp_tau_candidates = function() {
      tau_grid <- private$ccp_tau * private$ccp_tau_multipliers
      tau_grid <- pmin(
        pmax(tau_grid, private$ccp_tau_min),
        private$ccp_tau_max
      )
      unique(tau_grid[is.finite(tau_grid) & tau_grid > 0])
    },

    ccp_candidate_weights = function(tau) {
      private$ccp_candidate_count <- private$ccp_candidate_count + 1L
      private$ccp_ts$ccp_candidate_weights(
        private$ccp_log_w,
        private$ccp_direction,
        torch::jit_scalar(tau)
      )
    },

    ccp_try_candidate = function(ot, tau) {
      w <- private$ccp_candidate_weights(tau)
      ot$a <- w
      ot$potentials <- lapply(private$ccp_potentials, function(x) {
        x$detach()$clone()
      })
      eval <- private$ccp_eval_sinkhorn(ot)
      private$ccp_candidate_tau <- tau
      private$ccp_candidate_w <- w
      private$ccp_candidate_loss <- eval$loss$detach()$to(device = "cpu")$item()
      private$ccp_candidate_direction <- eval$direction$detach()
      private$ccp_candidate_potentials <- private$ccp_clone_potentials(ot)
    },

    ccp_update_tau = function(accepted, tau = private$ccp_tau) {
      if (accepted) {
        private$ccp_tau <- max(
          tau / private$ccp_tau_shrink,
          private$ccp_tau_min
        )
      } else {
        private$ccp_tau <- min(
          private$ccp_tau * private$ccp_tau_grow,
          private$ccp_tau_max
        )
      }
    },

    ccp_backtracking_search = function(ot) {
      accepted <- FALSE
      for (tau in private$ccp_tau_candidates()) {
        private$ccp_try_candidate(ot, tau)
        if (
          is.finite(private$ccp_candidate_loss) &&
            private$ccp_candidate_loss < private$ccp_loss
        ) {
          private$ccp_update_tau(TRUE, tau)
          accepted <- TRUE
          break
        }
      }
      if (!accepted) {
        private$ccp_update_tau(FALSE)
        private$ccp_candidate_tau <- private$ccp_tau
        private$ccp_candidate_w <- private$ccp_w
        private$ccp_candidate_loss <- private$ccp_loss
        private$ccp_candidate_potentials <- private$ccp_potentials
        private$ccp_candidate_direction <- private$ccp_direction
      }
    },

    ccp_accept_candidate = function(ot) {
      weights <- rlang::env()
      class(weights) <- c("weightEnv", class(weights))
      weights[[private$ccp_adapt_address]] <- private$ccp_candidate_w
      private$parameters_get_set(weights)
      private$weights <- weights
      ot$a <- private$ccp_candidate_w
      ot$potentials <- lapply(private$ccp_candidate_potentials, function(x) {
        x$detach()$clone()
      })
      private$ccp_w <- private$ccp_candidate_w$detach()
      private$ccp_log_w <- private$ccp_ts$ccp_log_weights(private$ccp_w)
      private$ccp_loss <- private$ccp_candidate_loss
      private$ccp_direction <- private$ccp_candidate_direction
      private$ccp_state_ready <- TRUE
    },

    ccp_step_single_tau = function(ot) {
      w <- private$ccp_candidate_weights(private$ccp_tau)
      ot$a <- w
      eval <- private$ccp_eval_sinkhorn(ot)
      private$ccp_candidate_tau <- private$ccp_tau
      private$ccp_candidate_w <- w
      private$ccp_candidate_loss <- eval$loss$detach()$to(device = "cpu")$item()
      private$ccp_candidate_potentials <- NULL
      private$ccp_candidate_direction <- eval$direction$detach()
      private$ccp_update_tau(TRUE, private$ccp_tau)
      weights <- rlang::env()
      class(weights) <- c("weightEnv", class(weights))
      weights[[private$ccp_adapt_address]] <- w
      private$parameters_get_set(weights)
      private$weights <- weights
      private$ccp_w <- w$detach()
      private$ccp_log_w <- private$ccp_ts$ccp_log_weights(private$ccp_w)
      private$ccp_loss <- private$ccp_candidate_loss
      private$ccp_direction <- private$ccp_candidate_direction
      private$ccp_state_ready <- TRUE
    },

    # Bregman proximal gradient with backtracking on relative smoothness:
    # accept tau if F(w+) <= F(w) + <g, w+ - w> + tau * KL(w+ || w), else
    # double tau. After an accepted step, halve tau so the step can grow
    # again, so the step size tracks the local curvature rather than lambda.
    ccp_step_adaptive = function(ot, max_backtracks = 30L) {
      if (private$ccp_iter <= 1L) {
        # first step moves log weights by at most one unit
        g_range <- (private$ccp_direction$max() -
          private$ccp_direction$min())$item()
        private$ccp_tau <- max(g_range, .Machine$double.eps)
      }
      eps <- if (identical(ot$a$dtype, torch::torch_float32())) {
        1.2e-7
      } else {
        .Machine$double.eps
      }
      tau <- private$ccp_tau
      accepted <- FALSE
      for (k in seq_len(max_backtracks)) {
        private$ccp_try_candidate(ot, tau)
        terms <- private$ccp_ts$ccp_model_terms(
          private$ccp_candidate_w,
          private$ccp_w,
          private$ccp_log_w,
          private$ccp_direction
        )
        bound <- private$ccp_loss +
          terms$lin$item() +
          tau * terms$kl$item()
        slack <- 10 * eps * max(abs(private$ccp_loss), 1)
        if (
          is.finite(private$ccp_candidate_loss) &&
            private$ccp_candidate_loss <= bound + slack
        ) {
          accepted <- TRUE
          break
        }
        tau <- tau * 2.0
      }
      if (accepted) {
        private$ccp_tau <- tau / 2.0
        private$ccp_accept_candidate(ot)
      } else {
        private$ccp_tau <- tau
        private$ccp_candidate_loss <- private$ccp_loss
      }
    },

    # FISTA-type momentum on the log weights with gradient restart
    # (O'Donoghue & Candes, 2015). The current state (ccp_w, ccp_direction,
    # ccp_loss) is the extrapolated point y_k where the gradient is evaluated;
    # x_{k+1} is the mirror step from y_k and y_{k+1} extrapolates x_{k+1}.
    # Uses only gradients, so it is not affected by noise in the loss values.
    ccp_step_accelerated = function(ot) {
      tau <- private$ccp_tau
      log_x <- (private$ccp_log_w + private$ccp_direction / tau)$log_softmax(1)
      x <- log_x$exp()
      t_old <- private$ccp_momentum_t
      if (!is.null(private$ccp_prev_x)) {
        # restart when the step from y_k goes uphill relative to x_k
        uphill <- (-private$ccp_direction$dot(x - private$ccp_prev_x))$item()
        if (uphill > 0) {
          t_old <- 1.0
        }
      }
      t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
      beta <- if (is.null(private$ccp_prev_log_x)) 0 else (t_old - 1) / t_new
      if (beta > 0) {
        y <- (log_x + beta * (log_x - private$ccp_prev_log_x))$log_softmax(1)$exp()
      } else {
        y <- x
      }
      private$ccp_prev_log_x <- log_x
      private$ccp_prev_x <- x
      private$ccp_momentum_t <- t_new

      private$ccp_candidate_count <- private$ccp_candidate_count + 1L
      ot$a <- y
      eval <- private$ccp_eval_sinkhorn(ot)
      private$ccp_candidate_tau <- tau
      private$ccp_candidate_w <- y
      private$ccp_candidate_loss <- eval$loss$detach()$to(device = "cpu")$item()
      weights <- rlang::env()
      class(weights) <- c("weightEnv", class(weights))
      weights[[private$ccp_adapt_address]] <- y
      private$parameters_get_set(weights)
      private$weights <- weights
      private$ccp_w <- y$detach()
      private$ccp_log_w <- private$ccp_ts$ccp_log_weights(private$ccp_w)
      private$ccp_loss <- private$ccp_candidate_loss
      private$ccp_direction <- eval$direction$detach()
      private$ccp_state_ready <- TRUE
    },

    ccp_should_run_diagnostics = function(iter, niter, converged = FALSE) {
      isTRUE(converged) ||
        iter >= niter ||
        iter %% private$ccp_diagnostics_every == 0L
    },

    ccp_compute_kkt = function(ot) {
      # relative Frank-Wolfe gap, (<g, w> - min g) / max(|loss|, scale), with
      # the loss corrected for the omitted y self term so it is the actual
      # divergence; the floor keeps the test attainable as the divergence -> 0
      fw_gap <- private$ccp_ts$ccp_kkt(
        ot$a$detach(),
        ot$potentials$f_xy$detach(),
        ot$potentials$f_xx$detach(),
        torch::jit_scalar(1e-10)
      )$fw_gap$item()
      private$phase1_cert <- list(fw_gap = fw_gap,
        divergence = private$ccp_loss - private$ccp_loss_offset)
      private$ccp_kkt_resid <- fw_gap /
        max(
          abs(private$ccp_loss - private$ccp_loss_offset),
          private$ccp_conv_scale
        )
    },

    iterate_over_lambda = function(niter, tol) {
      for (l in private$penalty_list$lambda) {
        private$set_lambda(l)

        res <- private$iterate_over_delta(niter, tol)

        private$weights_list[[as.character(l)]] <- res$weights
        private$iterations_run[[as.character(l)]] <- res$iter
        private$final_loss[[as.character(l)]] <- res$loss
        private$converged_run[[as.character(l)]] <- res$converged
        private$fit_status_run[[as.character(l)]] <- res$fit_status
      }
    },

    ccp_step = function(opt, osqp_args = NULL, tol) {
      ot <- private$ot_objects[[private$ccp_ot_name]]
      torch::with_no_grad({
        private$ccp_set_current_state(ot)
        if (identical(private$ccp_step_rule, "adaptive")) {
          private$ccp_step_adaptive(ot)
        } else if (identical(private$ccp_step_rule, "accelerated")) {
          private$ccp_step_accelerated(ot)
        } else if (isTRUE(private$ccp_backtracking)) {
          private$ccp_backtracking_search(ot)
          private$ccp_accept_candidate(ot)
        } else {
          private$ccp_step_single_tau(ot)
        }
        loss <- torch::torch_tensor(
          private$ccp_candidate_loss,
          dtype = ot$a$dtype,
          device = "cpu"
        )
      })
      loss
    },

    optimization_loop = function(niter, tol) {
      ot <- private$ot_objects[[private$ccp_ot_name]]
      # computed here, after set_lambda(), since it depends on the penalty
      private$ccp_loss_offset <- cot_ccp_y_self_term(ot)
      private$ccp_conv_scale <- cot_conv_scale(ot)
      torch::with_no_grad({
        private$ccp_set_current_state(ot)
      })
      loss_old <- private$ccp_loss
      check <- TRUE
      converged_flag <- FALSE
      private$phase1_cert <- list(fw_gap = NA_real_, divergence = NA_real_)

      for (i in 1:niter) {
        private$ccp_iter <- i
        loss <- private$optimization_step(
          private$opt,
          private$osqp_args,
          tol = tol
        )$detach()$to(device = "cpu")$item()
        private$ccp_obj_delta <- abs(loss - loss_old)
        check <- private$lr_reduce(loss)

        kkt_tol <- max(tol, 1e-6)
        # relative change of the divergence, not of the offset CCP loss
        loss_converged <- (i > 1) &&
          converged(
            loss - private$ccp_loss_offset,
            loss_old - private$ccp_loss_offset,
            tol,
            scale = private$ccp_conv_scale
          )
        if (private$ccp_should_run_diagnostics(i, niter, loss_converged)) {
          private$ccp_compute_kkt(private$ot_objects[[private$ccp_ot_name]])
        }
        # the debiased divergence is nonnegative, so it bounds f(w) - f(w*) as
        # well; it is the tighter bound when the optimum is near zero and the
        # Frank-Wolfe gap decays slowly
        kkt_ok <- (is.finite(private$ccp_kkt_resid) &&
          private$ccp_kkt_resid <= kkt_tol) ||
          (loss - private$ccp_loss_offset) <= kkt_tol * private$ccp_conv_scale
        if (check && loss_converged && kkt_ok) {
          converged_flag <- TRUE
          break
        }

        loss_old <- loss
      }

      private$torch_optim_reset()
      private$lbfgs_count <- private$lbfgs_reset <- 0L

      return(list(loss = loss, iter = i, converged = converged_flag,
        stopped_by_rule = converged_flag, fw_gap = private$phase1_cert$fw_gap,
        divergence = private$phase1_cert$divergence))
    }
  )
)
