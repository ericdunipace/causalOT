ccp_forward_code_tensorized <- "

 def ccp_log_weights(w: Tensor):
   return torch.clamp(torch.log(w), min=-100000.0)

 def ccp_candidate_weights(log_w: Tensor, direction: Tensor, tau: float):
   return (log_w + direction / tau).log_softmax(0).exp()

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
   return {'kkt_constant' : kappa, 'kkt_max_deviation' : max_deviation}

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
      private$iterate_over_lambda(niter, tol)

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
    ccp_diagnostics_every = "integer",
    ccp_iter = "integer",
    ccp_state_ready = "logical",
    ccp_candidate_count = "integer",
    ccp_w = "torch_tensor",
    ccp_log_w = "torch_tensor",
    ccp_direction = "torch_tensor",
    ccp_loss = "numeric",
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
      private$ccp_tau_min <- lambda / 10.0
      private$ccp_tau_max <- lambda * 10.0
      private$ccp_tau <- private$ccp_tau_min
      private$ccp_tau_shrink <- 1.5
      private$ccp_tau_grow <- 3.0
      private$ccp_tau_multipliers <- c(1.0, 3.0, 10.0)
      if (!is.logical(private$ccp_backtracking)) {
        private$ccp_backtracking <- FALSE
      }
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

    ccp_should_run_diagnostics = function(iter, niter, converged = FALSE) {
      isTRUE(converged) ||
        iter >= niter ||
        iter %% private$ccp_diagnostics_every == 0L
    },

    ccp_compute_kkt = function(ot) {
      private$ccp_kkt_resid <- private$ccp_ts$ccp_kkt(
        ot$a$detach(),
        ot$potentials$f_xy$detach(),
        ot$potentials$f_xx$detach(),
        torch::jit_scalar(1e-10)
      )$kkt_max_deviation$item()
    },

    iterate_over_lambda = function(niter, tol) {
      for (l in private$penalty_list$lambda) {
        private$set_lambda(l)

        res <- private$iterate_over_delta(niter, tol)

        private$weights_list[[as.character(l)]] <- res$weights
        private$iterations_run[[as.character(l)]] <- res$iter
        private$final_loss[[as.character(l)]] <- res$loss
        private$converged_run[[as.character(l)]] <- res$converged
      }
    },

    ccp_step = function(opt, osqp_args = NULL, tol) {
      ot <- private$ot_objects[[private$ccp_ot_name]]
      torch::with_no_grad({
        private$ccp_set_current_state(ot)
        if (isTRUE(private$ccp_backtracking)) {
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
      torch::with_no_grad({
        private$ccp_set_current_state(ot)
      })
      loss_old <- private$ccp_loss
      check <- TRUE
      converged_flag <- FALSE

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
        loss_converged <- (i > 1) && converged(loss, loss_old, tol)
        if (private$ccp_should_run_diagnostics(i, niter, loss_converged)) {
          private$ccp_compute_kkt(private$ot_objects[[private$ccp_ot_name]])
        }
        kkt_ok <- is.finite(private$ccp_kkt_resid) &&
          private$ccp_kkt_resid <= kkt_tol
        if (check && loss_converged && kkt_ok) {
          converged_flag <- TRUE
          break
        }

        loss_old <- loss
      }

      private$torch_optim_reset()
      private$lbfgs_count <- private$lbfgs_reset <- 0L

      return(list(loss = loss, iter = i, converged = converged_flag))
    }
  )
)
