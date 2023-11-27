setOldClass("summary_causalWeights")

#' Summary diagnostics for causalWeights
#'
#' @param object an object of class [causalWeights][causalOT::causalWeights-class]
#' @param r_eff The r_eff used in the PSIS calculation. See [PSIS_diag()][PSIS_diag()]
#' @param penalty The penalty parameter to use 
#' @param p The power of the Lp distance to use. Overridden by argument `cost.`
#' @param cost A user supplied cost function. Should take arguments `x1`, `x2`, `p`.
#' @param debias Should debiased optimal transport distances be used. TRUE or FALSE
#' @param online.cost Should the cost be calculated online? One of "auto","tensorized", or "online".
#' @param diameter the diameter of the covariate space. Default is NULL.
#' @param niter the number of iterations to run the optimal transport distances
#' @param tol the tolerance for convergence for the optimal transport distances
#' @param ... passed to [PSIS_diag()][PSIS_diag()]
#'
#' @return The summary method returns an object of class "summary_causalWeights".
#' @export
#'
#' @examples
#' if(torch::torch_is_installed()) {
#' n <- 2^6
#' p <- 6
#' overlap <- "high"
#' design <- "A"
#' estimand <- "ATE"
#' 
#' #### get simulation functions ####
#' original <- Hainmueller$new(n = n, p = p, 
#'                             design = design, overlap = overlap)
#' original$gen_data()
#' weights <- calc_weight(x = original, estimand = estimand, method = "Logistic")
#' s <- summary(weights)
#' plot(s)
#' }
summary.causalWeights <- function(object, r_eff = NULL, 
                                  penalty, p = 2, cost = NULL, 
                                  debias = TRUE, online.cost = "auto", 
                                  diameter = NULL, niter = 1000, 
                                  tol = 1e-07, ...) {
  
  stopifnot(inherits(object, "causalWeights"))
  if ( missing(penalty) || is.null(penalty) ) { 
    penalty <- NA_real_
  }
  
  # run diagnostic functions
  psis <- PSIS_diag(object, r_eff, ...)
  mb   <- mean_balance_diagnostic(object)
  ot   <- ot_distance(object, penalty = penalty, p = p,
                      cost = cost, debias = debias, online.cost = online.cost,
                      diameter = diameter, niter = niter, tol = tol)
  
  # pull out pareto info
  pareto_k                      <- sapply(psis, function(p) p[["pareto_k"]])
  n_eff                         <- sapply(psis, function(p) p[["n_eff"]])
  names(pareto_k)<-names(n_eff) <- c("controls", "treated")
  attributes(n_eff)             <- c(attributes(n_eff),     list(n_original = c(controls = length(object@w0), treated = length(object@w1))))
  
  # save results
  res   <- list(ot           = ot,
                pareto_k     = pareto_k,
                n_eff        = n_eff,
                mean.balance = mb,
                estimand     = object@estimand)
  
  class(res) <- "summary_causalWeights"
  return(res)
}

#' print.summary_causalWeights
#' 
#' @param x an object of class "summary_causalWeights" returned by the function [summary.causalWeights()][summary.causalWeights()]
#' @param ... Not used at this time.
#' @export
#' @method print summary_causalWeights
#' @describeIn summary.causalWeights print method
print.summary_causalWeights <- function(x,...) {
  
  object   <- x
  estimand <- object$estimand
  
  control.output <- treated.output <- NULL
  
  # res   <- list(ot           = ot,
  #               pareto_k     = pareto_k,
  #               n_eff        = n_eff,
  #               mean.balance = mb,
  #               estimand     = object@estimand)
  row.labels <- c("OT distance", "Pareto k", "N eff", "Avg. std. mean balance")
  
  if (estimand == "ATT" ) {
    control.output <- data.frame(pre  = c(object$ot$pre, NA_real_, 
                                          attr(object$n_eff,"n_original")[1], 
                                          mean(object$mean.balance$pre)),
                                 post = c(object$ot$post, object$pareto_k[1], 
                                          object$n_eff[1], 
                                          mean(object$mean.balance$post)),
                                 row.names = row.labels)
  } else if (estimand == "ATC") {
    treated.output <- data.frame(pre  = c(object$ot$pre, NA_real_, 
                                          attr(object$n_eff,"n_original")[2], 
                                          mean(object$mean.balance$pre)),
                                 post = c(object$ot$post, object$pareto_k[2], 
                                           object$n_eff[2], 
                                           mean(object$mean.balance$post)),
                                 row.names = row.labels)
  } else if (estimand == "ATE") {
    control.output <- data.frame(pre  = c(object$ot$pre[1], NA_real_, 
                                          attr(object$n_eff,"n_original")[1], 
                                          mean(object$mean.balance$pre$controls)),
                                post = c(object$ot$post[1], object$pareto_k[1], 
                                         object$n_eff[1], 
                                         mean(object$mean.balance$post$controls)),
                                row.names = row.labels)
    treated.output <- data.frame(pre  = c(object$ot$pre[2], NA_real_, 
                                          attr(object$n_eff,"n_original")[2], 
                                          mean(object$mean.balance$pre$treated)),
                                  post = c(object$ot$post[2], object$pareto_k[2], 
                                           object$n_eff[2], 
                                           mean(object$mean.balance$post$treated)),
                                  row.names = row.labels)
  }
  
  cat( paste0("Diagnostics for causalWeights for estimand ", estimand, "\n") )
  if(!is.null(control.output)) {
    cat( paste0("Control group\n") )
    print(control.output)
    cat("\n")
  }  
  if (!is.null(treated.output)) {
    cat( paste0("Treated group\n") )
    print(treated.output)
  }
}

# no roxygen needed
#' @keywords internal
setMethod("show", "summary_causalWeights", function(object){print(object)})

#' plot.summary_causalWeights
#' 
#' @param x an object of class "summary_causalWeights"
#' @param ... Not used
#'
#' @export
#' @importFrom rlang .data
#' @method plot summary_causalWeights
#' @describeIn summary.causalWeights plot method
plot.summary_causalWeights <- function(x, ...) {
  object   <- x
  estimand <- object$estimand
  
  plot_ot <- function(object, ...) {
    ot_dat <- if (object$estimand == "ATT") {
      data.frame(value = c(object$ot$pre, object$ot$post),
                 group = factor(rep("control",2), levels = c("control", "treated")),
                 period = factor(c("pre", "post"), levels = c("pre","post")))
    } else if (object$estimand == "ATC") {
      data.frame(value = c(object$ot$pre, object$ot$post),
                 group = factor(rep("treated",2), levels = c("control", "treated")),
                 period = factor(c("pre", "post"), levels = c("pre","post")))
    } else if (object$estimand == "ATE") {
      data.frame(value = unlist(object$ot),
                 group = factor(rep(c("control", "treated"), 2), levels = c("control", "treated")),
                 period = factor(rep(c("pre", "post"), each = 2), levels = c("pre","post")))
    }
    
    
    
    p <- ggplot2::ggplot(ot_dat, 
                         ggplot2::aes(x = .data$period, 
                                      y = .data$value, 
                                      group = .data$group, 
                                      color = .data$group)
                         ) +
      ggplot2::geom_line() + ggplot2::geom_point() + 
      ggplot2::scale_color_manual(values = c("black", "gray")) +
      ggplot2::ylab("Sinkhorn Distance") + ggplot2::xlab("") + 
      ggplot2::ggtitle("Sinkhorn Distances") +
      ggplot2::theme_bw()
    
    return(p)
  }
  
  plot_pareto_k <- function(object,  ...) {
    # list(pareto_k     = sapply(psis, function(p) p["pareto_k"]),
    #      n_eff        = sapply(psis, function(p) p["n_eff"]),
    #      mean.balance = mb,
    #      ot           = ot,
    #      estimand     = object@estimand)
    
    pareto_k_dat <- if(object$estimand == "ATT") {
      data.frame(value = c(object$pareto_k["controls"]),
                 group = factor("control", levels = c("control", "treated")))
    } else if (object$estimand == "ATC") {
      data.frame(value = c(object$pareto_k["treated"]),
                 group = factor("treated", levels = c("control", "treated")))
    } else if (object$estimand == "ATE") {
      data.frame(value = c(object$pareto_k),
                 group = factor(c("control", "treated"), levels = c("control", "treated")))
    }
    
    p <- ggplot2::ggplot(pareto_k_dat, 
                         ggplot2::aes(x = .data$group, 
                                      y = .data$value)) +
      ggplot2::geom_point() + 
      ggplot2::geom_hline(yintercept = 0.5, linetype = 2) +
      ggplot2::geom_hline(yintercept = 1, linetype = 1, color = "red") + 
      ggplot2::ylab("k") + ggplot2::xlab("") +
      ggplot2::ggtitle("Pareto k Statistics") +
      ggplot2::theme_bw()
    
    return(p)
  }
  
  plot_neff <- function(object, ...) {
    # list(pareto_k     = sapply(psis, function(p) p["pareto_k"]),
    #      n_eff        = sapply(psis, function(p) p["n_eff"]),
    #      mean.balance = mb,
    #      ot           = ot,
    #      estimand     = object@estimand)
    
    pareto_n_dat <- if(object$estimand == "ATT") {
      data.frame(value = c(attr(object$n_eff, "n_original")[1],
                           object$n_eff["controls"]),
                 group = factor("control", levels = c("control", "treated")),
                 period = factor(c("pre", "post"), levels = c("pre","post")))
    } else if (object$estimand == "ATC") {
      data.frame(value = c(attr(object$n_eff, "n_original")[2],
                           object$n_eff["treated"]),
                 group = factor("treated", levels = c("control", "treated")),
                 period = factor(c("pre", "post"), levels = c("pre","post")))
    } else if (object$estimand == "ATE") {
      data.frame(value = c(attr(object$n_eff, "n_original")[1],
                           object$n_eff[1],
                           attr(object$n_eff, "n_original")[2],
                           object$n_eff[2]
                           ),
                 group = factor(rep(c("control", "treated"), each = 2), levels = c("control", "treated")),
                 period = factor(rep(c("pre", "post"), 2), levels = c("pre","post")))
    }
    
    p <- ggplot2::ggplot(pareto_n_dat, 
                         ggplot2::aes(x = .data$group, 
                                      y = .data$value, 
                                      fill = .data$period,
                                      label = paste0("N_eff = ", round(.data$value, digits = 1)))) +
      ggplot2::geom_col( position = 'dodge' ) +
      ggplot2::scale_fill_manual(values = c("black","gray")) + 
      ggplot2::xlab("") + 
      ggplot2::scale_y_continuous("N effective", expand = c(0,0), 
                                  limits = c(0, max(pareto_n_dat$value) * 1.1)) +
      ggplot2::geom_text(position=ggplot2::position_dodge(width=0.9),vjust = -.5) +
      ggplot2::ggtitle("Effective sample size") + 
      ggplot2::theme_bw()
    
    return(p)
  }
  
  plot_mb <- function(object, ...) {
    # list(pareto_k     = sapply(psis, function(p) p["pareto_k"]),
    #      n_eff        = sapply(psis, function(p) p["n_eff"]),
    #      mean.balance = mb,
    #      ot           = ot,
    #      estimand     = object@estimand)
    
    d <- ifelse(object$estimand == "ATE", max(length(object$mean.balance$pre$control), 
                                              length(object$mean.balance$pre$treated)),
                length(object$mean.balance$pre))
    
    mb_dat <- if(object$estimand == "ATT") {
      data.frame(value = c(object$mean.balance$pre,
                           object$mean.balance$post),
                 covariate = rep(1:d, 2),
                 group = factor("control", levels = c("control", "treated")),
                 period = factor(rep(c("pre", "post"), each = d), levels = c("pre","post")))
    } else if (object$estimand == "ATC") {
      data.frame(value = c(object$mean.balance$pre,
                           object$mean.balance$post),
                 covariate = rep(1:d, 2),
                 group = factor("treated", levels = c("control", "treated")),
                 period = factor(rep(c("pre", "post"), each = d), levels = c("pre","post")))
    } else if (object$estimand == "ATE") {
      data.frame(value = c(unlist(object$mean.balance)),
                  covariate = rep(1:d, 4),
                  group = factor(rep(rep(c("control", "treated"), each = d), 2), 
                                 levels = c("control", "treated")),
                  period = factor(rep(c("pre", "post"), each = 2*d), levels = c("pre","post")))
    }
    
    p <- ggplot2::ggplot(mb_dat, ggplot2::aes(x = .data$period, y = .data$value, group = interaction(.data$covariate, .data$group), color = .data$group, size = .data$group)) +
      ggplot2::geom_hline(yintercept = 0.1, linetype = 2) +
      ggplot2::geom_hline(yintercept = 0.2, linetype = 1, color = "red") +
      ggplot2::geom_line(size = .5) + 
      ggplot2::geom_point() +
      ggplot2::scale_size_manual(values = c(1.5, 1)) +
      ggplot2::scale_color_manual(values = c("black", "grey")) +
      ggplot2::ylab("Differences") + ggplot2::xlab("") + 
      ggplot2::ggtitle("Absolute standardized mean difference") + 
      ggplot2::theme_bw()
    
    return(p)
  }
  
  plot_functions <- list(
    plot_ot,
    plot_pareto_k,
    plot_neff,
    plot_mb
  )
  
  cat( paste0("Plotting diagnostics for causalWeights for estimand ", estimand, "\n") )
  plots_output <- lapply(plot_functions, function(f) do.call(f, list(object = object,
                                                                     ...)))
  
  for(p in plots_output) {
    line <- readline(prompt="Press [enter] to continue")
    print(p)
  }
}

#' plot.causalWeights
#' 
#' @param x A [causalOT::causalWeights-class] object
#' @param r_eff The \eqn{r_\text{eff}} to use for the [causalOT::PSIS_diag()] function.
#' @param penalty The penalty of the optimal transport distance to use. If missing or NULL, the function will try to guess a suitable value depending if debias is TRUE or FALSE.
#' @param p \eqn{L_p} distance metric power
#' @param cost Supply your own cost function. Should take arguments `x1`, `x2`, and `p`.
#' @param debias TRUE or FALSE. Should the debiased optimal transport distances be used.
#' @param online.cost How to calculate the distance matrix. One of "auto", "tensorized", or "online".
#' @param diameter The diameter of the metric space, if known. Default is NULL.
#' @param niter The maximum number of iterations for the Sinkhorn updates
#' @param tol The tolerance for convergence
#' @param ... Not used at this time
#' 
#' @details The plot method first calls summary.causalWeights on the causalWeights object. Then plots the diagnostics from that summary object.
#' 
#' @seealso [causalOT::summary.causalWeights()]
#'
#' @return The plot method returns an invisible object of class summary_causalWeights.
#' @export
#' @method plot causalWeights
plot.causalWeights <- function(x, r_eff = NULL, penalty, p = 2, cost = NULL, 
                               debias = TRUE, online.cost = "auto", diameter = NULL, niter = 1000, 
                               tol = 1e-07, ...) {
  object   <- x
  mc <- match.call()
  mc[[1]] <- quote(summary.causalWeights)
  ml  <- as.list(mc)
  names(ml)[which(names(ml) == "x")] <- "object"
  summary.object <- eval(as.call(ml), envir = asNamespace("causalOT"))
  
  print(plot(summary.object, ...))
  return(invisible(summary.object))
}