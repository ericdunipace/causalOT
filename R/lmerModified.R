ot_lmer <- function(formula, data = NULL,
                    ot.penalty = NULL, ...) {
  
  browser()
  
  Terms <- terms(as.formula(formula), specials = "tx", data = data)
  txl <- attr(Terms, "specials")$tx
  tx <- attr(Terms, "variables")[[1 + txl]][[2]]
  Terms <- update(Terms, paste0(". ~ . - tx(",tx,")"))
  tx.form <- paste0(as.character(tx), " ~ 0")
  if (missing(weights)) weights <- NA_real_
  dH <- df2dataHolder(outcome.formula = Terms, treatment.formula = tx.form,
                      data = data, weights = weights)
  
  
  x0 <- get_x0(dH)
  x1 <- get_x1(dH)
  y0 <- get_y0(dH)
  y1 <- get_y1(dH)
  z <- get_z(dH)
  
  if(is.null(ot.penalty)) ot.penalty <- 1/get_n(dH)
  w <- get_w(dH)
  w0 <- w[z == 0]
  w1 <- w[z == 1]
  ot <- OT$new(x =x0,
         y = x1,
         a = w0,
         b = w1, debias = TRUE,
         penalty = ot.penalty,
         tensorized = "tensorized")
  
  ot$sinkhorn_opt()
  tmat <- ot$primal()
  pi_01 <- round_pi(as.matrix(tmat$xy), w0, w1)
  
  pi_00 <- round_pi(as.matrix(tmat$xx), w0, w0)
  pi_11 <- round_pi(as.matrix(tmat$yy), w1, w1)
    
  nonzero_01 <- which(pi_01 > 0, arr.ind = TRUE)
  nonzero_00 <- which(pi_00 > 0, arr.ind = TRUE)
  nonzero_11 <- which(pi_11 > 0, arr.ind = TRUE)
  
  # groups  <- factor(rep(1:nrow(nonzero), 2))
  groups  <- factor(c(nonzero_01[,2], nonzero_11[,1]))
  # groups  <- factor(c(nonzero_00[,1], nonzero_01[,1]))
  # groups  <- factor(c(nonzero_01[,2], 1:ncol(pi_01)))
  
  # x0_rep <- x0[nonzero[,1], , drop = FALSE]
  # x1_rep <- x1[nonzero[,2], , drop = FALSE]
  x0_rep <- x0[ nonzero_01[,1], , drop = FALSE]
  # x0_rep <- x0[ nonzero_00[,1], , drop = FALSE]
  x1_rep <- x1[ nonzero_11[,1], , drop = FALSE]
  # x1_rep <- x1[ nonzero_01[,2], , drop = FALSE]
  # x1_rep <- x1
  
  y0_rep <- y0[nonzero_01[,1], drop = FALSE]
  # y0_rep <- y0[nonzero_00[,1], drop = FALSE]
  # y1_rep <- y1[nonzero[,2], drop = FALSE]
  y1_rep <- y1[ nonzero_11[,1], drop = FALSE]
  # y1_rep <- y1[ nonzero_01[,2], drop = FALSE]
  # y1_rep <- y1
  
  z_rep <- c(rep(0, nrow(x0_rep)), 
             rep(1, nrow(x1_rep))
            )
  
  # weights_rep <- rep(c(weights[weights > 0]), 2)
  weights_rep <- c(pi_01[pi_01 > 0],
                   pi_11[pi_11 > 0])
  # weights_rep <- c(pi_00[pi_00>0],
  #                  pi_01[pi_01 > 0]
  #                  )
  # weights_rep <- c(pi_01[pi_01 > 0], 
  #                  colSums(pi_01))
  
  df <- data.frame(y = c(y0_rep, y1_rep),
                   z = z_rep,
                   groups = groups,
                   wt = weights_rep,
                   rbind(x0_rep, x1_rep))
  out_terms <- update(Terms, paste0("y ~ . -", tx))
  out_terms <- update(out_terms, paste0(". ~ . + z"))
  lmerTerms <- update(out_terms, paste0(". ~ . + (z | groups)"))
  fit <- lme4::lmer(lmerTerms, data = df,
                    weights = wt)
  
  # targeting just controls
  groups <- nonzero_00[,2]
  x0_rep <- x0[ nonzero_00[,1], , drop = FALSE]
  y0_rep <- y0[ nonzero_00[,1],   drop = FALSE]
  z_rep <- rep(0, nrow(x0_rep))
  weights_rep <- c(pi_00[pi_00>0]) 
  
  df <- data.frame(y = y0_rep,
                   z = z_rep,
                   groups = groups,
                   wt = weights_rep,
                   x0_rep)
  out_terms <- update(Terms, paste0("y ~ . -", tx))
  lmerTerms <- update(out_terms, paste0(". ~ . + (. | groups)"))
  fit <- lme4::lmer(lmerTerms, data = df,
                    weights = wt)
  
  return(fit)
  
}