sbw_jose <-  function(data, constraint,  estimate = c("ATT", "ATC","ATE","feasible"),
  solver = c("cples", "gurobi")) 
{
    x <- data$get_x()
    z <- data$get_z()
    ns <- get_n(data)
    
    df <- data.frame(z = z, x)
    bal_cov <- colnames(x)
    
    target <- ifelse(estimate == "ATT", "treated", "controls")
    
    sbw_out <- sbw::sbw(df, t_ind = "z", bal_covs = bal_cov, bal_tols = constraint, bal_tols_sd = TRUE,
                   target = target, l_norm = "l_2", w_min = 0, normalize = 1, solver = solver)
    output <- list(w0 = NULL, w1 = NULL, gamma = NULL)
    sol <- sbw_out$data_frame_weights$weights
    
    if(estimate == "ATT") {
      output$w0 <- renormalize(sol[ns["n1"] + 1:ns["n0"]])
      output$w1 <- rep.int(1/ns["n1"],ns["n1"])
    } else if (estimate == "ATC") {
      output$w0 <- rep.int(1/ns["n0"],ns["n0"])
      output$w1 <- renormalize(sol[1:ns["n1"]])
    }
    
    output$estimate <- estimate
    class(output) <- "causalWeights"
    return(output)
    
}