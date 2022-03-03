testthat::test_that("wassCGD runs", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("Rmosek"); testthat::skip_on_ci()
  causalOT:::skip_if_no_geomloss()
  set.seed(123123)
  data <- Hainmueller$new(n = 512, p = 6, design = "B", overlap = "high")
  data$gen_data()
  
  pd <- prep_data(data)
  z  <- pd$z
  sample_weight <- get_sample_weight(NULL, z)
  
  if (!is.null(pd$df$y)) {
    pd$df$y <- NULL
  }
  x  <- as.matrix(pd$df)
  x1 <- x[z == 1,, drop = FALSE]
  x0 <- x[z == 0,, drop = FALSE]
  
  n <- nrow(x)
  n0 <- nrow(x0)
  n1 <- nrow(x1)
  
  wass.dat <- cbind(z = z, x)
  
  cost <- NULL
  add.margins <- FALSE
  estimand <- "ATT"
  method <- "Wasserstein"
  joint.mapping <- TRUE
  add.joint <- TRUE
  p <- 1
  metric <- "mahalanobis"
  solver <- "mosek"
  
  cost <- wass_cost_default(x = x, z = z, estimand = estimand, 
                            metric = metric, method = method, p = p, 
                            rkhs.args = NULL, 
                            add.margins = add.margins)
  
  
  testthat::expect_silent(wo <- wassCGD$new(X1 = data$get_x0(), X2 = data$get_x1(), 
                    z = data$get_z(),
                    qp_constraint = list(joint = 1,
                                         margins = 1,
                                         penalty = 1),
                    qp_solver = "mosek",
                    qp_penalty = "none",
                    lambda = list(joint = 1,
                                  penalty = 1
                    ),
                    add.mapping = joint.mapping,
                    add.margins = add.margins,
                    penalty = "L2",
                    metric = "mahalanobis",
                    power = 4.,
                    niter = 1000,
                    tol = 1e-3,
                    sample_weight = NULL))
  testthat::expect_silent(cg(wo, verbose = FALSE))
  testthat::expect_equal(wo$return_cw()$args$cur_iter, 1)
})
