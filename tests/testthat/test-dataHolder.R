testthat::test_that("dataHolder object forms with matrix", {
  set.seed(20348)
  n <- 15
  d <- 4
  x <- matrix(stats::rnorm(n*d), n, d)
  z <- rbinom(n, 1, prob = 0.5)
  weights <- rep(1/n,n)
  
  testthat::expect_silent(dh <-causalOT::dataHolder(x = x, z = z))
  testthat::expect_silent(dh <-causalOT::dataHolder(x = x, z = z,
                                                    y = NULL))
  testthat::expect_silent(dh <-causalOT::dataHolder(x = x, z = z,
                                                    y = NULL,
                                                    weights = weights ))
  
  testthat::expect_error(dh <-causalOT::dataHolder(x = x))
  testthat::expect_error(dh <-causalOT::dataHolder(z = z))
  testthat::expect_error(dh <-causalOT::dataHolder(y = NULL))
  
  
})

testthat::test_that("dataHolder object forms with df2dataHolder", {
  set.seed(20348)
  n <- 15
  d <- 4
  x <- matrix(stats::rnorm(n*d), n, d)
  z <- stats::rbinom(n, 1, prob = 0.5)
  y <- stats::rnorm(n)
  weights <- rep(1/n,n)
  df <- data.frame(x, z, y)
  
  testthat::expect_silent(data <- causalOT:::df2dataHolder(treatment.formula = "z ~ .", outcome.formula = NA_character_, data = df))
  
  testthat::expect_true("y" %in% colnames(data@x))
  
  testthat::expect_true("terms" %in% names(attributes(data)))
  
  testthat::expect_silent(data <- causalOT:::df2dataHolder(treatment.formula = "z ~ .", outcome.formula = "y ~ .", data = df))
  
  testthat::expect_true(!("y" %in% colnames(data@x)))
  
  
})

testthat::test_that("dataHolder forms with DataSim objects", {
  
  set.seed(12312)
  hain <- causalOT:::Hainmueller$new(128)
  hain$gen_data()
  testthat::expect_silent(
    data <- dataHolder(hain)
  )
  testthat::expect_true(inherits(data, "dataHolder"))
  
  
  
})
