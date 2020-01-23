DataSim <- R6::R6Class("DataSim",
        public = list(get_x = function() { return( private$x)},
                      get_y = function() { return( private$y)},
                      get_z = function() { return( private$z)},
                      get_n = function() {return(c("n1" = private$n1, "n0" = private$n0))},
                      get_x1 = function() {
                        if(!is.character(private$x)) {
                          if(is.character(private$x1)) private$x1 <- private$x[private$z == 1,,drop=FALSE]
                          return(private$x1)
                        } else {
                          stop("x not initialized yet")
                        }
                      },
                      get_x0 = function(){
                        if(!is.character(private$x)) {
                          if(is.character(private$x0)) private$x0 <- private$x[private$z == 0,,drop=FALSE]
                          return(private$x0)
                        } else {
                          stop("x not initialized yet")
                        }
                      },
                      get_p = function() {
                        return(private$p)
                      }),
         private = list(n = "numeric",
                       p = "numeric",
                       x = "matrix",
                       y = "vector",
                       z = "vector",
                       param = "list",
                       n1 = "numeric",
                       n0 = "numeric",
                       x0 = "vector",
                       x1 = "vector",
                       check_data = function() {
                         complete <- all(is.matrix(private$x) & is.vector(private$z) & is.vector(private$y))
                         if(complete) {
                           private$n1 <- sum(private$z == 1)
                           private$n0 <- sum(private$z == 0)
                         }
                       }#,
                       # weights = "list",
                       # std_diff = "numeric",
                       # wass_param = "list",
                       # qp = "list",
                       # wass_dists = "list"
                       )#,
         # methods = list(outcome = function(estimate = c("ATT", "ATC", "ATE"),
         #                                  method = c("Hajek","HT"), 
         #                                  doubly.robust = c(TRUE, FALSE),
         #                                  weight.type = c("Wasserstein", "SBW", "Logistic")
         #                                  ) {
         #   est <- match.arg(estimate)
         #   meth <- match.arg(method)
         #   wt <- match.arg(weight.type)
         #   dr <- match.arg(doubly.roubst)
         #   
         # },
         # set_qp = function(method = c("Wasserstein", "SBW"),
         #                   estimate = c("ATC", "ATT", "ATE")) {
         #   est <- match.arg(estimate)
         #   meth <- match.arg(meth)
         #   if(meth== "SBW") {
         #     K <- rep(std_diff, p)
         #     qp$sbw <<- qp_sbw(x, z, K)
         #   } else if (meth == "Wasserstein") {
         #     K <- wass_param$max
         #     qp$wass <<- qp_w2(x, z, K, p = wass_param$power,)
         #   }
         # },
         # calc_weights = function(method = c("Wasserstein", "SBW","Logistic"),
         #                   estimate = c("ATC", "ATT", "ATE")) {
         #   est <- match.arg(estimate)
         #   meth <- match.arg(method)
         #   if(meth == "Logistic"){
         #     fit <- glm(z ~ x, family = "binomial")
         #     weights$logistic <<- predict(fit, type = "response")
         #   } else if (meth == "SBW") {
         #     fit <- ROI.plugin.cplex:::solve_OP(qp$sbw)
         #     weights$sbw <<- ROI::solution(fit)
         #   } else if (meth == "Wasserstein") {
         #     fit <- ROI.plugin.cplex:::solve_OP(qp$wass)
         #     weights$wasserstein <<- ROI::solution(fit)
         #   }
         # }
         # )
)
# setOldClass("DataSim")
# #### Generate Y functions ####
# setGeneric("gen_y")
