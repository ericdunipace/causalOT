# main cost function in the package
cost <- function(x, y, p = 2, tensorized = TRUE, cost_function = NULL) {
  if(missing(p) || is.null(p) || is.na(p)) p <- 2L
  cfm <- (missing(cost_function) || is.null(cost_function))
  if (inherits(cost_function, "character") ) {
    costObj <- costOnline$new(x, y, p = p, cost_function = cost_function)
    
  } else if (cfm && isTRUE(tensorized)) {
    costObj <- costTensor$new(x, y, p = p, cost_function = cost_function)
  } else if (isFALSE(tensorized) && cfm) {
    costObj <- costOnline$new(x, y, p = p, cost_function = cost_function)
  } else if(inherits(cost_function, "function") && isTRUE(tensorized) ) {
    costObj <- costTensor$new(x, y, p = p, cost_function = cost_function)
  } else if (is.function(cost_function) && isFALSE(tensorized) ){
    stop("cost_function must be either not given or a keops recognized character of functions.")
  } else {
    stop("cost options combination not accounted for! Please report this bug.")
  }
  return(costObj)
}
costParent <- R6::R6Class("cost",
            public = list(
              data = "torch_tensor",
              fun = "function",
              p   = "numeric",
              reduction = "function",
              algorithm = "character"
            ))
costTensor <- R6::R6Class("costTensor",
         inherit = costParent,
         public = list(
           initialize = function(x, y, p = 2, cost_function = NULL) {
             self$p = as.numeric(p)
             cfm <- (missing(cost_function) || is.null(cost_function))
             if (cfm) {
               self$fun <- function(x1, x2, p) {
                 if(!inherits(x1, "torch_tensor")) {
                   x1 <- torch::torch_tensor(x1, dtype = torch::torch_double())
                 }
                 if(!inherits(x2, "torch_tensor")) {
                   x2 <- torch::torch_tensor(x2, dtype = torch::torch_double())
                 }
                 ((1/p) * torch::torch_cdist(x1 = x1, 
                                             x2 = x2, 
                                             p = p)^p)$contiguous()
               }
               if( p  == 1) {
                 self$algorithm = "L1"
               } else if (p == 2) {
                 self$algorithm = "squared.euclidean"
               } else {
                 self$algorithm = "other"
               }
             } else if(inherits(cost_function, "function") ) {
              # self$fun <- function(x, y, p) {
              #   (1/p) * cost_function(x, y, p)^p
              # }
               self$fun <- cost_function                          
               self$algorithm = "user"
             } else {
               stop("cost function not found. please report this bug")
             }
             self$data =  self$fun(x, y, p)
             
             if(!inherits(self$data, "torch_tensor")) {
               self$data <-torch::torch_tensor(
                 self$data,
               dtype = torch::torch_double())$contiguous()
             }
           }
         ),
         active = list(
           to_device = function(device) {
             if(missing(device) || is.null(device) ){
               return(NULL)
             }
             self$data <- self$data$to(device = device)
             return(invisible(self))
           }
         )
)
costOnline <- R6::R6Class("costOnline",
         inherit = costParent,
         public = c(
           initialize = function(x, y, p = 2, cost_function = NULL) {
             self$p = p
             if (missing(cost_function) || is.null(cost_function) || is.na(cost_function)) {
               self$algorithm <- if (p == 2) {
                 "squared.euclidean"
               } else if (p == 1) {
                 "L1"
               } else {
                 "other"
               }
               cost_function <-  if (p == 2) {
                 "(SqDist(X,Y) / IntCst(2))"
               } else if (p == 1) {
                 "Sum(Abs(X-Y))"
               } else if (is.integer(p)) {
                 paste0("(Sum(Pow(Abs(X-Y),",p,")) /  IntCst(",p,"))")
               } else {
                 stop("'p' must be an integer for online cost functions.")
               }
             } else if (!is.character(cost_function)) {
               stop("cost_function must be not provided or a keops recognized character function.")
             } else {
               self$algorithm <-  "user"
             }
             
            self$data = list(x = x, 
                             y = y)
            self$fun = cost_function
            self$reduction = function(...){NULL}
           }
         ),
         active = list(
           to_device = function(device) {
             if (missing(device) || is.null(device)) {
               return(NULL)
             }
             self$data$x <- self$data$x$to(device = device)
             self$data$y <- self$data$y$to(device = device)
             return(invisible(self))
           }
      )
)


to_device <- function(cost, device) {UseMethod("to_device")}
# setGeneric("to_device", function(cost, device) standardGeneric("to_device"))

# setOldClass(c("costParent","R6"))
# setOldClass(c("costTensor","costParent"))
# setOldClass(c("costOnline", "costParent"))

to_device.costTensor <- function(cost, device) {
  function(cost, device) {
    cost$data <- cost$data$to(device = device)
    return(cost)
  }
}

# setMethod("to_device", signature(cost = "costTensor", device = "ANY"),
# function(cost, device) {
#   cost$data <- cost$data$to(device = device)
#   return(cost)
# }
# )

to_device.costOnline <- function(cost, device) {
  cost$data <- list(x = cost$data$x$to(device = device),
                    y = cost$data$y$to(device = device))
  return(cost)
}

# setMethod("to_device", signature(cost = "costOnline", device = "ANY"),
#           function(cost, device) {
#             cost$data <- list(x = cost$data$x$to(device = device),
#                               y = cost$data$y$to(device = device))
#             return(cost)
#           }
# )

update_cost <- function(cost, x, y) {UseMethod("update_cost")}
setGeneric("update_cost", function(cost, x, y) standardGeneric("update_cost"))

update_cost.costOnline <- function(cost, x, y) {
  n <- nrow(cost$data$x)
  m <- nrow(cost$data$y)
  stopifnot("data for cost rows has different number of rows" = (n == nrow(x)))
  stopifnot("data for cost columns has different number of rows" = (m == nrow(y)))
  stopifnot("data must have same number of columns" = ncol(x) == ncol(y))
  cost$data <- list(x = x, y = y)
}  
# setMethod("update_cost", signature(cost = "costOnline", x = "ANY", y = "ANY"),
# function(cost, x, y) {
#   n <- nrow(cost$data$x)
#   m <- nrow(cost$data$y)
#   stopifnot("data for cost rows has different number of rows" = (n == nrow(x)))
#   stopifnot("data for cost columns has different number of rows" = (m == nrow(y)))
#   stopifnot("data must have same number of columns" = ncol(x) == ncol(y))
#   cost$data <- list(x = x, y = y)
# }          
# )

update_cost.costTensor <- function(cost, x, y) {
  nm <- dim(cost$data)
  device <- cost$data$device
  dtype <- cost$data$dtype
  stopifnot("data for rows has different number of rows" = (nm[1] == nrow(x)))
  stopifnot("data for columns has different number of rows" = (nm[2] == nrow(y)))
  stopifnot("data must have same number of columns" = ncol(x) == ncol(y))
  cost$data <- cost$fun(x,y,cost$p)$to(device = device, dtype = dtype)
} 
# setMethod("update_cost", signature(cost = "costTensor", x = "ANY", y = "ANY"),
# function(cost, x, y) {
#   nm <- dim(cost$data)
#   device <- cost$data$device
#   dtype <- cost$data$dtype
#   stopifnot("data for rows has different number of rows" = (nm[1] == nrow(x)))
#   stopifnot("data for columns has different number of rows" = (nm[2] == nrow(y)))
#   stopifnot("data must have same number of columns" = ncol(x) == ncol(y))
#   cost$data <- cost$fun(x,y,cost$p)$to(device = device, dtype = dtype)
# }          
# )
