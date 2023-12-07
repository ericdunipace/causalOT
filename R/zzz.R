# 
# .onLoad <- function(libname, pkgname) {
#   softmin_jit <- torch::jit_load("data/softmin.zip")
#   
#   assign("softmin_jit", softmin_jit, envir = parent.env(environment()))
# }

# .onAttach <- function(libname, pkgname) {
#   where <- as.environment("package:causalOT")
#   clss <- list(
#     c("DataSim","R6"),
#     c("Hainmueller", "DataSim","R6"),
#     c("OT","R6"),
#     # "cotProblem",
#     c("COT","R6"),
#     c("SCM","R6"),
#     c("balanceFunction", "R6"),
#     c("EntropyBW", "balanceFunction", "R6"),
#     c("SBW", "balanceFunction", "R6")
#   )
#   sapply(clss, function(cls) {
#     try(setOldClass(cls, where = where))
#   })
# }