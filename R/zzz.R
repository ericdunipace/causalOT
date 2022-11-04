# 
# .onLoad <- function(libname, pkgname) {
#   softmin_jit <- torch::jit_load("data/softmin.zip")
#   
#   assign("softmin_jit", softmin_jit, envir = parent.env(environment()))
# }