ot_fun_path <- lp_python_path <- NULL


#' Python path
#'
#' @param libname default args. unused.
#' @param pkgname default args. unused.
#'
#' @details sets variables with paths to python functions after installation.
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # np <<- reticulate::import("numpy", delay_load = TRUE)
  # scipy <<- reticulate::import("scipy", delay_load = TRUE)
  # torch <<- reticulate::import("torch", delay_load = TRUE)
  ot_fun_path <- file.path(system.file(package="causalOT"), "Python","r_to_OT_imputer.py")
  pycot_path <- file.path(system.file(package="causalOT"), "Python")
  
  assign("ot_fun_path", ot_fun_path, envir = parent.env(environment()))
  assign("pycot_path", pycot_path, envir = parent.env(environment()))
}