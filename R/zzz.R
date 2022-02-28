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
  ot_fun_path <<- file.path(system.file(package="causalOT"), "Python/r_to_OT_imputer.py")
  lp_python_path <<- file.path(system.file(package="causalOT"), "Python/python_lp.py")
  # cot_stanmodels <- vector("list", 3L)
  # names(cot_stanmodels) <- c("barycenter_projection", "gp_hyper_dose", "gp_hyper")
  # assign("cot_stanmodels", cot_stanmodels, envir = topenv())
}