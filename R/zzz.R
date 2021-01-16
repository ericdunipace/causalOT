ot_fun_path <- torch <- NULL

.onLoad <- function(libname, pkgname) {
  # np <<- reticulate::import("numpy", delay_load = TRUE)
  # scipy <<- reticulate::import("scipy", delay_load = TRUE)
  torch <<- reticulate::import("torch", delay_load = TRUE)
  ot_fun_path <<- file.path(system.file(package="causalOT"), "Python/r_to_OT_imputer.py")
}