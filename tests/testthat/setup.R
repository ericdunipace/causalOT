# rkeops >= 2.0 features a python dependency, which, at this time,
# requires users to manage both a Python install and R install
if ( rlang::is_installed("rkeops") ) {
  two_debian_catch <- pkg_vers_number("2.0")
  vers_match   <- tryCatch(utils::packageVersion("rkeops") >= two_debian_catch,
                           error = function(e) FALSE)
  not_on_cran <-  !(!interactive() && !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false"))))
  if (not_on_cran && vers_match && rlang::is_installed("reticulate")) {
    if (is.character(find.package("reticulate"))) {
      rkeops_pyenv <- "rkeops_test_asdf"
      # setup.mess <- reticulate::py_capture_output({
      reticulate::virtualenv_create(envname = rkeops_pyenv)
      reticulate::use_virtualenv(virtualenv = rkeops_pyenv, required = TRUE)
      reticulate::py_config()
      reticulate::py_install("pykeops")
      # })
      mess <- reticulate::py_capture_output(
        testthat::capture_messages(rkeops_good <- rkeops::check_rkeops())
      )
      stopifnot("rkeops not working" = rkeops_good)
      withr::defer({
        reticulate::virtualenv_remove(envname = rkeops_pyenv)
      }, testthat::teardown_env())
    }
  }
}
