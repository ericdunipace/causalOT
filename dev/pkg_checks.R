
# setup cran comments
usethis::use_cran_comments()

# release issue
usethis::use_release_issue()

# build documents
devtools::document()

# url checks
urlchecker::url_check()

# command check local
devtools::check(remote = TRUE, manual = TRUE)

# command checks, remote
devtools::check_win_devel(quiet=TRUE)
devtools::check_win_release(quiet = TRUE)
devtools::check_win_oldrelease(quiet = TRUE)

# out <- rhub::check_for_cran(show_status = FALSE) #not work on windows with CRAN versions
# versions with rkeops
out <- rhub::check(platforms=c("ubuntu-gcc-release" ,"fedora-clang-devel","linux-x86_64-rocker-gcc-san"),
                   check_args = "--as-cran", 
                   env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "true", `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true"),
                   show_status = FALSE)

# don't force suggests since rkeops not available on windows
out_windows <- rhub::check(platforms= "windows-x86_64-devel",
                    check_args = "--as-cran", 
                    env_vars = c(`_R_CHECK_FORCE_SUGGESTS_` = "false", `_R_CHECK_CRAN_INCOMING_USE_ASPELL_` = "true"),
                    show_status = FALSE)

# reverse dependency
# run if no rev dep check devtools::install_github('r-lib/revdepcheck')
revdepcheck::revdep_reset()
revdepcheck::revdep_check(num_workers = 4,
                          timeout = as.difftime(240, units = "mins"), show_status = FALSE)


# to submit
devtools::submit_cran()

# once accepted, github release
usethis::use_github_release()

# increment version to dev on my machine
usethis::use_dev_version()
