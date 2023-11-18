
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
out <- rhub::check_for_cran(show_status = FALSE)

# reverse dependency
# run if no rev dep check devtools::install_github('r-lib/revdepcheck')
revdepcheck::revdep_reset()
revdepcheck::revdep_check(num_workers = 4,
                          timeout = as.difftime(240, units = "mins"))


# to submit
devtools::submit_cran()

# once accepted, github release
usethis::use_github_release()

# increment version to dev on my machine
usethis::use_dev_version()
