# runs vignette one time

if (interactive()) {
  old_wd <- getwd()
  setwd("vignettes/")
  knitr::knit("usage.Rmd.orig", output = "usage.Rmd")
  knitr::purl("usage.Rmd.orig", output = "usage.R")
  
  knitr::knit("oopCOT.Rmd.orig", output = "oopCOT.Rmd")
  knitr::purl("oopCOT.Rmd.orig", output = "oopCOT.R")
  setwd(old_wd)
}

