devtools::install_github("jjchern/lalonde")

library(lalonde)
library(dplyr)

lalonde_nsw <- lalonde::nsw
# usethis::use_data(lalonde_dat, overwrite = TRUE, internal = TRUE)

lalonde_full <- lalonde_nsw %>% filter(treat == 1) %>% 
  rbind(lalonde::cps_controls %>% select(!re74))


usethis::use_data(lalonde_full, lalonde_nsw, overwrite = TRUE, internal = TRUE)