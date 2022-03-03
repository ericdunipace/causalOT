#### Simulation array look-up ####
library(dplyr)

data <- c("KangSchafer",
          "Sonabend",
          "Hainmueller")

n <- 2^(9:13)
n <- 2^(9)

p    <- c(6, 4, 2)

overlap <- c("high", "medium","low")

design <- c("A", "B", "C")

metric <- c("sdLp")
penalty <- list(c("L2","entropy"))
formula <- list(c(NA_character_, "~ . + 0"))

nexperiments <- 1000
expernum <- 1:nexperiments

arrayset <- 1:4
arrayset.idx <- rep(arrayset,each = nexperiments/max(arrayset))

maxarray <- 1e4

df <- expand.grid(experiment.number = expernum,
                  data = data, n = n, p = p, overlap = overlap, 
                  design = design, 
                  penalty = penalty, 
                  metric = metric,
                  formula = formula
                  )
df$arrayset <-  rep(1:max(arrayset), length.out = nexperiments)


df.out <- df %>% filter((data != "Sonabend" & df$p != 2) |
                data == "Sonabend") %>% 
       filter(data != "Hainmueller" | (data == "Hainmueller" & p == 6)) %>% 
       filter((data != "FanEtAl" & design != "C") | (data == "FanEtAl")) %>%
       filter(data != "FanEtAl" | (data == "FanEtAl" & overlap == "high")) %>% 
       mutate(method = list(c("Logistic",
                         "SBW",
                         "NNM",
                         "Wasserstein",
                         "SCM")))


df.out <- df.out %>% filter(data != "FanEtAl") %>% filter(
  (data == "KangSchafer" & p != 6) |
    (data == "Sonabend" & !(p %in% c(4,6))) |
    data == "Hainmueller") %>%  
  # filter(design == "B") %>% 
  filter(arrayset == 1) %>% filter("Hainmueller" == data)


saveRDS(df.out, file = "code/original_sim/Hainmueller/sim_arraylookup.rds")

