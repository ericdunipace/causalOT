#### Hainmueller Simulations ####

#### load packages ####
library(causalOT)
library(dplyr)
library(ggplot2)
library(ggsci)
library(xtable)
library(cowplot)

#### Load Data ####
hain.file <- file.path("data","Hainmuller.rds")
if (!file.exists(hain.file)) {
  download.file(url = "https://dataverse.harvard.edu/api/access/datafile/6040527",
                destfile = hain.file)
}

hain <- readRDS(file = hain.file)

#### set true ATE ####
true <- 0

#### tables ####
haintab <- hain %>% 
  filter(estimate < 1e6 & estimate > -1e6) %>% #filter bad 4 estimates from SCM
  filter(model == "lm" & method != "NNM") %>% # not tackle NNM in the paper
  filter(delta == 0 | delta == 1e-4 | is.na(delta)) %>% # make sure proper delta's selected
  filter(penalty == "entropy" | penalty == "none" | is.na(penalty)) %>%  #only sinkhorn
  filter(add.divergence == "TRUE" | is.na(add.divergence)) %>% #only sinkhorn
  filter(wass_p == 2 | is.na(wass_p)) %>% # only squared cost
  filter(match == "FALSE" | is.na(match)) %>% # no barycenter
  
  # select methods
  mutate(method = factor(method, levels = c("Logistic",  
                                            "CBPS",
                                            "SBW",
                                            "SCM", 
                                            "Wasserstein"))) %>% 
  
  # minor data cleaning
  mutate(add.margins = as.logical(add.margins)) %>% 
  
  # translate constraint to easy text for table
  mutate(constraint = ifelse(!is.na(add.margins) & add.margins == TRUE, 1, 0) + 
           ifelse(!is.na(formula) & !grepl("I", formula), 2, 0) +
           ifelse(!is.na(formula) & grepl("I",formula), 5,0),
         p = as.numeric(as.character(p))
  ) %>% 
  mutate(constraint = ifelse(method == "Logistic", 0, constraint )) %>% 
  mutate(constraint = ordered(constraint, levels = c(0, 1,2,3,5),
                              labels = c("none","margins","means","m + m",
                                        "means + 2nd"))) %>%
  mutate(overlap = forcats::fct_relevel(overlap, "high", "medium", "low")) %>% 
  filter((constraint == "none" & method == "Logistic") |
           (constraint == "means" & method == "CBPS") |
           (constraint == "means" & method == "SBW") |
           (constraint == "none" & method == "SCM") |
           method == "Wasserstein") %>% 
  
  # rename methods
  mutate(method = forcats::fct_recode(method,
                                      "COT" = "Wasserstein",
                                      "GLM" = "Logistic"))

# get info for hajek esitmates
non.augment <- haintab %>% 
filter(model.augmentation == FALSE & split.model == TRUE & match == FALSE) %>% 
  group_by(design, overlap, method, constraint) %>% 
  summarise(Hajek = mean(estimate - true, na.rm = TRUE),
            RMSE = sqrt(mean((estimate - true)^2, na.rm = TRUE))) %>% 
  mutate(method = ifelse(duplicated(as.character(method)), "", as.character(method))) %>% 
  ungroup() %>% 
  group_by(design, overlap) %>% 
  mutate(hmin = abs(Hajek) <= min(abs(Hajek)),
         rmin = RMSE <= min(RMSE)) %>% 
  ungroup() %>% 
  mutate(Hajek = case_when(hmin ~ paste0("\\textbf{", format(round(Hajek, 2), width = 4,
                                                                                 justify = "right", nsmall=2), "}"),
                           TRUE ~ format(round(Hajek, digits = 2), nsmall = 2, width = 4, justify = "right")),
         RMSE = case_when(rmin ~ paste0("\\textbf{", format(round(RMSE, 2), justify = "right", nsmall = 2), "}"),
                          TRUE ~ format(round(RMSE, digits = 2), justify = "right", 
                                        nsmall = 2))) %>% 
  ungroup() %>% 
  group_by(design) %>% 
  mutate(overlap =
           ifelse(duplicated(overlap), "", as.character(overlap))) %>% 
  ungroup() %>% group_by(design) %>% 
  mutate(design =
           ifelse(duplicated(design), "", as.character(design))) %>% 
  select(!c(hmin,rmin)) 

#get info for augmented estimates
augment <- haintab %>% 
filter(model.augmentation == TRUE & split.model == TRUE & match == FALSE) %>% 
  group_by(design, overlap, method, constraint) %>% 
  summarise(DR = mean(estimate - true, na.rm = TRUE),
            RMSE = sqrt(mean((estimate - true)^2, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(design, overlap) %>% 
  mutate(dmin = abs(DR) <= min(abs(DR)),
         rmin = RMSE <= min(RMSE)) %>% 
  ungroup() %>% 
  mutate(DR = case_when(dmin ~ paste0("\\textbf{", format(round(DR, 2), width = 4,
                                                             justify = "right", nsmall=2), "}"),
                           TRUE ~ format(round(DR, digits = 2), nsmall = 2, width = 4, justify = "right")),
         RMSE = case_when(rmin ~ paste0("\\textbf{", format(round(RMSE, 2), justify = "right", nsmall = 2), "}"),
                          TRUE ~ format(round(RMSE, digits = 2), justify = "right", 
                                        nsmall = 2))) %>% 
  select(!c(dmin,rmin))

# get info for weighted least squares
wols <- haintab %>% 
  filter(split.model == FALSE & match == FALSE) %>% 
  group_by(design, overlap, method, constraint) %>% 
  summarise(WOLS = mean(estimate - true, na.rm = TRUE),
            RMSE = sqrt(mean((estimate - true)^2, na.rm = TRUE))) %>% 
  ungroup() %>% 
  group_by(design, overlap) %>% 
  mutate(wmin = abs(WOLS) <= min(abs(WOLS)),
         rmin = RMSE <= min(RMSE)) %>% 
  ungroup() %>% 
  mutate(WOLS = case_when(wmin ~ paste0("\\textbf{", format(round(WOLS, 2), width = 4,
                                                             justify = "right", nsmall=2), "}"),
                           TRUE ~ format(round(WOLS, digits = 2), nsmall = 2, width = 4, justify = "right")),
         RMSE = case_when(rmin ~ paste0("\\textbf{", format(round(RMSE, 2), justify = "right", nsmall = 2), "}"),
                          TRUE ~ format(round(RMSE, digits = 2), justify = "right", 
                                        nsmall = 2))) %>% 
  select(!c(wmin,rmin))

# setup table
tab <- cbind(non.augment %>% 
                  ungroup() %>% 
                  select(-c(RMSE)), 
              augment  %>% 
                  ungroup() %>% 
                  select(c("DR")),
              wols  %>% 
                ungroup() %>% 
                select(c("WOLS")),
              non.augment %>% 
                ungroup() %>% 
                select(-c(Hajek)) %>% 
                mutate(Hajek = RMSE) %>% 
                select(c("Hajek")), 
              augment  %>% 
                ungroup() %>% 
                select(-c(DR)) %>% 
                mutate(DR = RMSE) %>% 
                select(c("DR")),
              wols  %>% 
                ungroup() %>% 
                select(-c("WOLS")) %>% 
                mutate(WOLS = RMSE) %>% 
                select(c("WOLS"))
            )

addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- c("\\hline & & & & \\multicolumn{3}{c}{Bias} & \\multicolumn{3}{|c}{RMSE}\\\\")

ncond <- length(unique(non.augment$method))

#print table 
print(xtable(tab,
             align = "lllll|rrr|rrr",
             caption = "Performance of various weighting methods under the simulation settings of \\cite{Hainmueller2012}. Bold values are the values with the lowest bias or root mean-squared error (RMSE) of the methods under the same conditions. GLM refers to weighting by the inverse of the propensity score as calculated from a logistic regression model, CBPS is the covariate balancing propensity score, SBW is the stable balancing weights, SCM is the synthetic control method, and COT is the optimal transport formulation proposed in this paper. The estimators are Hajek weights (Hajek), doubly-robust augmented IPW (DR), and weighted least squares (WOLS). All weights are normalized to sum to 1. Constraints refer to balancing constraints and are one of ``none'' for no constraints or ``mean'' for mean constraints.",
             label = "tab:hain"), 
      hline.after = c(-1, 0, rep(ncond,ncond) + ncond * 0:(ncond-1), nrow(tab)),
      include.rownames = FALSE,
      add.to.row = addtorow,
      sanitize.text.function = function(x){x},
      table.placement = "!htb",
      size="\\fontsize{9pt}{10pt}\\selectfont",
      file = "tables/hainmueller.tex"
)
