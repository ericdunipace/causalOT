#### Hainmueller Simulations ####

#### load packages ####
library(causalOT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(xtable)
library(cowplot)
library(ggridges)

#### Load Data ####
hain.file <- file.path("data","Hainmueller.rds")
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
  filter(model == "lm" ) %>% 
  # filter(method != "NNM") %>% # not tackle NNM in the paper
  filter(delta == 0 | delta == 1e-4 | is.na(delta)) %>% # make sure proper delta's selected for GLM
  filter(penalty == "entropy" | penalty == "none" | is.na(penalty)) %>%  #only sinkhorn
  filter(add.divergence == "TRUE" | is.na(add.divergence)) %>% #only sinkhorn
  filter(wass_p == 2 | is.na(wass_p)) %>% # only squared cost
  filter(match == "FALSE" | is.na(match)) %>% # no barycenter
  filter(n == 512) %>% 
  filter(metric == "sdLp" | is.na(metric)) %>%
  # distinct(design, overlap, method, model, model.augmentation, match, split.model,
  #          solver, delta, add.margins, joint.mapping, penalty, formula,
  #          neg.weights, add.divergence, design, overlap, n, p, .keep_all = TRUE) %>% 
  
  # select methods
  mutate(method = factor(method, levels = c("Logistic",  
                                            "CBPS",
                                            "SBW",
                                            "SCM", 
                                            "NNM",
                                            "Wasserstein"))) %>% 
  
  filter(design == "B") %>% 
  
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
           (constraint == "none" & method == "NNM") |
           method == "Wasserstein") %>% 
  
  # rename methods
  mutate(method = forcats::fct_recode(method,
                                      "COT" = "Wasserstein",
                                      "GLM" = "Logistic"))

# get info for hajek esitmates
non.augment <- haintab %>% 
filter(model.augmentation == FALSE & split.model == TRUE & match == FALSE) %>% 
  # group_by(design, overlap, method, constraint) %>% 
  group_by(overlap, method, constraint) %>% 
  summarise(Hajek = mean(estimate - true, na.rm = TRUE),
            RMSE = sqrt(mean((estimate - true)^2, na.rm = TRUE))) %>% 
  mutate(method = ifelse(duplicated(as.character(method)), "", as.character(method))) %>% 
  ungroup() %>% 
  # group_by(design, overlap) %>% 
  group_by(overlap) %>% 
  mutate(hmin = abs(round(Hajek, digits = 2)) <= min(abs(round(Hajek, digits = 2))),
         rmin = round(RMSE, digits = 2) <= min(round(RMSE, digits = 2))) %>% 
  ungroup() %>% 
  mutate(Hajek = case_when(hmin ~ paste0("\\textbf{", format(round(Hajek, 2), width = 4,
                                                                                 justify = "right", nsmall=2), "}"),
                           TRUE ~ format(round(Hajek, digits = 2), nsmall = 2, width = 4, justify = "right")),
         RMSE = case_when(rmin ~ paste0("\\textbf{", format(round(RMSE, 2), justify = "right", nsmall = 2), "}"),
                          TRUE ~ format(round(RMSE, digits = 2), justify = "right", 
                                        nsmall = 2))) %>% 
  ungroup() %>% 
  # group_by(design) %>% 
  mutate(overlap =
           ifelse(duplicated(overlap), "", as.character(overlap))) %>% 
  ungroup() %>% 
  # group_by(design) %>% 
  # mutate(design =
  #          ifelse(duplicated(design), "", as.character(design))) %>% 
  select(!c(hmin,rmin)) 

#get info for augmented estimates
augment <- haintab %>% 
filter(model.augmentation == TRUE & split.model == TRUE & match == FALSE) %>% 
  # group_by(design, overlap, method, constraint) %>% 
  group_by(overlap, method, constraint) %>% 
  summarise(DR = mean(estimate - true, na.rm = TRUE),
            RMSE = sqrt(mean((estimate - true)^2, na.rm = TRUE))) %>% 
  ungroup() %>% 
  # group_by(design, overlap) %>% 
  group_by(overlap) %>% 
  mutate(dmin = abs(round(DR, digits = 2)) <= min(abs(round(DR, digits = 2))),
         rmin = round(RMSE, digits = 2) <= min(round(RMSE, digits = 2))) %>% 
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
  # group_by(design, overlap, method, constraint) %>% 
  group_by(overlap, method, constraint) %>% 
  summarise(WOLS = mean(estimate - true, na.rm = TRUE),
            RMSE = sqrt(mean((estimate - true)^2, na.rm = TRUE))) %>% 
  ungroup() %>% 
  # group_by(design, overlap) %>% 
  group_by(overlap) %>% 
  mutate(wmin = abs(round(WOLS, digits = 2)) <= min(abs(round(WOLS, digits = 2))),
         rmin = round(RMSE, digits = 2) <= min(round(RMSE, digits = 2))) %>%
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
# addtorow$command <- c("\\hline & & & & \\multicolumn{3}{c}{Bias} & \\multicolumn{3}{|c}{RMSE}\\\\")
addtorow$command <- c("\\hline & & & \\multicolumn{3}{c}{Bias} & \\multicolumn{3}{|c}{RMSE}\\\\")

ncond <- length(unique(non.augment$method))
lengthtabgrp <- nrow(tab)/ncond

#print table 
print(xtable(tab,
             align = "llll|rrr|rrr",
             # align = "lllll|rrr|rrr",
             caption = "Performance of various weighting methods under the simulation settings of \\cite{Hainmueller2012}. Bold values are the values with the lowest bias or root mean-squared error (RMSE) of the methods under the same conditions. GLM refers to weighting by the inverse of the propensity score as calculated from a logistic regression model, CBPS is the covariate balancing propensity score, SBW is the stable balancing weights, SCM is the synthetic control method, and COT is the optimal transport formulation proposed in this paper. The estimators are Hajek weights (Hajek), doubly-robust augmented IPW (DR), and weighted least squares (WOLS). All weights are normalized to sum to 1. Constraints refer to balancing constraints and are one of ``none'' for no constraints or ``mean'' for mean constraints.",
             label = "tab:hain"), 
      hline.after = c(-1, 0, rep(ncond,lengthtabgrp) + ncond * 0:(lengthtabgrp-1), nrow(tab)),
      include.rownames = FALSE,
      add.to.row = addtorow,
      sanitize.text.function = function(x){x},
      table.placement = "!htb",
      size="\\fontsize{9pt}{10pt}\\selectfont",
      file = "tables/hainmueller.tex"
)

#### Distributional balance ####
theme_cot <- function(base_size = 11, base_family = "", 
                      base_line_size = base_size/22, 
                      base_rect_size = base_size/22,
                      legend.position = 'bottom',
                      legend.box = "horizontal",
                      legend.justification = "center",
                      legend.margin = ggplot2::margin(0,0,0,0),
                      legend.box.margin = ggplot2::margin(-10,-10,0,-10)) { 
  ggplot2::`%+replace%`(ggplot2::theme_bw(base_size = base_size, base_family = "", 
                                          base_line_size = base_line_size, 
                                          base_rect_size = base_rect_size),
                        ggplot2::theme(
                          plot.title = ggplot2::element_text(hjust = 1),
                          panel.grid.minor = ggplot2::element_blank(),
                          panel.grid.major = ggplot2::element_blank(),
                          strip.background = ggplot2::element_blank(),
                          strip.text.x = ggplot2::element_text(face = "bold"),
                          strip.text.y = ggplot2::element_text(face = "bold"),
                          legend.position = legend.position,
                          legend.box = legend.box,
                          legend.justification = legend.justification,
                          legend.margin = legend.margin,
                          legend.box.margin = legend.box.margin
                          # change stuff here
                        ))
}
  
  haintab %>%
    group_by(overlap, method, constraint) %>% 
    summarise(m0 = mean(w0),
              v0 = var(w0),
              m1 = mean(w1),
              v1 = var(w1))
  
  haintab %>%
    ggplot(aes(x = method, y = w0, color = constraint)) +
    geom_point(position = "jitter") +
    facet_grid(rows = vars(overlap))

 #  pdf(file = "figures/hain_distribution.pdf", 
 #      width = 5.5, height = 2)  
  haintab %>%
    filter(model.augmentation == FALSE & split.model == TRUE & match == FALSE) %>%
    mutate(treated = w1, control = w0) %>%
    pivot_longer( cols  = treated:control,
                  names_to = "txgrp", values_to = "sink") %>%
    ggplot(aes(y = method, x = sink, fill = constraint)) +
    scale_fill_jama() +
    geom_density_ridges() +
    facet_grid(rows = vars(overlap), cols = vars(txgrp)) +
    scale_x_continuous(name = "2-Sinkhorn distance", limits = c(0,.15)) +
    theme_cot()
 # dev.off()
 
 pdf(file = "figures/hain_distribution.pdf", 
     width = 6.5, height = 4)  
 haintab %>%
   filter(model.augmentation == FALSE & split.model == TRUE & match == FALSE) %>%
   mutate(treated = w1, control = w0) %>% 
   pivot_longer( cols  = treated:control, 
                 names_to = "txgrp", values_to = "sink") %>% 
   ggplot(aes(x = method, y = sink, color = constraint)) +
   scale_fill_jama() + 
   scale_color_jama() +
   geom_boxplot() +
   facet_grid(rows = vars(overlap), cols = vars(txgrp)) + 
   scale_y_continuous(name = "2-Sinkhorn distance", expand = c(0,0)) +
   coord_cartesian(ylim = c(0,0.5)) +
   theme_cot()
 dev.off()
 
