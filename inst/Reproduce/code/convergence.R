#### convergence and ci ####

#### load packages ####
library(causalOT)
library(dplyr)
library(ggplot2)
library(ggsci)
library(xtable)
library(cowplot)
library(tidyr)

#expose plot theme
theme_cot <- causalOT:::theme_cot

#### load data ####
conv.file <- file.path("data","convergence.rds")
if (!file.exists(conv.file)) {
  download.file(url = "https://dataverse.harvard.edu/api/access/datafile/6040578",
                destfile = conv.file)
}

conv <- readRDS(file = conv.file)


#### set true ATE ####
true <- 0

#### graph by design, overlap ####
conv$which.wass <- NA_integer_
conv$which.wass[conv$n0 <= 2^9] <- c(rep(0,3), 1:6)
conv$which.wass[conv$n0 > 2^9 & conv$n0 <= 2^10] <- c(rep(0,3), 1,3:4,6)
conv$which.wass[conv$n0 > 2^10] <- c(rep(0,3), 1)

# get plots
full.conv.plot <- conv %>% 
  mutate(method.new = sapply(which.wass, function(ww) switch(ww + 1L,
                                                             NA,
                                                             "COT, divergence",
                                                             "COT, L2",
                                                             "COT, entropy",
                                                             "COT, divergence",
                                                             "COT, L2",
                                                             "COT, entropy")   )) %>% 
  mutate(method.new = ifelse(which.wass == 0, method, method.new)) %>%
  mutate(method.new = factor(method.new)) %>% 
  mutate(method.new = forcats::fct_recode(method.new, "GLM" = "Probit")) %>% 
  mutate(method.new = factor(method.new, labels = sort(levels(method.new)),
         levels = sort(levels(method.new)))) %>% 
  mutate(constraint = ifelse(which.wass < 4, "none", "means")) %>% 
  tidyr::pivot_longer(cols = c( "w_b_c", "w_b_t", 
                                "w_c"  , "w_t",
                                "l2_c" , "l2_t",
                                "and_c", "and_t")) %>% 
  mutate(constraint = factor(constraint, levels = c("none","means"))) %>% 
  mutate(name = forcats::fct_recode(name))


#plot for sinkhorn div
fconvpw <- full.conv.plot %>%
  filter(!(as.character(method.new) %in%  c("COT, entropy", "COT, L2"))) %>% 
  filter(!(name %in% c("l2_c", "l2_t","and_c", "and_t") )) %>% 
  mutate(name = forcats::fct_recode(name, "S[lambda](w[0], b) " = "w_b_c",
                                    "S[lambda](w[1], b) " = "w_b_t",
                                    "S[lambda](w[0], w^'*')" = "w_c",
                                    "S[lambda](w[1], w^'*')" = "w_t")) %>% 
  filter(constraint == "none") %>% 
  mutate(method = method.new) %>% 
  group_by(n, method, name) %>% 
  summarize(E = mean(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>% 
  ggplot(aes(y = E, x = n, color = method,
             fill = method)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = ggsci::pal_jama()(7)[c(1,4:7)]) +  
  scale_fill_manual(values = ggsci::pal_jama()(7)[c(1,4:7)]) +  
  xlab("N") + ylab("2-Sinkhorn Divergence") +
  facet_wrap(~name, nrow = 2, labeller = label_parsed) + 
  theme_cot() + 
  theme(
    strip.background = element_blank() #,
  ) + 
  scale_x_continuous(trans = "log2",
                     breaks = 2^(seq(6, 12, 2))) +
  scale_y_continuous(trans = "log",labels = scales::scientific,
                     breaks = 10^(-5:-1))

#plot for L_2 convergence
fconvpl <- full.conv.plot %>%
  filter(!(as.character(method.new) %in%  c("COT, entropy", "COT, L2"))) %>% 
  filter((name %in% c("l2_c", "l2_t") )) %>% 
  mutate(weight = forcats::fct_recode(name, "L2" = "l2_c",
                                      "L2" = "l2_t",
                                      "Anderson" = "and_c",
                                      "Anderson" = "and_t"
  )) %>% 
  mutate(name = forcats::fct_recode(name, "control" = "l2_c",
                                    "treated" = "l2_t",
                                    "control" = "and_c",
                                    "treated" = "and_t"
  )) %>% 
  filter(constraint == "none") %>% 
  mutate(method = method.new) %>% 
  group_by(n, method, name, weight) %>% 
  summarize(E = mean(value),
            lwr = quantile(value, 0.025),
            upr = quantile(value, 0.975)) %>% 
  ggplot(aes(y = E, x = n, color = method,
             fill = method)) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = ggsci::pal_jama()(7)[c(1,4:7)]) + 
  scale_fill_manual(values = ggsci::pal_jama()(7)[c(1,4:7)]) + 
  xlab("N") + ylab(expression("|"~ w^"*" - w ~"|"^2)) +
  facet_wrap(weight~name, nrow = 2
  ) + 
  facet_wrap(~name, nrow = 1
  ) + 
  theme_cot() + 
  theme(
    strip.background = element_blank() #,
  ) + 
  scale_x_continuous(trans = "log2",
                     breaks = 2^(seq(6, 12, 2))) +
  scale_y_continuous(trans = "log",labels = scales::scientific,
                     breaks = 10^(-6:-1))


#grayscale safe images
gray.labels <- c("COT", "GLM","NNM","SBW")
pdf("figures/gray_compare_meth_conv_l2.pdf",
    width = 7, height = 4)
print(fconvpl + geom_line(aes(linetype = method), size = 0.7) +
        scale_linetype_manual(values = 1:4,labels = gray.labels) +
        scale_color_manual(values = rje::cubeHelix(6)[1:4], labels = gray.labels) +
        scale_fill_manual(values = rje::cubeHelix(6)[1:4], labels = gray.labels))
dev.off()


pdf("figures/gray_compare_meth_conv_wass.pdf",
    width = 7, height = 4)
print(fconvpw + geom_line(aes(linetype = method), size = 0.7) +
        scale_linetype_manual(values = 1:4,labels = gray.labels) +
        scale_color_manual(values = rje::cubeHelix(6)[1:4], labels = gray.labels) +
        scale_fill_manual(values = rje::cubeHelix(6)[1:4], labels = gray.labels))
dev.off()

#### confidence interval ####
gfull.ci.plot <- conv %>% 
  mutate(method.new = sapply(which.wass, function(ww) switch(ww + 1L,
                                                             NA,
                                                             "COT, divergence",
                                                             "COT, L2",
                                                             "COT, entropy",
                                                             "COT, divergence",
                                                             "COT, L2",
                                                             "COT, entropy")   )) %>% 
  mutate(method.new = ifelse(which.wass == 0, method, method.new)) %>%
  mutate(method.new = factor(method.new)) %>% 
  mutate(method.new = forcats::fct_recode(method.new, "GLM" = "Probit")) %>% 
  mutate(constraint = ifelse(which.wass < 4, "none", "means")) %>% 
  mutate(constraint = factor(constraint, levels = c("none","means"))) %>% 
  group_by(n, method.new, constraint) %>% 
  summarize(cover.E = mean(ci.lwr <= mean(estimate) & ci.upr >= mean(estimate)),
            cover = mean(ci.lwr <= 0 & ci.upr >= 0),
            aug_cover.E = mean(ci.lwr.aug <= mean(estimate.aug) & ci.upr >= mean(estimate.aug)),
            aug_cover = mean(ci.lwr.aug <= 0 & ci.upr.aug >= 0),
            method = unique(method)) %>% 
  ungroup()


# get plot data
gci.plot.long <- gfull.ci.plot %>% 
  filter(method == "Wasserstein") %>% 
  filter(method.new == "COT, divergence") %>% 
  tidyr::pivot_longer(cols = c( "cover", "aug_cover")) %>%
  mutate(name = forcats::fct_recode( name, "Coverage of Hajek estimator" = "cover",
                                     "Coverage of augmented estimator" = "aug_cover"))

gcip_full <- gci.plot.long %>% 
  mutate(method = method.new) %>% 
  ggplot(aes(y = value, x = n,
             linetype = constraint)) +
  geom_hline(yintercept = 0.95) +
  geom_line(size = 1, color = "#374E55FF") +
  xlab("N") + ylab("Asymptotic C.I. coverage") +
  theme_cot() +
  theme(
    strip.background = element_blank() #,
  ) + 
  scale_x_continuous(trans = "log2",
                     breaks = 2^(seq(6, 12, 2))) +
  scale_y_continuous(breaks = seq(0.75, 1, 0.05),
                     minor_breaks = waiver(),
                     limits = c(0.75,1),
                     expand = c(0,0)) +
  facet_grid(cols = vars(name))

gcip_full_E <- gfull.ci.plot %>% 
  filter(method == "Wasserstein") %>% 
  filter(method.new == "COT, divergence") %>% 
  tidyr::pivot_longer(cols = c( "cover.E", "aug_cover.E")) %>%
  mutate(name = forcats::fct_recode( name, "Coverage of Hajek estimator" = "cover.E",
                                     "Coverage of augmented estimator" = "aug_cover.E")) %>% 
  mutate(method = method.new) %>% 
  ggplot(aes(y = value, x = n,
             linetype = constraint)) +
  geom_hline(yintercept = 0.95) +
  geom_line(size = 1, color = "#374E55FF") +
  xlab("N") + ylab("Asymptotic C.I. coverage") +
  theme_cot() +
  theme(
    strip.background = element_blank() #,
  ) + 
  scale_x_continuous(trans = "log2",
                     breaks = 2^(seq(6, 12, 2))) +
  scale_y_continuous(breaks = seq(0.75, 1, 0.05),
                     minor_breaks = waiver(),
                     limits = c(0.75,1),
                     expand = c(0,0)) +
  facet_grid(cols = vars(name))



pdf("figures/gray_ot_meth_ci_coverage_both.pdf",
    width = 7, height = 2.5)
print(gcip_full)
dev.off()

pdf("figures/gray_ot_meth_ci_coverage_both_expectation.pdf",
    width = 7, height = 2.5)
print(gcip_full_E)
dev.off()
