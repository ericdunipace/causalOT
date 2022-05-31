#### convergence and ci ####

#### load packages ####
library(causalOT)
library(dplyr)
library(ggplot2)
library(ggsci)
library(xtable)
library(cowplot)
library(tidyr)

#### theme function ####
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
  mutate(name = forcats::fct_recode(name, "S[lambda](w[0], a) " = "w_b_c",
                                    "S[lambda](w[1], a) " = "w_b_t",
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

# rates reported in table
l2_rates <- full.conv.plot %>%
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
  group_by(method, name, weight) %>% 
  summarize(rate = coef(lm(log(value) ~ log(n)))[2])

sink_rates <- full.conv.plot %>%
  filter(!(as.character(method.new) %in%  c("COT, entropy", "COT, L2"))) %>% 
  filter(!(name %in% c("l2_c", "l2_t","and_c", "and_t") )) %>% 
  mutate(tx_group = forcats::fct_recode(name, "control" = "w_b_c",
                                        "treated" = "w_b_t",
                                        "control" = "w_c",
                                        "treated" = "w_t")) %>% 
  mutate(name = forcats::fct_recode(name, "$S_\\lambda(\\boldw_0, \\bolda) $" = "w_b_c",
                                    "$S_\\lambda(\\boldw_1, \\bolda)$ " = "w_b_t",
                                    "$S_\\lambda(\\boldw_0, \\boldw^\\star)$" = "w_c",
                                    "$S_\\lambda(\\boldw_1, \\boldw^\\star)$" = "w_t")) %>% 
  filter(constraint == "none") %>% 
  mutate(method = method.new) %>% 
  group_by(method, name, tx_group) %>% 
  summarize(rate = coef(lm(log(value) ~ log(n)))[2])

rates_tab_control <- sink_rates %>% 
  filter(tx_group == "control") %>% 
  pivot_wider(id_cols = method, names_from = name, values_from = rate)
rates_tab_treated <- sink_rates %>% 
  filter(tx_group == "treated") %>% 
  pivot_wider(id_cols = method, names_from = name, values_from = rate)

rates_tab_control <- rates_tab_control %>% cbind(
  l2_rates %>% 
    filter(name == "control") %>% 
    mutate("$\\|\\boldw_0 - \\boldw^\\star\\|^2$" = rate) %>% 
    ungroup() %>% 
    select("$\\|\\boldw_0 - \\boldw^\\star\\|^2$")
)
rates_tab_treated <- rates_tab_treated %>% cbind(
  l2_rates %>% 
    filter(name == "treated") %>% 
    mutate("$\\|\\boldw_1 - \\boldw^\\star\\|^2$" = rate) %>% 
    ungroup() %>% 
    select("$\\|\\boldw_1 - \\boldw^\\star\\|^2$")
)
rates_tab <- cbind(rates_tab_control,
                   rates_tab_treated[,-1], .name_repair = "minimal")
rates_tab$method <- c("COT","GLM","NNM","SBW")
rates_tab <- rates_tab[c(2,4,3,1),]
add.to.row <- list(pos = list(-1),
                   command = c("  & \\multicolumn{3}{c}{control} & \\multicolumn{3}{|c}{treated}\\\\"))
rates_xtab <- xtable(rates_tab, 
                     caption = "Empirical rates of convergence of the listed methods under various metrics for treated and control groups. The numbers correspond to the power of $n$ at which the given metric decreases to zero---\\textit{e.g.} $n^\\text{rate}$. $S_\\lambda(\\boldw_1, \\bolda)$ is the 2-Sinkhorn divergence between the weights and the full sample, $S_\\lambda(\\boldw_z, \\boldw^\\star)$ is the 2-Sinkhorn divergence between the weights and the self-normalized importance sampling weights, and $\\|\\boldw_1 - \\boldw^\\star\\|^2$ is the $L_2$ norm between the weights and the self-normalized importance sampling weights. GLM is a probit regression (the true model), SBW is Stable Balancing Weights, NNM is nearest neighbor matching, and COT is Causal Optimal Transport.",
                     label = "tab:rates",
                     align = c("ll|ccc|ccc"))
print(rates_xtab, add.to.row = add.to.row, 
      sanitize.colnames.function = function(x){x},
      include.rownames = FALSE,
      file = "tables/conv_rates.tex")

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
