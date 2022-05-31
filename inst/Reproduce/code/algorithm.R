#### convergence and ci ####

#### load packages ####
library(causalOT)
library(dplyr)
library(ggplot2)
library(scales)

#### load data ####
alg.file <- file.path("data","algorithm.rds")
# if (!file.exists(conv.file)) {
#   download.file(url = "https://dataverse.harvard.edu/api/access/datafile/6040578",
#                 destfile = conv.file)
# }

alg <- readRDS(file = alg.file)


#### set true ATE ####
true <- 0

#### functions ####
n_labeller <- function (labels, multi_line = TRUE) 
{
  labels <- lapply(labels, function(x) paste0("n = ", x))
  if (multi_line) {
    labels
  }
  else {
    ggplot2:::collapse_labels_lines(labels)
  }
}

scientific_10 <- function(x) {
  text <- gsub("1e", "10^", scales::scientific_format()(x))
  text <- gsub("\\+","",text)
  return(parse(text=text))
}


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

#### graph of selection vs anderson darling statistic ####

algorithm <- alg %>% 
  mutate(id = factor(id)) %>% 
  mutate(lambda = lambda_0) %>% 
  group_by(id, lambda, n, bf ) %>% 
  summarize(ot0 = mean(ot0),
            ot1 = mean(ot1),
            msel0 = mean(sel0),
            msel1 = mean(sel1),
            sel0 = median(sel0),
            sel1 = median(sel1),
            E_Y0 = mean(E_Y0),
            E_Y1 = mean(E_Y1),
            E_Y0.aug = mean(E_Y0.aug),
            E_Y1.aug = mean(E_Y1.aug),
            u_and_c = quantile(and_c, 0.975),
            l_and_c = quantile(and_c, 0.025),
            u_and_t = quantile(and_t, 0.975),
            l_and_t = quantile(and_t, 0.025),
            and_c = mean(and_c),
            and_t = mean(and_t),
            w_0 = mean(w_0),
            w_1 = mean(w_1)
  ) 

selected.mean <- algorithm %>% ungroup() %>% 
  group_by(bf, n) %>%
  # filter(n == 256, bf == FALSE) %>%
  summarize(
    app_c = exp(approx(x = id, y = log(and_c), xout = msel0,
                       method = "linear")$y[1]),
    app_t = exp(approx(x = id, y = log(and_t), xout = msel1,
                       method = "linear")$y[1]),
    msel0 = msel0[1], 
    msel1 = msel1[1]
  )

sink <- algorithm %>%  filter(bf == FALSE)

max_andt <- max(sink$and_t)
max_andc <- max(sink$and_c)
min_andt <- min(sink$and_t)
min_andc <- min(sink$and_c)
lwr_t <- round(log(min_andt)/log(10)-1)
upr_t <- round(log(max_andt)/log(10)+1)
lwr_c <- round(log(min_andc)/log(10)-1)
upr_c <- round(log(max_andc)/log(10)+1)


treated <- alg %>% filter(bf == FALSE, id == sel1) %>% 
  ggplot() +
  geom_histogram(
    mapping = aes(x = sel1),
    stat = "count",
    fill = "#f0f0f0",
    color = "#949494") + 
  geom_line(data = sink  
            , mapping = aes(x = id, y = (log(and_t)-log(10^lwr_t))*50, group = factor(n), 
                            color = factor(n)),
            color = "#1f9aff",#"#1f9aff" blue,"#ff8b1f" orange 
            size = 1) +
  geom_ribbon(data = sink 
              , mapping = aes(x = id, y = (log(and_t)-log(10^lwr_t))*50,
                              ymin = (log(l_and_t)-log(10^lwr_t))*50, ymax = (log(u_and_t)-log(10^lwr_t))*50, 
                              fill = factor(n),
                              group = factor(n)),
              alpha = 0.3, color = NA, fill = "#0062ff"#"#1f9aff""#ff8b1f"
  ) + 
  scale_y_continuous(name= "Anderson-Darling Statistic",
                     sec.axis = sec_axis(trans = ~. , name = "Probability selected",
                                         breaks = seq(0, 1000, 250),
                                         labels = seq(0, 1, .25)
                     ),
                     limits = c(0, 1000),
                     
                     breaks = (log(10^seq(lwr_t, upr_t+3, 1)) - log(10^lwr_t))*50,
                     labels = scientific_10(10^seq(lwr_t, upr_t+3, 1)),
                     expand = c(0.0,0.0)
  )  +
  scale_x_discrete(name = expression("Penalty parameter,"~lambda),
                   labels = scientific_10(unique(sort(sink$lambda)))) +
  facet_wrap(facets=~n,
             labeller = n_labeller) +
  theme_cot()

control <- alg %>% filter(bf == FALSE, id == sel0) %>% 
  ggplot() +
  geom_histogram(
    mapping = aes(x = sel0),
    stat = "count",
    fill = "#f0f0f0",
    color = "#949494") + 
  geom_line(data = sink  
            , mapping = aes(x = id, y = (log(and_c)-log(10^lwr_c))*50, group = factor(n), 
                            color = factor(n)),
            color = "#1f9aff",#"#1f9aff" blue,"#ff8b1f" orange 
            size = 1) +
  geom_ribbon(data = sink 
              , mapping = aes(x = id, y = (log(and_c)-log(10^lwr_c))*50,
                              ymin = (log(l_and_c)-log(10^lwr_c))*50, ymax = (log(u_and_c)-log(10^lwr_c))*50, 
                              fill = factor(n),
                              group = factor(n)),
              alpha = 0.3, color = NA, fill = "#0062ff"#"#1f9aff""#ff8b1f"
  ) + 
  scale_y_continuous(name= "Anderson-Darling Statistic",
                     sec.axis = sec_axis(trans = ~. , name = "Probability selected",
                                         breaks = seq(0, 1000, 250),
                                         labels = seq(0, 1, .25)
                     ),
                     limits = c(0, 1000),
                     breaks = (log(10^seq(lwr_c, upr_c+5, 1)) - log(10^lwr_c))*50,
                     labels = scientific_10(10^seq(lwr_c, upr_c+5, 1))
                     , expand = c(0,0)
                     
  )  +
  scale_x_discrete(name = expression("Penalty parameter,"~lambda),
                   labels = scientific_10(unique(sort(sink$lambda)))) +
  facet_wrap(facets=~n,
             labeller = n_labeller) +
  theme_cot()

txs <- cbind(alg %>%
        filter(bf == FALSE) %>%
        filter(id == sel1) %>%
        select(n, E_Y1),
        alg %>%
        filter(bf == FALSE) %>%
        filter(id == sel0) %>%
        select(E_Y0)
) %>%  group_by(n) %>% summarize(tau = mean(E_Y1 - E_Y0),
                                 sd_tau = sd(E_Y1 - E_Y0),
                                 rmse = sqrt(mean( (E_Y1 - E_Y0)^2)))

#### Save ####
pdf("figures/algorithm_check_treated.pdf",
    width = 7, height = 4)
print(treated)
dev.off()

pdf("figures/algorithm_check_control.pdf",
    width = 7, height = 4)
print(control)
dev.off()
