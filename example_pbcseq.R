## Load required packages
library(survival); library(dplyr); library(ggplot2); library(ggpubr); 
library(FKmL); library(tidyr)

## -------------------------------------------------------------------------- ##
## Real data example
## -------------------------------------------------------------------------- ##

## pbcseq data
pbcseq %>% summary
pbcseq %>% dim
pbcseq$id %>% unique %>% length
pbcseq <- pbcseq %>% dplyr::select(id, day, bili, chol, albumin, alk.phos, ast, platelet, protime)

## Select variables of interest
df_pbc_refined <- pbcseq %>% dplyr::select(id, day, bili, chol, albumin, alk.phos, ast, platelet, protime)

## Remove missing values and keep only subjects with >= 3 time points
df_pbc_refined <- df_pbc_refined %>%
  na.omit %>% 
  group_by(id) %>% 
  mutate(n = n()) %>% 
  filter(n >= 3) %>% 
  select(-n) %>% 
  data.frame

df_pbc_refined %>% summary
df_pbc_refined %>% dim
df_pbc_refined$id %>% unique %>% length

## Log-transformation
df_pbc_refined <- df_pbc_refined %>% mutate(
  log_bili = log(bili),
  log_chol = log(chol),
  log_alk.phos = log(alk.phos),
  log_ast = log(ast),
  log_platelet = log(platelet)
)

df_pbc_refined %>% summary
df_pbc_refined %>% dim
df_pbc_refined$id %>% unique %>% length

## min-max scaling
df_pbc_refined <- df_pbc_refined %>%
  dplyr::select(id, day, log_bili, log_chol, albumin, 
                log_alk.phos, log_ast, log_platelet, protime) %>% 
  mutate(across(
    c("day", "log_bili", "log_chol", "albumin", 
      "log_alk.phos", "log_ast", "log_platelet", "protime"),
    ~ (. - min(.)) / (max(.) - min(.))
  ))

df_pbc_refined %>% summary

## -------------------------------------------------------------------------- ##
## Applying MFKmL
## -------------------------------------------------------------------------- ##

## `scales` vector
scales_pbc_refined <- rep(1, ncol(df_pbc_refined[, -1]))

## Apply MFKmL with equal weights
set.seed(11)
MFclust_pbc_refined_eqwght <- mfkml(dt = df_pbc_refined, clt_n = 3, 
                                    scales = scales_pbc_refined,
                                    weight = 1, maxIter = 50)

MFclust_pbc_refined_eqwght$Cluster %>% head()
df_pbc_refined$id %>% unique %>% head()
MFclust_pbc_refined_eqwght$Center %>% head()

p1 <- df_pbc_refined %>%
  ggplot(aes(x = day, y = log_bili, group = id)) +
  geom_line(alpha = 0.7, linewidth = 0.8, color = "gray") +
  geom_line(
    data = MFclust_pbc_refined_eqwght$Center,
    aes(x = day, y = log_bili, group = Cluster, color = factor(Cluster)),
    linewidth = 1.2
  ) +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray90")) +
  labs(
    x = "Day",
    y = "log(Bilirubin)",
    color = "Cluster"
  ) +
  theme(legend.position = "none")

p2 <- df_pbc_refined %>%
  ggplot(aes(x = day, y = log_chol, group = id)) +
  geom_line(alpha = 0.7, linewidth = 0.8, color = "gray") +
  geom_line(
    data = MFclust_pbc_refined_eqwght$Center,
    aes(x = day, y = log_chol, group = Cluster, color = factor(Cluster)),
    linewidth = 1.2
  ) +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray90")) +
  labs(
    x = "Day",
    y = "log(Cholesterol)",
    color = "Cluster"
  )

ggarrange(p1, p2, nrow = 1, ncol = 2, align = 'hv', common.legend = TRUE, legend = "top")

## Apply MFKmL with unequal weights (randomly assigned)
df_pbc_refined_wght <- data.frame(id = unique(df_pbc_refined$id))
set.seed(11)
df_pbc_refined_wght$weight <- sample(1:10, nrow(df_pbc_refined_wght), replace = T)

set.seed(11)
MFclust_pbc_refined_uneqwght <- mfkml(dt = df_pbc_refined, clt_n = 3,
                                      scales = scales_pbc_refined,
                                      weight = df_pbc_refined_wght, 
                                      maxIter = 50)
MFclust_pbc_refined_uneqwght$Cluster %>% head(10)
MFclust_pbc_refined_uneqwght$Center


## -------------------------------------------------------------------------- ##
## Applying SFKmL: Tuning the L1 norm constraint (s) only
## -------------------------------------------------------------------------- ##

## Compute Fréchet distance array
dist_ary_pbc_refined <- dist.array(dt = df_pbc_refined, 
                                   time_scale = scales_pbc_refined[1], 
                                   var_scales = scales_pbc_refined[-1])

## Candidate l1bound
l1b_pbc_refined <- seq(1.05, (sqrt(7) - 0.05), length.out = 10)

## Permutation-based gap statistic for l1bound selection
set.seed(11)
SFclust.permute_pbc_refined <- SFclust.permute(dist.ary = dist_ary_pbc_refined, 
                                               k = 3, nperms = 50, 
                                               l1b = l1b_pbc_refined)

SFclust.permute_pbc_refined

print(SFclust.permute_pbc_refined$gapplot.l1b)
print(SFclust.permute_pbc_refined$gapplot.nnz)

## -------------------------------------------------------------------------- ##
## Applying SFKmL: Joint tuning of time-scale (λ) and L1 constraint (s)
## -------------------------------------------------------------------------- ##

## Candidate time-scaling
const_c <- c(0.01, 0.1, 0.25, 0.5, 0.75, 1, seq(2, 5, by = 1))
lambda_pbc_refined <- scales_pbc_refined[1] * const_c

## Loop over λ values to compute gap statistics for each
dist_ary_list <- vector("list", length(const_c))
res_sfperm_list <- vector("list", length(const_c))
df_gaps_list <- vector("list", length(const_c))

for (j in seq_along(const_c)) {
  dist_ary_list[[j]] <- dist.array(dt = df_pbc_refined, 
                                   time_scale = scales_pbc_refined[1] * const_c[j], 
                                   var_scales = scales_pbc_refined[-1])
  
  set.seed(11)
  res_sfperm_list[[j]] <- SFclust.permute(dist.ary = dist_ary_list[[j]], k = 3, 
                                          nperms = 50, l1b = l1b_pbc_refined)
  
  df_gaps_list[[j]] <- data.frame(
    const_c = rep(const_c[j], length(res_sfperm_list[[j]]$l1bounds)),
    lam = rep(lambda_pbc_refined[j], length(res_sfperm_list[[j]]$l1bounds)),
    l1b = res_sfperm_list[[j]]$l1bounds,
    gap = res_sfperm_list[[j]]$gaps,
    sdgap = res_sfperm_list[[j]]$sdgaps
  )
}

## (1) Find the best (c, s) combination based on maximum gap
df_gaps <- bind_rows(df_gaps_list)
df_gaps[which.max(df_gaps$gap), ]
df_gaps_max_by_lam <- df_gaps %>% 
  group_by(const_c, lam) %>% 
  summarise(gap_max = max(gap), .groups = "drop")

c_optimal <- df_gaps$const_c[which.max(df_gaps$gap)]
l1b_optimal <- df_gaps$l1b[which.max(df_gaps$gap)]
c_opt_index <- which(const_c == c_optimal)

## Apply SFKmL with optimal (λ, s) based on maximum gap
set.seed(11)
SFclust_pbc_refined_optimal <- SFclust(k = 3, l1bound = l1b_optimal, 
                                       dist.ary = dist_ary_list[[c_opt_index]])
SFclust_pbc_refined_optimal$final.weight

## (2) Find the best (c, s) combination based on 1-SD rule
df_gaps_1sd <- df_gaps %>% 
  mutate(max_gap = max(gap),
         gap_upper_1sd = gap + sdgap) %>% 
  filter(gap_upper_1sd >= max_gap) %>%
  arrange(gap)
df_gaps_1sd

c_optimal_1sd <- df_gaps_1sd$const_c[which.min(df_gaps_1sd$gap)]
l1b_optimal_1sd <- df_gaps_1sd$l1b[which.min(df_gaps_1sd$gap)]
c_opt_1sd_index <- which(const_c == c_optimal_1sd)

## Apply SFKmL with optimal (λ, s) based on 1-SD rule
set.seed(11)
SFclust_pbc_refined_1sd <- SFclust(k = 3, l1bound = l1b_optimal_1sd, 
                                   dist.ary = dist_ary_list[[c_opt_1sd_index]])
SFclust_pbc_refined_1sd$final.weight

## Plotting: gap vs. const_c and variable weights vs. l1bound
weight_pbc_refined <- list()
for (j in 1:length(l1b_pbc_refined)) {
  dist_ary <- dist_ary_list[[c_opt_index]]
  set.seed(11)
  res_sfclust <- SFclust(k = 3, l1bound = l1b_pbc_refined[j], dist.ary = dist_ary)
  weight_pbc_refined[[j]] <- res_sfclust$final.weight
}
df_weight <- do.call(rbind, weight_pbc_refined) %>% data.frame

rownames(df_weight) <- l1b_pbc_refined
colnames(df_weight) <- colnames(df_pbc_refined)[-(1:2)]

df_weight_long <- data.frame(l1b_pbc_refined, df_weight) %>% 
  pivot_longer(
    cols = -l1b_pbc_refined,
    names_to = "variable",
    values_to = "value"
  )
df_weight_long$variable <- factor(
  df_weight_long$variable,
  levels = c("log_bili", "log_chol", "albumin", "log_alk.phos", "log_ast", "log_platelet", "protime")
)

p_gap_pbc <- df_gaps %>% ggplot(aes(x = const_c, y = gap)) +
  geom_point() +
  geom_line(data = df_gaps_max_by_lam, aes(x = const_c, y = gap_max)) +
  geom_vline(xintercept = c_optimal, size = 0.7, linetype = "dashed") +
  labs(x = "const_c (c)", y = "Gap") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_line(color = "gray90"))

p_wght_pbc <- ggplot(df_weight_long, aes(x = l1b_pbc_refined, y = value, color = variable)) +
  geom_line() +
  geom_vline(xintercept = l1b_pbc_refined, linetype = "dashed", color = "gray90") +
  geom_vline(xintercept = l1b_optimal, size = 0.7, linetype = "dashed", color = "gray10") +
  geom_vline(xintercept = l1b_optimal_1sd, size = 0.6, linetype = "dashed", color = "red", alpha = 0.8) +
  scale_x_continuous(breaks = sort(c(1, 1.5, 2, 2.5))) +
  labs(x = "l1bound (s)",
       y = "Weight",
       color = "Variable") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggarrange(p_gap_pbc, p_wght_pbc, nrow = 1, ncol = 2, align = 'hv')