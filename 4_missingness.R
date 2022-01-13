library(tidyverse)
library(glue)
library(broom)
library(broom.mixed)
library(summarytools)
library(mice)
library(miceadds)
library(splines)
library(ggeffects)
library(scales)
library(lme4)
library(Hmisc)

rm(list = ls())

# 1. Load Data ----
load("Data/hse_clean.Rdata")
load("Data/model_parameters.Rdata")

df <- df_clean %>% 
  filter(age >= !!age_low, 
         age <= !!age_high) %>%
  group_by(birth) %>%
  filter(min(age) <= !!age_high - !!fup + 1, 
         max(age) >= !!age_low + !!fup - 1) %>%
  ungroup() %>%
  mutate(miss_bmi = ifelse(is.na(bmi), 1, 0))

rm(df_clean)

# 2. Missingness Prevalence ----
df %>%
  descr() %>%
  tb() %>%
  select(variable, n.valid, pct.valid)

df %>%
  select(-edu) %>%
  add_count(name = "n_all") %>%
  drop_na(sex, bmi) %>%
  count(n_all, name = "n_obs") %>%
  mutate(pct = n_obs*100/n_all)

# 2. Missingness Patterns ----
# Year
lm(miss_bmi ~ factor(year), df) %>%
  tidy(conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(year = str_sub(term, -4) %>% as.integer()) %>%
  ggplot() +
  aes(x = year, y = estimate, ymin = conf.low, ymax = conf.high) +
  geom_hline(yintercept = 0) +
  geom_pointrange() +
  scale_x_continuous(breaks = seq(from = 1991, to = 2019, by = 4)) +
  theme_minimal() +
  labs(x = "Survey Year", y = "Difference in Missingness Rate")

# Age
lm(miss_bmi ~ ns(age, 3), df) %>%
  ggpredict(terms = "age") %>%
  as_tibble() %>%
  ggplot() +
  aes(x = x, y = predicted, 
      ymin = conf.low, ymax = conf.high) +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line() +
  theme_minimal() +
  labs(x = "Age",
       y = "Missingness Rate")

# Survey Year and Cohort
lm(miss_bmi ~ factor(birth) + factor(year), df) %>%
  tidy(conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(year = str_sub(term, -4) %>% as.integer(),
         type = ifelse(str_detect(term, "birth"), "Cohort", "Survey Year")) %>%
  ggplot() +
  aes(x = year, y = estimate, ymin = conf.low, ymax = conf.high) +
  facet_wrap(~ type, scales = "free") +
  geom_hline(yintercept = 0) +
  geom_pointrange() +
  theme_bw() +
  labs(x = "Survey Year", y = "Difference in Missingness Rate")


# 3. Multiple Imputation ----
df_mice <- df %>%
  select(bmi, age, birth, edu, wt_int) %>%
  mutate(age_s = rescale(age),
         obs = !is.na(bmi)) %>%
  select(bmi, age_s, age, birth, wt_int, obs)

set.seed(1e5)
fit <- lmer(bmi ~ age_s + (1 + age_s | birth), df_mice, weights = wt_int)

get_imp <- function(imp){
  set.seed(imp)
  
  fe_draw <- miceadds:::mice_multilevel_draw_rnorm1(mu = fixef(fit), Sigma = vcov(fit))
  
  re <- ranef(fit, condVar = TRUE)[[1]] 
  pv <- attr(re, "postVar")
  re <- as.matrix(re)
  re_names <- dimnames(re)[[1]] %>% as.integer()
  
  re_draw <- miceadds:::mice_multilevel_imputation_draw_random_effects(mu = re, Sigma = pv, ridge = 1e-06)
  
  cov_matrix <- cbind(1, df_mice$age_s)
  re_var <- df_mice$birth
  
  pred_fe <- cov_matrix %*% fe_draw %>% as.vector()
  pred_re <- rowSums(cov_matrix * re_draw[match(re_var, re_names), ])
  sigma <- attr(VarCorr(fit), "sc")
  full_draw <- miceadds:::mice_multilevel_imputation_draw_residuals(predicted = pred_fe + pred_re, 
                                                                    sigma = sigma)
  
  ry <- df_mice$obs
  predicted0 <- full_draw[!ry]
  predicted1 <- full_draw[ry]
  
  pred <- miceadds:::mice_multilevel_imputation_pmm5(y = df_mice$bmi, ry = ry, x = 1,
                                                     yhatobs = predicted1, yhatmis = predicted0, 
                                                     donors = 5, noise = 1e+05)
  
  df_mice$bmi[!ry] <- pred
  return(df_mice)
}

df_imp <- map(1:10, get_imp)


# 4. Pattern Mixture Models ----
df_spec <- df_mice %>%
  distinct(birth, age, age_s) %>%
  arrange(birth, age_s)

get_pred <- function(df_mod){
  mod <- lmer(bmi ~ age_s + (1 + age_s | birth), df_mod, weights = wt_int)
  
  fe <- fixef(mod)
  re <- ranef(mod, condVar = TRUE)[[1]] %>% as.matrix()
  re_names <- dimnames(re)[[1]] %>% as.integer()
  
  cov_matrix <- cbind(1, df_spec$age_s)
  re_var <- df_spec$birth
  
  
  pred_fe <- as.vector(cov_matrix %*% fe)
  pred_re <- rowSums(cov_matrix * re[match(re_var, re_names), ])
  
  df_spec %>%
    mutate(pred = pred_fe + pred_re) %>%
    left_join(re %>%
                as_tibble(rownames = "birth") %>%
                mutate(birth = as.integer(birth)) %>%
                rename(re_1 = 2, re_2 = 3),
              by = "birth") %>%
    select(birth, age, pred, re_1, re_2)
}

pool_lmer <- function(imps){
  map_dfr(imps, get_pred) %>%
    group_by(birth, age) %>%
    summarise(across(pred:re_2, mean),
              .groups = "drop")
}

get_pattern <- function(diff){
  map(df_imp, 
      ~ .x %>%
        mutate(bmi = ifelse(obs, bmi, bmi + diff))) %>%
    pool_lmer()
}

wtd.var(df_mice$bmi, df_mice$wt_int) %>%
  sqrt()

res_pattern <- tibble(diff = seq(from = -5, to = 5, by = 1)) %>%
  mutate(res = map(diff, get_pattern)) %>%
  unnest(res)


# 5. Pattern Mixture Plots ----
# Plot 1
res_pattern %>%
  mutate(sign = ifelse(diff < 0, "", "+ "),
         diff_clean = glue("{sign}{diff} kg/m2") %>%
           fct_reorder(diff)) %>%
  ggplot() +
  aes(x = age, y = pred, color = birth, group = birth) +
  facet_wrap(~ diff_clean, ncol = 3) +
  geom_line() +
  scale_color_viridis_c() +
  scale_x_continuous(breaks = c(2:6*10, 64)) +
  theme_bw() +
  labs(x = "Age", y = "Predicted BMI", color = "Cohort") +
  guides(color = guide_colorbar(title.position = 'top', 
                                title.hjust = .5,                                
                                barwidth = unit(20, 'lines'), 
                                barheight = unit(.5, 'lines'))) +
  theme(legend.position = "bottom")
ggsave("Images/mnar_fig1.png", width = 29.7, height = 21, units = "cm")

# Plot 2
res_pattern %>%
  distinct(diff, birth, re_1, re_2) %>%
  pivot_longer(-c(birth, diff)) %>%
  mutate(name = ifelse(name == "re_1",
                       "Random Intercept",
                       "Random Slope")) %>%
  ggplot() +
  aes(x = birth, y = value, color = diff, group = diff) +
  facet_wrap(~ name, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0) +
  geom_line() +
  scale_color_viridis_c() +
  scale_x_continuous(breaks = c(1931, 194:199*10, 1995)) +
  theme_bw() +
  labs(x = "Cohort", y = "Difference in BMI (kg/m2)",
       color = "Imputed Difference") +
  guides(color = guide_colorbar(title.position = 'top', 
                                title.hjust = .5,                                
                                barwidth = unit(20, 'lines'), 
                                barheight = unit(.5, 'lines'))) +
  theme(legend.position = "bottom")
ggsave("Images/mnar_fig2.png", width = 21, height = 16, units = "cm")
