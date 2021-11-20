library(tidyverse)
library(glue)
library(ggridges)
library(ggrepel)
library(Hmisc)
library(scales)
library(splines)
library(broom.mixed)
library(magrittr)
library(gamlss)
library(furrr)

rm(list = ls())

# 1. Load Data ----
load("Data/df_clean.Rdata")

age_low <- 20
age_high <- 64
fup <- 5

df_cohort <- df_clean %>%
  arrange(birth) %>%
  mutate(cohort = birth - birth %% 10) %>%
  filter(age >= age_low, age <= age_high) %>%
  group_by(cohort) %>%
  filter(max(age) >= age_low + !!fup, 
         min(age) <= age_high + !!fup) %>%
  ungroup() %>%
  filter(cohort > 1920)

sexes <- list(male = "Male", female = "Female", 
              all = c("Male", "Female"))

make_df <- function(cohort, sex){
  df_cohort %>%
    filter(cohort == !!cohort,
           sex %in% sexes[[!!sex]]) %>%
    group_by(year) %>%
    mutate(wt_int = wt_int*n()/sum(wt_int)) %>%
    ungroup()
}

# Sample Size
df_count <- df_clean %>%
  arrange(birth) %>%
  filter(age >= age_low, age <= age_high) %>%
  group_by(birth) %>%
  filter(max(age) >= age_low + !!fup, 
         min(age) <= age_high + !!fup) %>%
  ungroup() %>% 
  count(year, sort = TRUE)

summary(df_count$n)
df_count %>% slice(c(1, n()))


# 2. Run Models ----
get_linpred <- function(cohort, sex){
  data <- df_cohort %>%
    filter(cohort == !!cohort,
           sex %in% sexes[[!!sex]]) %>%
    group_by(year) %>%
    mutate(wt_int = wt_int*n()/sum(wt_int)) %>%
    ungroup()
  
  df_s <- data %>%
    distinct(age) %>%
    arrange(age) %>%
    mutate(ns(age, 2) %>%
             as_tibble() %>%
             mutate(across(everything(), as.double)) %>%
             rename_with(~ glue("age_ns{.x}")))
  
  df_mod <- data %>%
    left_join(df_s, by = "age")
  
  run_mod <- function(df){
    mod <- gamlss(bmi ~ age_ns1 + age_ns2,
                  sigma.formula = ~ age_ns1 + age_ns2,
                  nu.formula = ~ age_ns1 + age_ns2,
                  family = BCCG, 
                  data = df, weights = df$wt_int,
                  trace = FALSE)
    
    coefAll(mod) %>%
      map_dfr(as_tibble_row, .id = "param") %>%
      pivot_longer(-param, names_to = "term", values_to = "beta") %>%
      left_join(df_s %>%
                  mutate(`(Intercept)` = 1) %>%
                  pivot_longer(-age, names_to = "term", values_to = "coef"),
                by = "term") %>%
      group_by(param, age) %>%
      summarise(estimate = sum(beta * coef),
                .groups = "drop") %>%
      pivot_wider(names_from = param, values_from = estimate)
  }
  
  map_dfr(1:10, ~ sample_frac(df_mod, replace = TRUE) %>%
            run_mod(),
          .id = "boot") %>%
    mutate(boot = as.integer(boot)) %>%
    bind_rows(run_mod(df_mod) %>%
                mutate(boot = 0L), .)
}

res_gamlss <- df_cohort %>% 
  distinct(cohort) %>% 
  expand_grid(sex = names(sexes)) %>%
  mutate(res = map2(cohort, sex, get_linpred)) %>%
  unnest(res)

save(res_gamlss, file = "Data/gamlss_results.Rdata")


# 3. Predicted Values ----
get_ci <- function(x){
  quantile(x, probs = c(.5, .025, .975)) %>%
    set_names(c("beta", "lci", "uci")) %>%
    as_tibble_row()
}

# Parameters
res_param <- res_gamlss %>%
  mutate(sigma = exp(sigma)) %>%
  pivot_longer(c(mu,  sigma, nu), 
               names_to = "parameter", 
               values_to = "value") %>%
  mutate(cohort = glue("{cohort}-{cohort+9}") %>% ordered(),
         boot = ifelse(boot == 0, 0, 1)) %>%
  arrange(cohort, sex, age, parameter, boot) %>%
  group_by(cohort, sex, age, parameter, boot) %>%
  summarise(get_ci(value), .groups = "drop_last") %>%
  summarise(beta = nth(beta, 1), lci = nth(lci, 2), 
            uci = nth(uci, 2), .groups = "drop")


# Linear Predictions
res_linpred <- res_gamlss %>%
  uncount(99, .id = "centile") %>%
  mutate(cohort = glue("{cohort}-{cohort+9}") %>% ordered(),
         tau = centile/100,
         bmi = mu * (1 + nu * exp(sigma) * qnorm(tau))^(1/nu),
         boot = ifelse(boot == 0, 0, 1)) %>%
  arrange(cohort, sex, age, centile, boot) %>%
  group_by(cohort, sex, age, centile, boot) %>%
  summarise(get_ci(bmi), .groups = "drop_last") %>%
  summarise(beta = nth(beta, 1), lci = nth(lci, 2), 
            uci = nth(uci, 2), .groups = "drop")

# Difference Across Centiles
res_diff <- res_gamlss %>%
  uncount(99, .id = "centile") %>%
  filter(centile %in% c(25, 50, 75, 90)) %>%
  mutate(cohort = glue("{cohort}-{cohort+9}") %>% ordered(),
         tau = centile/100,
         bmi = mu * (1 + nu * exp(sigma) * qnorm(tau))^(1/nu),
         comp = ifelse(centile %in% c(25, 75), "75th vs 25th", "90th vs 10th")) %>%
  arrange(cohort, sex, comp, age, boot, comp, centile) %>%
  group_by(cohort, sex, age, comp, boot) %>%
  summarise(bmi = nth(bmi, 2) - nth(bmi, 1),
            .groups = "drop") %>%
  mutate(boot = ifelse(boot == 0, 0, 1)) %>%
  group_by(cohort, sex, age, comp, boot) %>%
  summarise(get_ci(bmi), .groups = "drop_last") %>%
  summarise(beta = nth(beta, 1), lci = nth(lci, 2), 
            uci = nth(uci, 2), .groups = "drop")



save(res_param, res_linpred, res_diff,
     file = "Data/gamlss_predictions.Rdata")


# 3. Plots ----
load("Data/gamlss_predictions.Rdata")


plot_1 <- function(sex){
  res_linpred %>%
    filter(sex == !!sex) %>%
    ggplot() +
    aes(x = centile, y = beta,
        ymin = lci, ymax = uci, 
        color = age, fill = age, group = age) +
    facet_wrap(~ cohort) +
    # geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    scale_color_viridis_c() +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "Centile", y = "Predicted BMI", 
         color = "Age", fill = "Age")
}

plot_2 <- function(sex){
  res_linpred %>%
    filter(sex == !!sex,
           centile %in% c(10, 25, 50, 75, 90)) %>%
    ggplot() +
    aes(x = age, y = beta,
        ymin = lci, ymax = uci,
        color = cohort, fill = cohort) +
    facet_wrap(~ centile, nrow = 1) +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "Age", y = "Predicted BMI",
         color = "Cohort", fill = "Cohort")
}

plot_3 <- function(sex){
  res_linpred %>%
    filter(age %% 5 == 0,
           sex == !!sex) %>%
    mutate(age_f = glue("Age {age}")) %>%
    ggplot() +
    aes(x = centile, y = beta,
        ymin = lci, ymax = uci, 
        color = cohort, fill = cohort) +
    facet_wrap(~ age_f) +
    geom_ribbon(color = NA, alpha = 0.2) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "Centile", y = "Predicted BMI", 
         color = "Cohort", fill=  "Cohort")
}


param_dict <- c(mu = "Median", sigma = "CoV", nu = "Skewness")

res_param %>%
  mutate(sex = str_to_title(sex),
         param_clean = factor(param_dict[parameter], param_dict)) %>%
  ggplot() +
  aes(x = age, y = beta, ymin = lci, ymax = uci,
      color = cohort, fill = cohort) +
  facet_grid(param_clean ~ sex, scales = "free_y", switch = "y") +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        strip.background.y = element_blank()) +
  labs(x = "Age", y = NULL, color = "Cohort", fill=  "Cohort")


res_diff %>%
  mutate(sex = str_to_title(sex)) %>%
  ggplot() +
  aes(x = age, y = beta, ymin = lci, ymax = uci,
      color = cohort, fill = cohort) +
  facet_grid(sex ~ comp, scales = "free_y", switch = "y") +
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        strip.background.y = element_blank()) +
  labs(x = "Age", y = NULL, color = "Cohort", fill=  "Cohort")
  
