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
library(tictoc)

rm(list = ls())

# 1. Load Data ----
load("Data/df_clean.Rdata")

age_low <- 20
age_high <- 64
fup <- 5
sexes <- list(male = "Male", female = "Female", 
              all = c("Male", "Female"))

save(age_low, age_high, fup, sexes,
     file = "Data/model_parameters.Rdata")

df_cohort <- df_clean %>%
  arrange(birth) %>%
  mutate(cohort = birth - birth %% 10) %>%
  filter(age >= age_low, age <= age_high) %>%
  group_by(cohort) %>%
  filter(max(age) >= age_low + !!fup - 1, 
         min(age) <= age_high + !!fup + 1) %>%
  ungroup() %>%
  filter(cohort > 1920)

# Sample Size
df_count <- df_clean %>%
  arrange(birth) %>%
  filter(age >= age_low, age <= age_high) %>%
  group_by(birth) %>%
  filter(max(age) >= age_low + !!fup - 1, 
         min(age) <= age_high + !!fup + 1) %>%
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
  
  map_dfr(1:200, ~ sample_frac(df_mod, replace = TRUE) %>%
            run_mod(),
          .id = "boot") %>%
    mutate(boot = as.integer(boot)) %>%
    bind_rows(run_mod(df_mod) %>%
                mutate(boot = 0L), .)
}

set.seed(1)
tic()
plan(multisession, workers = 4)
res_gamlss <- df_cohort %>% 
  distinct(cohort) %>% 
  expand_grid(sex = names(sexes)) %>%
  sample_frac() %>%
  mutate(res = future_map2(cohort, sex, get_linpred, .progress = TRUE)) %>%
  arrange(cohort, sex) %>%
  unnest(res)
future:::ClusterRegistry("stop")
toc()

save(res_gamlss, file = "Data/gamlss_results.Rdata")


# 3. Predicted Values ----
load("Data/gamlss_results.Rdata")

get_ci <- function(x){
  quantile(x, probs = c(.5, .025, .975), 
           na.rm = TRUE) %>%
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

# Plot 1
res_linpred %>%
  mutate(sex = str_to_title(sex)) %>%
  ggplot() +
  aes(x = centile, y = beta,
      ymin = lci, ymax = uci, 
      color = age, fill = age, group = age) +
  facet_grid(sex ~ cohort, switch = "y") +
  geom_line() +
  scale_color_viridis_c() +
  scale_fill_viridis_c() +   
  guides(color = guide_colorbar(title.position = 'top', 
                                title.hjust = .5,                                
                                barwidth = unit(20, 'lines'), 
                                barheight = unit(.5, 'lines'))) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        strip.background.y = element_blank()) +
  labs(x = "Centile", y = NULL, 
       color = "Age", fill = "Age")
ggsave("Images/gamlss_1.png", width = 29.7, height = 21, units = "cm")


# Plot 2
obese_cat <- tribble(
  ~centile, ~age, ~beta, ~label,
  "10th", -Inf, 18.5, "Normal",
  "10th", -Inf, 25, "Overweight",
  "10th", -Inf, 30, "Obese")

res_linpred %>%
  filter(centile %in% c(10, 25, 50, 75, 90)) %>%
  mutate(sex = str_to_title(sex),
         centile = glue("{centile}th")) %>%
  ggplot() +
  aes(x = age, y = beta) +
  facet_grid(sex ~ centile, switch = "y") +
  geom_hline(yintercept = c(18.5, 25, 30), linetype = "dashed", color = "grey60") + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = cohort), color = NA, alpha = 0.2) +
  geom_line(aes(color = cohort)) + 
  geom_text(data = obese_cat, aes(label = label),
            vjust = -0.5, hjust = -0.1, color = "grey50") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        strip.background.y = element_blank()) +
  labs(x = "Age", y = NULL, color = "Cohort", fill = "Cohort")
ggsave("Images/gamlss_2.png", width = 29.7, height = 21, units = "cm")

# Plot 3a
obese_cat <- tribble(
  ~centile, ~ age_f, ~age, ~beta, ~label,
  -Inf, "Age 25", -Inf, 18.5, "Normal",
  -Inf, "Age 25", -Inf, 25, "Overweight",
  -Inf, "Age 25", -Inf, 30, "Obese")

res_linpred %>%
  filter(age %in% c(25, 35, 45, 55),
         sex == "all") %>%
  mutate(age_f = glue("Age {age}")) %>%
  ggplot() +
  aes(x = centile, y = beta) +
  facet_grid( ~ age_f, switch = "y") +
  geom_hline(yintercept = c(18.5, 25, 30), linetype = "dashed", color = "grey60") + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = cohort), color = NA, alpha = 0.2) +
  geom_line(aes(color = cohort)) +
  geom_text(data = obese_cat, aes(label = label),
            vjust = -0.5, hjust = -0.1, color = "grey50") +
  coord_cartesian(ylim = c(NA, 45)) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        strip.background.y = element_blank()) +
  labs(x = "Centile", y = NULL, color = "Cohort", fill=  "Cohort")
ggsave("Images/gamlss_3a.png", width = 29.7, height = 21, units = "cm")


# Plot 3
res_linpred %>%
  filter(age %in% c(25, 35, 45, 55),
         sex != "all") %>%
  mutate(age_f = glue("Age {age}"),
         sex = str_to_title(sex)) %>%
  ggplot() +
  aes(x = centile, y = beta) +
  facet_grid(sex ~ age_f, switch = "y") +
  geom_hline(yintercept = c(18.5, 25, 30), linetype = "dashed", color = "grey60") + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = cohort), color = NA, alpha = 0.2) +
  geom_line(aes(color = cohort)) +
  geom_text(data = obese_cat, aes(label = label),
            vjust = -0.5, hjust = -0.1, color = "grey50") +
  coord_cartesian(ylim = c(NA, 50)) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        strip.background.y = element_blank()) +
  labs(x = "Centile", y = NULL, color = "Cohort", fill=  "Cohort")
ggsave("Images/gamlss_3b.png", width = 29.7, height = 21, units = "cm")

# Plot 4
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
ggsave("Images/gamlss_4.png", width = 29.7, height = 21, units = "cm")


# Plot 5
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
ggsave("Images/gamlss_5.png", width = 21, height = 16, units = "cm")
