library(tidyverse)
library(glue)
library(haven)
library(labelled)
library(ggridges)
library(ggrepel)
library(Hmisc)
library(lme4)
library(scales)
library(splines)
library(quantreg)
library(broom.mixed)
library(magrittr)

rm(list = ls())

# 1. Load Data ----
load("Data/df_clean.Rdata")

df_clean <- df_clean %>% 
  mutate(obese = ifelse(bmi >= 30, 1, 0))


# 2. Run lmer ----
make_df <- function(age_low, age_high, fup, cohort_y, last_year){
  df_clean %>%
    filter(between(bmi, 13, 70),
           age >= !!age_low, 
           age <= !!age_high,
           year <= !!last_year) %>%
    mutate(birth = birth - (birth %% !!cohort_y),
           age_s = rescale(age)) %>%
    group_by(birth) %>%
    filter(min(age) <= !!age_high - !!fup, 
           max(age) >= !!age_low + !!fup) %>%
    group_by(year) %>%
    mutate(wt_int = wt_int*n()/sum(wt_int)) %>%
    ungroup() %>%
    group_by(birth) %>%
    mutate(wt_int = wt_int*n()/sum(wt_int)) %>% # Carle (2009) Weighting Scaling
    ungroup()
}

run_mod <- function(age_low, age_high, fup, cohort_y, last_year, outcome){
  
  df_mod <- make_df(age_low, age_high, fup,
                    cohort_y, last_year)
  
  if (outcome == "bmi"){
    mod <- lmer(bmi ~ 1 + age_s + (1 + age_s | birth), 
                df_mod, REML = TRUE, weights = wt_int)
  } else{
    mod <- glmer(obese ~ 1 + age_s + (1 + age_s | birth), 
                 df_mod, binomial, weights = wt_int)
  }
  
  df_birth <- df_mod %>%
    group_by(birth) %>%
    summarise(min_age = min(age),
              max_age = max(age))
  
  df_age <- df_mod %>%
    distinct(age, age_s) %>%
    arrange(age) %>%
    mutate(`(Intercept)` = 1) %>%
    pivot_longer(2:last_col(), 
                 names_to = "term", 
                 values_to = "value")
  
  tibble(mod = list(mod), df_age = list(df_age), df_birth = list(df_birth))
}

res_lmer <- expand_grid(age_low = c(20, 44),
                        age_high = c(44, 64),
                        fup = 4,
                        cohort_y = 1:4,
                        last_year = c(2014, 2019),
                        outcome = c("bmi", "obese")) %>%
  filter(age_low < age_high, outcome == "bmi") %>%
  rowwise() %>%
  mutate(res = run_mod(age_low, age_high, fup, cohort_y, last_year, outcome)) %>%
  ungroup() %>% 
  tidyr::unpack(res) %>%
  mutate(singular = map_lgl(mod, isSingular),
         cor = map_dbl(mod, ~ tidy(.x) %>%
                         filter(str_detect(term, "^cor")) %>%
                         pull(4)))

# 4. Make Plots ----
plot_pred <- function(mod, df_age, df_birth, outcome){
  df_term <- left_join(
    fixef(mod) %>%
      enframe(name = "term", value = "fixef"),
    ranef(mod)[[1]] %>%
      as_tibble(rownames = "birth") %>%
      pivot_longer(-birth, names_to = "term", values_to = "ranef"),
    by = "term"
  )
  
  df_pred <- df_term %>%
    mutate(coef = fixef + ranef) %>%
    uncount(64 - 20 + 1, .id = "age") %>%
    mutate(age = age + 19,
           birth = as.double(birth)) %>%
    left_join(df_age, by = c("term", "age")) %>%
    mutate(beta = coef*value) %>%
    group_by(birth, age) %>%
    summarise(estimate = sum(beta),
              .groups = "drop") %>%
    left_join(df_birth, by = "birth") %>%
    filter(age >= min_age, age <= max_age)
  
  if (outcome == "obese"){
    df_pred <- df_pred %>%
      mutate(estimate = exp(estimate)/(1+(exp(estimate))))
  }
  
  ggplot(df_pred) +
    aes(x = age, y = estimate, color = birth, group = birth) +
    geom_line() +
    theme_bw() +
    scale_color_viridis_c() +
    labs(x = "Age", y = "Predicted BMI", color = "Birth Year")
}

plot_re <- function(mod){
  ran_effects <- ranef(mod)$birth %>%
    as_tibble(rownames = "birth")
  
  attributes(ran_effects)$postVar %>%
    plyr::alply(.margins = 3) %>%
    map_dfr(~ diag(.x) %>% 
              sqrt() %>%
              set_names(c(names(ran_effects)[-1])) %>%
              as_tibble_row(),
            .id = "row") %>%
    pivot_longer(-row, names_to = "term", values_to = "se") %>%
    mutate(row = as.integer(row)) %>%
    left_join(ran_effects %>%
                mutate(row = row_number()) %>%
                pivot_longer(-c(birth, row), names_to = "term", values_to = "ranef"),
              by = c("row", "term")) %>%
    mutate(lci = qnorm(.025, ranef, se),
           uci = qnorm(.975, ranef, se),
           birth = as.integer(birth),
           term_clean = ifelse(term == "(Intercept)", 
                               "Random Intercept", 
                               "Random Slope")) %>%
    ggplot() +
    aes(x = birth, y = ranef, ymin = lci, ymax = uci) +
    facet_wrap(~ term_clean, scales = "free", ncol = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange() +
    theme_bw() +
    labs(x = "Cohort", y = NULL)
}

res_lmer %>%
  slice(10) %$%
  plot_pred(mod[[1]], df_age[[1]], df_birth[[1]], outcome)

plot_re(res_lmer$mod[[10]])
