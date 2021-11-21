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
library(furrr)
library(tictoc)

rm(list = ls())

# 1. Load Data ----
load("Data/df_clean.Rdata")
load("Data/model_parameters.Rdata")

df_clean <- df_clean %>% 
  mutate(obese = ifelse(bmi >= 30, 1, 0)) %>%
  filter(between(bmi, 13, 70),
         age >= !!age_low, 
         age <= !!age_high)


# 2. Run lmer ----
make_df <- function(sex, cohort_y, last_year){
  df_clean %>%
    filter(year <= !!last_year,
           sex %in% sexes[[!!sex]]) %>%
    mutate(birth = birth - (birth %% !!cohort_y),
           age_s = rescale(age)) %>%
    group_by(birth) %>%
    filter(min(age) <= !!age_high - !!fup, 
           max(age) >= !!age_low + !!fup) %>%
    group_by(year) %>%
    mutate(wt_int = wt_int*n()/sum(wt_int)) %>%
    ungroup() # %>%
  # group_by(birth) %>%
  # mutate(wt_int = wt_int*n()/sum(wt_int)) %>% # Carle (2009) Weighting Scaling
  # ungroup()
}

run_mod <- function(sex, cohort_y, last_year, outcome){
  
  df_mod <- make_df(sex, cohort_y, last_year)
  
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

set.seed(1)
tic()
plan(multisession, workers = 4)
res_lmer <- expand_grid(sex = names(sexes),
                        cohort_y = 1:4,
                        last_year = c(2014, 2019),
                        outcome = c("bmi", "obese")) %>%
  sample_frac() %>%
  mutate(res = future_pmap(list(sex, cohort_y, 
                                last_year, outcome), 
                           run_mod, .progress = TRUE)) %>%
  arrange(sex, cohort_y, last_year, outcome) %>%
  unnest(res) %>%
  mutate(singular = map_lgl(mod, isSingular),
         cor = map_dbl(mod, ~ tidy(.x) %>%
                         filter(str_detect(term, "^cor")) %>%
                         pull(4)),
         sex_clean = str_to_title(sex),
         outcome_clean = ifelse(outcome == "bmi", 
                                "Predicted BMI",
                                "Probability of Obesity"))
future:::ClusterRegistry("stop")
toc()


# 3. Add Predicted Values ----
make_pred <- function(mod, df_age, df_birth, outcome){
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
  
  list(df_pred)
}

res_lmer <- res_lmer %>%
  rowwise() %>%
  mutate(pred = make_pred(mod, df_age, df_birth, outcome)) %>%
  ungroup()

make_re <- function(mod){
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
    select(birth, term_clean, ranef, lci, uci)
}

res_lmer <- res_lmer %>%
  mutate(re = map(mod, make_re))


res_lmer2 <- res_lmer %>%
  select(-mod)
save(res_lmer2, file = "Data/lmer_results.Rdata")

# 4. Make Plots ----
# Predictions
plot_pred <- function(cohort_y = 1, last_year = 2019){
  res_lmer %>%
    filter(cohort_y == !!cohort_y, last_year == !!last_year) %>%
    select(sex_clean, outcome_clean, pred) %>%
    unnest(pred) %>%
    ggplot() +
    aes(x = age, y = estimate, color = birth, group = birth) +
    facet_grid(outcome_clean ~ sex_clean, scales = "free_y", switch = "y") +
    geom_line() +
    theme_bw() +
    scale_color_viridis_c() +
    labs(x = "Age", y = "Predicted BMI", color = "Birth Year") +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0),
          strip.background.y = element_blank()) +   
    guides(color = guide_colorbar(title.position = 'top', 
                                  title.hjust = .5,                                
                                  barwidth = unit(20, 'lines'), 
                                  barheight = unit(.5, 'lines'))) +
    labs(x = "Centile", y = NULL, color = "Cohort")
}


# Random Effects
plot_re <- function(outcome, sexes = "all", cohort_y = 1, last_year = 2019){
  res_lmer %>%
    filter(cohort_y == !!cohort_y, 
           last_year == !!last_year,
           sex %in% !!sexes, 
           outcome == !!outcome) %>%
    unnest(re) %>%
    ggplot() +
    aes(x = birth, y = ranef, 
        ymin = lci, ymax = uci,
        color = sex_clean) +
    facet_wrap(~ term_clean, scales = "free", ncol = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "Cohort", y = NULL, color = NULL)
}

