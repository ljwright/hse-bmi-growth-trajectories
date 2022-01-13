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
load("Data/hse_clean.Rdata")
load("Data/model_parameters.Rdata")

df <- df_clean %>% 
  select(-edu) %>%
  mutate(obese = ifelse(bmi >= 30, 1, 0)) %>%
  filter(age >= !!age_low, 
         age <= !!age_high)

rm(df_clean)


# 2. Model Objects ----
make_df <- function(sex, cohort_y, last_year){
  df %>%
    filter(year <= !!last_year,
           sex %in% sexes[[!!sex]]) %>%
    mutate(birth = birth - (birth %% !!cohort_y),
           age_s = rescale(age)) %>%
    group_by(birth) %>%
    filter(min(age) <= !!age_high - !!fup + 1, 
           max(age) >= !!age_low + !!fup - 1) %>%
    group_by(year) %>%
    mutate(wt_int = wt_int*n()/sum(wt_int)) %>%
    ungroup() # %>%
  # group_by(birth) %>%
  # mutate(wt_int = wt_int*n()/sum(wt_int)) %>% # Carle (2009) Weighting Scaling
  # ungroup()
}

mod_specs <- expand_grid(sex = names(sexes),
                         cohort_y = c(1, 5),
                         last_year = c(2014, 2019),
                         outcome = c("bmi", "obese"),
                         ran_pars = c("1 + age_s", "-1 + age_s")) %>%
  filter(cohort_y == 1 | last_year == 2019) %>%
  mutate(spec_id = row_number(), .before = 1)

run_mod <- function(spec_id){
  spec <- slice(mod_specs, !!spec_id)
  
  df_mod <- make_df(spec$sex, spec$cohort_y, spec$last_year)
  
  mod_form <- glue("{spec$outcome} ~ 1 + age_s + ({spec$ran_pars} | birth)") %>%
    as.formula()
  
  if (spec$outcome == "bmi"){
    mod <- lmer(mod_form, df_mod, REML = TRUE, weights = wt_int)
  } else{
    mod <- glmer(mod_form, df_mod, binomial, weights = wt_int)
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
  
  tibble(mod = list(mod), 
         df_age = list(df_age), 
         df_birth = list(df_birth))
}


# 3. Run Models ----
set.seed(1)
tic()
plan(multisession, workers = 4)
res_lmer <- mod_specs %>%
  sample_frac() %>%
  filter(spec_id == 30) %>%
  mutate(res = future_map(spec_id, run_mod, .progress = TRUE)) %>%
  arrange(spec_id) %>%
  unnest(res) %>%
  mutate(singular = map_lgl(mod, isSingular),
         tidy = map(mod, tidy),
         glance = map(mod, glance),
         cor = map(tidy, ~ .x %>%
                     filter(str_detect(term, "^cor")) %>%
                     pull(4)),
         cor = map_dbl(cor, ~ ifelse(is.null(.x), NA, .x)),
         sex_clean = str_to_title(sex),
         outcome_clean = ifelse(outcome == "bmi", 
                                "Predicted BMI",
                                "Probability of\nObesity"))
future:::ClusterRegistry("stop")
toc()


# 4. Add Predicted Values and Random Effects ----
make_pred <- function(mod, df_age, df_birth, outcome){
  # re <- ranef(mod)[[1]]
  # if (str_detect(dimnames(re)[[2]], "Intercept") == FALSE){
  #   re <- cbind(`(Intercept)` = 0, re)
  # }
  # 
  # df_term <- left_join(
  #   fixef(mod) %>%
  #     enframe(name = "term", value = "fixef"),
  #   re %>%
  #     as_tibble(rownames = "birth") %>%
  #     pivot_longer(-birth, names_to = "term", values_to = "ranef"),
  #   by = "term"
  # )
  # 
  # df_pred <- df_term %>%
  #   mutate(coef = fixef + ranef) %>%
  #   uncount(!!age_high - !!age_low + 1, .id = "age") %>%
  #   mutate(age = age + !!age_low + 1,
  #          birth = as.double(birth)) %>%
  #   left_join(df_age, by = c("term", "age")) %>%
  #   mutate(beta = coef*value) %>%
  #   group_by(birth, age) %>%
  #   summarise(estimate = sum(beta),
  #             .groups = "drop") %>%
  #   left_join(df_birth, by = "birth") %>%
  #   filter(age >= min_age, age <= max_age)
  
  df_pred <- fixef(mod) %>%
    enframe(name = "term", value = "fixef") %>%
    expand_grid(df_birth) %>%
    full_join(ranef(mod)[[1]] %>%
                as_tibble(rownames = "birth") %>%
                mutate(birth = as.double(birth)) %>%
                pivot_longer(-birth, names_to = "term", values_to = "ranef"), 
              by = c("birth", "term")) %>%
    mutate(across(c(fixef, ranef),
                  ~ ifelse(is.na(.x), 0, .x))) %>%
    mutate(coef = fixef + ranef) %>%
    uncount(max_age - min_age + 1, .id = "age") %>%
    mutate(age = age + min_age - 1) %>%
    left_join(df_age, by = c("term", "age")) %>%
    mutate(beta = coef*value) %>%
    group_by(birth, age) %>%
    summarise(estimate = sum(beta),
              .groups = "drop")
  
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
    map_dfr(~ as.matrix(.x) %>%
              diag() %>% 
              sqrt() %>%
              set_names(names(ran_effects)[-1]) %>%
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


# 5. Make Plots ----
load("Data/lmer_results.Rdata")

# Predictions
plot_pred <- function(cohort_y = 1, last_year = 2019, save_p = FALSE){
  p <- res_lmer2 %>%
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
    labs(x = "Age", y = NULL, color = "Cohort")
  
  if (save_p == TRUE){
    glue("Images/lmer_{cohort_y}_{last_year}.png") %>%
      ggsave(plot = p, width = 29.7, height = 21, units = "cm")
  }
  
  return(p)
}

res_lmer2 %>%
  distinct(cohort_y, last_year) %$%
  map2(cohort_y, last_year, plot_pred, TRUE)


# Random Effects
plot_re <- function(outcome, sexes = "all", cohort_y = 1, 
                    last_year = 2019, save_p = FALSE){
  
  if (outcome == "bmi"){
    y_lab <- bquote('Difference in BMI (kg/'~m^2*')')
  } else{
    y_lab <- "Difference in (Log-Odds) Obesity"
  }
  
  p <- res_lmer2 %>%
    filter(cohort_y == !!cohort_y, 
           last_year == !!last_year,
           sex %in% !!sexes, 
           outcome == !!outcome) %>%
    unnest(re) %>%
    ggplot() +
    aes(x = birth, y = ranef, 
        ymin = lci, ymax = uci,
        color = sex_clean, shape = sex_clean) +
    facet_wrap(~ term_clean, scales = "free", ncol = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(position = position_dodge(0.7), size = 0.3) +
    scale_color_brewer(palette = "Set1") +
    scale_shape_manual(values = 16:15) +
    scale_x_continuous(breaks = seq(from = 1930, to = 1995, by = 5)) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(x = "Cohort", y = y_lab, color = NULL, shape = NULL)
  
  if (identical(sexes, "all")){
    p <- p + theme(legend.position = "none")
  }
  
  sexes_stub <- ifelse(identical(sexes, "all"), "all", "sex")
  
  if (save_p == TRUE){
    glue("Images/re_{outcome}_{sexes_stub}_{cohort_y}_{last_year}.png") %>%
      ggsave(plot = p, width = 21, height = 16, units = "cm")
  }
  
  return(p)
}

res_lmer2 %>%
  distinct(outcome, cohort_y, last_year) %>%
  expand_grid(sexes = list("all", c("male", "female"))) %$%
  pmap(list(outcome, sexes, cohort_y, last_year), plot_re, TRUE)
