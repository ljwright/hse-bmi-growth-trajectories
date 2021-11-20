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

rm(list = ls())

# 1. Load Data ----
load("Data/df_clean.Rdata")

run_mod <- function(age_low, age_high, fup, ){
  df_s <- df_clean %>%
    filter(age >= age_low, age <= age_high) %>%
    distinct(age) %>%
    arrange(age) %>%
    mutate(age_s = rescale(age),
           # ns(age_s, 2) %>% 
             poly(age_s, 2) %>%
             as_tibble() %>%
             mutate(across(everything(), as.double)) %>%
             rename_with(~ glue("age_ns{.x}")))
  
  df_mod <- df_clean %>%
    filter(between(bmi, 13, 70)) %>%
    right_join(df_s, by = "age") %>%
    group_by(year) %>%
    mutate(wt_int = wt_int*n()/sum(wt_int),
           year_f = factor(year)) %>%
    group_by(birth) %>%
    filter(min(age) <= !!age_high - 5, 
           max(age) >= !!age_low + 5) %>%
    ungroup()
  
  mod <- lmer(bmi ~ age_ns1 + age_ns2 +
                (1 + age_ns1 + age_ns2 | birth), 
              df_mod, REML = TRUE)
  
  list(mod = mod, df_mod = df_mod, df_s = df_s)
}

mod <- run_mod(20, 64)
summary(mod$mod)


# Main
df_s <- df_clean %>%
  filter(age >= 20, age <= 64) %>%
  distinct(age) %>%
  arrange(age) %>%
  mutate(age_s = rescale(age),
         ns(age_s, 2) %>%
           as_tibble() %>%
           mutate(across(everything(), as.double)) %>%
           rename_with(~ glue("age_ns{.x}"))) # %>%
  # mutate(age_ns1 = age_s) %>%
  # select(-age_ns2)

df_mod <- df_clean %>%
  filter(between(bmi, 13, 70)) %>%
  right_join(df_s, by = "age") %>%
  group_by(year) %>%
  mutate(wt_int = wt_int*n()/sum(wt_int),
         year_f = factor(year)) %>%
  group_by(birth) %>%
  filter(min(age) <= 50, max(age) >= 34) %>%
  ungroup()

# write_dta(df_mod, "test.dta")

mod <- lmer(bmi ~ 1 + age_ns1 + age_ns2 +
              (1 + age_ns1 + age_ns2 | birth), 
            df_mod, REML = TRUE)
summary(mod)

mod <- lmer(bmi ~ 1 + age_ns1 +
              (1 + age_ns1 | birth), 
            df_mod, REML = TRUE)
summary(mod)

with(mod@optinfo$derivs,
     solve(Hessian, gradient)) %>%
  abs() %>%
  max()


# Predicted Effects
df_term <- left_join(
  fixef(mod) %>%
    enframe(name = "term", value = "fixef"),
  ranef(mod)[[1]] %>%
    as_tibble(rownames = "birth") %>%
    pivot_longer(-birth, names_to = "term", values_to = "ranef"),
  by = "term"
)

# df_term <- left_join(
#   read_dta("Data/regsave.dta") %>%
#     filter(str_detect(var, "bmi")) %>%
#     mutate(term = str_replace(var, "bmi\\:", "")) %>%
#     select(term, fixef = coef),
#   read_dta("Data/reffects.dta") %>%
#     pivot_longer(-birth, names_to = "term", values_to = "ranef") %>%
#     mutate(term = str_replace(term, "term_", "")),
#   by = "term"
# ) %>%
#   mutate(term = str_replace(term, "_cons", "(Intercept)"))

df_term %>%
  mutate(coef = fixef + ranef) %>%
  uncount(64 - 20 + 1, .id = "age") %>%
  mutate(age = age + 19,
         birth = as.double(birth)) %>%
  left_join(
    df_s %>%
      dplyr::select(-age_s) %>%
      mutate(`(Intercept)` = 1) %>%
      pivot_longer(2:last_col(), names_to = "term", values_to = "value"),
    by = c("term", "age")
  ) %>%
  mutate(beta = coef*value) %>%
  group_by(birth, age) %>%
  summarise(estimate = sum(beta),
            .groups = "drop") %>%
  left_join(
    df_mod %>%
      group_by(birth) %>%
      summarise(min_age = min(age),
                max_age = max(age)),
    by = "birth"
  ) %>%
  filter(age >= min_age, age <= max_age) %>%
  ggplot() +
  aes(x = age, y = estimate, color = birth, group = birth) +
  geom_line() +
  theme_bw() +
  scale_color_viridis_c() +
  labs(x = "Age", y = "Predicted BMI", color = "Birth Year")

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
         uci = qnorm(.975, ranef, se)) %>%
  # mutate(birth = as.integer(birth)) %>%
  ggplot() +
  aes(x = birth, y = ranef, ymin = lci, ymax = uci) +
  facet_wrap(~ term, scales = "free_x") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange() +
  coord_flip() +
  theme_bw()
  

df_g <- df_clean %>%
  mutate(birth_lb = birth - birth %% 5,
         birth_ub = birth_lb + 4,
         age_lb = year - birth_ub,
         age_ub = year - birth_lb) %>%
  filter(between(bmi, 9, 90),
         age_lb >= 20,
         age_lb %% 5 == 0,
         age_ub <= 59) %>%
  mutate(birth_range = glue("{birth_lb}-{birth_ub}") %>%
           fct_reorder(birth_lb) %>%
           ordered(),
         age_range = glue("{age_lb}-{age_ub}") %>%
           fct_reorder(age_lb))




library(gamlss)
library(broom.mixed)

get_gamlss <- function(df){
  df_gamlss <- df %>%
    mutate(birth_range = fct_drop(birth_range))
  
  mod_g <- gamlss(bmi ~ birth_range,
                  sigma.formula = ~ birth_range,
                  nu.formula = ~ birth_range,
                  family = BCCGo, 
                  data = df_gamlss,
                  trace = FALSE)
  
  birth_levels <- levels(df_gamlss$birth_range)
  
  coefAll(mod_g) %>%
    map_dfr(as_tibble_row, .id = "param") %>%
    rename(intercept = `(Intercept)`) %>%
    pivot_longer(matches("^birth_range"), names_to = "term", values_to = "beta") %>%
    mutate(term = str_replace(term, "birth_range", "")) %>%
    group_by(param) %>%
    uncount(ifelse(row_number() == 1, 2, 1), .id = "row") %>%
    ungroup() %>%
    mutate(term = ifelse(row == 2, birth_levels[1], term),
           beta = ifelse(row == 2, intercept, beta + intercept)) %>%
    pivot_wider(term, names_from = param, values_from = beta)
}

gamlss_res <- df_q %>%
  mutate(mod = map(data, get_gamlss)) %>%
  dplyr::select(-data) %>%
  unnest(mod) %>%
  uncount(99, .id = "centile") %>%
  mutate(tau = centile/100,
         bmi = exp(mu) * (1 + nu * exp(sigma) * qnorm(tau))^(1/nu),
         term = ordered(term)) # %>%

ggplot(gamlss_res) +
  facet_wrap(~ age_range, nrow = 2) +
  aes(x = centile, y = bmi, colour = term, label = term) +
  geom_line() +
  geom_text_repel(data = gamlss_res %>%
                    filter(centile == 99)) +
  labs(x = "Centile", y = "BMI") +
  guides(color = "none") +
  theme_bw()
ggsave("Images/gamlss_bmi.png", width = 29.7, height = 21, units = "cm")


# Non-Linear Models  ----
df_mod2 <- df_mod %>%
  mutate(bmi_cat = ifelse(bmi > 30, 1, 0))

# 2. Linear Models
mod <- glmer(bmi_cat ~ age_ns1 + age_ns2 +
               (1 + age_ns1 + age_ns2 | birth), 
             df_mod2, binomial)
summary(mod)

df_term <- left_join(
  fixef(mod) %>%
    enframe(name = "term", value = "fixef"),
  ranef(mod)[[1]] %>%
    as_tibble(rownames = "birth") %>%
    pivot_longer(-birth, names_to = "term", values_to = "ranef"),
  by = "term"
)

df_term %>%
  mutate(coef = fixef + ranef) %>%
  uncount(64 - 18 + 1, .id = "age") %>%
  mutate(age = age + 17,
         birth = as.double(birth)) %>%
  left_join(
    df_s %>%
      dplyr::select(-age_s) %>%
      mutate(`(Intercept)` = 1) %>%
      pivot_longer(2:4, names_to = "term", values_to = "value"),
    by = c("term", "age")
  ) %>%
  mutate(beta = coef*value) %>%
  group_by(birth, age) %>%
  summarise(estimate = sum(beta),
            .groups = "drop") %>%
  left_join(
    df_mod %>%
      group_by(birth) %>%
      summarise(min_age = min(age),
                max_age = max(age)),
    by = "birth"
  ) %>%
  mutate(estimate = exp(estimate)/1+(exp(estimate))) %>%
  filter(age >= min_age, age <= max_age) %>%
  ggplot() +
  aes(x = age, y = estimate, color = birth, group = birth) +
  geom_line() +
  theme_bw() +
  scale_color_viridis_c() +
  labs(x = "Age", y = "Predicted BMI", color = "Birth Year")
#https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-9-49
# Alternatively, group together survey years (1992-1997, 2012-2017) and ages (20-29, etc.) and compare density plots!





df_mod <- df_clean %>%
  mutate(birth_decade = birth - birth %% 5) %>%
  filter(age >= 20, age <= 64) %>%
  right_join(df_s, by = "age") %>%
  filter(birth_decade >= 1930, birth_decade <= 1985) %>%
  nest(data = -birth_decade)

library(gamlss)

get_johnson <- function(df){
  mod_g <- gamlss(bmi ~ age_ns1 + age_ns2,
                  sigma.formula = ~ age_ns1 + age_ns2,
                  nu.formula = ~ age_ns1 + age_ns2,
                  family = BCCGo, 
                  data = df,
                  trace = FALSE)
  
  coefAll(mod_g) %>%
    map_dfr(as_tibble_row, .id = "param") %>%
    rename(age_ns0 = `(Intercept)`) %>%
    pivot_longer(-param, names_to = "term", values_to = "beta") %>%
    left_join(df_s %>%
                dplyr::select(-age_s) %>%
                mutate(age_ns0 = 1) %>%
                pivot_longer(-age, names_to = "term", values_to = "value"),
              by = "term") %>%
    group_by(param, age) %>%
    summarise(estimate = sum(beta*value),
              .groups = "drop") %>%
    pivot_wider(age, names_from = param, values_from = estimate) %>%
    filter(age >= min(df$age),
           age <= max(df$age))
}

df_mod %>%
  mutate(mod = map(data, get_johnson)) %>%
  unnest(mod) %>%
  uncount(99, .id = "centile") %>%
  mutate(birth_decade = ordered(birth_decade),
         tau = centile/100,
         bmi = exp(mu) * (1 + nu * exp(sigma) * qnorm(tau))^(1/nu)) %>%
  filter(centile %in% c(2, 5, 10, 50, 90, 95, 98)) %>%
  ggplot() +
  aes(x = age, y = bmi, color = birth_decade) +
  facet_wrap(~ centile) +
  geom_line()
