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
library(summarytools)

rm(list = ls())

# 1. Load Data ----
load("Data/hse_clean.Rdata")
load("Data/model_parameters.Rdata")

# 2. Make Dataset ----
df_cohort <- df_clean %>%
  dplyr::select(-edu) %>%
  arrange(birth) %>%
  mutate(cohort = age - age %% 5) %>%
  filter(age >= age_low, age <= age_high) %>%
  drop_na()

# 3. Run Models ----
get_linpred <- function(cohort, sex){
  data <- df_cohort %>%
    filter(cohort == !!cohort,
           sex %in% sexes[[!!sex]]) %>%
    group_by(year) %>%
    mutate(wt_int = wt_int*n()/sum(wt_int)) %>%
    ungroup()
  
  df_s <- data %>%
    distinct(year) %>%
    arrange(year) %>%
    mutate(ns(year, 3) %>%
             as_tibble() %>%
             mutate(across(everything(), as.double)) %>%
             rename_with(~ glue("year_ns{.x}")))
  
  df_mod <- data %>%
    left_join(df_s, by = "year")
  
  run_mod <- function(df){
    mod <- gamlss(bmi ~ year_ns1 + year_ns2 + year_ns3,
                  sigma.formula = ~ year_ns1 + year_ns2 + year_ns3,
                  nu.formula = ~ year_ns1 + year_ns2 + year_ns3,
                  family = BCCG, 
                  data = df, weights = df$wt_int,
                  trace = FALSE)
    
    coefAll(mod) %>%
      map_dfr(as_tibble_row, .id = "param") %>%
      pivot_longer(-param, names_to = "term", values_to = "beta") %>%
      left_join(df_s %>%
                  mutate(`(Intercept)` = 1) %>%
                  pivot_longer(-year, names_to = "term", values_to = "coef"),
                by = "term") %>%
      group_by(param, year) %>%
      summarise(estimate = sum(beta * coef),
                .groups = "drop") %>%
      pivot_wider(names_from = param, values_from = estimate)
  }
  
  set.seed(1)
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
res_age <- df_cohort %>% 
  distinct(cohort) %>% 
  expand_grid(sex = names(sexes)) %>%
  sample_frac() %>%
  mutate(res = future_map2(cohort, sex, get_linpred, .progress = TRUE)) %>%
  arrange(cohort, sex) %>%
  unnest(res)
future:::ClusterRegistry("stop")
toc()

save(res_age, file = "Data/gamlss_age_results.Rdata")


# 4. Predicted Values ----
load("Data/gamlss_age_results.Rdata")

get_ci <- function(x){
  quantile(x, probs = c(.5, .025, .975), 
           na.rm = TRUE) %>%
    set_names(c("beta", "lci", "uci")) %>%
    as_tibble_row()
}

# Linear Predictions
age_linpred <- res_age %>%
  uncount(99, .id = "centile") %>%
  mutate(cohort = glue("{cohort}-{cohort+4}") %>% ordered(),
         tau = centile/100,
         bmi = mu * (1 + nu * exp(sigma) * qnorm(tau))^(1/nu),
         boot = ifelse(boot == 0, 0, 1)) %>%
  arrange(cohort, sex, year, centile, boot) %>%
  group_by(cohort, sex, year, centile, boot) %>%
  summarise(get_ci(bmi), .groups = "drop_last") %>%
  summarise(beta = nth(beta, 1), lci = nth(lci, 2), 
            uci = nth(uci, 2), .groups = "drop")

save(age_linpred, file = "Data/gamlss_age_predictions.Rdata")


# 5. Plots ----
load("Data/gamlss_age_predictions.Rdata")

# Plot 1
obese_cat <- tribble(
  ~centile, ~year, ~beta, ~label,
  "10th", -Inf, 18.5, "Normal",
  "10th", -Inf, 25, "Overweight",
  "10th", -Inf, 30, "Obese")

age_linpred %>%
  filter(centile %in% c(10, 25, 50, 75, 90),
         sex == "all") %>%
  mutate(centile = glue("{centile}th")) %>%
  ggplot() +
  aes(x = year, y = beta) +
  facet_grid(~ centile, switch = "y") +
  geom_hline(yintercept = c(18.5, 25, 30), linetype = "dashed", color = "grey60") + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = cohort), color = NA, alpha = 0.2) +
  geom_line(aes(color = cohort)) + 
  geom_text(data = obese_cat, aes(label = label),
            vjust = -0.5, hjust = -0.1, color = "grey50") +
  scale_x_continuous(breaks = c(1991, 2000, 2010, 2019)) +
  scale_color_viridis_d(begin = 0.2, end = 0.8, option = "B") +
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "B") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.spacing = unit(1, "lines"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        strip.background.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Year", y = NULL, color = "Age Group", fill = "Age Group") +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
ggsave("Images/gamlss_age_1.png", width = 29.7, height = 21, units = "cm")

# Plot 2
age_linpred %>%
  filter(centile %in% c(10, 25, 50, 75, 90),
         sex != "all") %>%
  mutate(sex = str_to_title(sex),
         centile = glue("{centile}th")) %>%
  ggplot() +
  aes(x = year, y = beta) +
  facet_grid(sex ~ centile, switch = "y") +
  geom_hline(yintercept = c(18.5, 25, 30), linetype = "dashed", color = "grey60") + 
  geom_ribbon(aes(ymin = lci, ymax = uci, fill = cohort), color = NA, alpha = 0.2) +
  geom_line(aes(color = cohort)) + 
  geom_text(data = obese_cat, aes(label = label),
            vjust = -0.5, hjust = -0.1, color = "grey50") +
  scale_x_continuous(breaks = c(1991, 2000, 2010, 2019)) +
  scale_color_viridis_d(begin = 0.2, end = 0.8, option = "B") +
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "B") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.spacing = unit(1, "lines"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0),
        strip.background.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Year", y = NULL, color = "Age Group", fill = "Age Group")+
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
ggsave("Images/gamlss_age_2.png", width = 29.7, height = 21, units = "cm")

# 6. Rename Files ----
file.copy("Images/lmer_1_2019.png", "Paper Items/figure_1.png", overwrite = TRUE)
file.copy("Images/gamlss_2a.png", "Paper Items/figure_2.png", overwrite = TRUE)

file.copy("Images/re_bmi_all_1_2019.png", "Paper Items/figure_s1.png", overwrite = TRUE)
file.copy("Images/re_obese_all_1_2019.png", "Paper Items/figure_s2.png", overwrite = TRUE)
file.copy("Images/re_bmi_sex_1_2019.png", "Paper Items/figure_s3.png", overwrite = TRUE)
file.copy("Images/re_obese_sex_1_2019.png", "Paper Items/figure_s4.png", overwrite = TRUE)
file.copy("Images/gamlss_3a.png", "Paper Items/figure_s5.png", overwrite = TRUE)
file.copy("Images/gamlss_4a.png", "Paper Items/figure_s6.png", overwrite = TRUE)
file.copy("Images/gamlss_2b.png", "Paper Items/figure_s7.png", overwrite = TRUE)
file.copy("Images/gamlss_3b.png", "Paper Items/figure_s8.png", overwrite = TRUE)
file.copy("Images/gamlss_4b.png", "Paper Items/figure_s9.png", overwrite = TRUE)
file.copy("Images/gamlss_age_1.png", "Paper Items/figure_s10.png", overwrite = TRUE)
file.copy("Images/gamlss_age_2.png", "Paper Items/figure_s11.png", overwrite = TRUE)