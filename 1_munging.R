library(tidyverse)
library(glue)
library(haven)
library(labelled)

rm(list = ls())

# 1. Create Lookup
hse_fld <- "D:/HSE"
# hse_fld <- "/Volumes/CLS DATA/HSE"

get_label <- function(var){
  lbl <- attr(var, "label", exact = TRUE)
  lbl <- ifelse(is.null(lbl), "", lbl)
  
  tibble(label = lbl)
}

get_lookup <- function(file){
  read_dta(file, n_max = 1) %>%
    map_dfr(get_label, .id = "name")
}

df_files <- tibble(fld = list.files(hse_fld)) %>%
  filter(str_detect(fld, "^[1-9]")) %>%
  mutate(full_fld = glue("{hse_fld}/{fld}"),
         fls = map(full_fld, list.files, full.names = TRUE)) %>%
  unnest(fls) %>%
  filter((fld %in% c("1991-1992", "1993") |
            str_detect(fls, "[0-9](ai|gp3|gpa)") |
            str_detect(fls, "eul"))) %>%
  mutate(size = map_dbl(fls, file.size)  / (1024^2),
         lookup = map(fls, get_lookup)) %>%
  arrange(fld, desc(size)) %>%
  group_by(fld) %>%
  mutate(fld_N = n(),
         fld_n = row_number()) %>%
  ungroup() %>%
  unnest(lookup)

save(df_files, file = "Data/var_lookup.Rdata")
write_csv(df_files, "Data/var_lookup.csv")


# 2. Open Files ---
load("Data/var_lookup.Rdata")

any_list <- c("year", "age", "Age", "Age90", "Age35g",
              "bmi", "BMI", "bmiok", "BMIok", "BMIOK",
              "bmival", "BMIval", "samptype",
              "SampType", "sex", "Sex", "wt_int")
grep_list <- "(ethni|ethin|ethcind|origin|Origin|topqual|TopQual|Topqual)"

edu_lbl <- c("NVQ 4/5", "HE Below Degree", "NVQ3", "NVQ2",
             "NVQ1", "Foreign/Other", "None", "Current Student")


df <- df_files %>%
  distinct(fld, fls) %>%
  mutate(dta = map(fls, 
                   ~ read_dta(.x, col_select = c(any_of(any_list), 
                                                 matches(grep_list))) %>%
                     rename_all(str_to_lower)))

save(df, file = "Data/df_raw.Rdata")

# 3. Clean Files ----
load("Data/df_raw.Rdata")

set.seed(1)
clean_df <- function(fld, dta){
  fld_year <- str_sub(fld, 1, 4) %>% as.numeric()
  
  if (fld_year == 1991){ # ADD YEAR
    dta <- dta %>%
      mutate(year = 1900 + year)
  } else{
    dta <- dta %>%
      mutate(year = !!fld_year)
  }
  
  dta <- dta %>%
    mutate(bmi = ifelse(bmiok == 1, bmi, NA))
  
  if (fld_year <= 1997){ # ETHNICITY
    dta <- dta %>%
      filter(ethnic == 1)
  } else if (fld_year == 1998){
    dta <- dta %>%
      filter(nethnic == 1)
  } else if (between(fld_year, 1999, 2003)){
    dta <- dta %>%
      filter(ethnici == 1)
  } else if (fld_year == 2004){
    dta <- dta %>%
      filter(ethcind == 1)
  } else if (between(fld_year, 2005, 2007)){
    dta <- dta %>%
      filter(ethinda == 1)
  } else if (between(fld_year, 2008, 2010)){
    dta <- dta %>%
      filter(between(origin, 1, 3))
  } else if (between(fld_year, 2011, 2013)){
    dta <- dta %>%
      filter(between(origin, 1, 4))
  } else if (fld_year %in% 2014:2019){
    dta <- dta %>%
      filter(origin2 == 1)
  }
  
  sample_age <- function(age_range){
    lbls <- names(val_labels(age_range))[val_labels(age_range) %in% 1:21]
    
    sample_range <- function(age_cat){
      x <- as.double(age_cat)
      
      if (is.na(x)){
        smp <- NA
      } else{
        sample(x[1]:x[2], 1)
      }
    }
    
    factor(age_range, levels = 1:21, labels = lbls) %>%
      as.character() %>%
      str_split("\\-") %>%
      map_int(sample_range)
    
  }
  
  if (fld_year == 2014){ # Age
    dta$age <- ifelse(dta$age90 == 90, NA, dta$age90)
  } else if (fld_year %in% 2015:2019){
    dta$age <- sample_age(dta$age35g)
  }
  
  
  if (fld_year %in% c(2000, 2002)){ # CORRECT SAMPLE
    dta <- dta %>%
      filter(samptype == 2)
  } else if (fld_year %in% c(2005, 2007, 2008, 2009, 2010, 2015)){
    dta <- dta %>%
      filter(samptype == 1)
  } else if (fld_year == 2006){
    dta <- dta %>%
      filter(between(samptype, 1, 2))
  }
  
  if (fld_year %in% 1995:1996){
    dta$edu <- factor(rep(NA, nrow(dta)), levels = 1:8, labels = edu_lbl)
  } else{
    dta$edu <- factor(dta$topqual2, levels = 1:8, labels = edu_lbl)
  }
  
  
  if (is.null(dta$wt_int)){ # INTERVIEW WEIGHT
    dta$wt_int <- 1
  }
  
  dta <- dta %>%
    mutate(sex = ifelse(sex == 1, "Male", "Female") %>%
             factor()) %>%
    select(year, sex, age, bmi, edu, wt_int)
  
  return(dta)
}

df_clean <- df %>%
  mutate(clean = map2(fld, dta, clean_df)) %>%
  select(clean) %>%
  unnest(clean) %>%
  zap_labels() %>%
  zap_label() %>%
  zap_formats() %>%
  mutate(age = ifelse(age <= 0, NA, age),
         birth = year - age,
         bmi = ifelse(between(bmi, 13, 70), bmi, NA))

save(df_clean, file = "Data/hse_clean.Rdata")
write_dta(df_clean, "Data/hse_clean.dta")


