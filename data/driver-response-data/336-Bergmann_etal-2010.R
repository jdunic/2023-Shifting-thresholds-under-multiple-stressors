library(tidyverse)
library(lubridate)
library(readxl)
library(fuzzyjoin)
library(IRanges)
#BiocManager::install("IRanges")

bergmann1 <- 
  read_xlsx('data/raw-data/driver-response-data/336-Bergmann_etal-2010.xlsx', sheet = 1) %>% 
  mutate(start_date = ymd(start_date), measurement_date = ymd(measurement_date)) %>% 
  #mutate(date_interval = interval(ymd(start_date), ymd(measurement_date))) %>% 
  mutate(numeric_start = as.numeric(start_date), numeric_end = as.numeric(measurement_date)) %>% 
  unite(key, site, treatment, remove = FALSE, sep = "_")

bergmann2 <- 
  read_xlsx('data/raw-data/driver-response-data/336-Bergmann_etal-2010.xlsx', sheet = 2) %>% 
  mutate(numeric_date = as.numeric(ymd(date))) %>% 
  mutate(numeric_start = numeric_date, numeric_end = numeric_date) %>% 
  unite(key, site, treatment, remove = FALSE, sep = "_") %>% 
  mutate(dli = `light hours` * `light µmol photons/m2/s` * 3600 / 1e6) %>% 
  select(key, numeric_start, numeric_end, dli, temperature)

dat <- 
  fuzzyjoin::genome_left_join(bergmann1, bergmann2,
                     by = c('key' = 'key', 'numeric_start' = 'numeric_start', 'numeric_end' = 'numeric_end')) %>%  
    dplyr::rename("key" = "key.x") %>% 
    group_by(site, treatment, start_date, measurement_date) %>%
    summarise(seasonal_mean_dli = mean(dli), 
              seasonal_mean_temp = mean(temperature)) %>%
    ungroup()

names(dat)

bergmann <- 
  left_join(bergmann1, dat) %>% 
  select(-key)

write_csv(bergmann, 'data_outputs/seasonal_driver_data/336-Bergmann_etal_2010_seasonal.csv')