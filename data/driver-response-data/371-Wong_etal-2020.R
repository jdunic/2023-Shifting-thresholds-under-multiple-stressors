library(tidyverse)
library(lubridate)
library(readxl)
library(fuzzyjoin)
library(IRanges)
#BiocManager::install("IRanges")

wong1 <- 
  read_xlsx('data/raw-data/driver-response-data/371-Wong_etal-2020.xlsx', sheet = 1) %>% 
  mutate(start_date = ymd(start_date), measurement_date = ymd(measurement_date)) %>% 
  mutate(numeric_start = as.numeric(start_date), numeric_end = as.numeric(measurement_date))

wong2 <- 
  read_xlsx('data/raw-data/driver-response-data/371-Wong_etal-2020.xlsx', sheet = 2)

dates <- seq(ymd('2017-06-15'),ymd('2017-11-29'), by = '1 day')
week_days <- seq(ymd('2017-06-15'),ymd('2017-11-30'), by = '4 weeks')

paste(week_days, collapse = ",")

# wong3 <- 
#   wong2 %>% 
#     group_by(treatment) %>% 
#     dplyr::slice(rep(1:n(), 12)) %>% 
#     mutate(date = dates) %>% 
#     ungroup()
# write_csv(wong3, 'data/raw-data/driver-response-data/371-Wong_etal-2020_seasonal-light_sheet3.csv')

wong3 <- 
  read_xlsx('data/raw-data/driver-response-data/371-Wong_etal-2020.xlsx', sheet = "Sheet3") %>% 
  mutate(numeric_start = as.numeric(ymd(date)), numeric_end = as.numeric(ymd(date))) %>% 
  mutate(dli = `light hours` * `light µmol photons/m2/s` * 3600 / 1e6)

dat <- 
  fuzzyjoin::genome_left_join(wong1, wong3,
                     by = c('treatment' = 'treatment', 'numeric_start' = 'numeric_start', 'numeric_end' = 'numeric_end')) %>%  
    dplyr::rename("treatment" = "treatment.x") %>% 
    group_by(treatment, start_date, measurement_date) %>%
    summarise(seasonal_mean_dli = mean(dli)) %>%
    ungroup()

wong <- left_join(wong1, dat)

write_csv(wong, 'data_outputs/seasonal_driver_data/371-Wong_etal_2020_seasonal.csv')