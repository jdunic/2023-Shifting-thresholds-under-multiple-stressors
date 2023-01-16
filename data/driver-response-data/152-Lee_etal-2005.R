library(tidyverse)
library(lubridate)
library(readxl)
library(fuzzyjoin)
library(IRanges)
#BiocManager::install("IRanges")

lee1 <- 
  read_xlsx('data/raw-data/driver-response-data/152-Lee_etal-2005.xlsx', sheet = 1) %>% 
  mutate(start_date = ymd(start_date), measurement_date = ymd(measurement_date)) %>% 
  #mutate(date_interval = interval(ymd(start_date), ymd(measurement_date))) %>% 
  mutate(numeric_start = as.numeric(start_date), numeric_end = as.numeric(measurement_date))
# Sheet 3 has the imputed values using webplot digitizer to extract along the line
lee2 <- 
  read_xlsx('data/raw-data/driver-response-data/152-Lee_etal-2005.xlsx', sheet = 3) %>% 
  mutate(light_date = ymd(light_date), temp_date = ymd(temp_date)) %>% 
  mutate(numeric_light_date = as.numeric(light_date), numeric_temp_date = as.numeric(temp_date))

light_values <- 
  lee2 %>% 
  select(treatment, light_date, dli) %>% 
  drop_na() %>% 
  mutate(date_end = light_date) %>% 
  mutate(numeric_start = as.numeric(light_date), numeric_end = as.numeric(date_end))

lee_light <- 
  fuzzyjoin::genome_left_join(lee1, light_values,
                   by = c('treatment' = 'treatment', 'numeric_start' = 'numeric_start', 'numeric_end' = 'numeric_end')) %>%  
  dplyr::rename("treatment" = "treatment.x") %>% 
  group_by(treatment, start_date, measurement_date) %>%
  summarise(seasonal_mean_dli = mean(dli)) %>%
  ungroup()

temp_values <- 
  lee2 %>% 
  select(treatment, temp_date, temperature) %>% 
  drop_na() %>% 
  mutate(date_end = temp_date) %>% 
  mutate(numeric_start = as.numeric(temp_date), numeric_end = as.numeric(date_end))

lee_temp <- 
  genome_left_join(lee1, temp_values,
                   by = c('treatment' = 'treatment', 'numeric_start' = 'numeric_start', 'numeric_end' = 'numeric_end')) %>%  
  dplyr::rename("treatment" = "treatment.x") %>% 
  group_by(treatment, start_date, measurement_date) %>%
  summarise(seasonal_mean_temp = mean(temperature)) %>%
  ungroup()

lee <- 
  left_join(lee1, lee_light) %>%
  left_join(lee_temp)

#lee %>% select(-notes) %>% view

write_csv(lee, 'data_outputs/seasonal_driver_data/152-Lee_etal-2005_seasonal.csv')