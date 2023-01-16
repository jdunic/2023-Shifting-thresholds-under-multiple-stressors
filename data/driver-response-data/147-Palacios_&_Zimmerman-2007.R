library(tidyverse)
library(lubridate)
library(readxl)
library(fuzzyjoin)
library(IRanges)

#BiocManager::install("IRanges")

pal1 <- 
  read_xlsx('data/raw-data/driver-response-data/147-Palacios_&_Zimmerman-2007.xlsx', sheet = 1) %>% 
  mutate(start_date = ymd(start_date), measurement_date = ymd(measurement_date)) %>% 
  mutate(date_interval = interval(ymd(start_date), ymd(measurement_date)))
pal2 <- 
  read_xlsx('data/raw-data/driver-response-data/147-Palacios_&_Zimmerman-2007.xlsx', sheet = 2) %>% 
  mutate(light_date = ymd(light_date), temp_date = ymd(temp_date))

light_values <- 
  pal2 %>% 
  select(light_date, dli) %>% 
  mutate(date_end = light_date)

temp_values <- 
  pal2 %>% 
  drop_na() %>% 
  select(temperature, temp_date) %>% 
  mutate(date_end = temp_date)

pal_light <- 
  interval_left_join(pal1, light_values, by = c('start_date' = 'light_date', 'measurement_date' = 'date_end')) %>% 
  group_by(treatment, date_interval) %>%
  summarise(seasonal_mean_dli = mean(dli)) %>%
  ungroup()

pal_temp <- 
  interval_left_join(pal1, temp_values, by = c('start_date' = 'temp_date', 'measurement_date' = 'date_end')) %>% 
  group_by(treatment, date_interval) %>%
  summarise(seasonal_mean_temp = mean(temperature)) %>%
  ungroup()

pal <- 
  left_join(pal1, pal_light) %>%
  left_join(pal_temp)

write_csv(pal, 'data_outputs/seasonal_driver_data/147-Palacios_&_Zimmerman-2007_seasonal.csv')