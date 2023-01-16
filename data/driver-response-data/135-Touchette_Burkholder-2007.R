library(imputeTS)

# Relies on study_data because we need to get light hours
source('R/driver-data-prep.R')

touch1 <- 
  read_xlsx('data/raw-data/driver-response-data/135-Touchette_etal_2007.xlsx', sheet = 1) %>% 
  mutate(start_date = ymd(start_date), measurement_date = ymd(measurement_date)) %>% 
  mutate(date_interval = interval(ymd(start_date), ymd(measurement_date))) %>% 
  mutate(numeric_start = as.numeric(start_date), numeric_end = as.numeric(measurement_date))

touch2 <- 
  read_xlsx('data/raw-data/driver-response-data/135-Touchette_etal_2007.xlsx', sheet = 2) %>% 
  pivot_longer(control:`HT + HN`, names_to = "treatment", values_to = "value") %>% 
  separate(value, into = c('value', 'se'), sep = "±") %>% 
  # NEED TO FIX and parse <0.01 if using nutrient values in the future
  mutate(date = ymd(date), value = as.numeric(value), se = as.numeric(se))


# Visualisations of seasonal data/imputation
# ------------------------------------------------------------------------------
touch2 %>%
ggplot(data = ., aes(x = date, y = value, colour = treatment)) +
  geom_point() + 
  geom_line() + 
  facet_wrap(~ parameter, scales = 'free_y')

touch2 %>%
ggplot(data = ., aes(x = date, y = value, colour = treatment)) +
  geom_point() + 
  stat_smooth() +
  facet_wrap(~ parameter, scales = 'free_y')

# ------------------------------------------------------------------------------
# Impute/aggregate seasonal data
# ------------------------------------------------------------------------------

touch_light <- 
  touch2 %>% 
    filter(parameter == "Light") %>% 
    arrange(treatment, date) %>%
    select(treatment, date, value) %>% 
    rename(light = value)

imputed_light <- 
  touch_light %>% 
  complete(treatment, date = full_seq(date, period = 1)) %>% 
  na_interpolation %>% 
  mutate(numeric_start = as.numeric(date), numeric_end = as.numeric(date))

light_hours <- 
  # Get lat / lon
  driver_data %>% 
  select(studyid, site, treatment, lat, lon) %>% 
  filter(., studyid == 135) %>% 
  left_join(., imputed_light) %>%
  rowwise() %>% 
  mutate(x = list(getSunlightTimes(date = date, lat = lat, lon = lon, keep = c('sunrise', 'sunset')))) %>% 
  unnest() %>%
  mutate(hours = as.numeric(sunset - sunrise)) %>% 
  mutate(sun_dli = light * hours * 3600 / 1e6)

seasonal_light <- 
  fuzzyjoin::genome_left_join(touch1, light_hours,
                   by = c('treatment' = 'treatment', 'numeric_start' = 'numeric_start', 'numeric_end' = 'numeric_end')) %>%  
  dplyr::rename("treatment" = "treatment.x") %>% 
  group_by(treatment, start_date, measurement_date) %>%
  summarise(seasonal_mean_dli = mean(sun_dli)) %>%
  ungroup()

# Get seasonal temperature means

touch_temp <- 
  touch2 %>% 
    filter(parameter == "Temperature") %>% 
    arrange(treatment, date) %>%
    select(treatment, date, value) %>% 
    rename(temperature = value)

imputed_temp <- 
  touch_temp %>% 
  complete(treatment, date = full_seq(date, period = 1)) %>% 
  na_interpolation %>% 
  mutate(numeric_start = as.numeric(date), numeric_end = as.numeric(date))

ggplot(data = imputed_temp, aes(x = date, y = temperature, colour = treatment)) + 
  geom_line()

seasonal_temp <- 
  fuzzyjoin::genome_left_join(touch1, imputed_temp,
                   by = c('treatment' = 'treatment', 'numeric_start' = 'numeric_start', 'numeric_end' = 'numeric_end')) %>%  
  dplyr::rename("treatment" = "treatment.x") %>% 
  group_by(treatment, start_date, measurement_date) %>%
  summarise(seasonal_mean_temp = mean(temperature)) %>%
  ungroup()

touch <- 
  left_join(touch1, seasonal_light) %>%
  left_join(seasonal_temp) %>% 
  mutate(`exposure time` = "calculate")

write_csv(touch, 'data_outputs/seasonal_driver_data/135-Touchette_Burkholder-2007.csv')

