library(tidyverse)
library(readxl)

# source(here('R', '01_functions.R'))

study_data <- 
  read_xlsx(here('data', 'study-data.xlsx'), na = c("", "NA", "np")) %>% 
  mutate(clean_lat = dms2dec(lat), clean_lon = dms2dec(lon)) %>% 
  # For now just grab first value, as of Oct 4, I have cleaned the spreadsheet so this works
  # This may not hold in the future. 
  mutate(clean_lat2 = if_num_parse_num(lat), clean_lon2 = if_num_parse_num(lon)) %>% 
  mutate(clean_lat = as.numeric(ifelse(is.na(clean_lat), clean_lat2, clean_lat)), 
         clean_lon = as.numeric(ifelse(is.na(clean_lon), clean_lon2, clean_lon))) %>% 
  select(-clean_lat2, -clean_lon2) %>% 
  rename('pub_year' = 'year')

driver_data <- 
  read_xlsx(here('data', 'drivers.xlsx'), na = c("", "NA")) %>%
  # Clean driver data
  # parsing extracts first number of string (e.g., 144 +/- 10 --> 144)
  dplyr::rename(surface_irradiance =`surface irradiance (µmol photons m-2 s-1)`) %>% 
  mutate(pen_light = if_num_parse_num(`penetrating irradiance`), 
         surface_light = if_num_parse_num(surface_irradiance), 
         light_hours = if_num_parse_num(`light hours`), 
         mean_dli = if_num_parse_num(`mean DLI (mol m-2day-1)`)) %>%
  left_join(., select(study_data, studyid, author, pub_year, location, clean_lat, clean_lon), 
            by = c('studyid', 'site' = 'location')) %>% 
  rename(lat = clean_lat, lon = clean_lon) %>%
  mutate_at(c("pen_light", "surface_light", "light_hours", "mean_dli"), as.numeric) %>% 
  # Deal with percent surface irradiance - for now use the reduced light as the pen light
  mutate(pen_light = ifelse((is.na(pen_light) & !is.na(surface_light) & !is.na(`%SI`)), 
                            (surface_light * as.numeric(`%SI`) / 100), pen_light)) %>% 
  # Use penetrating irradiance preferentially, if 'np', use surface light
  mutate(DLI = case_when(!is.na(pen_light) ~ pen_light * light_hours * 3600 / 1e6, 
                         (is.na(pen_light) & !is.na(mean_dli)) ~ mean_dli, 
                         (`penetrating irradiance` == "np" & is.na(mean_dli) & !is.na(surface_light)) ~ 
                           surface_light * light_hours * 3600 / 1e6)) %>%  # Units = mol / m2 / day
  mutate(clean_temp = parse_number(temperature_c)) %>% 
  mutate(parsed_collection_date = str_extract(collection_date, "(^\\d{4}-\\d{2})")) %>% 
  mutate(parsed_collection_date = parse_date(parsed_collection_date, "%Y-%m"))

write_csv(driver_data, here('data_outputs', 'clean_driver_data.csv'))

# Load Seasonal Data
#-------------------------------------------------------------------------------
seasonal_files <- 
  tibble(folder = here('data', 'seasonal_driver_data'), 
         files = list.files(here('data', 'seasonal_driver_data'))
  ) %>% 
  mutate(path = here(folder, files)) %>% 
  select(path) %>% 
  filter(., str_detect(path, "Icon\\\r", negate = TRUE))

seasonal <- 
  apply(seasonal_files, 1, load_seasonal_data) %>% 
  bind_rows()
