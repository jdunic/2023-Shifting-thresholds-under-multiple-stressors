library(tidyverse)

#source(here('R', '01_functions.R'))

# Omit responses where there were other covariates
inclusion_strings <- readxl::read_xlsx(here('data', 'light-temperature_inclusion_strings.xlsx'))

# Factor levels
experiment_levels <- c("lab mesocosm", 
                       "greenhouse mesocosm", "outdoor mesocosm", 
                       "field manipulation", "field observation")

clean_pop_growth_df <- 
  read_csv(here('data', 'clean_pop_growth_df.csv')) %>% 
  mutate(studyid = factor(studyid), 
         treatment = factor(treatment),
         study_treatment = factor(paste(studyid, treatment, sep = "_")), 
         study_site = factor(paste(studyid, site, sep = "_")), 
         study_site_treatment = factor(paste(studyid, site, treatment, sep = "_"))) %>% 
  drop_na(DLI, clean_temp) %>% 
  mutate(field_light = ifelse(str_detect(`light hours`, 'field'), 'suncalc light hours', 'provided light hours')) %>% 
  mutate(sign = case_when(pop_growth_rate < 0 ~ 'loss', 
                                            pop_growth_rate > 0 ~ 'gain', 
                                            pop_growth_rate == 0 ~ 'no change')) %>% 
  left_join(., inclusion_strings) %>% 
  filter(include == 'yes') %>% 
  mutate(lab_outdoors = ifelse(str_detect(`exp setup`, 'lab'), 'lab', 'outdoors')) %>% 
  mutate(surface_pen_light = ifelse(is.na(pen_light), 'likely_not_pen', 'pen')) 

# Manually truncate observed values to improve visualisation 
lower_tail <- -0.1
upper_tail <-  0.1

pop_growth_df <-  
  clean_pop_growth_df %>% 
  filter(studyid != 135) %>%  # This study light values would need to be corrected
  mutate(pop_growth_rate_truncated = case_when(pop_growth_rate < lower_tail ~ lower_tail, 
                                               pop_growth_rate > upper_tail  ~ upper_tail, 
                                               TRUE ~ pop_growth_rate))
