library(tidyverse)

brakel <- 
  read_delim('data/raw-data/driver-response-data/312-Brakel_etal-2017_raw/datasets/z_marina_L_zosterae_growth.tab', 
             delim = '\t', skip = 17) %>% 
  mutate(treatment = case_when((`T:nutrients` == 6 & Treat == "no") ~ 'N- L-', 
                               (`T:nutrients` == 6 & Treat == "yes") ~ 'N- L+',
                               (`T:nutrients` == 57 & Treat == "no") ~ 'N+ L-', 
                               (`T:nutrients` == 57 & Treat == "yes") ~ 'N+ L+'))

brakel %>% 
  group_by(treatment) %>% 
  summarise(mean_shoot_production = mean(`Shoots [#]`), se = sd(`Shoots [#]`) / n()) 
