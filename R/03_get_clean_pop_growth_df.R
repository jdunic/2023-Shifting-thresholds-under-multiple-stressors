library(tidyverse)
library(lubridate)
library(suncalc)

# source(here('R', '01_functions.R'))
# source(here('R', '02_driver-data-prep.R'))

# ------------------------------------------------------------------------------
# Clean mortality data
# ----------------------------------------------------------------------------
get_mortality_data <- function() {

mortality_data <- 
  load_response_data(response = 'mortality') %>%
  # survival is inverse of mortality
  # NEED TO DEAL WITH MULTIPLE RESPONSES AND STANDARDISING
  mutate(clean_response = case_when(str_detect(response, "survival|mortality") ~ "shoot mortality", 
                                    str_detect(response, "production") ~ "shoot production", 
                                    str_detect(response, "density") ~ "shoot density")) %>%
  # Remove shoot densities, they are handled separately
  filter(clean_response != "shoot production") %>%
  mutate_at(c("studyid", "value", "upper se", "se"), as.numeric) %>%
  # standardise by exposure time to days
  mutate(`exposure time` = ifelse(`exposure time` == "calculate", 
                                  paste0(as.character(parse_date(measurement_date) - parse_date(start_date)), " days"), 
                                  `exposure time`)) %>% 
  mutate(std_time = duration(`exposure time`) %>% seconds_to_period %>% day) %>% 
  # standardise to mortality proportion
  mutate(clean_value = case_when((response == "shoot survival" & units == "proportion") ~ (1 - value),
                                 response != "shoot survival" & units == "proportion" ~ value, 
                                 units == "%" ~ (value / 100), 
                                 units == "% d-1" ~ (value / 100) * std_time)
  ) %>% 
  mutate(mortality_rate = case_when(units == "% d-1" ~ clean_value, 
                                    units != "% d-1" ~ clean_value / std_time)) %>% 
  mutate(lifestage = ifelse(str_detect(response, "seedling"), "seedling", "adult")) %>% 
  # Add driver data after responses are clean
  mutate(studyid = as.numeric(studyid)) %>%  # studyid type needs to be the same
  left_join(., driver_data, by = c('studyid', 'site', 'treatment')) %>% 
  # Use lambda --> Nt / N0; 
  mutate(lambda = 1 - clean_value)

assign("mortality_data", mortality_data, envir = .GlobalEnv)
}

# ------------------------------------------------------------------------------
# Clean shoot production
# ------------------------------------------------------------------------------
get_shoot_production_data <- function() {

shoot_prod <- 
  load_response_data(response = ('change in shoot density|lateral shoot production|shoot density')) %>% 
  # Get only shoot production values 
  filter(str_detect(response, 'change in shoot density|lateral shoot production|shoot density')) %>%
  mutate_at(c("studyid", "value", "upper se", "se"), as.numeric) %>%
  mutate(lifestage = ifelse(str_detect(response, "seedling"), "seedling", "adult")) %>% 
  mutate(clean_response = 
    case_when(str_detect(response, 'change in shoot density|shoot density') ~ 'shoot density change', 
              str_detect(response, 'lateral shoot production') ~ 'lateral shoot production')
  ) %>% 
  drop_na(studyid) %>%  # drop these missing values
  filter(studyid != "118") %>%  # drop this study for now because the data were a mess
  #filter(response %in% c('change in shoot density', 'lateral shoot production', 'shoot density')) %>% 
  mutate(`exposure time` = ifelse(`exposure time` == "calculate", 
                                  paste0(as.character(parse_date(measurement_date) - parse_date(start_date)), " days"), 
                                  `exposure time`)) %>% 
  mutate(std_time = duration(`exposure time`) %>% seconds_to_period %>% day) %>% 
  # Add driver data after responses are clean
  left_join(., driver_data, by = c('studyid', 'site', 'treatment'))

  assign("shoot_prod", shoot_prod, envir = .GlobalEnv)

# Let's compare everything to initial shoot density because that is the most 
# common baseline we will have in the dataset
# Deal with uncertainty later... or not...
min_shoots <- 
  shoot_prod %>%
    filter(response == "shoot density") %>% 
    # Do not include studies with seasonal data, these require diffs between rows
    filter(!(studyid %in% unique(seasonal$studyid))) %>%
    group_by(studyid, site, treatment) %>% 
    slice(which.min(std_time)) %>% 
    select(studyid, site, treatment, value) %>% 
    rename(initial_density = "value") %>% 
    ungroup()

assign("min_shoots", min_shoots, envir = .GlobalEnv)

#### *****
#### NEED to check that lateral shoot / shoot makes sense to translate directly to % change???
# Yes it is, because it is the per capita rate of growth
#### e.g., studyid 24
shoot_change_data <- 
  #left_join(shoot_prod %>% filter(response == "shoot density"), min_shoots) %>% 
  left_join(shoot_prod, min_shoots) %>% 
  mutate(prop_change = (value - initial_density) / initial_density) %>% 
  mutate(clean_value = case_when(is.na(prop_change) & units == "%" ~ value / 100, 
                                 is.na(prop_change) & units == "lateral shoot / shoot" ~ value, 
                                 is.na(prop_change) & units == "shoot / shoot / day" ~ value * std_time,
                                 studyid == 135 ~ value / 675,  # Deal with studyid == 135 manually
                                 !is.na(prop_change) ~ prop_change)) %>% 
  # Studies with seasonal data will require diffs between rows, handle separately
  filter(response != "shoot density" & studyid %in% unique(seasonal$studyid) |
         !(studyid %in% unique(seasonal$studyid))) %>% 
  group_by(studyid, site, treatment) %>% 
  mutate(lambda = case_when(str_detect(units, "shoot \\/ shoot$") ~ value, 
                            str_detect(units, "shoot \\/ shoot \\/ day$") ~ value * std_time,
                            (str_detect(clean_response, "shoot density change") & units != "%") ~ value / initial_density,
                            (str_detect(clean_response, "shoot density change") & units == "%") ~ value + 1,
                           TRUE ~ 999)) %>% 
  ungroup()

study_135_shoot_change_data <- filter(shoot_change_data, studyid == 135)
study_135_lambda <- 
  study_135_shoot_change_data %>% 
  select(studyid, site, treatment, start_date, measurement_date, initial_density, response, units, value, clean_response, std_time, clean_value) %>% 
  mutate(initial_density = 675) %>% 
  group_by(site, treatment) %>%
  mutate(t = row_number()) %>%
  mutate(time_delta = (std_time - lag(std_time))) %>% 
  mutate(Nt = ifelse(is.na(time_delta), value * std_time + 675, 999)) %>% 
  mutate(Nt = ifelse(!is.na(time_delta), value * time_delta + 675, Nt)) %>% 
  ungroup() %>% 
  mutate(lambda = ifelse(!is.na(time_delta), Nt / lag(Nt), Nt / 675)) %>% 
  select(studyid, site, treatment, lambda)
study_135_shoot_change_data <- left_join(study_135_shoot_change_data %>% select(-lambda), study_135_lambda)
shoot_change_data <- 
  shoot_change_data %>% 
  filter(studyid != 135) %>%
  bind_rows(study_135_shoot_change_data)


assign("shoot_change_data", shoot_change_data, envir = .GlobalEnv)

# Get shoot change data for studies that have accompanying seasonal environmental 

# data, e.g., studyid 147, also 57 does have but I have not eshoot_prodtracted the seasonal data
seasonal_shoot_change_data <- 
  shoot_prod %>%
    filter(response == "shoot density") %>% 
    filter(studyid %in% unique(seasonal$studyid)) %>%
    arrange(studyid, site, treatment, measurement_date) %>% 
    group_by(studyid, site, treatment) %>% 
    mutate(prop_change = (value - lag(value)) / lag(value)) %>% 
    mutate(clean_value = case_when(is.na(prop_change) & units == "%" ~ value / 100, 
                                 is.na(prop_change) & units == "lateral shoot / shoot" ~ value, 
                                 is.na(prop_change) & units == "shoot / day" ~ value * std_time, 
                                 !is.na(prop_change) ~ prop_change)) %>% 
    ungroup() %>% 
    mutate(lambda = value / lag(value))

assign("seasonal_shoot_change_data", seasonal_shoot_change_data, envir = .GlobalEnv)

}

# ------------------------------------------------------------------------------
# Create Population Growth Dataframe
# ------------------------------------------------------------------------------

get_pop_growth_df <- function() {

pop_growth_df_before_net <-  
  bind_rows(mortality_data, shoot_change_data, seasonal_shoot_change_data) %>%
  # Set mortality value as negative growth
  mutate(pop_growth = ifelse(clean_response == "shoot mortality", -1 * clean_value, clean_value))

assign("pop_growth_df_before_net", pop_growth_df_before_net, envir = .GlobalEnv)

# Get Daylight Hours
# ------------------------------------------------------------------------------
# Intermission - get daylight hours - may want to tweak this at some point to 
# account for low sun angle in early/late hours
need_light_hours_studies <- 
  pop_growth_df_before_net %>%
    filter(`light hours` == "field conditions") %>% 
    distinct(studyid) %>% 
    .[['studyid']]

assign("need_light_hours_studies", need_light_hours_studies, envir = .GlobalEnv)

light_hours <- 
  get_study_dates(need_light_hours_studies) %>% 
    left_join(., driver_data) %>%
    filter(`light hours` == "field conditions") %>%
    distinct(studyid, site, treatment, start_date, measurement_date, response, `exposure time`, lat, lon) %>% 
    rowwise() %>% 
    mutate(x = list(getSunlightTimes(date = seq(from = ymd(start_date), to = ymd(measurement_date), by = "1 day"), 
               lat  = lat, 
               lon = lon,
               keep = c("sunrise", "sunset")))) %>%
    unnest(cols = c(x), names_repair = tidyr_legacy) %>%  
    select(-c(lat1, lon1)) %>%
    mutate(hours = sunset - sunrise) %>% 
    group_by(studyid, site, treatment, response, `exposure time`) %>% 
    summarise(suncalc_light_hours = as.numeric(mean(hours))) %>%
    ungroup() %>% 
    left_join(., driver_data) %>% 
    select(studyid, site, treatment, response, `exposure time`, lat, lon, 
       `penetrating irradiance`, pen_light, surface_light, mean_dli, suncalc_light_hours) %>%
    mutate(light_hours = suncalc_light_hours) %>% 
    mutate(suncalc_DLI = case_when(!is.na(pen_light) ~ pen_light * light_hours * 3600 / 1e6, 
                          (is.na(pen_light) & !is.na(mean_dli)) ~ mean_dli, 
                          (`penetrating irradiance` == "np" & is.na(mean_dli) & !is.na(surface_light)) ~ 
                           surface_light * light_hours * 3600 / 1e6)) %>% 
    select(-light_hours)

assign("light_hours", light_hours, envir = .GlobalEnv)

# See Brock 1981
# I4 = I5 / (2 * (1 + cos(2 * P1 * (T - 12) / L1))
# I4 = solar radiation per unit area per hour

# Use seasonal data environmental data
# ------------------------------------------------------------------------------
# Specifically this can/should(?) be used when plots are resampled over time 
# and continuous environmental data are provided and change seasonally
# e.g., studyid 147

seasonal_df <- 
  seasonal %>%
  select(studyid, site, treatment, start_date, measurement_date, seasonal_mean_dli, seasonal_mean_temp) %>% 
  mutate(start_date = strftime(start_date, "%Y-%m-%d", tz = 'UTC'), 
         measurement_date = strftime(measurement_date, "%Y-%m-%d", tz = 'UTC')) %>% 
  # Need to correct seasonal_mean_dli for %SI
  mutate(seasonal_mean_dli = case_when(studyid == 147 & str_detect(treatment, 'ambient') ~ seasonal_mean_dli * 0.33, 
                                            studyid == 147 & str_detect(treatment, 'low') ~ seasonal_mean_dli * 0.05, 
                                            studyid != 147 ~ seasonal_mean_dli))

assign("seasonal_df", seasonal_df, envir = .GlobalEnv)
# A quick check for what studies have 'duplicate' mortality/production values
# Whatever study shows up more than 1x here needs to be investigated
check_mort_prod_duplication <- 
  pop_growth_df_before_net %>%
    distinct(studyid, clean_response) %>% 
    arrange(studyid) %>% 
    group_by(studyid) %>% 
    summarise(n = n()) %>% 
    filter(n > 1)

assign("check_mort_prod_duplication", check_mort_prod_duplication, envir = .GlobalEnv)

# So far one study (and the way I collected data another one) has mortality and
# lateral shoot production values... so I need to get that into a net change
birth_and_death_provided <- 
  pop_growth_df_before_net %>%
    filter(studyid %in% check_mort_prod_duplication$studyid) %>% 
    select(studyid, site, treatment, clean_response, std_time, initial_density, `planting density (shoots / area)`, prop_change, pop_growth) %>% 
      group_by(studyid, site, treatment, std_time) %>% 
      pivot_wider(names_from = clean_response, values_from = pop_growth) %>%
      ungroup() %>% 
      mutate(net_growth = 
        case_when((!is.na(`shoot mortality`) & !is.na(`lateral shoot production`)) ~ `shoot mortality` + `lateral shoot production`, 
                  (!is.na(`shoot mortality`) &  is.na(`lateral shoot production`)) ~ `shoot mortality`, 
                  (is.na(`shoot mortality`) &  !is.na(`lateral shoot production`)) ~ `lateral shoot production`)
      ) %>%
      mutate(lambda = case_when(
        (!is.na(`shoot mortality`) & !is.na(`lateral shoot production`)) ~ ((1 + `shoot mortality`) * `lateral shoot production`), 
        (!is.na(`shoot mortality`) &  is.na(`lateral shoot production`)) ~ (1 + `shoot mortality`), 
        (is.na(`shoot mortality`) &  !is.na(`lateral shoot production`)) ~ `lateral shoot production`)
      ) %>%
      select(studyid, site, treatment, std_time, net_growth, lambda)

assign("birth_and_death_provided", birth_and_death_provided, envir = .GlobalEnv)

# Net growth after accounting for studies that may measure both mortality and 
# lateral shoot production
pop_growth_df <- 
  # Deal with light hours from field conditions
  pop_growth_df_before_net %>% 
  left_join(., light_hours) %>% 
  mutate(light_hours = ifelse(is.na(suncalc_light_hours), light_hours, suncalc_light_hours),
         DLI = ifelse(is.na(suncalc_DLI), DLI, suncalc_DLI)) %>% 
  # Deal with seasonal environmental data
  left_join(., seasonal_df) %>%
  mutate(DLI = ifelse(is.na(seasonal_mean_dli), DLI, seasonal_mean_dli), 
         clean_temp = ifelse(is.na(seasonal_mean_temp), clean_temp, seasonal_mean_temp)) %>%
  # Exclude initial values if they exist
  filter(std_time != 0) %>% 
  mutate(temp_bin = cut_interval(clean_temp, n = 3)) %>% 
  mutate(dli_bin  = case_when(DLI <= 5 ~ "DLI <=5", 
                              (DLI > 5 & DLI <= 10) ~ ">5 DLI <=10", 
                              DLI > 10 ~ "DLI > 10")) %>% 
  mutate(dli_bin = factor(dli_bin, levels = c("DLI <=5", ">5 DLI <=10", "DLI > 10"))) %>% 
  left_join(., birth_and_death_provided) %>% 
  mutate(pop_growth = ifelse(is.na(net_growth), pop_growth, net_growth)) %>%
  distinct(studyid, site, treatment, std_time, .keep_all = TRUE) %>% 
  mutate(pop_growth_rate = pop_growth / std_time) %>% 
  mutate(season = case_when(month(parsed_collection_date) %in% c(1:4, 11:12) ~ "winter", 
                            month(parsed_collection_date) %in% c(5:10) ~ "summer",
                            str_detect(collection_date, 'use summer') ~ "summer", 
                            str_detect(collection_date, 'seasonal data') ~ as.character(month(parse_date(start_date, "%Y-%m-%d")))
  )) %>% 
  mutate(season = case_when((str_detect(season, "winter|summer") | is.na(season)) ~ season, 
                             as.numeric(season) %in% c(1:4, 11:12) ~ "winter", 
                             as.numeric(season) %in% c(5:10) ~ "summer"))

assign("pop_growth_df", pop_growth_df, envir = .GlobalEnv)

}

get_mortality_data()

get_shoot_production_data()

get_pop_growth_df()

pop_growth_df %>% filter(is.na(studyid) | is.na(site) | is.na(treatment))

distinct(pop_growth_df, studyid) %>% arrange(studyid) %>% as.data.frame

write_csv(pop_growth_df, here('data_outputs', 'clean_pop_growth_df.csv'))