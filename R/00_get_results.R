library(here)

# Load functions used for data cleaning and analysis
source(here('R', '01_functions.R'))

# Load and clean the driver data
source(here('R', '02_driver-data-prep.R'))

# Combine and clean driver and response data to get clean data used in analysis
source(here('R', '03_get_clean_pop_growth_df.R'))

# A final clean of the data in preparation to run the GAMM
source(here('R', '04_prep_gam_data.R'))

# Fit the GAMM
source(here('R', '05_fit_gamm.R'))

# Generate main text and supplemental figures and get the critical values
source(here('R', '06_get_critical_values.R'))