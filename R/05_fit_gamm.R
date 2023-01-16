library(tidyverse)

#source(here('R', '01_functions.R'))
#source(here('R', '04_prep_gam_data.R'))

# Fit GAMM
# ------------------------------------------------------------------------------
# Cubic spline
fit_cr <- mgcv::gam(pop_growth_rate ~ te(DLI, clean_temp, bs="cr") +
            s(study_site_treatment, bs = "re"),
            data = pop_growth_df,
            method = "REML")
#gam.check(fit_cr)
summary(fit_cr)
draw(fit_cr) + ggtitle('pop_growth_rate ~ te(DLI, clean_temp, bs="cr"')

gam_fit = fit_cr
exclude_re = "s(study_site_treatment)"
saveRDS(gam_fit, file = here("data_outputs", "gam_fit.RDS"))


# gratia::appraise(gam_fit)
# gratia::draw(gam_fit)
# **Note:** fit of the GAMM is not great, with overdispersion of residuals at high positive values; high observed values are being underpredicted.
# The mean trend is likely on point but there are other sources of variability 
# (i.e., see discussion) that could be accounted for to improve model fit for 
# more applied situations. E.g., seasonality, different splines, number of knots, 
# random effects structures, ambient condition information, etc.,