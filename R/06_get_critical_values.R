library(viridis)
library(patchwork)
library(metR)  # does the contour labelling
library(ggrepel)
library(janitor)
library(flextable)
library(beepr)
library(ggtext)


# source(here('R', '05_fit_gamm.R'))

# Step 1: get simulated draws for combinations along DLI and Temp
# Step 2: Figure out how to get the mode
# Step 3: Get derivative of GAM
# Step 4: Combine x-intercept and sign of derivative at x value

# ------------------------------------------------------------------------------

# Set key parameters
# ------------------
gam_fit <- readRDS(gam_fit, file = here("data_outputs", "gam_fit.RDS"))

terms <- smooths(gam_fit)[[1]]  #"te(DLI,clean_temp)"
exclude_re <- smooths(gam_fit)[[2]] #"s(study_site_treatment)"
type <- "central"

nsim = 100
level = 0.95  # for 95% CIs if used (likely won't use them)
h <- 1e-7

temp  <- seq(5, 30, by = 0.1)
light <- seq(0, 30, by = 0.1) 

# Setup fitted values and posterior draws
pred_data <- 
  crossing(clean_temp = temp, DLI = light) %>% 
  mutate(study_site_treatment = factor(0))


# Get fitted values and posterior draws
# -------------------------------------
# Get dataframe of fitted values
fits <- fitted_values(data = pred_data, object = gam_fit, exclude = exclude_re, unconditional = TRUE)
beepr::beep()

# Get simulated draws over pred_data to get a visual idea of uncertainty
post_fit <- fitted_samples(gam_fit, seed = 1, n = nsim, newdata = pred_data, exclude = exclude_re)
beepr::beep()
fit_draws <- 
  pred_data %>% 
  mutate(row = row_number()) %>% 
  left_join(., post_fit)


# ---------------
# Get Derivatives
# ---------------
lpc_deriv <- get_fderiv(model = gam_fit, newdata = pred_data, terms = smooths(gam_fit)[[1]], 
                   term_columns = c('clean_temp', 'DLI'), x_axis = DLI, h = h, nsim = nsim, 
                   exclude_re = exclude_re)
beepr::beep()
tpc_deriv  <- get_fderiv(model = gam_fit, newdata = pred_data, terms = smooths(gam_fit)[[1]], 
                   term_columns = c('clean_temp', 'DLI'), x_axis = clean_temp, h = h, nsim = nsim, 
                   exclude_re = exclude_re)
beepr::beep()

# Visual check that derivatives are correct
lpc_deriv %>% 
  filter(clean_temp %in% c(5, 10, 15, 20, 25, 30)) %>% 
ggplot(data = .) + 
  geom_line(aes(x = DLI, y = derivative, colour = clean_temp, group = clean_temp)) +
  geom_line(aes(x = DLI, y = fit / 10, group = clean_temp), colour = "red") +
  geom_hline(yintercept = 0) + 
  mytheme() + 
  facet_wrap(~ clean_temp)

tpc_deriv %>% 
  filter(DLI %in% c(0, 5, 10, 15, 20, 25)) %>% 
ggplot(data = .) + 
  geom_line(aes(x = clean_temp, y = derivative, colour = DLI, group = DLI)) +
  geom_line(aes(x = clean_temp, y = fit / 10, group = DLI), colour = "red") +
  geom_hline(yintercept = 0) + 
  mytheme() + 
  facet_wrap(~ DLI)

# ------------------------------------------------------------------------------

# Plot Figure 1 - surface plot
# ----------------------------
# Reminder of what the truncated extreme tails are for plotting Figure 1
c(lower_tail, upper_tail)


# Figure 1a
# ----------------------------
gam_fit_plot <- 
  fits %>% 
    ggplot(data = .) + 
      geom_raster(aes(x = DLI, y = clean_temp, fill = fitted)) + 
      geom_contour(aes(x = DLI, y = clean_temp, z = fitted), colour = "grey75") + 
      geom_jitter(data = pop_growth_df, width = 0.5, height = 0.5, alpha = 0.8,
                 aes(x = DLI, y = clean_temp, fill = pop_growth_rate_truncated, shape = surface_pen_light,
                     text = paste('Studyid: ', studyid, '</br>', 
                                  'Treat: ', treatment, '</br>',
                                  'Rate: ', round(pop_growth_rate, digits = 3), '</br>',
                                  'Time: ', std_time)), 
                 shape = 21, 
                 size = 4, colour = "black") + 
      geom_contour(aes(x = DLI, y = clean_temp, z = fitted, colour = ..level..), breaks = 0, colour = "black", size = 1) + 
      scale_shape_manual(values = c(21, 24)) + 
      colorspace::scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, limits = c(-0.1, 0.1), 
        breaks = c(-0.1, -0.05, 0, 0.05, 0.1),
        labels = c("≤ -0.1", "-0.05", "0", "0.05", "≥ 0.1")) +
      colorspace::scale_colour_continuous_divergingx(palette = 'RdBu', mid = 0, limits = c(-0.1, 0.1), 
        breaks = c(-0.1, -0.05, 0, 0.05, 0.1),
        labels = c("≤ -0.1", "-0.05", "0", "0.05", "≥ 0.1")) + 
      geom_text_contour(aes(x = DLI, y = clean_temp, z = fitted), skip = 1, colour = "grey25",
        stroke = 0.1, stroke.colour = "white", check_overlap = FALSE, label.placer = label_placer_flattest(ref_angle = 0)) + 
      coord_equal() + 
      xlab(expression(Daily~light~integral~(mol~m^-2~day^-1))) + 
      labs(y = "Temperature (˚C)", 
           fill = "r (day-1)") + 
      geom_rug(data = pop_growth_df, aes(x = DLI, y = clean_temp), position = "jitter",
               outside = FALSE, length = unit(0.01, "npc"), alpha = 0.5) + 
      scale_x_continuous(breaks = scales::pretty_breaks(n = 6), oob = scales::squish) + 
      scale_y_continuous(breaks = scales::pretty_breaks(n = 6), oob = scales::squish) + 
      mytheme() + 
      guides(colour = "none") + 
      theme(legend.position = "bottom") + 
      guides(fill = guide_colorbar(title = expression(r~(day^-1)),
                               label.position = "bottom", label.vjust = -0.5,
                               title.position = "left", title.vjust = 1.25, 
                               barwidth = 15,
                               barheight = 0.8))
set.seed(18)  # set seed to to improve have a jitter that look decent
gam_fit_plot
  #geom_point(data = pop_growth_df, aes(x = DLI, y = clean_temp, shape = surface_pen_light)) +
  #geom_text(data = pop_growth_df, aes(x = DLI, y = clean_temp, label = studyid))

# Without trimming tails for colour scale (not shown in MS)
gam_fit_plot + 
colorspace::scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, limits = c(-0.1, 0.1)) +
    colorspace::scale_colour_continuous_divergingx(palette = 'RdBu', mid = 0, limits = c(-0.1, 0.1)) 

#plotly::ggplotly(gam_fit_plot)
set.seed(18)
gam_fit_plot

ggsave('figures/Figure1a_response-surface.png', width = 16.9 / cm(1), height = 16.9 / cm(1))

# ----------------------------
# Figure 1b + 1c
# ----------------------------
dev.new(width = 16.9 / cm(1), height = 22 / cm(1))

ipcc_df <- tibble(DLI = 4.25, clean_temp = 19, pop_growth_rate = 0, text_x = 0)
max_dli <- 15
m_pink <- "#FF0076"
d_pink <- "#CC005E"
m_blue <- "#15ACE7"
m_green <- "#9ade16"
m_yellow <- "#d8b717"
fig1b <- 
  fits %>% 
    filter(DLI <= max_dli) %>% 
    ggplot(data = .) + 
      geom_raster(aes(x = DLI, y = clean_temp, fill = fitted)) + 
      geom_contour(aes(x = DLI, y = clean_temp, z = fitted), colour = "grey75") + 
      geom_contour(aes(x = DLI, y = clean_temp, z = fitted, colour = ..level..), breaks = 0, colour = "black", size = 1) + 
      colorspace::scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, limits = c(-0.1, 0.1), 
        breaks = c(-0.1, -0.05, 0, 0.05, 0.1),
        labels = c("≤ -0.1", "-0.05", "0", "0.05", "≥ 0.1")) +
      colorspace::scale_colour_continuous_divergingx(palette = 'RdBu', mid = 0, limits = c(-0.1, 0.1), 
        breaks = c(-0.1, -0.05, 0, 0.05, 0.1),
        labels = c("≤ -0.1", "-0.05", "0", "0.05", "≥ 0.1")) +
      mytheme() + 
      # 1.5 shift y-direction
      geom_segment(data = ipcc_df, aes(x = 1, xend = DLI + 1.55, y = clean_temp + 1.5, yend = clean_temp + 1.5), 
                   colour = m_pink, size = 1) +  # 1.5 horizontal shift
      geom_segment(data = ipcc_df, aes(x = (DLI / 2) + 0.75,  xend = (DLI / 2) + 0.75, y = clean_temp, yend = clean_temp + 1.5), 
                   colour = m_pink, linetype = "11", size = 1) +  # 1.5 vertical arrow
      geom_segment(data = ipcc_df, aes(x = (DLI / 2) + 0.75,  xend = (DLI / 2) + 0.75, y = clean_temp + 1, yend = clean_temp + 1.5), 
                   colour = m_pink, size = 0.5, arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +  # 1.5 vertical arrow
      # 1.5 shift x-direction
      geom_segment(data = ipcc_df, aes(x = DLI + 1.55,  xend = DLI + 1.55, y = -Inf, yend = clean_temp + 1.5), 
                   colour = m_pink, linetype = "solid", size = 1) +  # 1.5 vertical shift
      geom_segment(data = ipcc_df, aes(x = DLI, xend = DLI + 1.55, y = (clean_temp / 2) + 3, yend = (clean_temp / 2) + 3), 
                   colour = m_pink, size = 1, linetype = "11") +  # horizontal dotted
      geom_segment(data = ipcc_df, aes(x = DLI + 1.5, xend = DLI + 1.55, y = (clean_temp / 2) + 3, yend = (clean_temp / 2) + 3), 
                   colour = m_pink, size = 1, arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +  # horizontal arrowhead
      # 4.8 shift y-direction
      geom_segment(data = ipcc_df, aes(x = 1, xend = DLI + 4.7, y = clean_temp + 4.8, yend = clean_temp + 4.8), 
                   colour = d_pink, size = 1) +
      geom_segment(data = ipcc_df, aes(x = (DLI / 2) + 1.5,  xend = (DLI / 2) + 1.5, y = clean_temp, yend = clean_temp + 4.8), 
                   colour = d_pink, linetype = "21", size = 1) +  # 4.8 vertical
      geom_segment(data = ipcc_df, aes(x = (DLI / 2) + 1.5,  xend = (DLI / 2) + 1.5, y = clean_temp + 4.3, yend = clean_temp + 4.8), 
                   colour = d_pink, size = 1, arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +  # 4.8 vertical arrowhead
      # 4.8 shift x-direction
      geom_segment(data = ipcc_df, aes(x = DLI + 4.7,  xend = DLI + 4.7, y = -Inf, yend = clean_temp + 4.8), 
                   colour = d_pink, size = 1) +  # 4.8 vertical shift
      geom_segment(data = ipcc_df, aes(x = DLI, xend = DLI + 4.7, y = (clean_temp / 2) + 4, yend = (clean_temp / 2) + 4), 
                   colour = d_pink, size = 1, linetype = "21",
                   arrow = arrow(length = unit(0.03, "npc"))) +   # arrow line
      geom_segment(data = ipcc_df, aes(x = DLI + 4.5, xend = DLI + 4.7, y = (clean_temp / 2) + 4, yend = (clean_temp / 2) + 4), 
                   colour = d_pink, size = 1, arrow = arrow(length = unit(0.03, "npc"), type = "closed")) +  # arrowhead
      xlab(expression(Daily~light~integral~(mol~m^-2~day^-1))) + 
      add_ipcc_lines_manual(df = ipcc_df, colour = m_blue, text_x = 0, annotation_size = 3.5) +
      labs(y = "Temperature (˚C)") + 
      coord_cartesian(xlim = c(0, max_dli), clip = 'off', expand = 0) + 
      guides(fill = guide_colorbar(title = expression(r~(day^-1)),
                           label.position = "right", label.vjust = -0.5,
                           title.position = "top", title.vjust = 5, 
                           barwidth = 0.8,
                           barheight = 10), 
             colour = "none") + 
      theme(legend.position = "right") #+ 
fig1b
    # geom_point(data = ipcc_df, aes(x = DLI, y = clean_temp + 1.5), colour = "black", fill = "white", size = 3.5, shape = 21) + 
    # geom_point(data = ipcc_df, aes(x = DLI, y = clean_temp + 4.8), colour = "black", fill = "white", size = 3.5, shape = 21) + 
    # geom_point(data = ipcc_df, aes(x = DLI + 1, y = clean_temp + 1.5), colour = "black", fill = "white", size = 3.5, shape = 21) + 
    # geom_point(data = ipcc_df, aes(x = DLI + 3.5, y = clean_temp + 4.8), colour = "black", fill = "white", size = 3.5, shape = 21)

fig1c <- 
  gam_fit_plot + 
  guides(fill = "none") + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) + 
  coord_cartesian(xlim = c(0, 30), clip = 'off')

layout <- "
  AA
  BC
"

(fig1c + fig1b + guide_area()) + plot_layout(design = layout, heights = c(1, 0.6), widths = c(1, 0.4), guides = "collect") + 
plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold', vjust = -3, hjust = 0), legend.position = "right")


ggsave(here('figures', 'Figure1_surface-with-ipcc-management-shifts.png'), width = 16.9 / cm(1), height = 22 / cm(1))

# ------------------------------------------------------------------------------
# Get critical values (Figure 2 + Figure 3)
# ------------------------------------------------------------------------------
# Set useful bins
temp_bin <- c(10, 14, 16, 20, 21, 23, 25, 27)
dli_bin  <- c(3, 4, 5, 6, 10, 15, 20, 25)

# ------------------------------------------------------------------------------
# R = R_max
# ------------------------------------------------------------------------------

# Get mean critical value for DLI_opt
# --------------------------
opt_dli <- 
  lpc_deriv %>% 
    group_by(clean_temp) %>% 
    arrange(clean_temp) %>% 
    mutate(forward = lead(fit), backward = lag(fit), max_fit = max(fit)) %>% 
    filter(abs(derivative) < 1e-4) %>% 
    filter(clean_temp >= 8) %>% 
    filter(DLI <= 26) %>% 
    filter(fit > 0) %>% 
    filter(fit == max_fit) %>% 
    arrange(DLI, clean_temp) %>% 
    summarise(mean_dli_opt = mean(DLI))

opt_dli_draws <- 
  fit_draws %>%
    left_join(., lpc_deriv) %>%  
      group_by(clean_temp, draw) %>% 
      mutate(max_fit = max(fitted)) %>% 
      filter(clean_temp >= 8) %>% 
      filter(DLI <= 26) %>% 
      filter(fit > 0) %>% 
      filter(fitted == max_fit) %>% 
      summarise(mean_dli_opt = mean(DLI)) %>%  # get mean because sometimes multiple values
      ungroup()
beepr::beep()

opt_dli_se <- 
  opt_dli_draws %>% 
  group_by(clean_temp) %>% 
  summarise(mode_dli_opt = get_mode(mean_dli_opt), 
            mean_mean_dli_opt = mean(mean_dli_opt),
            lower = quantile(mean_dli_opt, prob = 0.025), 
            upper = quantile(mean_dli_opt, prob = 0.975))

opt_dli_plot <- 
  opt_dli %>%
  ggplot(data = .) + 
    geom_ribbon(data = opt_dli_se, aes(x = clean_temp, ymin = lower, ymax = upper), fill = "grey", alpha = 0.5) +
    geom_path(aes(x = clean_temp, y = mean_dli_opt, colour = clean_temp), size = 1) + 
    mytheme() + 
    theme(axis.text.y = element_markdown()) +
    labs(x = expression(Temperature~(degree*C)), y = expression(DLI[r==rmax]~(mol~m^-2~day^-1))) + 
    scale_colour_viridis(aesthetics = c("colour"), lim = c(8, 28), option = "plasma") + 
    lims(x = c(8, 28), y = c(0, 23)) +
    guides(colour = "none")

T_bin_DLI_opt  <- c(8, 10, 12, 14, 16, 18, 20, 22)
shifting_opt_dli <- 
  fit_draws %>% 
    filter(clean_temp %in% T_bin_DLI_opt) %>% 
  ggplot(data = .) + 
    mytheme() + 
    geom_ribbon(data = fits %>% filter(clean_temp %in% T_bin_DLI_opt), aes(x = DLI, ymin = lower, ymax = upper), fill = "grey", alpha = 0.3) +
    geom_ribbon(data = fits %>% filter(clean_temp %in% T_bin_DLI_opt), aes(x = DLI, ymin = lower, ymax = upper, fill = clean_temp), alpha = 0.2) +
    geom_line(data = fits %>% filter(clean_temp %in% T_bin_DLI_opt), aes(x = DLI, y = fitted, colour = clean_temp), size = 1) +
    geom_hline(yintercept = 0, colour = "black") + 
    geom_vline(data = opt_dli %>% filter(clean_temp %in% T_bin_DLI_opt), aes(xintercept = mean_dli_opt), linetype = 'dashed') + 
    facet_wrap(~ clean_temp, ncol = 4) + 
    labs(y = expression(Population~growth~rate~(d^-1)), x = expression(Daily~light~integral~(mol~m^-2~day^-1))) + 
    ylim(c(-0.065, 0.065)) + 
    scale_colour_viridis(aesthetics = c("fill", "colour"), lim = c(8, 28), option = "plasma") + 
    guides(colour = "none", fill = "none")
shifting_opt_dli

((shifting_opt_dli / plot_spacer()) / opt_dli_plot) + plot_layout(heights = c(2.1, 0.1, 1.1)) +
plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

dev.new(width = 16.9 / cm(1), height = 13 / cm(1))
shifting_opt_dli
ggsave(filename = here('figures', 'SOM', 'FigureS2_shifting-opt-dli-temp>=8-dli<=26.png'), width = 16.9 / cm(1), height = 13 / cm(1))


# --------------------------
# Get mean critical value for minimum light
# --------------------------
min_light <- 
  lpc_deriv %>% 
    group_by(clean_temp) %>% 
    mutate(forward = lead(fit), backward = lag(fit)) %>%
    filter(clean_temp >= 8) %>% 
    filter(DLI <= 26) %>% 
    filter(derivative > 0) %>% 
    filter(forward * backward <= 0 & forward >= backward) %>% 
    arrange(DLI) %>%  
    summarise(mean_dli_min = mean(DLI))

min_light_draws <- 
  fit_draws %>%
    left_join(., lpc_deriv) %>%  
      group_by(clean_temp, draw) %>% 
      mutate(forward = lead(fitted), backward = lag(fitted)) %>% 
      filter(clean_temp >= 8) %>% 
      filter(DLI <= 26) %>% 
      filter(derivative > 0) %>%   # using mean derivative reduces messed up mean(mean_dli_min); but mode is fine regardless
      filter((forward * backward <= 0) & (forward >= backward)) %>% 
      summarise(mean_dli_min = mean(DLI)) %>%  # get mean because sometimes multiple values
      ungroup()
beepr::beep()

min_light_se <- 
  min_light_draws %>% 
  group_by(clean_temp) %>% 
  summarise(mode_dli_min = get_mode(mean_dli_min),  # reminder to self that these were checks to see how they compared to the mean fit, just using the mean estimated performance curve was good for the critical value. These draw shenanigans were ultimately to get some kind of uncertainty estimate.
            mean_mean_dli_min = mean(mean_dli_min),
            lower = quantile(mean_dli_min, prob = 0.025), 
            upper = quantile(mean_dli_min, prob = 0.975))

min_light_temp_plot <-   
  min_light %>% 
  ggplot(data = .) + 
    geom_ribbon(data = min_light_se, aes(x = clean_temp, ymin = lower, ymax = upper), fill = "grey", alpha = 0.5) +
    geom_path(aes(x = clean_temp, y = mean_dli_min, colour = clean_temp), size = 1) + 
    mytheme() + 
    labs(x = expression(Temperature~(degree*C)), y = expression(DLI[r==0]~(mol~m^-2~day^-1)), colour = expression(Temp~(degree*C))) +
    scale_colour_viridis(option = "plasma") + 
    guides(colour = "none") + 
    lims(x = c(8, 28), y = c(0, 23)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))
    #guides(colour = guide_colourbar(barwidth = 0.8, title.vjust = 3, title.hjust = 0))

shifting_min_light <- 
  fit_draws %>% 
    filter(clean_temp %in% temp_bin) %>% 
  ggplot(data = .) + 
    mytheme() + 
    geom_ribbon(data = fits %>% filter(clean_temp %in% temp_bin), aes(x = DLI, ymin = lower, ymax = upper), fill = "grey", alpha = 0.6) +
    geom_ribbon(data = fits %>% filter(clean_temp %in% temp_bin), aes(x = DLI, ymin = lower, ymax = upper, fill = clean_temp), alpha = 0.1) +
    geom_line(data = fits %>% filter(clean_temp %in% temp_bin), aes(x = DLI, y = fitted, colour = clean_temp), size = 1) +
    geom_hline(yintercept = 0, colour = "black") + 
    geom_vline(data = min_light %>% filter(clean_temp %in% temp_bin), aes(xintercept = mean_dli_min), linetype = 'dashed') + 
    facet_wrap(~ clean_temp, ncol = 4) + 
    labs(y = expression(Population~growth~rate~(d^-1)), x = expression(Daily~light~integral~(mol~m^-2~day^-1))) + 
    ylim(c(-0.075, 0.075)) + 
    scale_colour_viridis(option = "plasma", aesthetics = c("fill", "colour")) + 
    guides(colour = "none", fill = "none")

(shifting_min_light / plot_spacer() / (min_light_temp_plot | opt_dli_plot)) + plot_layout(heights = c(2, 0.1, 1.1)) + 
plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here('figures', 'Figure2_min-light.png', width = 16.9 / cm(1), height = 18 / cm(1))


# T_crits
# ------------------------------------------------------------------------------

# --------------------------
# Get mean critical value for T_opt
# --------------------------
opt_temp <- 
  tpc_deriv %>% 
    group_by(DLI) %>% 
    arrange(DLI) %>% 
    mutate(forward = lead(fit), backward = lag(fit), max_fit = max(fit)) %>% 
    filter(abs(derivative) < 1e-4) %>% 
    filter(clean_temp >= 8) %>% 
    filter(DLI <= 26) %>% 
    filter(fit > 0) %>% 
    filter(fit == max_fit) %>% 
    arrange(clean_temp, DLI) %>% 
    summarise(mean_temp_opt = mean(clean_temp), 
              median_temp_opt = median(clean_temp))

opt_temp_draws <- 
  fit_draws %>%
    left_join(., tpc_deriv) %>%  
      group_by(DLI, draw) %>% 
      mutate(max_fit = max(fitted)) %>% 
      filter(clean_temp >= 8) %>% 
      filter(DLI <= 26) %>% 
      filter(fit > 0) %>% 
      filter(fitted == max_fit) %>% 
      summarise(mean_temp_opt = mean(clean_temp), 
                median_temp_opt = median(clean_temp)) %>%  # get mean because sometimes multiple values
      ungroup()
beepr::beep()

opt_temp_se <- 
  opt_temp_draws %>% 
  group_by(DLI) %>% 
  summarise(mode_temp_opt = get_mode(mean_temp_opt), 
            mean_mean_temp_opt = mean(mean_temp_opt),
            mean_median_temp_opt = mean(median_temp_opt),
            lower = quantile(mean_temp_opt, prob = 0.025), 
            upper = quantile(mean_temp_opt, prob = 0.975),
            lower_median = quantile(median_temp_opt, prob = 0.025), 
            upper_median = quantile(median_temp_opt, prob = 0.975))

opt_temp_plot <- 
  opt_temp %>%
  ggplot(data = .) + 
    geom_ribbon(data = opt_temp_se, aes(x = DLI, ymin = lower, ymax = upper), fill = "grey", alpha = 0.5) +
    geom_path(aes(x = DLI, y = mean_temp_opt, colour = DLI), size = 1) + 
    mytheme() + 
    labs(x = expression(Daily~light~integral~(mol~m^-2~day^-1)), y = expression(T[r==rmax]~(degree*C))) + 
    scale_colour_viridis(aesthetics = c("colour"), lim = c(0, 30)) + 
    xlim(c(3.7, 26)) +
    ylim(c(8, 30)) +
    guides(colour = "none")

dli_bin_T_opt  <- c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26)
shifting_opt_temp <- 
  fit_draws %>% 
    filter(DLI %in% dli_bin_T_opt) %>% 
  ggplot(data = .) + 
    mytheme() + 
    geom_ribbon(data = fits %>% filter(DLI %in% dli_bin_T_opt), aes(x = clean_temp, ymin = lower, ymax = upper), fill = "grey", alpha = 0.3) +
    geom_ribbon(data = fits %>% filter(DLI %in% dli_bin_T_opt), aes(x = clean_temp, ymin = lower, ymax = upper, fill = DLI), alpha = 0.2) +
    geom_line(data = fits %>% filter(DLI %in% dli_bin_T_opt), aes(x = clean_temp, y = fitted, colour = DLI), size = 1) +
    geom_hline(yintercept = 0, colour = "black") + 
    geom_vline(data = opt_temp %>% filter(DLI %in% dli_bin_T_opt), aes(xintercept = mean_temp_opt), linetype = 'dashed') + 
    facet_wrap(~ DLI, ncol = 4) + 
    labs(y = expression(Population~growth~rate~(d^-1)), x = "Temperature (˚C)") + 
    ylim(c(-0.075, 0.075)) + 
    scale_colour_viridis(aesthetics = c("fill", "colour"), lim = c(0, 30)) + 
    guides(colour = "none", fill = "none")
shifting_opt_temp


dev.new(width = 16.9 / cm(1), height = 13 / cm(1))
shifting_opt_temp
ggsave(filename = here('figures', 'SOM', 'FigureS3_shifting-opt-temp-temp>=8-dli<=26.png'), width = 16.9 / cm(1), height = 13 / cm(1))


# --------------------------
# Get mean critical value for max temperature
# --------------------------
max_temp <- 
  tpc_deriv %>% 
    group_by(DLI) %>% 
    mutate(forward = lead(fit), backward = lag(fit)) %>% 
    filter(derivative < 0) %>% 
    filter((forward * backward <= 0) & (forward <= backward)) %>% 
    filter(clean_temp >= 8) %>% 
    filter(DLI <= 26) %>% 
    arrange(clean_temp, DLI) %>% 
    summarise(mean_temp_max = mean(clean_temp)) %>% 
    filter()

max_temp_draws <- 
  fit_draws %>%
    left_join(., tpc_deriv) %>%  
      group_by(DLI, draw) %>% 
      mutate(forward = lead(fitted), backward = lag(fitted)) %>% 
      filter(derivative < 0) %>%   # using mean derivative reduces messed up mean(mean_dli_min); but mode is fine regardless
      filter((forward * backward <= 0) & (forward <= backward)) %>% 
      filter(clean_temp >= 8) %>% 
      filter(DLI <= 26) %>% 
      summarise(mean_temp_max = mean(clean_temp)) %>%  # get mean because sometimes multiple values
      ungroup()
beepr::beep()

max_temp_se <- 
  max_temp_draws %>% 
  group_by(DLI) %>% 
  summarise(mode_temp_max = get_mode(mean_temp_max), 
            mean_mean_temp_max = mean(mean_temp_max),
            lower = quantile(mean_temp_max, prob = 0.025), 
            upper = quantile(mean_temp_max, prob = 0.975))

max_temp_plot <- 
  max_temp %>%  
  ggplot(data = .) + 
    geom_ribbon(data = max_temp_se %>% filter(DLI > 3.7), aes(x = DLI, ymin = lower, ymax = upper), fill = "grey", alpha = 0.5) +
    geom_path(aes(x = DLI, y = mean_temp_max, colour = DLI), size = 1.5) + 
    mytheme() + 
    labs(x = expression(Daily~light~integral~(mol~m^-2~day^-1)), y = expression(T[r==0]~(degree*C))) + 
    scale_colour_viridis(aesthetics = c("colour"), lim = c(0, 30)) +
    guides(colour = "none") + 
    ylim(c(8, 30))

dli_bin  <- c(3, 5, 7, 9, 12, 16, 20, 24)
shifting_max_temp <- 
  fit_draws %>% 
    filter(DLI %in% dli_bin) %>% 
  ggplot(data = .) + 
    mytheme() + 
    geom_ribbon(data = fits %>% filter(DLI %in% dli_bin), aes(x = clean_temp, ymin = lower, ymax = upper, fill = DLI), alpha = 0.2) +
    geom_line(data = fits %>% filter(DLI %in% dli_bin), aes(x = clean_temp, y = fitted, colour = DLI), size = 1) +
    geom_hline(yintercept = 0, colour = "black") + 
    geom_vline(data = max_temp %>% filter(DLI %in% dli_bin), aes(xintercept = mean_temp_max), linetype = 'dashed') + 
    facet_wrap(~ DLI, ncol = 4) + 
    labs(y = expression(Population~growth~rate~(d^-1)), x = expression(Temperature~(degree*C))) + 
    ylim(c(-0.065, 0.065)) + 
    scale_colour_viridis(aesthetics = c("fill", "colour"), lim = c(0, 30)) + 
    guides(colour = "none", fill = "none")

(shifting_max_temp / plot_spacer() / (max_temp_plot | opt_temp_plot)) + plot_layout(heights = c(2.1, 0.1, 1.1)) + 
plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(face = 'bold'))

ggsave(filename = here('figures', 'Figure3_max-temp.png'), width = 16.9 / cm(1), height = 20 / cm(1))




# Threshold values for contour plot
# ------------------------------------------------------------------------------
# Get range of DLI_0
message(cat("DLI_0 range:", min(min_light$mean_dli_min), "-", max(min_light$mean_dli_min)))

# Get range of T_max
message(cat("T_max range:", min(max_temp$mean_temp_max), "-", max(max_temp$mean_temp_max)))

# Get optimum growth rate combo:
message(cat("r_opt temp:", slice(fits, which.max(fitted))[["clean_temp"]], "DLI:", slice(fits, which.max(fitted))[["DLI"]]))

critical_values <- list(DLI_0_min = min(min_light$mean_dli_min), T_max_max = max(max_temp$mean_temp_max), 
                        DLI_opt = slice(fits, which.max(fitted))[["DLI"]], 
                        T_opt   = slice(fits, which.max(fitted))[["clean_temp"]])
saveRDS(critical_values, file = here('data_outputs', 'critical_values.RDS'))


# Markdown to print out these values
# population growth occurred (Quick reference for adding to text; autogenerated values referenced in the text will update here.): 

#    + above light levels of `r critical_values$DLI_0_min` mol m^-2^ day^-1^ 
#    + below temperatures of `r critical_values$T_max_max`˚C. 

# Maximum shoot production occurred at:

#    + light levels of `r critical_values$DLI_opt` mol m^-2^ day^-1^ 
#    + temperatures of `r critical_values$T_opt` ˚C. 

# Model used:    

#    + `r gam_fit_formula`
#    + `r gam_fit_data`
#    + `r gam_fit_method` 

# ------------------------------------------------------------------------------



# Do some visual investigations
# --------------------------
# Zoom in on low temp curve back on DLI minimum value;
# It looks like it's just because of the tails of the data;
# Just seems constrained to some choke points causing some funky intercepts
lpc_deriv %>% 
  filter(clean_temp %in% 5:15) %>% 
ggplot(data = .) + 
  geom_line(aes(x = DLI, y = fit, colour = clean_temp, group = clean_temp)) + xlim(c(0, 7)) + ylim(c(-0.01, 0.01)) + 
  geom_hline(yintercept = 0)

lpc_deriv %>% 
  filter(clean_temp %in% 5:10) %>% 
ggplot(data = .) + 
  geom_line(aes(x = DLI, y = fit, colour = clean_temp, group = clean_temp)) + 
  geom_hline(yintercept = 0)

lpc_deriv %>% 
  filter(clean_temp %in% 25:30) %>% 
ggplot(data = .) + 
  geom_line(aes(x = DLI, y = fit, colour = clean_temp, group = clean_temp)) + 
  geom_hline(yintercept = 0)
# --------------------------
