# Reference cleaning functions
# ------------------------------------------------------------------------------
# Keywords plus are a pain in my butt. Do my own search to get Z. marina only
# Not used in the codebase here, buut this is how the keywords plus limitation was done
rem_kw_plus <- function(df) {
  df %>% 
  filter(str_detect(tolower(TI), c('Eelgrass|eelgrass|Zostera marina|zostera marina|Z\\. marina|z\\. marina')) | 
         str_detect(tolower(AB), c('Eelgrass|eelgrass|Zostera marina|zostera marina|Z\\. marina|z\\. marina')) | 
         str_detect(tolower(DE), c('Eelgrass|eelgrass|Zostera marina|zostera marina|Z\\. marina|z\\. marina')))
}

# Cleaning functions
# ------------------------------------------------------------------------------
clean_latlon <- function(df) {
  df %>% 
    mutate(lat_dms = str_detect(lat, '\\u00B0'),
           lon_dms = str_detect(lon, '\\u00B0')) %>% 
    mutate(lat = ifelse(lat_dms, dms2dec(lat), lat), 
           lon = ifelse(lon_dms, dms2dec(lon), lon)) %>% 
    select(-lat_dms, -lon_dms)
  }

# Parsing numbers fails if there is no number; This way you can also keep 
# non-number but non-NA values like 'np'.
if_num_parse_num <- function(x) {
  ifelse(grepl("\\d", x), parse_number(x), x)
} 

# Conversions
# ------------------------------------------------------------------------------
# Good discussion on str_extract/str_match and positive look behind
# https://stackoverflow.com/questions/35804379/stringr-str-extract-how-to-do-positive-lookbehind
dms2dec <- function(string) {
  string <- 
    str_squish(string) %>% 
    str_replace_all(., "\\s", "")

  values <- str_match(string, "(\\d*)[\\u00B0|˚]\\s?(\\d*\\.?\\d*?)('|′|’)\\s?(\\d*?\\.?\\d*?)[\"|″]?\\s?(N|n|E|e|S|s|W|w)")

  original <- values[, 1]
  degrees <- values[, 2]
  minutes <- (values[, 3])
  seconds <- values[, 5]
  direction <- values[, 6]

  degrees <- as.numeric(degrees)
  minutes <- ifelse(!is.na(as.numeric(minutes)), as.numeric(minutes), 0)
  seconds <- ifelse(!is.na(as.numeric(seconds)), as.numeric(seconds), 0)

  dms <- degrees + minutes / 60 + seconds / 3600

  dms <- ifelse(direction %in% c("N", "n", "E", "e"), 1 * dms, -1 * dms)

  return(dms)
}

# Convert photon flux density (m^-2 s-1) to daily light integral
ppfd_2_dli <- function(light_photon_m2_s1, light_hours) {
  light_photon_m2_s1 * light_hours * 3600 / 1e6
}

# ------------------------------------------------------------------------------
# Load response data
# ------------------------------------------------------------------------------

read_response_xlsx <- function(path, na = c("", "NA")) {
  response_df <- tryCatch(
    expr = {
      read_xlsx(path = path, na = na, col_types = "text")
    }, 
    error = function(e) {
      message("\nFile: ", path, " failed to load")
      message("\n    ", e)
      name_peak <- as_tibble(names(read_excel(path, n_max = 0)))
      return(name_peak)
    }
  )
  return(response_df)
}

load_response_data <- function(response = character()) {
  study_data <- read_xlsx(here('data', 'study-data.xlsx'), na = c("", "NA", "np"))
  raw_data_files <-
      tibble(folder = here('data', 'raw-data', 'driver-response-data'), 
             files = list.files(here('data', 'raw-data', 'driver-response-data'))
      ) %>% 
      unite(path, folder, files, sep = "", remove = FALSE)

  studies <- 
    study_data %>% 
      filter(str_detect(`raw data extracted`, {{response}})) %>% 
      .$studyid %>% 
      paste(., collapse = '|') %>% 
      paste0("^(", ., ")-")

  response_files <- 
    raw_data_files %>% 
      filter(., str_detect(files, "~\\$", negate = TRUE)) %>%
      filter(., str_detect(files, studies)) %>% 
      filter(., str_detect(files, ".R$|_raw$|Icon\\\r|^Icon", negate = TRUE)) %>% 
      distinct(path)
  
  response_data_list <- apply(response_files, 1, read_response_xlsx)

# Probably should omit failed files from read_response_xlsx at some point
  response_data_dfs <- 
    response_data_list %>% 
    bind_rows() %>% 
    rename(response = 'response (formal name)') %>% 
    mutate(studyid = as.numeric(studyid))

  return(response_data_dfs)
}

load_seasonal_data <- function(path = character()) {
  seasonal_data <- 
    read_csv(file = path, na = c("", "NA", "np")) %>% 
    mutate(`upper se` = as.character(`upper se`), 
            se = as.character(se))
  return(seasonal_data)
}


# Get study dates for calculating daylight hours
get_study_dates <- function(studyid) {
  studies <- 
    paste(studyid, collapse = '|') %>% 
    paste0("^(", ., ")-")
  raw_data_files <-
      tibble(folder = here('data', 'driver-response-data'), 
             files = list.files((here('data', 'driver-response-data')))
      ) %>% 
      unite(path, folder, files, sep = "", remove = FALSE) %>% 
      filter(., str_detect(files, "~\\$", negate = TRUE)) %>%
      filter(., str_detect(files, studies)) %>% 
      filter(., str_detect(files, ".R$|_raw$", negate = TRUE)) %>% 
      distinct(path)

  response_data <- 
    apply(raw_data_files, 1, read_response_xlsx) %>% 
    bind_rows() %>% 
    rename(response = 'response (formal name)') %>%
    distinct(studyid, site, treatment, start_date, measurement_date, response, `exposure time`) %>% 
    mutate(studyid = as.numeric(studyid))
  return(response_data)
}

# ------------------------------------------------------------------------------
# Plotting functions
# ------------------------------------------------------------------------------
mytheme <- function(base_size = 12, base_family="helvetica", axis_text_adj = 2, ...) {
  theme_minimal(base_size = base_size) + 
  theme(axis.line = element_line(colour="black"), 
        panel.grid = element_blank(), 
        axis.text.x = element_text(size = base_size - axis_text_adj), 
        axis.text.y = element_text(size = base_size - axis_text_adj), 
        axis.title.y = element_text(angle = 90, vjust = 0, margin = margin(r = 10, l = 0.5)),
        axis.title.x = element_text(vjust = -0.2, margin = margin(t = 2, b = 0.5)), 
        plot.margin = unit(c(2, 2, 2, 1), "mm"), 
        strip.background = element_blank()
        #strip.background = element_rect(colour = "white", fill = "transparent")
      ) + 
  theme(...)
}

# https://stackoverflow.com/questions/44628130/ggplot2-dealing-with-extremes-values-by-setting-a-continuous-color-scale
trim_tails <- function(range = c(-Inf, Inf), extended_breaks = 7) scales::trans_new("trim_tails", 
                transform = function(x) {
                  force(range)
                  desired_breaks <- scales::extended_breaks(n = extended_breaks)(x[x >= range[1] & x <= range[2]])
                  break_increment <- diff(desired_breaks)[1]
                  x[x < range[1]] <- range[1] - break_increment
                  x[x > range[2]] <- range[2] + break_increment
                  x
                },
                inverse = function(x) x,
#
                breaks = function(x) {
                  force(range)
                  scales::extended_breaks(n = extended_breaks)(x)
                },
                format = function(x) {
                  force(range)
                  x[1] <- paste("<", range[1])
                  x[length(x)] <- paste(">", range[2])
                  x
                })

# Add lines for IPCC warming conditions to surface plot
add_ipcc_lines <- function(df, colour = "goldenrod") {
  list(
    coord_cartesian(xlim = c(0, 30), clip = 'off'),
    scale_y_continuous(sec.axis = sec_axis(~ . * 1, name = "")), 
    theme(plot.margin = unit(c(0, 15, 0, 1), "mm")),
    geom_segment(data = df, aes(x = DLI,  xend = DLI, y = -Inf, yend = clean_temp + 5), 
                 colour = colour, size = 2, alpha = 0.8), 
    geom_segment(data = df, aes(x = 33, xend = DLI, y = clean_temp, yend = clean_temp), 
                 colour = colour, size = 2, alpha = 0.8), 
    geom_segment(data = df, aes(x = 33, xend = DLI, y = clean_temp + 1.5, yend = clean_temp + 1.5), 
                 colour = colour, size = 2, alpha = 0.8), 
    geom_segment(data = df, aes(x = 33, xend = DLI, y = clean_temp + 4.8, yend = clean_temp + 4.8), 
                 colour = colour, size = 2, alpha = 0.8), 
    geom_text(data = df, aes(x = 33.5, y = clean_temp), label = "Current", hjust = 0), 
    geom_text(data = df, aes(x = 33.5, y = clean_temp + 1.5), label = expression(1.5*degree*C), hjust = 0), 
    geom_text(data = df, aes(x = 33.5, y = clean_temp + 4.8), label = expression(4.8*degree*C), hjust = 0), 
    geom_point(data = df, aes(x = DLI, y = clean_temp, colour = pop_growth_rate), size = 4),
    geom_point(data = df, aes(x = DLI, y = clean_temp), colour = "black", size = 5, shape = 21)
  )
}

# Add lines for IPCC warming conditions to surface plot
add_ipcc_lines_manual <- function(df, colour = "goldenrod", annotation_size = 4, text_x = 0) {
  list(
    coord_cartesian(xlim = c(0, 30), clip = 'off'),
    geom_segment(data = df, aes(x = DLI,  xend = DLI, y = -Inf, yend = clean_temp), 
                 colour = colour, size = 1),  # Current vertical
    geom_segment(data = df, aes(x = 0, xend = DLI, y = clean_temp, yend = clean_temp), 
                 colour = colour, size = 1),  # Current horizontal
    geom_label(data = df, aes(x = text_x, y = clean_temp), label = "Current",   # Current label
               fill = "white", label.size = 0, label.padding = unit(0.15, "lines"), hjust = 0, size = annotation_size), 
    geom_label(data = df, aes(x = text_x, y = clean_temp + 1.5), label = expression(+1.5~degree*C),  # +1.5 label
               fill = "white", label.size = 0, label.padding = unit(0.15, "lines"), hjust = 0, size = annotation_size), 
    geom_label(data = df, aes(x = text_x, y = clean_temp + 4.8), label = expression(+4.8~degree*C),  # +4.8 label
               fill = "white", label.size = 0, label.padding = unit(0.15, "lines"), hjust = 0, size = annotation_size), 
    geom_point(data = df, aes(x = DLI, y = clean_temp), colour = "black", fill = "white", size = 3.5, shape = 21)
  )
}

# ------------------------------------------------------------------------------
# Statistical functions
# ------------------------------------------------------------------------------
rmvn <- function(n, mu, sig) { ## MVN random deviates
    L <- mroot(sig)
    m <- ncol(L)
    t(mu + L %*% matrix(rnorm(m*n), m, n))
}

# Check out this: https://stackoverflow.com/questions/58785930/r-find-maximum-of-density-plot
# for dealing with densities not matching raw mode calculations. 
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Get GAM fits with simultaneous CIs
get_gam_fits <- function(gam_fit, new_data, nsim = 1000, exclude = NULL) {
  # Extracts the Bayesian posterior covariance matrix of the parameters
  # `freq` = FALSE, `unconditional` = TRUE: Bayesian smoothing parameter uncertainty 
  # corrected covariance matrix is returned, if available.
  Vb <- vcov(gam_fit, freq = FALSE, unconditional = TRUE)   

  # Get predictions
  # "The standard errors produced by predict.gam are based on the Bayesian 
  #  posterior covariance matrix of the parameters Vp in the fitted gam object."
  pred <- predict(gam_fit, new_data, se.fit = TRUE, unconditional = TRUE, exclude = exclude_re)
  se.fit <- pred$se.fit

  # Next, we want nsim draws from [β̂ − β, û − u], which is approximately 
  # distributed multivariate normal with mean vector 0 and covariance matrix Vb
  BUdiff <- rmvn(nsim, mu = rep(0, nrow(Vb)), sig = Vb)


  # Now we calculate f̂ (x)−f(x), which is given by Cg[β̂ −βû −u] evaluated at the grid of x values
  # Cg is the evaluation of the basis functions at the locations g, and the thing in square brackets is the bias in the estimated model coefficients, which we assume to be mean 0 and follows, approximately, a multivariate normal distribution with mean vector 0 and covariance matrix Vb
  Cg <- predict(gam_fit, new_data, type = "lpmatrix", exclude = exclude)  # evaluates the basis function at g
  simDev <- Cg %*% t(BUdiff)  # computes the deviations between the fitted and true parameters

  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))  # find the absolute values of the standardized deviations from the true model

  masd <- apply(absDev, 2L, max)  # maximum of the absolute standardized deviations at the grid of x values for each simulation

  crit <- quantile(masd, prob = 0.95, type = 8)

  # I forgot that this contains point and simultaneous CIs to illustrate difference
  pred <- transform(cbind(data.frame(pred), new_data),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit)) %>% 
          as_tibble()

 return(pred)
}

get_GAM_draws <- function(gam_fit, new_data, nsim = 1000, exclude = NULL) {
  Vb <- vcov(gam_fit, freq = FALSE, unconditional = TRUE)  # Bayesian posterior CV
  sims <- rmvn(n = nsim, mu = coef(gam_fit), sig = Vb)  # mean vector given by the estimated model coefficients
  Cg <- predict(gam_fit, new_data, type = "lpmatrix", exclude = exclude)
  fits <- Cg %*% t(sims)  # contains `nsims` draws from the posterior
  
  return(fits)
}

stack_fits <- function(sim_fits, nrnd) {
  rnd <- sample(ncol(sim_fits), nrnd)

  stackFits <- stack(as.data.frame(sim_fits[, rnd]))
  stackFits <- 
    transform(stackFits,
              dli = rep(new_data$DLI, length(rnd)), 
              clean_temp = rep(new_data$clean_temp, length(rnd))) %>% 
    tibble()

  return(stackFits)
}

add_re_placeholders <- function(df) {
  df %>% 
  mutate(studyid = "127",
         site = "Seden Strand, Odense Fjiord, Denmark", 
         treatment = "T: 15",
         season = "summer") %>% 
  mutate(study_treatment = paste(studyid, treatment, sep = "_"), 
         study_site_treatment = factor(paste(studyid, site, treatment, sep = "_")))
}


# First derivative
# ------------------------------------------------------------------------------
shift_values <- function(df, column, h, FUN = `+`) {
      i <- grep(column, names(df))

      FUN <- match.fun(FUN)

      result <- df
      result[, i] <- FUN(result[, i], h)

      return(result)
}

central_diff1 <- function(model, newdata, x_axis, h) {
  ndf <- shift_values(df = newdata, column = x_axis, h = h / 2, FUN = `+`)
  ndb <- shift_values(df = newdata, column = x_axis, h = h / 2, FUN = `-`)
  
  Xf <- predict(gam_fit, ndf, type = "lpmatrix")#, exclude = exclude_re)
  Xb <- predict(gam_fit, ndb, type = "lpmatrix")#, exclude = exclude_re)
  Xdiff <- (Xf - Xb) / h

  return(Xdiff)
}

get_fderiv <- function(model, newdata, terms, term_columns, x_axis, h, nsim, exclude_re) {
  x_axis <- deparse(substitute(x_axis))

  Xdiff <- central_diff1(model, newdata, x_axis, h)

  Vb <- vcov(gam_fit, freq = FALSE, unconditional = TRUE)  # Bayesian posterior CV
  betas <- coef(gam_fit)  

  lpmat_col_ids <- grep(terms, colnames(Xdiff), fixed = TRUE)

  Xi <- Xdiff * 0  # zero out the Xp matrix

  Xi[, lpmat_col_ids] <- Xdiff[, lpmat_col_ids]  # copy bits of Xp we need

  fit <- predict(gam_fit, newdata, exclude = exclude_re)

  d <- drop(Xi %*% betas)  # drop some NULL dimension that shows up and get derivative
  se <- rowSums(Xi %*% Vb * Xi)^0.5

  d1_tbl <- tibble(smooth = rep(terms, length(d)), 
                newdata = newdata[term_columns],
                fit = fit,
                derivative = d, 
                se = se) %>% 
           unnest(c('newdata'))

  buDiff <- mvnfast::rmvn(n = nsim, mu = rep(0, nrow(Vb)), sigma = Vb)
  simDev <- tcrossprod(Xi, buDiff) # Xi %*% t(bu) # simulate deviations from expected
  absDev <- abs(sweep(simDev, 1L, d1_tbl[["se"]], FUN = "/")) # absolute deviations
  masd <- apply(absDev, 2L, max)  # & max abs deviation per sim
  ## simultaneous interval critical value
  crit <- quantile(masd, prob = level, type = 8)  # type 8 is recommended by Hyndman and Fan (1996)
  adj <- (crit * d1_tbl[["se"]])
  derivative <- add_column(d1_tbl,
                           crit  = rep(crit, nrow(d1_tbl)),
                           lower = d1_tbl[["derivative"]] - adj,
                           upper = d1_tbl[["derivative"]] + adj)#, 
                           #fit = predict(gam_fit, new_data, exclude = exclude_re))
  return(derivative)
}