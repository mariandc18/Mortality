library(mgcv)
library(tidyverse)
library(jsonlite)

#### Set the seed for this script ####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(42)

#### Load in data ####
# Registros individuales de muertes
load("../Generated_Data/registros.Rda")

# Población desde JSON
pop_json <- fromJSON("../Generated_Data/population.json")
pop_data <- bind_rows(lapply(names(pop_json), function(yr) {
  data.frame(
    year       = as.integer(yr),
    sex        = c("hombres", "mujeres"),
    population = c(pop_json[[yr]]$hombres, pop_json[[yr]]$mujeres)
  )
}))

#### Agregar registros individuales a conteos mensuales (Opción A) ####
exp_obs_bymonthyear <- registros %>%
  group_by(year, month, sex, age_group) %>%
  summarise(observed = n(), .groups = "drop") %>%
  complete(year, month, sex, age_group, fill = list(observed = 0)) %>%
  filter(year < 2020) %>%
  left_join(pop_data, by = c("year", "sex")) %>%
  mutate(observed = floor(observed)) %>%
  arrange(year, month, sex, age_group)

#### Create data frame to store expected monthly mortality in 2020-2023 ####
pred_years       <- 2020:2023
pred_months      <- 1:12
sex_levels       <- unique(exp_obs_bymonthyear$sex)
age_group_levels <- unique(exp_obs_bymonthyear$age_group)

pred_grid <- expand.grid(
  year      = pred_years,
  month     = pred_months,
  sex       = sex_levels,
  age_group = age_group_levels,
  KEEP.OUT.ATTRS = FALSE
) %>%
  left_join(pop_data, by = c("year", "sex")) %>%
  mutate(
    expected_acm        = NA_real_,
    expected_acm_se     = NA_real_,
    expected_log_acm    = NA_real_,
    expected_log_acm_se = NA_real_,
    gamma_E             = NA_real_,
    gamma_delta         = NA_real_,
    gamma_E_nb          = NA_real_,
    gamma_delta_nb      = NA_real_
  )

exp_obs_bymonthyear <- exp_obs_bymonthyear %>%
  mutate(
    expected_acm        = NA_real_,
    expected_acm_se     = NA_real_,
    expected_log_acm    = NA_real_,
    expected_log_acm_se = NA_real_,
    gamma_E             = NA_real_,
    gamma_delta         = NA_real_
  )

num_samples <- 10000
n_pred      <- nrow(pred_grid)

#### Fit model ####
temp <- exp_obs_bymonthyear

if (length(unique(temp$year)) < 3) {
  annual_model <- gam(
    observed ~ year +
      s(month, bs = "cc", k = length(unique(temp$month))) +
      sex + age_group +
      offset(log(population)),
    data = temp, family = nb(theta = NULL, link = "log")
  )
} else {
  annual_model <- gam(
    observed ~ s(year, k = length(unique(temp$year))) +
      s(month, bs = "cc", k = length(unique(temp$month))) +
      sex + age_group +
      offset(log(population)),
    data = temp, family = nb(theta = NULL, link = "log")
  )
}
overd <- exp(annual_model$family$getTheta())

#### Get predictions for 2020-2023 ####
pred <- predict(annual_model,
                se.fit  = TRUE,
                type    = "response",
                newdata = pred_grid)
pred_grid$expected_acm    <- pred$fit
pred_grid$expected_acm_se <- pred$se.fit

pred_log <- predict(annual_model,
                    se.fit  = TRUE,
                    newdata = pred_grid)
pred_grid$expected_log_acm    <- pred_log$fit
pred_grid$expected_log_acm_se <- pred_log$se.fit

# Gamma parameters for each prediction row
gamma_E        <- numeric(n_pred)
gamma_delta    <- numeric(n_pred)
gamma_E_nb     <- numeric(n_pred)
gamma_delta_nb <- numeric(n_pred)

for (j in seq_len(n_pred)) {
  samples <- exp(rnorm(num_samples,
                       mean = pred_log$fit[j],
                       sd   = pred_log$se.fit[j]))
  gamma_E[j]        <- mean(samples)
  gamma_delta[j]    <- (gamma_E[j]^2) / var(samples)

  samples_nb         <- rnbinom(num_samples, size = overd, mu = samples)
  gamma_E_nb[j]      <- mean(samples_nb)
  gamma_delta_nb[j]  <- (gamma_E_nb[j]^2) / var(samples_nb)
}
pred_grid$gamma_E        <- gamma_E
pred_grid$gamma_delta    <- gamma_delta
pred_grid$gamma_E_nb     <- gamma_E_nb
pred_grid$gamma_delta_nb <- gamma_delta_nb

#### Fitted values for historical periods ####
pred_hist <- predict(annual_model, se.fit = TRUE, type = "response")
exp_obs_bymonthyear$expected_acm    <- pred_hist$fit
exp_obs_bymonthyear$expected_acm_se <- pred_hist$se.fit

pred_log_hist <- predict(annual_model, se.fit = TRUE)
exp_obs_bymonthyear$expected_log_acm    <- pred_log_hist$fit
exp_obs_bymonthyear$expected_log_acm_se <- pred_log_hist$se.fit

num_hist         <- nrow(exp_obs_bymonthyear)
gamma_E_hist     <- numeric(num_hist)
gamma_delta_hist <- numeric(num_hist)

for (j in seq_len(num_hist)) {
  samples              <- exp(rnorm(num_samples,
                                   mean = pred_log_hist$fit[j],
                                   sd   = pred_log_hist$se.fit[j]))
  gamma_E_hist[j]      <- mean(samples)
  gamma_delta_hist[j]  <- (gamma_E_hist[j]^2) / var(samples)
}
exp_obs_bymonthyear$gamma_E     <- gamma_E_hist
exp_obs_bymonthyear$gamma_delta <- gamma_delta_hist

#### Save all expecteds (historical + 2020-2023) ####
acm_monthly_predictions_hist <- exp_obs_bymonthyear %>%
  select(year, month, sex, age_group,
         expected_acm, expected_acm_se,
         expected_log_acm, expected_log_acm_se,
         gamma_E, gamma_delta) %>%
  bind_rows(
    pred_grid %>%
      select(year, month, sex, age_group,
             expected_acm, expected_acm_se,
             expected_log_acm, expected_log_acm_se,
             gamma_E, gamma_delta)
  ) %>%
  arrange(year, month, sex, age_group)

save(acm_monthly_predictions_hist,
     file = "../Generated_Data/acm_monthly_predictions_tier1_hist.RData")

#### Save expecteds for 2020-2023 only ####
acm_monthly_predictions_tier1 <- pred_grid %>%
  select(year, month, sex, age_group,
         expected_acm, expected_acm_se,
         expected_log_acm, expected_log_acm_se,
         gamma_E, gamma_delta, gamma_E_nb, gamma_delta_nb) %>%
  mutate(
    gamma_sd    = sqrt((gamma_E^2)    / gamma_delta),
    gamma_sd_nb = sqrt((gamma_E_nb^2) / gamma_delta_nb)
  )

save(acm_monthly_predictions_tier1,
     file = "../Generated_Data/acm_monthly_predictions_tier1.RData")