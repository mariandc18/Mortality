library(mgcv)
library(tidyverse)
library(jsonlite)
library(readr)

# Set seed 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(42)

# Load data for death registration and population
death_records <- read_csv("./mortalidad_filtrada.csv")
pop_json <- fromJSON("./population_by_gender.json")
pop_data <- bind_rows(lapply(names(pop_json), function(yr) {
  bind_rows(lapply(c("hombres", "mujeres"), function(sx) {
    data.frame(
      year       = as.integer(yr),
      sex        = sx,
      age_group  = 1:6,         
      population = unlist(unname(pop_json[[yr]][[sx]]))
    )
  }))
}))

pop_data <- pop_data %>%
  mutate(age_group = factor(age_group, levels = 1:6))

print(head(pop_data))

codes <- c(
  'J13X', 'J15', 'J154', 'J158', 'J159', 'J17', 'J170', 'J178', 'J18', 'J182', 'J188', 'J189', 'J440', 'J850', 'J851',
  'G001', 'G009', 'A403', 'A409', 'A419')

cause_columns <- paste0("CAUSA", 1:10)

# Filter records
filtered_records <- death_records %>%
  filter(
    if_any(all_of(cause_columns), ~ . %in% codes)
  ) %>%
  mutate(
    sex = case_when(
      SEXO == 1 ~ "hombres",
      SEXO == 2 ~ "mujeres",
      TRUE ~ NA_character_
    ),
    year  = ANO,
    month = MES,
    age_group = case_when(
      ECANT >= 60 & ECANT <= 64 ~ 1,
      ECANT >= 65 & ECANT <= 69 ~ 2,
      ECANT >= 70 & ECANT <= 74 ~ 3,
      ECANT >= 75 & ECANT <= 79 ~ 4,
      ECANT >= 80 & ECANT <= 84 ~ 5,
      ECANT >= 85               ~ 6,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(sex), !is.na(age_group))

print(head(filtered_records %>% select(year, month, sex, age_group)))
print(paste("Registros de muertes filtrados:", nrow(registros_filtrados)))

# Add monthly counts
exp_obs_bymonthyear <- filtered_records %>%
  group_by(year, month, sex, age_group) %>%
  summarise(observed = n(), .groups = "drop") %>%
  complete(year, month, sex, age_group, fill = list(observed = 0)) %>%
  filter(year < 2020) %>%
  mutate(
    sex       = as.factor(sex),
    age_group = factor(age_group, levels = 1:6) 
  ) %>%
  left_join(pop_data, by = c("year", "sex", "age_group")) %>%
  mutate(observed = floor(observed)) %>%
  arrange(year, month, sex, age_group)

print(head(exp_obs_bymonthyear))
print(dim(exp_obs_bymonthyear))
print(summary(exp_obs_bymonthyear$observed))
print(any(is.na(exp_obs_bymonthyear$population)))
print(unique(pop_data$year))

# Create data frame to store expected monthly mortality in 2020-2023
pred_years       <- 2020:2023
pred_months      <- 1:12
sex_levels       <- levels(exp_obs_bymonthyear$sex)
age_group_levels <- levels(exp_obs_bymonthyear$age_group)

sex_levels       <- c("hombres", "mujeres") 
age_group_levels <- 1:6                       

pred_grid <- expand.grid(
  year      = pred_years,
  month     = pred_months,
  sex       = sex_levels,
  age_group = age_group_levels,
  KEEP.OUT.ATTRS = FALSE
) %>%
  mutate(
    sex       = factor(sex,       levels = sex_levels),
    age_group = factor(age_group, levels = 1:6)
  ) %>%
  left_join(pop_data, by = c("year", "sex", "age_group")) %>%
  mutate(
    expected_deaths        = NA_real_,
    expected_deaths_se     = NA_real_,
    expected_log_deaths    = NA_real_,
    expected_log_deaths_se = NA_real_,
    gamma_E             = NA_real_,
    gamma_delta         = NA_real_,
    gamma_E_nb          = NA_real_,
    gamma_delta_nb      = NA_real_
  )

print(any(is.na(pred_grid$population))) 
print(head(pred_grid))


exp_obs_bymonthyear <- exp_obs_bymonthyear %>%
  mutate(sex = factor(sex, levels = sex_levels))

exp_obs_bymonthyear <- exp_obs_bymonthyear %>%
  mutate(
    expected_deaths        = NA_real_,
    expected_deaths_se     = NA_real_,
    expected_log_deaths    = NA_real_,
    expected_log_deaths_se = NA_real_,
    gamma_E             = NA_real_,
    gamma_delta         = NA_real_
  )

num_samples <- 10000
n_pred      <- nrow(pred_grid)


# Fit model
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

# Get predictions for 2020-2023 
pred <- predict(annual_model,
                se.fit  = TRUE,
                type    = "response",
                newdata = pred_grid)
pred_grid$expected_deaths    <- pred$fit
pred_grid$expected_deaths_se <- pred$se.fit

pred_log <- predict(annual_model,
                    se.fit  = TRUE,
                    newdata = pred_grid)
pred_grid$expected_log_deaths    <- pred_log$fit
pred_grid$expected_log_deaths_se <- pred_log$se.fit

# Gamma parameters for each prediction row
gamma_E        <- numeric(n_pred)
gamma_delta    <- numeric(n_pred)
gamma_E_nb     <- numeric(n_pred)
gamma_delta_nb <- numeric(n_pred)

for (j in seq_len(n_pred)) {
  samples <- exp(rnorm(num_samples,
                       mean = pred_log$fit[j],
                       sd   = pred_log$se.fit[j]))
  gamma_E[j]       <- mean(samples)
  gamma_delta[j]   <- (gamma_E[j]^2) / var(samples)

  samples_nb        <- rnbinom(num_samples, size = overd, mu = samples)
  gamma_E_nb[j]     <- mean(samples_nb)
  gamma_delta_nb[j] <- (gamma_E_nb[j]^2) / var(samples_nb)
}
pred_grid$gamma_E        <- gamma_E
pred_grid$gamma_delta    <- gamma_delta
pred_grid$gamma_E_nb     <- gamma_E_nb
pred_grid$gamma_delta_nb <- gamma_delta_nb

# Fitted values for historical periods 
pred_hist <- predict(annual_model, se.fit = TRUE, type = "response")
exp_obs_bymonthyear$expected_deaths    <- pred_hist$fit
exp_obs_bymonthyear$expected_deaths_se <- pred_hist$se.fit

pred_log_hist <- predict(annual_model, se.fit = TRUE)
exp_obs_bymonthyear$expected_log_deaths    <- pred_log_hist$fit
exp_obs_bymonthyear$expected_log_deaths_se <- pred_log_hist$se.fit

num_hist         <- nrow(exp_obs_bymonthyear)
gamma_E_hist     <- numeric(num_hist)
gamma_delta_hist <- numeric(num_hist)

for (j in seq_len(num_hist)) {
  samples             <- exp(rnorm(num_samples,
                                  mean = pred_log_hist$fit[j],
                                  sd   = pred_log_hist$se.fit[j]))
  gamma_E_hist[j]     <- mean(samples)
  gamma_delta_hist[j] <- (gamma_E_hist[j]^2) / var(samples)
}
exp_obs_bymonthyear$gamma_E     <- gamma_E_hist
exp_obs_bymonthyear$gamma_delta <- gamma_delta_hist

# Save all expected (historical and 2020-2023) 
deaths_monthly_predictions_hist <- exp_obs_bymonthyear %>%
  select(year, month, sex, age_group,
         expected_deaths, expected_deaths_se,
         expected_log_deaths, expected_log_deaths_se,
         gamma_E, gamma_delta) %>%
  bind_rows(
    pred_grid %>%
      select(year, month, sex, age_group,
             expected_deaths, expected_deaths_se,
             expected_log_deaths, expected_log_deaths_se,
             gamma_E, gamma_delta)
  ) %>%
  arrange(year, month, sex, age_group)

save(deaths_monthly_predictions_hist,
     file = "./deaths_monthly_predictions_tier1_hist.RData")

# Save expected only for 2020-2023 
deaths_monthly_predictions_tier1 <- pred_grid %>%
  select(year, month, sex, age_group,
         expected_deaths, expected_deaths_se,
         expected_log_deaths, expected_log_deaths_se,
         gamma_E, gamma_delta, gamma_E_nb, gamma_delta_nb) %>%
  mutate(
    gamma_sd    = sqrt((gamma_E^2)    / gamma_delta),
    gamma_sd_nb = sqrt((gamma_E_nb^2) / gamma_delta_nb)
  )

save(deaths_monthly_predictions_tier1,
     file = "./deaths_monthly_predictions_tier1.RData")
