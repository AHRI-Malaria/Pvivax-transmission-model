
### Author: Legesse Alamerie Ejigu

########## final code beta binomial Bayesian joint measurement error model####

library(rstan)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(bayesplot)
library(patchwork)
library(knitr)
library(data.table)
library(tidybayes)
library(cmdstanr)
library(posterior)
library(loo)
library(viridis)
library(ggpubr)
library(scales)

#reading data 
df=read.csv("data/pvivax_transmission_data.csv")

# Data Preparation

df <- df %>%
  mutate(
    logG_obs = log10(gam),
    logP_obs = log10(par)
  )

df$study <- factor(df$study)
df$study_id <- as.numeric(df$study)

df <- df %>% mutate(study_indiv = paste(study_id, id, sep = "_")) %>% 
  mutate(indiv_id = as.numeric(as.factor(study_indiv)))

N <- nrow(df)
S <- length(unique(df$study_id))
I <- length(unique(df$indiv_id))

study_of_indiv <- df %>% group_by(indiv_id) %>% summarise(study_id_for_indiv = first(study_id)) %>% 
  arrange(indiv_id) %>% pull(study_id_for_indiv)


summary(log10(df$gam))

df$fever <- ifelse(df$fever == "Yes", 1, 0)


table(df$fever)
# Generate list of data for stan

stan_data <- list(
  N = nrow(df),
  J = length(unique(df$study)),
  m = df$number_dissect,
  y = df$number_inf,
  age = scale(df$age, center = TRUE, scale = FALSE)[,1],
  fever = as.integer(df$fever),  
  logP_obs = log10(df$par+0.1),
  logG_obs = log10(df$gam+0.1),
  study = as.integer(factor(df$study)),
  age_mean = mean(df$age),
  age_sd = sd(df$age),
  N_pred = 100,
  logG_pred_grid = seq(0, 6.6, length.out = 100),
  logP_pred_grid = seq(0, 6.9, length.out = 100)  #
)


### distribution of gam , par and age 
hist(stan_data$logP_obs)
hist(stan_data$logG_obs)
hist(stan_data$age)
hist(df$age)
hist(df$number_inf)
hist(df$number_dissect)


##reading the stan file 

mod <- cmdstan_model("stan/Joint_model_stan.stan")


### Fit the model 
fit <- mod$sample(
  data = stan_data,
  chains = 4,
  seed = 123,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 8000,
)


# Posterior predictive check

y_rep <- fit$draws("y_rep") |> as_draws_matrix()
pp_c=ppc_dens_overlay(stan_data$y, y_rep[1:100, ])  



ggsave(
  filename = "plots/S7_Fig.tiff",
  plot = pp_c,
  device = "tiff",
  width = 6,      
  height = 5,
  dpi = 600,
  units = "in",
  compression = "lzw"
)

## model accuarcy 
# LOO
log_lik <- fit$draws("log_lik", format = "matrix")
loo_result <- loo(log_lik)
loo_result

# trace plot to check convergence
trace_p=mcmc_trace(fit$draws(c("mu_P", "mu_G", "sigma_P", "sigma_G", "gamma_P_age",   
                               "gamma_G_age",   
                               "gamma_G_P",  "beta_fever",                  
                               "beta_study", "alpha",         
                               "beta_G","beta_P", "beta_age","phi","total_gametocyte_effect")))


ggsave(
  filename = "plots/S8_Fig.tiff",
  plot = trace_p,
  device = "tiff",
  width = 9,      
  height = 6,
  dpi = 600,
  units = "in",
  compression = "lzw"
)



mcmc_areas(fit$draws(c("mu_P", "mu_G", "sigma_P", "sigma_G", "gamma_P_age",   
                       "gamma_G_age",   
                       "gamma_G_P", "beta_fever",                   
                       "alpha", "beta_study",          
                       "beta_G","beta_P", "beta_age","phi","total_gametocyte_effect")), prob = 0.89)

# Density plots for all parameters
dens_p=mcmc_dens(
  fit$draws(variables = c("mu_P", "mu_G", "sigma_P", "sigma_G", 
                          "gamma_P_age", "gamma_G_age", "gamma_G_P", 
                          "alpha", "beta_fever", "beta_G","beta_P", "beta_age","beta_study","phi","total_gametocyte_effect")),
  facet_args = list(ncol = 4)
) +
  ggtitle("")


ggsave(
  filename = "plots/S9_Fig.tiff",
  plot = dens_p,
  device = "tiff",
  width = 9,      
  height = 6,
  dpi = 600,
  units = "in",
  compression = "lzw"
)
### autocorrelation plot 

outcome_pars <- c(
  "alpha",
  "beta_G",
  "beta_P",
  "beta_age",
  "beta_fever"
  
)
library(bayesplot)

acf_p <- mcmc_acf(
  fit$draws(outcome_pars),
  lags = 20
)

acf_p


ggsave(
  filename = "plots/S10_Fig.tiff",
  plot = acf_p,
  device = "tiff",
  width = 9,      
  height = 6,
  dpi = 600,
  units = "in",
  compression = "lzw"
)

####

outcome_pars <- c(
  "beta_G",
  "beta_P",
  "beta_age",
  "beta_fever"
)

library(bayesplot)

pairs_p <- mcmc_pairs(
  fit$draws(outcome_pars),
  diag_fun = "dens",
  off_diag_fun = "hex"
)

pairs_p

ggsave(
  filename = "plots/S11_Fig.tiff",
  plot = pairs_p,
  device = "tiff",
  width = 11,      
  height = 7,
  dpi = 600,
  units = "in",
  compression = "lzw"
)

###  model summary
fit$summary(variables = c("mu_P", "mu_G", "sigma_P", "sigma_G",  "gamma_P_age",   
                          "gamma_G_age",   
                          "gamma_G_P", "beta_fever",                    
                          "beta_study",          
                          "beta_G","beta_P","beta_age","phi","IE_age_gam","TE_age_gam","IE_age_inf","TE_age_inf",
                          "IE_par_inf","TE_par_inf"))


# To get summary for all parameters 
param_summary <- fit$summary(
  variables = c("mu_P", "mu_G", "sigma_P", "sigma_G",
                "gamma_P_age", "gamma_G_age", "gamma_G_P",
                "alpha",  "beta_G","beta_P","beta_fever","beta_age","beta_study"),
  mean, sd, ~quantile(.x, probs = c(0.025, 0.5, 0.975)), rhat, ess_bulk
)


kable(param_summary, caption = "Parameter Estimates with 95% Credible Intervals")




### prediction plot 
# Extract posterior draws 
posterior_draws <- fit$draws(variables = "pred_prob", format = "draws_matrix")

# pred_prob 
# Reshape into a 3D array
pred_prob <- array(posterior_draws, 
                   dim = c(nrow(posterior_draws), stan_data$N_pred, stan_data$J))
library(tidyverse)

pred_df <- map_df(1:stan_data$J, function(j) {
  data.frame(
    logG = stan_data$logG_pred_grid,
    study = j,
    median = apply(pred_prob[,,j], 2, median),
    lower = apply(pred_prob[,,j], 2, function(x) quantile(x, 0.025)),
    upper = apply(pred_prob[,,j], 2, function(x) quantile(x, 0.975))
  )
})
observed_data <- data.frame(
  logG_obs = stan_data$logG_obs,
  study = factor(stan_data$study),
  infected = stan_data$y / stan_data$m,  # Observed proportion infected
  size = stan_data$m
)



pred_df$study=as.factor(pred_df$study)


study_colors <- c("study 1" = "#1f77b4", "study 2" = "#ff7f0e", "study 3" = "#2ca02c")

infection_plot <- ggplot() +
  geom_ribbon(
    data = pred_df,
    aes(x = logG, ymin = lower, ymax = upper, fill = study),
    alpha = 0.2
  ) +
  geom_line(
    data = pred_df,
    aes(x = logG, y = median, color = study),
    linewidth = 1
  ) +
  geom_point(
    data = observed_data,
    aes(x = logG_obs, y = infected, color = study),
    alpha = 0.6
  ) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  labs(
    x = expression(paste("Log"[10], " Pvs25 transcript/μL")),
    y = "Proportion of infected mosquitoes, %"
  ) +
  theme_pubr() +  scale_x_continuous(limits = c(0,6.7),breaks = c(0,1,2,3,4,5,6))+
  
  theme(legend.position = "top")

infection_plot




# Scale predicted values
pred_df <- pred_df %>%
  mutate(
    median = median * 100,
    lower = lower * 100,
    upper = upper * 100
  )

# Scale observed values
observed_data <- observed_data %>%
  mutate(
    infected = infected * 100
  )

infection_plot <- ggplot() +
  geom_ribbon(
    data = pred_df,
    aes(x = logG, ymin = lower, ymax = upper, fill = study),
    alpha = 0.2
  ) +
  geom_line(
    data = pred_df,
    aes(x = logG, y = median, color = study),
    linewidth = 1
  ) +
  geom_point(
    data = observed_data,
    aes(x = logG_obs, y = infected, color = study),
    alpha = 0.6
  ) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  labs(
    x = expression(paste("Log"[10], " Pvs25 transcript/μL")),
    y = "Proportion of infected mosquitoes (%)"
  ) +
  theme_pubr() +
  scale_x_continuous(limits = c(0, 6.7), breaks = 0:6) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  theme(legend.position = "top")

infection_plot



######

ggsave(
  filename = "plots/modelfit_study22.tiff",
  plot = infection_plot,
  device = "tiff",
  width = 6.5,      
  height = 5.5,
  dpi = 600,
  units = "in",
  compression = "lzw"
)
##########

combined_pred <- pred_df %>%
  group_by(logG) %>%
  summarize(
    median = mean(median),
    lower = mean(lower),
    upper = mean(upper)
  )

all_pred <- ggplot() +
  geom_ribbon(
    data = combined_pred,
    aes(x = logG, ymin = lower, ymax = upper),
    fill = "gray70", alpha = 0.4
  ) +
  geom_line(
    data = combined_pred,
    aes(x = logG, y = median),
    color = "black", linewidth = 1.2
  ) +
  geom_point(
    data = observed_data,
    aes(x = logG_obs, y = infected, color = study),
    size = 2, alpha = 0.7
  ) +
  
  labs(
    x = expression(paste("Log"[10], "Pvs25 transcript/μL")),
    y = "Proportion of infected mosquitoes, %"
  ) +
  theme_pubr() +  scale_x_continuous(limits = c(0,6.7),breaks = c(0,1,2,3,4,5,6))+
  
  
  scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c")) +
  theme_pubr()

all_pred


ggsave(
  filename = "plots/Fig_5.tiff",
  plot = all_pred,
  device = "tiff",
  width = 6.5,      
  height = 5.5,
  dpi = 600,
  units = "in",
  compression = "lzw"
)

### combined effect of parasite and gametocyte density
############


#  Extract posterior draws

posterior_draws <- fit$draws(variables = c("pred_prob","alpha","beta_G","beta_P","beta_age","beta_fever","beta_study"), format = "draws_matrix")



#  Parasite + Gametocyte combined effect grid

logP_seq <- seq(min(stan_data$logP_obs) - 0.1, max(stan_data$logP_obs) + 0.1, length.out = 100)
logG_seq <- seq(min(stan_data$logG_obs) - 0.1, max(stan_data$logG_obs) + 0.1, length.out = 100)
pred_grid <- expand.grid(logP = logP_seq, logG = logG_seq)

# Posterior medians
draws_df <- as_draws_df(fit$draws(variables = c("alpha","beta_G","beta_P","beta_age","beta_fever","beta_study"), format = "draws_df"))
posterior_medians <- summarise_all(draws_df, median)

# Compute predicted infectivity (fixed age, study1, afebrile)
fixed_age <- mean(stan_data$age)
study_effect <- 0  # study1 reference
pred_grid$logit_p <- posterior_medians$alpha +
  posterior_medians$beta_G * pred_grid$logG +
  posterior_medians$beta_P * pred_grid$logP +
  posterior_medians$beta_age * fixed_age +
  posterior_medians$beta_fever * 0 +
  study_effect

pred_grid$infectivity <- plogis(pred_grid$logit_p)

# Observed data points
observed_df <- data.frame(
  logP = stan_data$logP_obs,
  logG = stan_data$logG_obs,
  infected = ifelse(stan_data$y > 0, "Yes", "No")
)


#  Combined effect plot

comp_effect <- ggplot(pred_grid, aes(x = logP, y = logG)) +
  geom_tile(aes(fill = infectivity)) +
  scale_fill_viridis(option = "C", limits = c(0,1), name = "Predicted\ninfectivity", direction = -1) +
  geom_point(data = observed_df, aes(x = logP, y = logG, color = infected), size = 2.2, alpha = 0.75) +
  scale_color_manual(values = c("Yes" = "#A63603", "No" = "#BABABA"), name = "Observed\nInfection") +
  labs(x = expression(paste("Log"[10], " Pv18S copies /μL")),
       y = expression(paste("Log"[10], " Pvs25 transcript /μL"))) +
  scale_x_continuous(breaks = seq(floor(min(logP_seq)), ceiling(max(logP_seq)), 1)) +
  scale_y_continuous(breaks = seq(floor(min(logG_seq)), ceiling(max(logG_seq)), 1)) +
  theme_pubr(legend = "right") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 9),
        panel.grid = element_blank())

ggsave("plots/Fig_6.tiff", comp_effect, width = 6.5, height = 5.5, dpi = 600, units = "in", compression = "lzw")



#### density thresholed

# Extract posterior draws

library(posterior)

draws <- as_draws_df(fit)

alpha  <- draws$alpha
beta_G <- draws$beta_G
beta_P <- draws$beta_P
phi    <- draws$phi


###
logP_ref <- mean(colMeans(as_draws_matrix(fit, variables = "logP")))

#Define beta–binomial Pr(Y ≥ 1)
prob_at_least_one <- function(p, m, phi) {
  1 - exp(
    lbeta(p * phi, m + (1 - p) * phi) -
      lbeta(p * phi, (1 - p) * phi)
  )
}

#
m_ref <- 30
target_prob <- 0.5

logG_threshold_onset <- rep(NA_real_, length(alpha))

#beta–binomial probability function

prob_at_least_one <- function(p, m, phi) {
  1 - exp(
    lbeta(p * phi, m + (1 - p) * phi) -
      lbeta(p * phi, (1 - p) * phi)
  )
}


#Solve for the onset threshold

#Pr(Y≥1)=0.05with m=30

logG_threshold_onset <- rep(NA_real_, length(alpha))

for (i in seq_along(alpha)) {
  
  f <- function(logG) {
    
    lp <- alpha[i] +
      beta_G[i] * logG +
      beta_P[i] * logP_ref
    
    p <- plogis(lp)
    
    if (!is.finite(p)) return(NA_real_)
    
    prob_at_least_one(p, m_ref, phi[i]) - target_prob
  }
  
  lower <- -10
  upper <-  10
  
  f_lower <- f(lower)
  f_upper <- f(upper)
  
  if (is.finite(f_lower) && is.finite(f_upper) && f_lower * f_upper < 0) {
    logG_threshold_onset[i] <- uniroot(
      f,
      lower = lower,
      upper = upper,
      tol = 1e-8
    )$root
  }
}

##

quantile(
  logG_threshold_onset,
  probs = c(0.025, 0.5, 0.975),
  na.rm = TRUE
)

threshold_onset_summary <- quantile(
  logG_threshold_onset,
  probs = c(0.025, 0.5, 0.975),
  na.rm = TRUE
)

threshold_onset_summary

10^quantile(logG_threshold_onset, c(0.025, 0.5, 0.975), na.rm = TRUE) - 0.1

G_threshold_onset <- 10^threshold_onset_summary - 0.1

G_threshold_onset


mean(is.na(logG_threshold_onset))



#####both estimate 

prob_at_least_one <- function(p, m, phi) {
  1 - exp(
    lbeta(p * phi, m + (1 - p) * phi) -
      lbeta(p * phi, (1 - p) * phi)
  )
}


### Function to compute thresholds

compute_threshold <- function(alpha, beta_G, beta_P, phi, logP_ref, m_ref = 30, target_prob = 0.05) {
  logG_threshold <- rep(NA_real_, length(alpha))
  
  for (i in seq_along(alpha)) {
    
    f <- function(logG) {
      lp <- alpha[i] + beta_G[i] * logG + beta_P[i] * logP_ref
      p <- plogis(lp)
      if (!is.finite(p)) return(NA_real_)
      prob_at_least_one(p, m_ref, phi[i]) - target_prob
    }
    
    lower <- -10
    upper <- 10
    f_lower <- f(lower)
    f_upper <- f(upper)
    
    if (is.finite(f_lower) && is.finite(f_upper) && f_lower * f_upper < 0) {
      logG_threshold[i] <- uniroot(f, lower = lower, upper = upper, tol = 1e-8)$root
    }
  }
  return(logG_threshold)
}


###

draws <- as_draws_df(fit)

alpha  <- draws$alpha
beta_G <- draws$beta_G
beta_P <- draws$beta_P
phi    <- draws$phi

# Parasite density reference (mean over posterior)
logP_ref <- mean(draws$logP)


##

# Onset threshold (5% probability)
logG_onset <- compute_threshold(alpha, beta_G, beta_P, phi, logP_ref, m_ref = 30, target_prob = 0.05)

# Likely transmission threshold (50% probability)
logG_50 <- compute_threshold(alpha, beta_G, beta_P, phi, logP_ref, m_ref = 30, target_prob = 0.5)


#

summary_onset <- quantile(logG_onset, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
summary_50 <- quantile(logG_50, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

# Back-transform to gametocytes/µL
G_onset <- 10^summary_onset - 0.1
G_50    <- 10^summary_50 - 0.1

G_onset
G_50

###

#Plot posterior curves with thresholds
library(ggplot2)

# Grid of logG for plotting probability curves
logG_grid <- seq(-2, 4, length.out = 200)

# Compute posterior median probability of infecting ≥1 mosquito
prob_grid <- sapply(logG_grid, function(lg) {
  median(sapply(1:length(alpha), function(i) {
    lp <- alpha[i] + beta_G[i]*lg + beta_P[i]*logP_ref
    p  <- plogis(lp)
    prob_at_least_one(p, m = 30, phi = phi[i])  # call with explicit argument names
  }))
})


# Grid of logG for plotting probability curves
logG_grid <- seq(0, 6.9, length.out = 200)

# Compute posterior median probability of infecting ≥1 mosquito
prob_grid <- sapply(logG_grid, function(lg) {
  median(sapply(1:length(alpha), function(i) {
    lp <- alpha[i] + beta_G[i]*lg + beta_P[i]*logP_ref
    p  <- plogis(lp)
    prob_at_least_one(p, m = 30, phi = phi[i])  # call with explicit argument names
  }))
})


