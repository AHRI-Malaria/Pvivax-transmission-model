functions {
  real beta_binomial_mean_precision_lpmf(int y, int m, real p, real phi) {
    real pp = fmin(1 - 1e-10, fmax(1e-10, p));
    return beta_binomial_lpmf(y | m, pp * phi, (1 - pp) * phi);
  }
}

data {
  int<lower=1> N;
  int<lower=1> J;
  array[N] int<lower=1> m;
  array[N] int<lower=0> y;
  vector[N] age10;
  array[N] int<lower=0, upper=1> fever;
  array[N] int<lower=1, upper=J> study;

  array[N] int<lower=0, upper=1> P_positive;
  array[N] int<lower=0, upper=1> G_positive;
  vector[N] logP_obs;
  vector[N] logG_obs;
  vector[J] log_lod_P;
  vector[J] log_lod_G;

  real<lower=1e-6> sigma_me_P;
  real<lower=1e-6> sigma_me_G;
  real center_P;
  real center_G;

  int<lower=0, upper=1> include_beta_P;
  int<lower=0, upper=1> prior_only;
  array[N] int<lower=0, upper=1> use_y;

  // All Normal priors use mean and standard deviation.
  real prior_mu_P_mean;
  real<lower=0> prior_mu_P_sd;
  real prior_mu_G_mean;
  real<lower=0> prior_mu_G_sd;
  real<lower=0> prior_latent_sd_scale;
  real prior_gamma_GP_mean;
  real<lower=0> prior_gamma_GP_sd;
  real<lower=0> prior_age_density_sd;
  real<lower=0> prior_density_study_sd;

  real prior_alpha_mean;
  real<lower=0> prior_alpha_sd;
  real<lower=0> prior_beta_density_sd;
  real<lower=0> prior_beta_age_sd;
  real<lower=0> prior_beta_fever_sd;
  real<lower=0> prior_beta_study_sd;
  real prior_log_phi_mean;
  real<lower=0> prior_log_phi_sd;
}

parameters {
  // Latent asexual-stage density process.
  real mu_P;
  real<lower=0> sigma_P;
  real gamma_P_age;
  vector[J - 1] delta_P_study_raw;
  vector[N] z_P;

  // Latent gametocyte density process conditional on asexual-stage density.
  real mu_G;
  real<lower=0> sigma_G;
  real gamma_G_P;
  real gamma_G_age;
  vector[J - 1] delta_G_study_raw;
  vector[N] z_G;

  // Mosquito infection process.
  real alpha;
  real beta_G;
  real beta_P_raw;
  real beta_age;
  real beta_fever;
  vector[J - 1] beta_study_raw;
  real log_phi;
}

transformed parameters {
  vector[J] delta_P_study;
  vector[J] delta_G_study;
  vector[J] beta_study;
  vector[N] logP;
  vector[N] logG;
  real beta_P;
  real<lower=0> phi;

  delta_P_study[1] = 0;
  delta_G_study[1] = 0;
  beta_study[1] = 0;
  if (J > 1) {
    delta_P_study[2:J] = delta_P_study_raw;
    delta_G_study[2:J] = delta_G_study_raw;
    beta_study[2:J] = beta_study_raw;
  }

  beta_P = include_beta_P * beta_P_raw;
  phi = exp(log_phi);

  for (i in 1:N) {
    logP[i] = mu_P
              + delta_P_study[study[i]]
              + gamma_P_age * age10[i]
              + sigma_P * z_P[i];

    logG[i] = mu_G
              + delta_G_study[study[i]]
              + gamma_G_P * (logP[i] - center_P)
              + gamma_G_age * age10[i]
              + sigma_G * z_G[i];
  }
}

model {
  // Latent-process priors.
  mu_P ~ normal(prior_mu_P_mean, prior_mu_P_sd);
  mu_G ~ normal(prior_mu_G_mean, prior_mu_G_sd);
  sigma_P ~ normal(0, prior_latent_sd_scale);
  sigma_G ~ normal(0, prior_latent_sd_scale);
  gamma_G_P ~ normal(prior_gamma_GP_mean, prior_gamma_GP_sd);
  gamma_P_age ~ normal(0, prior_age_density_sd);
  gamma_G_age ~ normal(0, prior_age_density_sd);
  delta_P_study_raw ~ normal(0, prior_density_study_sd);
  delta_G_study_raw ~ normal(0, prior_density_study_sd);
  z_P ~ std_normal();
  z_G ~ std_normal();

  // Weakly informative, sceptical infectivity priors.
  alpha ~ normal(prior_alpha_mean, prior_alpha_sd);
  beta_G ~ normal(0, prior_beta_density_sd);
  beta_P_raw ~ normal(0, prior_beta_density_sd);
  beta_age ~ normal(0, prior_beta_age_sd);
  beta_fever ~ normal(0, prior_beta_fever_sd);
  beta_study_raw ~ normal(0, prior_beta_study_sd);
  log_phi ~ normal(prior_log_phi_mean, prior_log_phi_sd);

  if (prior_only == 0) {
    for (i in 1:N) {
      // Positive microscopy values contribute a Gaussian measurement density.
      // A recorded zero contributes a left-censoring probability rather than
      // being discarded or replaced by an arbitrary constant.
      if (P_positive[i] == 1) {
        target += normal_lpdf(logP_obs[i] | logP[i], sigma_me_P);
      } else {
        target += normal_lcdf(log_lod_P[study[i]] | logP[i], sigma_me_P);
      }

      if (G_positive[i] == 1) {
        target += normal_lpdf(logG_obs[i] | logG[i], sigma_me_G);
      } else {
        target += normal_lcdf(log_lod_G[study[i]] | logG[i], sigma_me_G);
      }

      if (use_y[i] == 1) {
        real eta = alpha
                   + beta_study[study[i]]
                   + beta_G * (logG[i] - center_G)
                   + beta_P * (logP[i] - center_P)
                   + beta_age * age10[i]
                   + beta_fever * fever[i];
        target += beta_binomial_mean_precision_lpmf(
          y[i] | m[i], inv_logit(eta), phi
        );
      }
    }
  }
}

generated quantities {
  array[N] int y_rep;
  vector[N] p_infect;
  vector[N] p_any;
  vector[N] log_lik;

  // Model-implied association components on the log-odds scale.
  real assoc_P_via_G = gamma_G_P * beta_G;
  real assoc_P_residual = beta_P;
  real assoc_P_total = beta_P + gamma_G_P * beta_G;

  real assoc_age_via_P_residual = gamma_P_age * beta_P;
  real assoc_age_via_P_G = gamma_P_age * gamma_G_P * beta_G;
  real assoc_age_via_G = gamma_G_age * beta_G;
  real assoc_age_residual = beta_age;
  real assoc_age_total = beta_age
                         + gamma_P_age * beta_P
                         + gamma_P_age * gamma_G_P * beta_G
                         + gamma_G_age * beta_G;

  for (i in 1:N) {
    real eta = alpha
               + beta_study[study[i]]
               + beta_G * (logG[i] - center_G)
               + beta_P * (logP[i] - center_P)
               + beta_age * age10[i]
               + beta_fever * fever[i];
    real p = fmin(1 - 1e-10, fmax(1e-10, inv_logit(eta)));
    real a = p * phi;
    real b = (1 - p) * phi;
    real theta = beta_rng(a, b);
    real p_zero = exp(beta_binomial_lpmf(0 | m[i], a, b));

    p_infect[i] = p;
    p_any[i] = 1 - p_zero;
    y_rep[i] = binomial_rng(m[i], theta);
    log_lik[i] = beta_binomial_mean_precision_lpmf(y[i] | m[i], p, phi);
  }
}
