data {
  int<lower=1> N;                         // Number of individuals
  int<lower=1> J;                         // Number of studies
  array[N] int<lower=0> m;                // Mosquitoes fed
  array[N] int<lower=0> y;                // Infected mosquitoes
  vector[N] age;                          // Age already centered in R
  vector[N] logP_obs;                     // Observed parasite density (log10(par+0.1))
  vector[N] logG_obs;                     // Observed gametocyte density (log10(gam+0.1))
  array[N] int<lower=1, upper=J> study;   // Study ID
  array[N] int<lower=0, upper=1> fever;   // Fever 0/1

  // Prediction grid
  int<lower=1> N_pred;
  vector[N_pred] logG_pred_grid;
  vector[N_pred] logP_pred_grid;
}

parameters {
  // Parasite latent model
  real mu_P;
  real<lower=0> sigma_P;
  real gamma_P_age;
  vector[N] logP_raw;

  // Gametocyte latent model (depends on centered logP)
  real mu_G;
  real<lower=0> sigma_G;
  real gamma_G_P;
  real gamma_G_age;
  vector[N] logG_raw;

  // Infection model
  real alpha;
  vector[J - 1] beta_study;
  real beta_G;
  real beta_P;
  real beta_age;
  real beta_fever;

  // Overdispersion for Beta-Binomial
  real<lower=0> phi;
}

transformed parameters {
  vector[N] logP;         // latent true parasite density
  vector[N] logP_c;       // centered parasite
  real mean_logP;         // mean over the sample

  vector[N] logG;         // latent true gametocyte density

  // 1. Parasite latent process
  for (i in 1:N) {
    logP[i] = mu_P +
              gamma_P_age * age[i] +
              sigma_P * logP_raw[i];
  }

  // 2. Center parasite density
  mean_logP = mean(logP);

  for (i in 1:N) {
    logP_c[i] = logP[i] - mean_logP;  // centered version for gametocyte model
  }

  // 3. Gametocyte latent process (uses centered logP_c)
  for (i in 1:N) {
    logG[i] = mu_G +
              gamma_G_P * logP_c[i] +
              gamma_G_age * age[i] +
              sigma_G * logG_raw[i];
  }
}

model {
  // ---- PRIORS ----
  // Parasite priors
  mu_P ~ normal(3.5, 0.5);
  sigma_P ~ exponential(2);
  gamma_P_age ~ student_t(3, 0, 2);//normal(0, 1);
  logP_raw ~ std_normal();

  // Gametocyte priors
  mu_G ~student_t(3, 4, 2);//  mu_G ~ normal(2.5, 0.5);
  sigma_G ~ exponential(2);
  gamma_G_P ~ normal(0.3, 0.5);
  gamma_G_age ~ student_t(3, 0, 2);//normal(0.05, 0.2);
  logG_raw ~ std_normal();

  // Infection model priors
  alpha ~ normal(-3, 1);
  beta_study ~ normal(0, 0.5);
  beta_G ~ normal(1, 0.05);
  beta_P ~ normal(0.75, 0.05);
  beta_age ~student_t(3, 0, 2);// normal(0, 1);
  beta_fever ~ normal(-0.3, 0.5);
  phi ~ exponential(2);

  // ---- MEASUREMENT ERROR MODELS ----
  logP_obs ~ normal(logP, 0.5);
  logG_obs ~ normal(logG, 0.5);

  // ---- INFECTION MODEL ----
  for (i in 1:N) {
    real study_effect = (study[i] == 1) ? 0 : beta_study[study[i] - 1];
    real logit_p = alpha +
                   study_effect +
                   beta_G * logG[i] +
                   beta_P * logP[i] +
                   beta_age * age[i] +
                   beta_fever * fever[i];
    real p = inv_logit(logit_p);

    target += lgamma(m[i] + 1)
              - lgamma(y[i] + 1)
              - lgamma(m[i] - y[i] + 1)
              + lbeta(y[i] + p * phi,
                      m[i] - y[i] + (1 - p) * phi)
              - lbeta(p * phi, (1 - p) * phi);
  }
}

generated quantities {
  array[N] int y_rep;
  vector[N] log_lik;
  matrix[N_pred, J] pred_prob;

  // Prediction grid assumes fever = 0, age = mean, centered parasite
  for (i in 1:N_pred) {
    for (j in 1:J) {
      real study_effect = (j == 1) ? 0 : beta_study[j - 1];
      real logit_pred = alpha +
                        study_effect +
                        beta_G * logG_pred_grid[i] +
                        beta_P * logP_pred_grid[i];
      pred_prob[i, j] = inv_logit(logit_pred);
    }
  }

  // Posterior predictive checks
  for (i in 1:N) {
    real study_effect = (study[i] == 1) ? 0 : beta_study[study[i] - 1];
    real logit_p = alpha +
                   study_effect +
                   beta_G * logG[i] +
                   beta_P * logP[i] +
                   beta_age * age[i] +
                   beta_fever * fever[i];

    real p = inv_logit(logit_p);

    real a = p * phi;
    real b = (1 - p) * phi;

    real theta = beta_rng(a, b);
    y_rep[i] = binomial_rng(m[i], theta);

    log_lik[i] = lgamma(m[i] + 1)
                 - lgamma(y[i] + 1)
                 - lgamma(m[i] - y[i] + 1)
                 + lbeta(y[i] + a,
                         m[i] - y[i] + b)
                 - lbeta(a, b);
  }

  // mediation  effect
  real IE_age_gam = gamma_P_age * gamma_G_P;
  real TE_age_gam = gamma_G_age + IE_age_gam;
  real IE_age_inf = gamma_P_age * beta_P + beta_G * (gamma_G_age + gamma_P_age * gamma_G_P);
  real TE_age_inf = beta_age + IE_age_inf;
  real IE_par_inf = gamma_G_P * beta_G;
  real TE_par_inf = beta_P + IE_par_inf;
}

