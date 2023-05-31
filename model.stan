data {
  int<lower=1> J;
  array[J] int<lower=0> y;  // viral read counts
  vector[J] n;              // total read counts
  vector[J] x;              // estimated predictor (prevalence or incidence)
}
transformed data {
  vector[J] mu = log(x) - mean(log(x));
  real<lower=0> sigma = 0.5;
  vector[J] exposure = log(n) - mean(log(n));
  real log_mean_y = log(mean(y));
}
parameters {
  real b_std;               // P2RA coefficient (on standardized scale)
  real<lower=0> phi;        // inverse overdispersion
  vector[J] theta_std;      // standardized true predictor for each sample
}
model {
  b_std ~ normal(0, 2);
  phi ~ gamma(2, 2);
  theta_std ~ normal(mu, sigma);
  y ~ neg_binomial_2_log(b_std + theta_std + exposure + log_mean_y, phi);
}
generated quantities {
  // posterior predictive viral read counts
  array[J] int<lower=0> y_tilde
    = neg_binomial_2_log_rng(b_std + theta_std + exposure + log_mean_y, phi);
  // posterior true prevalence for each sample
  vector[J] theta = theta_std + mean(log(x));
  // posterior P2RA coefficient
  real b = b_std - mean(log(x)) - mean(log(n)) + log_mean_y;
}
