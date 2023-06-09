data {
  int<lower=1> J;
  array[J] int<lower=0> y;  // viral read counts
  array[J] int<lower=0> n;  // total read counts
  vector[J] x;              // estimated predictor (prevalence or incidence)
}
transformed data {
  vector[J] mu = log(x) - mean(log(x));
  real log_mean_y = log(mean(y));
  real log_mean_n = log(mean(n));
}
parameters {
  real<lower=0> sigma;      // standard deviation of true predictors
  real b_std;               // P2RA coefficient (on standardized scale)
  vector[J] theta_std;      // standardized true predictor for each sample
}
model {
  sigma ~ exponential(1);
  b_std ~ normal(0, 2);
  theta_std ~ normal(mu, sigma);
  y ~ binomial_logit(n, b_std + theta_std + log_mean_y - log_mean_n);
}
generated quantities {
  // posterior predictive viral read counts
  array[J] int<lower=0> y_tilde
    = binomial_rng(n, inv_logit(b_std + theta_std + log_mean_y - log_mean_n));
  // posterior true prevalence for each sample
  vector[J] theta = theta_std + mean(log(x));
  // posterior P2RA coefficient
  real b = b_std - mean(log(x)) + log_mean_y - log_mean_n;
}
