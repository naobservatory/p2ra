data {
  int<lower=1> J;
  array[J] int<lower=0> y;
  vector[J] n;
  vector[J] x;
}
transformed data {
  vector[J] mu = log(x) - mean(log(x));
  real<lower=0> sigma = 0.5;
  vector[J] exposure = log(n) - mean(log(n));
  real log_mean_y = log(mean(y));
}
parameters {
  real b_std;
  real<lower=0> phi;      // inverse overdispersion
  vector[J] theta_std;
}
model {
  b_std ~ normal(0, 2);
  phi ~ gamma(2, 2);

  theta_std ~ normal(mu, sigma);
  y ~ neg_binomial_2_log(b_std + theta_std + exposure + log_mean_y, phi);
}
generated quantities {
  array[J] int<lower=0> y_tilde
    = neg_binomial_2_log_rng(b_std + theta_std + exposure + log_mean_y, phi);
  vector[J] theta = theta_std + mean(log(x));
  real b = b_std - mean(log(x)) - mean(log(n)) + log_mean_y;
}
