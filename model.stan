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
}
parameters {
  real b;                 // log conversion factor
  real<lower=0> phi;      // inverse overdispersion
  vector[J] theta_std;
}
model {
  b ~ normal(0, 2);
  phi ~ gamma(2, 2);

  theta_std ~ normal(mu, sigma);
  y ~ neg_binomial_2_log(b + theta_std + exposure, phi);
}
generated quantities {
  array[J] int<lower=0> y_tilde
    = neg_binomial_2_log_rng(b + theta_std + exposure, phi);
  vector[J] theta = theta_std + mean(log(x));
}
