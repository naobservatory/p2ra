data {
  int<lower=1> J;
  array[J] int<lower=0> y;
  vector[J] n;
  vector[J] x;
}
transformed data {
  vector[J] mu = log(x);
  real<lower=0> sigma = 0.5;
}
parameters {
  real b;                 // log conversion factor
  real<lower=0> phi;      // inverse overdispersion
  vector[J] theta;        // true log incidence
}
model {
  b ~ normal(0, 10);
  phi ~ gamma(2, 2);

  theta ~ normal(mu, sigma);
  y ~ neg_binomial_2_log(b + theta + log(n), phi);
}
generated quantities {
  array[J] int<lower=0> y_tilde
    = neg_binomial_2_log_rng(b + theta + log(n), phi);
}
