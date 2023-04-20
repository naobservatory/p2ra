import stan  # type: ignore

stan_code = """
data {
  int<lower=1> J;         // number of samples
  int<lower=0> y[J];      // reads mapped to virus
  vector[J] n;            // total reads
  real mu;                // mean log prevalence
  real<lower=0>  sigma;   // std log prevalence
}
parameters {
  real b;                 // log conversion factor
  real<lower=0> phi;      // inverse overdispersion
  real theta;             // true log prevalence
}
model {
  b ~ normal(-20, 10);
  phi ~ gamma(2, 2);

  theta ~ normal(mu, sigma);
  y ~ neg_binomial_2(exp(b + theta) * n, phi);
}
"""

data = {
    "J": 3,
    "y": [5, 10, 15],
    "n": [50, 110, 140],
    "mu": 0.1,
    "sigma": 0.25,
}
seed = 1
posterior = stan.build(stan_code, data=data, random_seed=seed)
fit = posterior.sample(num_chains=4, num_samples=1000)
print(fit["b"])
