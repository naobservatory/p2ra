from pathlib import Path
from typing import IO

import numpy as np
import numpy.typing as npt
import stan  # type: ignore


def naive_relative_abundance(
    virus_counts: npt.ArrayLike,
    all_counts: npt.ArrayLike,
    prev_per_100k: np.float_,
) -> float:
    total_virus = np.sum(virus_counts)
    total_counts = np.sum(all_counts)
    return total_virus / total_counts / prev_per_100k


stan_code = """
data {
  int<lower=1> J;               // number of samples
  array[J] int<lower=0> y;      // reads mapped to virus
  vector[J] n;                  // total reads
  vector[J] mu;                 // mean log prevalence
}
transformed data {
  real<lower=0> sigma = 0.5;         // std log prevalence
}
parameters {
  real b;                 // log conversion factor
  real<lower=0> phi;      // inverse overdispersion
  vector[J] theta;        // true log prevalence
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
"""


# TODO: Log stan output rather than writing to stderr
def fit_model(
    num_samples: int,
    viral_read_counts: npt.ArrayLike,
    total_read_counts: npt.ArrayLike,
    mean_log_prevalence: npt.ArrayLike,
    random_seed: int,
) -> stan.fit.Fit:
    if isinstance(mean_log_prevalence, float):
        mu = np.ones(num_samples) * mean_log_prevalence
    else:
        mu = np.array(mean_log_prevalence)
    data = {
        "J": num_samples,
        "y": viral_read_counts,
        "n": total_read_counts,
        "mu": mu,
    }
    model = stan.build(stan_code, data=data, random_seed=random_seed)
    fit = model.sample(num_chains=4, num_samples=1000)
    return fit
