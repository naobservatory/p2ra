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
  int<lower=1> J;
  array[J] int<lower=0> viral_reads;
  vector[J] total_reads;
  vector[J] prevalence_per100k;
}
transformed data {
  vector[J] mu = log(prevalence_per100k);
  real<lower=0> sigma = 0.5;
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
  viral_reads ~ neg_binomial_2_log(b + theta + log(total_reads), phi);
}
generated quantities {
  array[J] int<lower=0> y_tilde
    = neg_binomial_2_log_rng(b + theta + log(total_reads), phi);
}
"""


# TODO: Log stan output rather than writing to stderr
def fit_model(
    data: dict,
    random_seed: int,
) -> stan.fit.Fit:
    stan_data = {
        k: v
        for k, v in data.items()
        if k in ["viral_reads", "total_reads", "prevalence_per100k"]
    }
    stan_data["J"] = len(stan_data["viral_reads"])
    model = stan.build(stan_code, data=stan_data, random_seed=random_seed)
    fit = model.sample(num_chains=4, num_samples=1000)
    return fit
