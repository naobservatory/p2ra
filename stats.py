from typing import TypeVar

import numpy as np  # type: ignore
import stan  # type: ignore

T = TypeVar("T")
# TODO: Specify the dtype of ndarray
ArrayLike = list[T] | np.ndarray


def naive_relative_abundance(
    virus_counts: ArrayLike[int],
    all_counts: ArrayLike[int],
    prev_per_100k: float,
) -> float:
    total_virus = sum(virus_counts)
    total_counts = sum(all_counts)
    return total_virus / total_counts / prev_per_100k


stan_code = """
data {
  int<lower=1> J;               // number of samples
  array[J] int<lower=0> y;      // reads mapped to virus
  vector[J] n;                  // total reads
  vector[J] mu;                 // mean log prevalence
  real<lower=0>  sigma;         // std log prevalence
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
  y ~ neg_binomial_2(exp(b + theta) .* n, phi);
}
"""


# TODO: Log stan output rather than writing to stderr
def fit_model(
    num_samples: int,
    viral_read_counts: ArrayLike[int],
    total_read_counts: ArrayLike[int],
    mean_log_prevalence: float | ArrayLike[float],
    std_log_prevalence: float,
    random_seed: int,
) -> stan.fit.Fit:
    mu: ArrayLike[float]
    if isinstance(mean_log_prevalence, float):
        mu = num_samples * [mean_log_prevalence]
    else:
        mu = mean_log_prevalence
    data = {
        "J": num_samples,
        "y": viral_read_counts,
        "n": total_read_counts,
        "mu": mu,
        "sigma": std_log_prevalence,
    }
    model = stan.build(stan_code, data=data, random_seed=random_seed)
    fit = model.sample(num_chains=4, num_samples=1000)
    return fit
