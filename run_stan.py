import numpy as np  # type: ignore
import stan  # type: ignore


def naive_ra_at_one_percent(
    virus_counts: list[int], all_counts: list[int], prev_per_100k: float
) -> float:
    total_virus = sum(virus_counts)
    total_counts = sum(all_counts)
    return 1e3 * total_virus / total_counts / prev_per_100k


stan_code = """
data {
  int<lower=1> J;               // number of samples
  array[J] int<lower=0> y;      // reads mapped to virus
  vector[J] n;                  // total reads
  real mu;                      // mean log prevalence
  real<lower=0>  sigma;         // std log prevalence
}
parameters {
  real b;                 // log conversion factor
  real<lower=0> phi;      // inverse overdispersion
  real theta;             // true log prevalence
}
model {
  b ~ normal(0, 10);
  phi ~ gamma(2, 2);

  theta ~ normal(mu, sigma);
  y ~ neg_binomial_2(exp(b + theta) * n, phi);
}
"""

# Norovirus data from Rothman
all_reads_millions = [40, 166, 6.7, 79, 42, 103, 97, 41, 86]
all_reads = [int(x * 1e6) for x in all_reads_millions]
n_studies = len(all_reads)
relative_abundance = [
    2.5e-8,
    2.9e-6,
    0.0,
    2.5e-8,
    1.9e-7,
    1.1e-6,
    7.1e-7,
    4.4e-6,
    6.7e-7,
]
virus_reads = [
    round(tot * ra) for tot, ra in zip(all_reads, relative_abundance)
]
assert len(virus_reads) == n_studies
norovirus_prevalence_per100k = 0.8
mu = np.log(norovirus_prevalence_per100k)
sigma = 1.0

naive = naive_ra_at_one_percent(
    virus_reads, all_reads, norovirus_prevalence_per100k
)


data = {
    "J": n_studies,
    "y": virus_reads,
    "n": all_reads,
    "mu": mu,
    "sigma": sigma,
}
seed = 1
posterior = stan.build(stan_code, data=data, random_seed=seed)
fit = posterior.sample(num_chains=4, num_samples=1000)
ra_at_one_percent = np.exp(fit["b"]) * 1e3
percentiles = [5, 25, 50, 75, 95]
sep = "\t\t"
title = "Relative abundance of norovirus at 1% prevalence"
print("-" * len(title))
print(title)
print("-" * len(title))
print("Naive estimate:")
print(f"{naive:0.5f}")
print("Posterior quantiles:")
print(*(f"{p}%" for p in percentiles), sep=sep)
print(
    *(f"{x:0.5f}" for x in np.percentile(ra_at_one_percent, percentiles)),
    sep=sep,
)
