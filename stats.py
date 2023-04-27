import stan  # type: ignore


def naive_ra(
    virus_counts: list[int], all_counts: list[int], prev_per_100k: float
) -> float:
    total_virus = sum(virus_counts)
    total_counts = sum(all_counts)
    return total_virus / total_counts / prev_per_100k


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


# TODO: Allow mean and log prevalences to be lists
# TODO: Figure out the actual required datatypes for list-like objects
# TODO: Log stan output rather than writing to stderr
def fit_model(
    num_studies: int,
    viral_read_counts: list[int],
    total_read_counts: list[int],
    mean_log_prevalence: float,
    std_log_prevalence: float,
    random_seed: int,
) -> stan.fit.Fit:
    data = {
        "J": num_studies,
        "y": viral_read_counts,
        "n": total_read_counts,
        "mu": mean_log_prevalence,
        "sigma": std_log_prevalence,
    }
    model = stan.build(stan_code, data=data, random_seed=random_seed)
    fit = model.sample(num_chains=4, num_samples=1000)
    return fit
