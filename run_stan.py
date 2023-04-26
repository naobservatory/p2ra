#!/usr/bin/env python3

from collections import Counter
import numpy as np  # type: ignore
import stan  # type: ignore
import mgs

from pathogens import pathogens


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

if __name__ == "__main__":
    repo = mgs.GitHubRepo(
        user="naobservatory", repo="mgs-pipeline", branch="main"
    )
    bp_data = mgs.load_bioprojects(repo)
    sample_data = mgs.load_sample_attributes(repo)
    counts = mgs.load_sample_counts(repo)
    taxtree = mgs.load_tax_tree(repo)

    bioproject = mgs.BioProject("PRJNA729801")  # Rothman
    samples = bp_data[bioproject]
    sample_attribs = {s: sample_data[s] for s in samples}
    all_reads = Counter()
    for attribs in sample_attribs.values():
        all_reads[attribs.fine_location] += attribs.reads
    studies = all_reads.keys()

    pathogen = "norovirus"
    taxid = pathogens[pathogen].pathogen_chars.taxid
    subtree = taxtree[taxid]
    assert subtree is not None

    samples_by_study = {
        study: [s for s in samples if sample_attribs[s].fine_location == study]
        for study in studies
    }

    virus_counts = mgs.count_reads(subtree, counts)
    virus_reads = {
        study: sum(virus_counts[s] for s in samples_by_study[study])
        for study in studies
    }

    estimates = pathogens[pathogen].estimate_prevalences()
    prevalence_per100k = estimates[0].infections_per_100k
    mu = np.log(prevalence_per100k)
    sigma = 1.0

    data = {
        "J": len(studies),
        "y": [virus_reads[s] for s in studies],
        "n": [all_reads[s] for s in studies],
        "mu": mu,
        "sigma": sigma,
    }
    naive = naive_ra_at_one_percent(data["y"], data["n"], np.exp(data["mu"]))
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
