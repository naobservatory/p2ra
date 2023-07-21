## Statistical model

The `stats` module uses the [Stan](https://mc-stan.org/) statistical programming language to fit a Bayesian model of viral prevalences and metagenomic sequencing read counts.
Our objective is to estimate a coefficient that connects a viral clade's estimated prevalence to the number of metagenomic reads assigned to that clade across a set of samples.
For each virus and study, we estimate both an overall coefficient and a coefficient that is specific to each stampling location in the study.

The stan code defining the model is contained in `model.stan`.
In the following sections, we will show and explain each block of the code in turn.

### Data

The `data` block defines the data that we provide to the model:

```stan
data {
  int<lower=1> J;           // number of samples
  array[J] int<lower=0> y;  // viral read counts
  array[J] int<lower=0> n;  // total read counts
  vector[J] x;              // estimated predictor (prevalence or incidence)
  int<lower=1> L;           // number of sampling locations
  array[J] int<lower=1, upper=L> ll;  // sampling locations
}
```

In particular, we must provide:

- The number metagenomic samples, `J`
- The number of reads assigned to our focal virus in each sample, $y_j$ with $j = 1, \ldots, J$. 
- The total number of reads in each sample, $n_j$.
- The estimated public health predictor (prevalence or incidence) of the virus in the population contributing to each sample, $\mu_j$.
- The number of sampling locations in the study, `L`
- The location of each sample `ll`. These are provided as integer indexes ranging from 1 to `L`.

### Transformed data

The `transformed data` block defines quantities derived from the data that are used in the `model` block:

```stan
transformed data {
  vector[J] x_std = log(x) - mean(log(x));
  real log_mean_y = 0;
  if (sum(y) > 0)           // can't normalize by this if there are no viral reads
    log_mean_y = log(mean(y));
  real log_mean_n = log(mean(n));
}
```

Because our logistic regression model works with the log predictors, we log-transform and mean-center the predictors to get `x_std`.

In order to center the coefficient, we also need a central value for the read counts.
We use the log of the average viral and total read counts, because the viral read counts have a lot of zeros.

### Parameters

When we run the program, Stan uses a [Hamiltonian Monte Carlo](https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo) algorithm to sample from the posterior distribution of the model parameters, defined in the `parameters` block: 

```stan
parameters {
  vector[J] theta_std;      // standardized true predictor for each sample
  real<lower=0> sigma;      // standard deviation of true predictors
  real mu;                  // mean P2RA coefficient (on standardized scale)
  real<lower=0> tau;        // std of P2RA coefficients per location
  vector[L] b_l;            // P2RA coefficient per location
}
```

### Model

The `model` block defines the model to be fit:

```stan
model {
  sigma ~ gamma($sigma_alpha, $sigma_beta);
  theta_std ~ normal(x_std, sigma);
  mu ~ normal(0, $mu_sigma);
  tau ~ gamma($tau_alpha, $tau_beta);
  b_l ~ normal(mu, tau);
  for (j in 1:J){
    y[j] ~ binomial_logit(n[j], b_l[ll[j]] + theta_std[j] + log_mean_y - log_mean_n);
  }
}
```

The first five lines give prior distributions for the parameters:

* The variance parameteters `sigma` and `tau` are given gamma priors with hyperparameters supplied at runtime.
* The true standardized predictor for each sample `theta_std` is given a normal prior centered on the estimated value
* The coefficients linking predictors to relative abundance are given a hierarchical model, where the overall coefficient `mu` has a prior centered at zero (because of the mean-centering) and the location-specific coefficients are centered at `mu`.

The last two lines define the likelihood.
We assume that the read counts follow a binomial distribuion and are independent, conditional on the parameters.
The expected relative abundance in each sample is given by the inverse logit of the sum of:

* `b_l[ll[j]]`, the coefficient for the location of the sample
* `theta_std[j]`, the true value of the (standardized) public health predictor
* `log_mean_y - log_mean_n`, included to center the coefficients near zero.

### Generated quantities

```stan
generated quantities {
  // posterior predictive viral read counts
  array[J] int<lower=0> y_tilde;
  for (j in 1:J){
    y_tilde[j] =
      binomial_rng(
        n[j],
        inv_logit(b_l[ll[j]] + theta_std[j] + log_mean_y - log_mean_n)
      );
  }
  // posterior true prevalence for each sample
  vector[J] theta = theta_std + mean(log(x));
  // for convenience, a single vector with the location coefficients and
  // overall coefficient in the final position
  vector[L + 1] b;
  b[:L] = b_l;
  b[L + 1] = mu;
  // location-specific expected relative abundance
  // last element is the overall coefficient
  // Converting from 1:100K to 1:1K means multiplying by 100
  vector[L + 1] ra_at_1in1000 = inv_logit(
    b - mean(log(x)) + log_mean_y - log_mean_n + log(100)
  );
}
```

### Limitations and future directions

Our biggest current limitation is that we model the true prevalences as independent of one another.
It would be more realistic for them to have shared as well as independent components of variation.
For example, our prevalence estimates may be systematically biased in one direction or another and that would induce a correlated effect on all of the samples.
The amount of shared variation between samples might also depend on features they have in common such as a shared sampling location or similar sampling times.
(Issue [#59](https://github.com/naobservatory/p2ra/issues/59) outlines a simple way forward here.)

A related issue is that we will usually have prevalence estimates for multiple times and places that don't line up exactly with the sampling times and locations.
In this case, we may want to have a stochastic process model (e.g., a Gaussian process) of the true prevalence that can be updated with our prevalence estimates and then emit the parameters for the negative binomial regressions.

As mentioned in the [Model](#model) section, it would also be useful to extend the negative binomial regression to include predictors other than sample identity.
By including sample processing methods, or viral characteristics as predictors, we could combine information across studies and viruses to develop a general model of MGS abundance. 

Finally, we currently do not do any evaluation of model fit or sensitivity to prior assumptions.
We will develop a workflow to check that the model is appropriate for the data and to improve it if not.
This will include steps like examining the mixing of the Markov chains during the HMC runs and comparing the posterior predictive distributions of read counts to the observed counts.
