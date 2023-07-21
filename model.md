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

```stan
transformed data {
  vector[J] x_std = log(x) - mean(log(x));
  real log_mean_y = 0;
  if (sum(y) > 0)           // can't normalize by this if there are no viral reads
    log_mean_y = log(mean(y));
  real log_mean_n = log(mean(n));
}
```

### Parameters

```stan
parameters {
  vector[J] theta_std;      // standardized true predictor for each sample
  real<lower=0> sigma;      // standard deviation of true predictors
  real mu;                  // mean P2RA coefficient (on standardized scale)
  real<lower=0> tau;        // std of P2RA coefficients pre location
  vector[L] b_l;            // P2RA coefficient per location
}
```

When we run the program, Stan uses a [Hamiltonian Monte Carlo](https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo) algorithm to sample from the posterior distribution of the model parameters. 
They are:

- The true log-prevalence of the focal clade in the population contributing to each sample, $\theta_j$.
- The coefficient linking the log-prevalence to the expected number of viral reads, $b$. Note that this is a scalar quantity, assumed to be the same across all samples.
- An inverse overdispersion parameter $\phi$, representing the extra variation in read counts beyond a Poisson distribution. Larger $\phi$ is more Poisson-like, smaller $\phi$ is more overdispersed.

### Model

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

We assume that the read counts follow a negative binomial distribution and are independent (conditional on the parameters):

$$
y_i \sim NB(n_i \exp(b + \theta_i), \phi),
$$

using the (mean, reciprocal overdispersion) parameterization.
Note that expected number of reads is proportional to the total number of reads times the exponentiated log-prevalence.
The constant of proportionality is $e^b$.
In other words, $b$ is the intercept term of a negative binomial regression with a log link function.
If we wanted to include other predictors of the abundance, $x$, we would multiply the mean by $e^{\beta x}$. 

We give our true log-prevalence parameters independent Gaussian priors centered on our log-prevalence estimates:

$$
\theta_i \sim Normal(\mu_i, \sigma).
$$

This is equivalent to assuming a flat prior on the log-prevalence and updating with Gaussian-likelihood to our prevalence data.

We put a weakly informative Gaussian prior on $b$.
We can center the prior at zero by choosing the right scale for prevalence inputs.
That is, we can convert prevalence to cases per $n$ individuals and choose $n$ so that when $\theta = \log prevalence = 0$, we expect to see $O(1)$ reads.

Finally, we give $\phi$ a Gamma prior with a mode at $\phi = 1$.
(This is largely arbitrary at this stage, but enforces $\phi > 0$.)

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
