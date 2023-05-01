## Statistical model

The `stats` module uses the [Stan](https://mc-stan.org/) statistical programming language to define a Bayesian model of viral prevalences and metagenomic sequencing read counts.
Our objective is to estimate a coefficient that connects estimates of the prevalence of a viral clade to the number of metagenomic reads assigned to that clade across a set of samples.

The model is a work in progress.
The goal so far has been to have something that is a reasonable representation of the problem that we can feed our various data sources into.
Now that we have data and working glue code, we will work to assess and improve the model.

Here is the Stan code.
In the following section, we will explain each component.
If you are not familiar with Stan, you may skip to the Data section.

```stan
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
```

### Data

The input to the model is data on read counts and viral prevalence.
In particular, we must provide:

- The number metagenomic samples, $J$
- The number of reads assigned to our focal viral clade in each sample, $y_j$ with $j = 1, \ldots, J$. 
- The total number of reads in each sample, $n_j$.
- The expected log-prevalence of the viral clade in the population contributing to each sample, $\mu_j$.
- An estimate of the uncertainty in the log-prevalence, represented as a standard deviation $\sigma$, assumed to be the same across samples.

### Parameters

When we run the program, Stan uses a [Hamiltonian Monte Carlo](https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo) algorithm to sample from the posterior distribution of the model parameters. 
They are:

- The true log-prevalence in the population contributing to each sample, $\theta_j$.
- The coefficient linking the log-prevalence to the expected number of viral reads, $b$. Note that this is a scalar quantity, assumed to be the same across all samples.
- An inverse overdispersion parameter $\phi$, representing the extra variation in read counts beyond a Poisson distribution. Larger $\phi$ is more Poisson-like, smaller $\phi$ is more overdispersed.

### Model

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

This is equivalent to assuming a flat prior on the log-abundance and updating with Gaussian-likelihood to our prevalence data.

We put a weakly informative Gaussian prior on $b$.
(We can center this at zero by normalizing abundance inputs.)

Finally, we give $\phi$ a Gamma prior with a mode at $\phi = 1$.
(This is largely arbitrary at this stage, but enforces $\phi > 0$.)

### Limitations and future directions

- Priors
- Independence of estimates
- Partial pooling of prevalence estimates
- Assessing model fit
- Multiple studies/methodologies
