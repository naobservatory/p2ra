## Statistical model

The `stats` module uses the [Stan](https://mc-stan.org/) statistical programming language to define a Bayesian model of viral prevalences and metagenomic sequencing read counts.
Our objective is to estimate a coefficient that connects estimates of the prevalence of a viral clade to the number of metagenomic reads assigned to that clade across a set of samples.

The model is a work in progress.
The goal so far has been to have something that is a reasonable representation of the problem that we can feed our various data sources into.
Now that we have data and working glue code, we will work to assess and improve the model.

Here is the Stan code.
In the following sections, we will explain each component.
If you are not familiar with Stan, you may skip to the [Data](#data) section.

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

Our biggest current limitation is that we model the true prevalences as independent of one another.
It would be more realistic for them to have shared as well as independent components of variation.
For example, our prevalence estimates may be systematically biased in one direction or another and that would induce a correlated effect on all of the samples.
The amount of shared variation between samples might also depend on features they have in common such as a shared sampling location or similar sampling times.
(Issue [#59](https://github.com/naobservatory/p2ra/issues/59) outlines a simple way forward here.)

A related issue is that we may have prevalence data that is more relevant for some samples than for others, but is not perfectly aligned with the sampling times and locations.
In this case, we may want to have a stochastic process model (e.g., a Gaussian process) of the true prevalence that can be updated with our prevalence estimates and then emit the parameters for the negative binomial regressions.

As mentioned in the [Model](#model) section, it would also be useful to extend the negative binomial regression to include predictors other than sample identity.
By including sample processing methods, or viral characteristics as predictors, we could combine information across studies and viruses to develop a general model of MGS abundance. 

Finally, we currently do not do any evaluation of model fit or sensitivity to prior assumptions.
We will develop a workflow to check that the model is appropriate for the data and to improve it if not.
