## Statistical model

The `stats` module uses the [Stan](https://mc-stan.org/) statistical programming language to define a Bayesian model of viral prevalences and metagenomic sequencing read counts.
Our objective is to estimate a coefficient that 

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

### Parameters

### Model

### Limitations and future directions

- Partial pooling of prevalence estimates
- Assessing model fit
- Multiple studies/methodologies
