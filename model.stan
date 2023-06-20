data {
  int<lower=1> J;           // number of samples
  array[J] int<lower=0> y;  // viral read counts
  array[J] int<lower=0> n;  // total read counts
  vector[J] x;              // estimated predictor (prevalence or incidence)
  int<lower=1> L;           // number of sampling locations
  array[J] int<lower=1, upper=L> ll;  // sampling locations
}
transformed data {
  vector[J] x_std = log(x) - mean(log(x));
  real log_mean_y = 0;
  if (sum(y) > 0)           // can't normalize by this if there are no viral reads
    log_mean_y = log(mean(y));
  real log_mean_n = log(mean(n));
}
parameters {
  vector[J] theta_std;      // standardized true predictor for each sample
  real<lower=0> sigma;      // standard deviation of true predictors
  real mu;                  // mean P2RA coefficient (on standardized scale)
  real<lower=0> tau;        // std of P2RA coefficients pre location
  vector[L] b_l;            // P2RA coefficient per location
}
model {
  sigma ~ gamma(2, 1);
  theta_std ~ normal(x_std, sigma);
  mu ~ normal(0, 4);
  tau ~ gamma(2, 1);
  b_l ~ normal(mu, tau);
  for (j in 1:J){
    y[j] ~ binomial_logit(n[j], b_l[ll[j]] + theta_std[j] + log_mean_y - log_mean_n);
  }
}
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
  real b = mu - mean(log(x)) + log_mean_y - log_mean_n;
  vector[L] b_loc = b_l - mean(log(x)) + log_mean_y - log_mean_n;
  // posterior P2RA coefficient
  real ra_at_1in1000 = inv_logit(b + log(100));
}
