data {
  int<lower=1> J;
  vector[J] y;
  vector<lower=0>[J] sigma;
  real<lower=0> tau_mu;
  real<lower=0> tau_sigma;
}
parameters {
  vector[J] theta;
  real mu;
  real<lower=0> tau;
}
model {
  y ~ normal(theta, sigma);
  theta ~ normal(mu, tau);

  mu ~ normal(0, 10);
  tau ~ lognormal(tau_mu, tau_sigma);
}
