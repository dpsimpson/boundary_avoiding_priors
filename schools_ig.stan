data {
  int<lower=1> J;
  vector[J] y;
  vector<lower=0>[J] sigma;
  real<lower=0> tau_shape;
  real<lower=0> tau_rate;
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
  tau ~ inv_gamma(tau_shape, tau_rate);
}
