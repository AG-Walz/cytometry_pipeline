// Gaussian response with a fixed-effects design matrix X (instrument dummies and protected
// covariates, K may be 0), unit and day random intercepts. Layer A and layer BG intensity.
data {
  int<lower=1> N;
  int<lower=1> P;              // units (patients)
  int<lower=1> D;              // acquisition days
  int<lower=0> K;              // fixed-effect columns
  array[N] int<lower=1, upper=P> unit;
  array[N] int<lower=1, upper=D> day;
  matrix[N, K] X;
  vector[N] y;
}
parameters {
  real alpha;
  vector[K] beta;
  vector[P] z_p;
  vector[D] z_d;
  real<lower=0> tau_p;
  real<lower=0> tau_d;
  real<lower=0> sigma;
}
transformed parameters {
  vector[P] a_p = tau_p * z_p;
  vector[D] c_d = tau_d * z_d;
}
model {
  vector[N] mu = alpha + a_p[unit] + c_d[day];
  if (K > 0) mu += X * beta;
  alpha ~ normal(0, 5);
  beta ~ normal(0, 1);
  z_p ~ std_normal();
  z_d ~ std_normal();
  tau_p ~ normal(0, 1);
  tau_d ~ normal(0, 1);
  sigma ~ normal(0, 1);
  y ~ normal(mu, sigma);
}
