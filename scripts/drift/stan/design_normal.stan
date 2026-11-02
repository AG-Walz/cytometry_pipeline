// Gaussian response with a fixed-effects design matrix X (instrument dummies and protected
// covariates, K may be 0), unit and day random intercepts and an optional visit (sample) random
// intercept. Layer A and layer BG intensity.
data {
  int<lower=1> N;
  int<lower=1> P;              // units (patients)
  int<lower=1> D;              // acquisition days
  int<lower=0> K;              // fixed-effect columns
  int<lower=1> V;              // visits (1 when has_visit = 0)
  int<lower=0, upper=1> has_visit;
  array[N] int<lower=1, upper=P> unit;
  array[N] int<lower=1, upper=D> day;
  array[N] int<lower=1, upper=V> visit;   // all 1 when has_visit = 0
  matrix[N, K] X;
  vector[N] y;
}
parameters {
  real alpha;
  vector[K] beta;
  vector[P] z_p;
  vector[D] z_d;
  vector[V] z_v;
  real<lower=0> tau_p;
  real<lower=0> tau_d;
  real<lower=0> tau_v;
  real<lower=0> sigma;
}
transformed parameters {
  vector[P] a_p = tau_p * z_p;
  vector[D] c_d = tau_d * z_d;
  vector[V] e_v = has_visit * tau_v * z_v;
}
model {
  vector[N] mu = alpha + a_p[unit] + c_d[day] + e_v[visit];
  if (K > 0) mu += X * beta;
  alpha ~ normal(0, 5);
  beta ~ normal(0, 1);
  z_p ~ std_normal();
  z_d ~ std_normal();
  z_v ~ std_normal();
  tau_p ~ normal(0, 1);
  tau_d ~ normal(0, 1);
  tau_v ~ normal(0, 1);
  sigma ~ normal(0, 1);
  y ~ normal(mu, sigma);
}
