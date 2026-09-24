// Counts within a parent population, beta-binomial, with a fixed-effects design matrix X
// (instrument dummies and protected covariates, K may be 0), unit random intercept and
// per-unit time slope, day random intercept and an optional stimulation random intercept.
// Layer B (has_stim = 1) and layer BG counts (has_stim = 0).
data {
  int<lower=1> N;
  int<lower=1> P; int<lower=1> D; int<lower=1> C;
  int<lower=0> K;
  int<lower=0, upper=1> has_stim;
  array[N] int<lower=1, upper=P> unit;
  array[N] int<lower=1, upper=D> day;
  array[N] int<lower=1, upper=C> stim;   // all 1 when has_stim = 0
  matrix[N, K] X;
  vector[N] week_c;            // centred, in years
  array[N] int<lower=0> n_parent;
  array[N] int<lower=0> k;
}
parameters {
  real alpha;
  vector[K] beta;
  real mu_week;
  vector[P] z_p; vector[P] z_w; vector[D] z_d; vector[C] z_c;
  real<lower=0> tau_p; real<lower=0> tau_w; real<lower=0> tau_d; real<lower=0> tau_c;
  real<lower=0> phi;           // beta-binomial precision
}
transformed parameters {
  vector[P] a_p = tau_p * z_p;
  vector[P] w_p = mu_week + tau_w * z_w;
  vector[D] c_d = tau_d * z_d;
  vector[C] s_c = has_stim * tau_c * z_c;
}
model {
  vector[N] eta = alpha + a_p[unit] + w_p[unit] .* week_c + c_d[day] + s_c[stim];
  if (K > 0) eta += X * beta;
  alpha ~ normal(-5, 2.5);
  beta ~ normal(0, 1);
  mu_week ~ normal(0, 1);
  z_p ~ std_normal(); z_w ~ std_normal(); z_d ~ std_normal(); z_c ~ std_normal();
  tau_p ~ normal(0, 1); tau_w ~ normal(0, 1); tau_d ~ normal(0, 1); tau_c ~ normal(0, 1);
  phi ~ gamma(2, 0.02);
  {
    vector[N] m = inv_logit(eta);
    k ~ beta_binomial(n_parent, m * phi, (1 - m) * phi);
  }
}
