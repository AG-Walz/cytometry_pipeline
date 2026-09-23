// Layer B: cytokine-positive counts within a parent population, beta-binomial.
// Protected covariates: prior vaccination (fixed), week (patient-specific slope), cocktail (random).
data {
  int<lower=1> N;
  int<lower=1> P; int<lower=1> D; int<lower=1> C;
  array[N] int<lower=1, upper=P> patient;
  array[N] int<lower=1, upper=D> day;
  array[N] int<lower=1, upper=C> cocktail;
  vector[N] instr;
  vector[N] prior_vacc;        // 1 = previously vaccinated
  vector[N] week_c;            // centred, in years
  array[N] int<lower=0> n_parent;
  array[N] int<lower=0> k;
}
parameters {
  real alpha;
  real b_instr;
  real b_vacc;
  real mu_week;
  vector[P] z_p; vector[P] z_w; vector[D] z_d; vector[C] z_c;
  real<lower=0> tau_p; real<lower=0> tau_w; real<lower=0> tau_d; real<lower=0> tau_c;
  real<lower=0> phi;           // beta-binomial precision
}
transformed parameters {
  vector[P] a_p = tau_p * z_p;
  vector[P] w_p = mu_week + tau_w * z_w;
  vector[D] c_d = tau_d * z_d;
  vector[C] s_c = tau_c * z_c;
}
model {
  alpha ~ normal(-5, 2.5);
  b_instr ~ normal(0, 1);
  b_vacc ~ normal(0, 1);
  mu_week ~ normal(0, 1);
  z_p ~ std_normal(); z_w ~ std_normal(); z_d ~ std_normal(); z_c ~ std_normal();
  tau_p ~ normal(0, 1); tau_w ~ normal(0, 1); tau_d ~ normal(0, 1); tau_c ~ normal(0, 1);
  phi ~ gamma(2, 0.02);
  for (n in 1:N) {
    real m = inv_logit(alpha + a_p[patient[n]] + w_p[patient[n]] * week_c[n] + c_d[day[n]]
                       + s_c[cocktail[n]] + b_instr * instr[n] + b_vacc * prior_vacc[n]);
    k[n] ~ beta_binomial(n_parent[n], m * phi, (1 - m) * phi);
  }
}
