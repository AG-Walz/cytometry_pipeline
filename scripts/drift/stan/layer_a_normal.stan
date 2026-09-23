// Layer A: one neg-tube statistic, patient and day random intercepts, instrument fixed effect.
data {
  int<lower=1> N;
  int<lower=1> P;              // patients
  int<lower=1> D;              // acquisition days
  array[N] int<lower=1, upper=P> patient;
  array[N] int<lower=1, upper=D> day;
  vector[N] instr;             // 1 = non-reference instrument, 0 = reference
  vector[N] y;
}
parameters {
  real alpha;
  real b_instr;
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
  alpha ~ normal(0, 5);
  b_instr ~ normal(0, 1);
  z_p ~ std_normal();
  z_d ~ std_normal();
  tau_p ~ normal(0, 1);
  tau_d ~ normal(0, 1);
  sigma ~ normal(0, 1);
  y ~ normal(alpha + a_p[patient] + c_d[day] + b_instr * instr, sigma);
}
