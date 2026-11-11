// Stage 2: hierarchical model of latent usual household iron intake (per 1000 kcal of need),
// given household-visit log intake y (already corrected for the measurement model in stage 1).
// Outputs posterior prevalence of inadequate absorbed iron among surveyed WRA (probability approach).
data {
  int<lower=1> HV; int<lower=1> H; int<lower=1> R;
  array[HV] int<lower=1, upper=H> hh; // household of each visit
  vector[HV] y;                       // log( mg Fe/day per 1000 kcal need ), corrected
  array[H] int<lower=1, upper=R> reg;
  vector[H] z;                        // standardised log income
  // women of reproductive age
  int<lower=1> N;
  array[N] int<lower=1, upper=H> w_hh;
  array[N] int<lower=1, upper=R> w_reg;
  vector[N] w_need;                   // kcal/day
  vector[N] w_med; vector[N] w_sdlog; // requirement distribution (absorbed iron mg/day)
  real<lower=0> absorption;
  // informative priors on what the data cannot identify
  real<lower=0> share_sd_mean; real<lower=0> share_sd_sd; // sharing/bulk-buying noise (log sd)
  real<lower=0> alloc_sd_mean; real<lower=0> alloc_sd_sd; // within-household allocation noise (log sd)
}
parameters {
  vector[R] mu;          // regional mean of log usual intake
  real beta;             // income slope
  real<lower=0> tau;     // between-household sd of CONSUMED intake
  real<lower=0> sigma_w; // within-household (week-to-week + recall) sd
  real<lower=0> sigma_s; // sharing noise sd (affects purchases, not consumption)
  vector[H] c_raw;
  vector[H] s_raw;
}
transformed parameters {
  vector[H] logx;
  for (h in 1:H) logx[h] = mu[reg[h]] + beta * z[h] + tau * c_raw[h];
}
model {
  mu ~ normal(2, 2); beta ~ normal(0, 0.5);
  tau ~ normal(0, 1); sigma_w ~ normal(0, 1);
  sigma_s ~ normal(share_sd_mean, share_sd_sd);
  c_raw ~ std_normal(); s_raw ~ std_normal();
  for (i in 1:HV) y[i] ~ normal(logx[hh[i]] + sigma_s * s_raw[hh[i]], sigma_w);
}
generated quantities {
  real prev_all; vector[R] prev_reg; vector[N] p_i; vector[N] intake_i;
  {
    real sa = abs(normal_rng(alloc_sd_mean, alloc_sd_sd));
    vector[R] num = rep_vector(0, R); vector[R] den = rep_vector(0, R);
    for (i in 1:N) {
      real a = exp(normal_rng(-0.5 * sa^2, sa));
      real it = (w_need[i] / 1000) * exp(logx[w_hh[i]]) * a * absorption;
      intake_i[i] = it;
      p_i[i] = 1 - normal_cdf(log(it) | log(w_med[i]), w_sdlog[i]);
      num[w_reg[i]] += p_i[i]; den[w_reg[i]] += 1;
    }
    prev_all = sum(p_i) / N;
    for (r in 1:R) prev_reg[r] = den[r] > 0 ? num[r] / den[r] : not_a_number();
  }
}
