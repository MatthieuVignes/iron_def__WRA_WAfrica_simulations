# Which error source costs how much accuracy?  (~2 min to run)
# Switch sources off one at a time (naive estimator A, p_other = 5%)
# Observation is the last step of the simulator, so for a given seed the TRUTH 
#  is identical across observation-side variants (clean paired comparison)
for (f in c("R/00_params.R", "R/01_simulate.R", "R/02_estimators.R")) source(f)
base <- default_params()
no_share <- list(p_cluster = 0, share_out = 0, nonres_out = 0, bulk_buyer_prob = 0)
variants <- list(
  "Everything on (baseline)"          = list(),
  "No omissions"                      = list(omit_rate = 0),
  "No recall error"                   = list(qty_err_sd = 0),
  "All quantities reported in kg"     = list(p_unit_reported = 0),
  "No survey error at all"            = list(omit_rate = 0, qty_err_sd = 0, p_unit_reported = 0),
  "No sharing / outflow"              = no_share,
  "No error AND no sharing"           = c(no_share, list(omit_rate = 0, qty_err_sd = 0, p_unit_reported = 0)))
nrep <- 12
res <- do.call(rbind, lapply(names(variants), function(v) {
  do.call(rbind, lapply(seq_len(nrep), function(r) {
    s <- simulate_survey(modifyList(base, variants[[v]]), seed = 500 + r, table_seed = 500 + r)
    e <- est_naive(s)
    data.frame(variant = v, rep = r, truth = s$truth$prev_all, est = e$prev_all,
               iron_ratio = mean(e$wra$intake) / s$truth$mean_fe_abs)
  }))
}))
saveRDS(res, "results/bias_budget.rds")
print(aggregate(cbind(bias = est - truth, iron_ratio) ~ variant, res, mean))
