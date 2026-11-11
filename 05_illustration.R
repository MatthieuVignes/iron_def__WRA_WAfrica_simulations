# One worked example (p_other = 20%) used for figures (in the report; ~1 min to run)
for (f in c("R/00_params.R", "R/01_simulate.R", "R/02_estimators.R", "R/03_bayes.R")) source(f)
par <- modifyList(default_params(), list(p_other = 0.20))
s <- simulate_survey(par, seed = 7, table_seed = 7)
m <- get_stan_model()
e_naive <- est_naive(s); e_eb <- est_eb_other(s)
b0 <- est_bayes(s, use_other = FALSE, J = 4, model = m, seed = 7)
b1 <- est_bayes(s, use_other = TRUE,  J = 4, model = m, seed = 7)
saveRDS(list(par = par, wra = s$wra[, c("hh", "region", "status", "fe_abs", "req", "req_median", "req_sdlog", "inadequate", "gap")],
             truth = s$truth[c("prev_all", "prev_region", "frac_iron_other", "mean_fe_abs")],
             naive = e_naive, eb = e_eb, b0 = b0[c("prev_all", "lo", "hi", "draws", "wra")],
             b1 = b1[c("prev_all", "lo", "hi", "draws", "wra", "prev_region", "region_lo", "region_hi", "diag")]),
        "results/illustration.rds")
cat("truth", s$truth$prev_all, "| naive", e_naive$prev_all, "| B1", e_eb$prev_all,
    "| C0", b0$prev_all, "| C1", b1$prev_all, "[", b1$lo, b1$hi, "]\n")
