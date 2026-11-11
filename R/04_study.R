# Scenario definitions and the replication
library(future)
library(furrr)

# Scenario library (overrides default_params)
accuracy_levels <- list(
  good = list(omit_rate = 0.03, qty_err_sd = 0.10, spend_err_sd = 0.08),
  base = list(omit_rate = 0.08, qty_err_sd = 0.25, spend_err_sd = 0.15),
  poor = list(omit_rate = 0.20, qty_err_sd = 0.50, spend_err_sd = 0.30))
sharing_levels <- list(
  none   = list(p_cluster = 0,   share_out = 0,    nonres_out = 0,    bulk_buyer_prob = 0),
  base   = list(p_cluster = 0.2, share_out = 0.15, nonres_out = 0.07, bulk_buyer_prob = 0.5),
  strong = list(p_cluster = 0.5, share_out = 0.35, nonres_out = 0.15, bulk_buyer_prob = 0.8))

make_scenarios <- function() {
  S <- list()
  add <- function(name, factor, level, ov, bayes = FALSE) S[[name]] <<- list(name = name, factor = factor, level = level, over = ov, bayes = bayes)
  for (p in c(0, .05, .10, .20, .30)) add(sprintf("p_other=%.2f", p), "p_other", p, list(p_other = p, p_other_sd = min(.02, p / 2 + 1e-9)),
                                         bayes = p %in% c(.05, .20))
  for (a in c("good", "poor")) add(paste0("accuracy=", a), "accuracy", a, accuracy_levels[[a]], bayes = FALSE)
  add("accuracy=poor,p=.20", "accuracy x p_other", "poor,0.20", c(accuracy_levels$poor, list(p_other = .20)), bayes = TRUE)
  for (h in c("none", "strong")) add(paste0("sharing=", h), "sharing", h, sharing_levels[[h]])
  for (n in c(100, 200, 800)) add(paste0("n_hh=", n), "n_hh", n, list(n_hh = n))
  S # baseline = "p_other=0.05" (all other settings at their defaults)
}

# deliberately wrong analyst belief about "other": staple-heavy centre (like the household diet)
misspec_center <- c(.25, .25, .10, .20, .05, .02, .02, .02, .02, .02, .02, .03)

# One replicate: simulate, run every estimator, return a tidy table
tidy_est <- function(e, name, truth, par) {
  R <- par$n_regions
  data.frame(estimator = name, region = c(0, 1:R),
             truth = c(truth$prev_all, as.numeric(truth$prev_region)),
             est = c(e$prev_all, as.numeric(e$prev_region)),
             lo = if (!is.null(e$lo)) c(e$lo, e$region_lo) else NA_real_,
             hi = if (!is.null(e$hi)) c(e$hi, e$region_hi) else NA_real_)
}

run_replicate <- function(scn, rep, base_seed = 2026, bayes = TRUE, J = 4, iter = 500, warmup = 250,
                          stan_file = "stan/intake_hier.stan") {
  par <- modifyList(default_params(), scn$over)
  seed <- base_seed + rep * 7919                      # same across scenarios (common random numbers)
  s <- simulate_survey(par, seed = seed, table_seed = base_seed + rep)
  asm <- default_assumed(par); wrong <- asm; wrong$center <- misspec_center / sum(misspec_center)
  out <- list(
    tidy_est(est_naive(s),                          "A: naive (ignore other)", s$truth, par),
    tidy_est(est_other(s),                          "A+: add other (plug-in)", s$truth, par),
    tidy_est(est_eb_noother(s),                     "B0: usual-intake, ignore other", s$truth, par),
    tidy_est(est_eb_other(s),                       "B1: usual-intake + other", s$truth, par),
    tidy_est(est_eb_other(s, assumed = wrong),      "B1w: usual-intake + other (wrong centre)", s$truth, par))
  frac_other <- s$truth$frac_iron_other
  if (bayes && scn$bayes) {
    m <- readRDS(sub("\\.stan$", ".rds", stan_file))
    pri <- default_prior(par); priw <- pri; priw$center <- wrong$center
    out <- c(out, list(
      tidy_est(est_bayes(s, pri,  use_other = FALSE, J = J, iter = iter, warmup = warmup, model = m, seed = seed %% 1e6), "C0: Bayes, ignore other", s$truth, par),
      tidy_est(est_bayes(s, pri,  use_other = TRUE,  J = J, iter = iter, warmup = warmup, model = m, seed = seed %% 1e6), "C1: Bayes + other", s$truth, par),
      tidy_est(est_bayes(s, priw, use_other = TRUE,  J = J, iter = iter, warmup = warmup, model = m, seed = seed %% 1e6), "C1w: Bayes + other (wrong centre)", s$truth, par)))
  }
  res <- do.call(rbind, out)
  res$scenario <- scn$name; res$factor <- scn$factor; res$level <- as.character(scn$level)
  res$rep <- rep; res$frac_iron_other <- frac_other; res$n_wra <- nrow(s$wra)
  res
}

# Parallel study
run_study <- function(scenarios, nrep = 20, bayes_nrep = nrep, workers = max(1, parallel::detectCores()),
                      bayes = TRUE, J = 4, out_file = "results/study.rds", verbose = TRUE,
                      cache_dir = "results/cache", max_minutes = NULL) {
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)     # resumable: finished replicates are cached
  deadline <- if (is.null(max_minutes)) NULL else Sys.time() + max_minutes * 60
  tasks <- do.call(rbind, lapply(scenarios, function(sc)
    data.frame(scn = sc$name, rep = seq_len(if (sc$bayes && bayes) bayes_nrep else nrep))))
  plan(multisession, workers = workers)
  t0 <- Sys.time()
  res <- future_map(seq_len(nrow(tasks)), function(i) {
    for (f in c("R/00_params.R", "R/01_simulate.R", "R/02_estimators.R", "R/03_bayes.R")) source(f)
    sc <- scenarios[[tasks$scn[i]]]
    cf <- file.path(cache_dir, sprintf("%s__%03d.rds", gsub("[^A-Za-z0-9=._,-]", "_", tasks$scn[i]), tasks$rep[i]))
    if (file.exists(cf)) return(readRDS(cf))
    if (!is.null(deadline) && Sys.time() > deadline) return(NULL)      # budget reached: leave for the next call
    r <- tryCatch(run_replicate(sc, tasks$rep[i], bayes = bayes, J = J), error = function(e) { message("FAIL ", tasks$scn[i], " ", tasks$rep[i], ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) saveRDS(r, cf)
    r
  }, .options = furrr_options(seed = TRUE, scheduling = 1), .progress = verbose && interactive())
  plan(sequential)
  complete <- !any(vapply(res, is.null, TRUE))
  res <- do.call(rbind, res)
  if (!complete) { message(sprintf("Partial: %d of %d replicates cached; call again to continue.", length(list.files(cache_dir)), nrow(tasks))); return(invisible(NULL)) }
  attr(res, "elapsed_min") <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  saveRDS(res, out_file)
  if (verbose) message(sprintf("Done: %d rows, %.1f min", nrow(res), attr(res, "elapsed_min")))
  res
}
