# Usage:  Rscript run_study.R quick   (~2-4 min: sanity check)
#         Rscript run_study.R full    (target: <= 1 h on 4 cores)
args <- commandArgs(trailingOnly = TRUE); mode <- if (length(args)) args[1] else "quick"
for (f in c("R/00_params.R", "R/01_simulate.R", "R/02_estimators.R", "R/03_bayes.R", "R/04_study.R")) source(f)
invisible(get_stan_model())                          # compile once (cached in stan/intake_hier.rds)
sc <- make_scenarios()
cfg <- switch(mode,
  quick = list(nrep = 3,  bayes_nrep = 2,  J = 3, workers = parallel::detectCores(), keep = c("p_other=0.05", "p_other=0.20")),
  full  = list(nrep = 12, bayes_nrep = 6,  J = 4, workers = parallel::detectCores(), keep = names(sc)))
res <- run_study(sc[cfg$keep], nrep = cfg$nrep, bayes_nrep = cfg$bayes_nrep, workers = cfg$workers, J = cfg$J,
                 out_file = sprintf("results/study_%s.rds", mode), cache_dir = sprintf("results/cache_%s", mode),
                 max_minutes = if (length(args) > 1) as.numeric(args[2]) else NULL) # e.g. Rscript run_study.R full 8
if (!is.null(res)) print(aggregate(cbind(err = est - truth) ~ estimator + scenario, subset(res, region == 0), mean))
