# Modular Bayesian estimator (rstan).
# Stage 1: draw measurement-model uncertainty from priors (unit conversions, iron composition,
#   composition of "other", omission, outflow) and rebuild household-visit iron intake J times
# Stage 2: for each draw, fit the Stan hierarchical model (latent usual intake, within-household
#   variance, sharing noise) and compute posterior prevalence. Pool draws across the J fits

library(rstan)
rstan_options(auto_write = TRUE)

get_stan_model <- function(stan_file = "stan/intake_hier.stan") {
  rds <- sub("\\.stan$", ".rds", stan_file)
  if (file.exists(rds) && file.mtime(rds) >= file.mtime(stan_file)) readRDS(rds) else stan_model(stan_file)
}

beta_ms <- function(m, s) { nu <- max(m * (1 - m) / s^2 - 1, 2); rbeta(1, m * nu, (1 - m) * nu) }

# Analyst's prior beliefs (extends default_assumed)
default_prior <- function(par) c(default_assumed(par), list(
  other_conc = par$other_conc, markup_sd = 0.25, omit_sd = 0.02, out_sd = 0.03,
  share_sd_sd = 0.05, alloc_sd_sd = 0.08, unit_prior_sd = par$unit_prior_sd))

# Matrices of reported quantities by household-visit (built once)
build_mats <- function(s, pr) {
  K <- nrow(s$par$ft); hv <- pr$hv; o <- pr$o
  row <- match(paste(o$hh, o$visit), paste(hv$hh, hv$visit))
  kgm <- matrix(0, nrow(hv), K); um <- matrix(0, nrow(hv), K)
  isk <- o$type == "kg"
  kgm[cbind(row[isk], o$food[isk])] <- o$value[isk]
  um[cbind(row[!isk], o$food[!isk])] <- o$value[!isk]
  list(kgm = kgm, um = um)
}

# One stage-1 draw -> corrected log intake y (length = n household-visits)
stage1_draw <- function(s, pr, mats, prior, use_other) {
  ft <- s$par$ft; K <- nrow(ft); R <- s$par$n_regions; hv <- pr$hv
  fe <- ft$fe_100g * exp(rnorm(K, -log(1 + ft$fe_cv^2) / 2, sqrt(log(1 + ft$fe_cv^2))))
  U <- outer(rep(1, R), ft$unit_kg) * exp(matrix(rnorm(R * K, 0, prior$unit_prior_sd), R, K))
  kg <- mats$kgm + mats$um * U[hv$region, , drop = FALSE]
  fe_known <- as.numeric(kg %*% (fe * 10))
  fe_other <- 0
  if (use_other) {
    g <- rgamma(K, pmax(prior$other_conc * prior$center, 1e-6)); w <- g / sum(g)
    mk <- prior$markup * exp(rnorm(1, 0, prior$markup_sd))
    lam <- sapply(1:R, function(r) sum(w * fe * 10 / (pr$price_hat[r, ] * mk)))
    fe_other <- prior$other_mult * hv$other_spend * lam[hv$region]
  }
  omit <- beta_ms(prior$omit, prior$omit_sd); out <- beta_ms(prior$out, prior$out_sd)
  M <- (fe_known + fe_other) / 7 / (hv$need / 1000) * (1 - out) / (1 - omit)
  log(pmax(M, 1e-3))
}

est_bayes <- function(s, prior = default_prior(s$par), use_other = TRUE, J = 8,
                      iter = 600, warmup = 300, model = get_stan_model(), seed = 1, adapt_delta = 0.9) {
  set.seed(seed)
  pr <- prep_survey(s, prior); mats <- build_mats(s, pr); hv <- pr$hv; par <- s$par
  W <- s$wra
  dat0 <- list(HV = nrow(hv), H = nrow(s$hh), R = par$n_regions, hh = hv$hh, reg = s$hh$region, z = s$hh$z_inc,
               N = nrow(W), w_hh = W$hh, w_reg = W$region, w_need = W$need,
               w_med = W$req_median, w_sdlog = W$req_sdlog, absorption = par$absorption,
               share_sd_mean = prior$share_sd, share_sd_sd = prior$share_sd_sd,
               alloc_sd_mean = prior$alloc_sd, alloc_sd_sd = prior$alloc_sd_sd)
  pa <- vector("list", J); pr_reg <- vector("list", J); pi_mean <- 0; it_med <- 0; diag <- NULL
  for (j in seq_len(J)) {
    dat <- c(dat0, list(y = stage1_draw(s, pr, mats, prior, use_other)))
    fit <- sampling(model, data = dat, chains = 1, iter = iter, warmup = warmup, seed = seed * 100 + j,
                    refresh = 0, control = list(adapt_delta = adapt_delta),
                    pars = c("prev_all", "prev_reg", "p_i", "intake_i"))
    ex <- rstan::extract(fit)
    pa[[j]] <- ex$prev_all; pr_reg[[j]] <- ex$prev_reg
    pi_mean <- pi_mean + colMeans(ex$p_i) / J; it_med <- it_med + apply(ex$intake_i, 2, median) / J
    sm <- summary(fit, pars = "prev_all")$summary
    diag <- rbind(diag, data.frame(j = j, rhat = sm[, "Rhat"], ess = sm[, "n_eff"],
                                   ndiv = sum(get_sampler_params(fit, inc_warmup = FALSE)[[1]][, "divergent__"])))
  }
  pall <- unlist(pa); preg <- do.call(rbind, pr_reg)
  list(prev_all = mean(pall), lo = unname(quantile(pall, .025)), hi = unname(quantile(pall, .975)),
       draws = pall, prev_region = colMeans(preg),
       region_lo = apply(preg, 2, quantile, .025), region_hi = apply(preg, 2, quantile, .975),
       wra = data.frame(hh = W$hh, region = W$region, status = W$status, intake = it_med, p = pi_mean),
       diag = diag)
}
