# Analyst-side estimators: only use the SURVEY data (households, composition, 
#   reported quantities/spend) + published tables
# Prevalence of inadequacy uses the probability approach: mean_i P(R_i > intake_i).

req_prob <- function(intake, med, sdlog) 1 - plnorm(intake, log(med), sdlog)

# What the analyst assumes (can be mis-specified on purpose for sensitivity analyses)
default_assumed <- function(par) list(
  omit = par$omit_rate, # assumed omission rate
  out = par$nonres_out, # assumed share of purchases leaving the household
  share_sd = 0.15,      # assumed sd (log) of between-household sharing noise
  alloc_sd = 0.25,      # assumed sd (log) of within-household allocation noise
  markup = par$other_markup, # price markup of prepared food
  center = par$other_center / sum(par$other_center), # assumed composition of "other" (spend shares)
  other_mult = 1        # global multiplier on the imputed iron of "other" (sensitivity)
)

# Shared data preparation
prep_survey <- function(s, assumed = default_assumed(s$par)) {
  ft <- s$par$ft; K <- nrow(ft); R <- s$par$n_regions
  o <- s$obs
  o$kg <- ifelse(o$type == "kg", o$value, o$value * ft$unit_kg[o$food])      # nominal conversion
  o$fe_mg <- o$kg * ft$fe_100g[o$food] * 10                                   # nominal composition
  # regional prices estimated from kg-reported lines (spend / kg)
  kgl <- o[o$type == "kg" & o$kg > 0, ]
  ph <- matrix(NA_real_, R, K)
  for (r in 1:R) for (k in 1:K) {
    v <- kgl$spend[kgl$region == r & kgl$food == k] / kgl$kg[kgl$region == r & kgl$food == k]
    ph[r, k] <- if (length(v) >= 5) median(v) else median(kgl$spend[kgl$food == k] / kgl$kg[kgl$food == k])
  }
  # iron per GHS of prepared food, by region, under the assumed composition
  lam <- sapply(1:R, function(r) sum(assumed$center * ft$fe_100g * 10 / (ph[r, ] * assumed$markup)))
  hv <- s$other
  key <- paste(o$hh, o$visit); hvkey <- paste(hv$hh, hv$visit)
  hv$fe_known <- as.numeric(rowsum(o$fe_mg, key)[hvkey, 1]); hv$fe_known[is.na(hv$fe_known)] <- 0
  hv$fe_other <- hv$other_spend * lam[hv$region]
  hv$need <- s$hh$need[hv$hh]
  list(hv = hv, price_hat = ph, lambda = lam, o = o)
}

# Core estimator
# use_other: add imputed iron of "other"; correct: undo assumed omission / outflow;
# eb: usual-intake adjustment (shrink household means, account for within-hh and sharing variance)
est_core <- function(s, assumed = default_assumed(s$par), use_other = FALSE, correct = FALSE,
                     eb = FALSE, L = 100, seed = 1) {
  pr <- prep_survey(s, assumed); hv <- pr$hv; par <- s$par
  hv$fe_day <- (hv$fe_known + (use_other * assumed$other_mult) * hv$fe_other) / 7     # mg/day, household
  if (correct) hv$fe_day <- hv$fe_day * (1 - assumed$out) / (1 - assumed$omit)
  W <- s$wra; W$need_hh <- s$hh$need[W$hh]
  if (!eb) {
    hhm <- tapply(hv$fe_day, hv$hh, mean)
    W$intake <- W$need / W$need_hh * hhm[as.character(W$hh)] * par$absorption
    W$p <- req_prob(W$intake, W$req_median, W$req_sdlog)
  } else {
    set.seed(seed)
    hv$y <- log(pmax(hv$fe_day, 1e-3) / (hv$need / 1000))
    n_h <- c(tapply(hv$y, hv$hh, length)); m_h <- c(tapply(hv$y, hv$hh, mean))
    ss <- c(tapply(hv$y, hv$hh, function(z) sum((z - mean(z))^2)))
    sw2 <- sum(ss) / sum(n_h - 1) # within-household variance
    reg <- s$hh$region[as.integer(names(m_h))]
    mu_r <- as.numeric(c(tapply(m_h, reg, mean))[as.character(reg)])
    sb2_tot <- max(var(m_h - mu_r) - mean(sw2 / n_h), 0.01) # between-hh variance of logs
    sb2 <- max(sb2_tot - assumed$share_sd^2, 0.01) # remove sharing noise: consumption variance
    Bh <- sb2 / (sb2 + sw2 / n_h + assumed$share_sd^2)
    post_m <- mu_r + Bh * (m_h - mu_r); post_v <- (1 - Bh) * sb2
    hi <- match(W$hh, as.integer(names(m_h)))
    z <- matrix(rnorm(nrow(W) * L), nrow(W), L); a <- matrix(rnorm(nrow(W) * L), nrow(W), L)
    x <- exp(post_m[hi] + sqrt(post_v[hi]) * z) # mg/day per 1000 kcal need
    al <- exp(assumed$alloc_sd * a - assumed$alloc_sd^2 / 2)
    inta <- (W$need / 1000) * x * par$absorption * al
    W$p <- rowMeans(req_prob(inta, W$req_median, W$req_sdlog)); W$intake <- apply(inta, 1, median)
  }
  list(prev_all = mean(W$p), prev_region = tapply(W$p, W$region, mean), wra = W[, c("hh", "region", "status", "intake", "p")])
}

# Named estimators we use in the study
est_naive  <- function(s, ...) est_core(s, use_other = FALSE, correct = FALSE, eb = FALSE, ...)  # A: ignore "other"
est_other  <- function(s, ...) est_core(s, use_other = TRUE,  correct = FALSE, eb = FALSE, ...)  # A+: add "other"
est_eb_noother <- function(s, ...) est_core(s, use_other = FALSE, correct = TRUE, eb = TRUE, ...) # B0: corrected, usual intake
est_eb_other   <- function(s, ...) est_core(s, use_other = TRUE,  correct = TRUE, eb = TRUE, ...) # B1: + "other"
