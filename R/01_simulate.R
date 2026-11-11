# Generative model: households -> true purchases -> sharing -> (...)
# (...) true individual iron intake -> noisy survey 
# Returns survey data (what the analyst sees) and hidden truth (what the simulation knows: oracle)

rdirichlet1 <- function(alpha) { g <- rgamma(length(alpha), pmax(alpha, 1e-8)); g / sum(g) }
rbeta_ms <- function(n, m, s) { # Beta with given mean and sd
  if (m <= 0) return(rep(0, n))
  nu <- max(m * (1 - m) / s^2 - 1, 2); rbeta(n, m * nu, (1 - m) * nu)
}
lnorm_unbiased <- function(n, sd) exp(rnorm(n, -sd^2 / 2, sd)) # mean-1 lognormal noise
sdlog_from_cv <- function(cv) sqrt(log(1 + cv^2))

simulate_survey <- function(par = default_params(), seed = NULL, table_seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ft <- par$ft; K <- nrow(ft); R <- par$n_regions; H <- par$n_hh; W <- par$n_weeks_truth
  stopifnot(length(par$other_center) == K)
  center <- par$other_center / sum(par$other_center)

  # 1. households and people
  region <- sample.int(R, H, replace = TRUE)
  hh_n   <- pmin(1 + rpois(H, par$hh_size_mean - 1), 14)
  people <- do.call(rbind, lapply(seq_len(H), function(h) {
    n <- hh_n[h]
    cls <- sample(c("child", "adult", "older"), n, TRUE, prob = c(.38, .47, .15))
    age <- ifelse(cls == "child", sample(0:14, n, TRUE),
           ifelse(cls == "adult", sample(15:49, n, TRUE), sample(50:80, n, TRUE)))
    sex <- sample(c("F", "M"), n, TRUE)
    data.frame(hh = h, region = region[h], sex = sex, age = age, stringsAsFactors = FALSE)
  }))
  people$wra <- people$sex == "F" & people$age >= 15 & people$age <= 49
  u <- runif(nrow(people))
  people$status <- ifelse(!people$wra, NA_character_,
                   ifelse(u < par$p_pregnant, "pregnant",
                   ifelse(u < par$p_pregnant + par$p_lactating, "lactating", "other")))
  base_need <- with(people, ifelse(age < 15, 1000 + 65 * age,
                     ifelse(age <= 49, ifelse(sex == "F", 2100, 2700), ifelse(sex == "F", 1900, 2300))))
  people$need <- base_need + ifelse(is.na(people$status), 0,
                   ifelse(people$status == "pregnant", 340, ifelse(people$status == "lactating", 500, 0)))
  hh_need <- as.numeric(tapply(people$need, people$hh, sum))           # kcal/day

  # income (log per capita), lower in North
  inc_mu <- c(North = 4.0, East = 4.4, South = 4.8, West = 4.6)[par$region_names][region]
  log_inc <- rnorm(H, inc_mu, 0.6); z_inc <- as.numeric(scale(log_inc))

  # 2. sharing network
  cluster <- rep(NA_integer_, H); cl_id <- 0L
  in_cl <- sample(H, round(par$p_cluster * H))
  # clusters formed within region (neighbours)
  for (r in seq_len(R)) {
    pool <- in_cl[region[in_cl] == r]; pool <- pool[sample.int(length(pool))]   # (not sample(): breaks for length-1 pools)
    while (length(pool) >= par$cluster_size[1]) {
      sz <- min(sample(par$cluster_size[1]:par$cluster_size[2], 1), length(pool))
      cl_id <- cl_id + 1L; cluster[pool[seq_len(sz)]] <- cl_id; pool <- pool[-seq_len(sz)]
    }
  }
  out_frac <- pmin(par$nonres_out * lnorm_unbiased(H, 0.5), 0.5)       # share of purchases leaving to non-residents

  # 3. true composition (unknown to the analyst)
  # the "true" composition table, prices and unit conversions can use their own seed (table_seed) so that
  # they stay identical across scenarios (common random numbers); main RNG restored afterwards
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  if (!is.null(table_seed)) set.seed(table_seed)
  fe_true <- ft$fe_100g * exp(rnorm(K, -sdlog_from_cv(ft$fe_cv)^2 / 2, sdlog_from_cv(ft$fe_cv)))
  price_r <- outer(rep(1, R), ft$price) * exp(matrix(rnorm(R * K, 0, par$price_sd_region), R, K))
  unit_true <- outer(rep(1, R), ft$unit_kg) * exp(matrix(rnorm(R * K, 0, rep(ft$unit_sd, each = R)), R, K))
  if (!is.null(table_seed) && !is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())

  # 4. latent weekly demand
  reg_tilt <- as.matrix(ft[, paste0("reg_", par$region_names)]) # K x R
  pi_h <- t(sapply(seq_len(H), function(h) {
    lw <- log(ft$base) + reg_tilt[, region[h]] + ft$inc_slope * z_inc[h]
    rdirichlet1(par$diet_conc * exp(lw) / sum(exp(lw)))
  })) # H x K mean energy shares
  Dk <- array(0, c(H, W, K)); Do <- array(0, c(H, W, K)); Ospend <- matrix(0, H, W)
  for (h in seq_len(H)) for (w in seq_len(W)) {
    E <- hh_need[h] * 7 * par$purchase_ratio_mean * exp(0.1 * z_inc[h]) * lnorm_unbiased(1, par$energy_ratio_sd)
    s <- rbeta_ms(1, par$p_other, par$p_other_sd)
    pi_hw <- rdirichlet1(par$visit_diet_conc * pi_h[h, ])
    kg <- E * (1 - s) * pi_hw / (ft$kcal_100g * 10)
    kg[kg < par$min_purchase_kg] <- 0 # negligible amounts = item not bought that week
    pr <- price_r[region[h], ] * lnorm_unbiased(K, 0.1)
    Ospend[h, w] <- s / (1 - s) * sum(kg * pr)
    Dk[h, w, ] <- kg
    Do[h, w, ] <- Ospend[h, w] * rdirichlet1(par$other_conc * center) / (pr * par$other_markup)
  }

  # 5. sharing: pooled consumption and who is recorded as buyer
  cons <- Dk; Qbuy <- Dk
  for (cl in seq_len(cl_id)) {
    m <- which(cluster == cl); wt <- hh_need[m] / sum(hh_need[m])
    for (w in seq_len(W)) for (k in seq_len(K)) {
      d <- Dk[m, w, k]; pooled <- sum(d) * wt
      if (runif(1) < par$bulk_buyer_prob * ft$bulk_prone[k]) {         # bulk purchase recorded at one buyer
        b <- sample.int(length(m), 1); Qbuy[m, w, k] <- 0; Qbuy[m[b], w, k] <- sum(d)
        cons[m, w, k] <- pooled
      } else {
        cons[m, w, k] <- (1 - par$share_out) * d + par$share_out * pooled
      }
    }
  }
  Qbuy <- Qbuy / (1 - out_frac) # purchases also cover non-resident relatives (recorded, not eaten here)
  Ospend_rec <- Ospend / (1 - out_frac)
  cons_all <- cons + Do # prepared food is eaten in the household

  # 6. true individual iron intake and requirement
  fe_hw <- apply(cons_all, c(1, 2), function(x) sum(x * fe_true * 10)) # mg/week, whole household
  fe_known_hw <- apply(cons, c(1, 2), function(x) sum(x * fe_true * 10))
  fe_hh_day <- rowMeans(fe_hw) / 7 # usual intake, mg/day
  frac_other <- 1 - sum(fe_known_hw) / sum(fe_hw)
  people$share <- NA_real_
  for (h in seq_len(H)) {
    ix <- which(people$hh == h); n <- length(ix)
    wgt <- people$need[ix] * exp(par$women_bias * (people$wra[ix]))
    people$share[ix] <- if (n == 1) 1 else rdirichlet1(par$alloc_conc * wgt / sum(wgt))
  }
  absr <- par$absorption * lnorm_unbiased(nrow(people), par$absorption_sd)
  people$fe_abs <- people$share * fe_hh_day[people$hh] * absr # absorbed mg/day
  med <- c(pregnant = par$req_pregnant_median, lactating = par$req_lact_median, other = par$req_other_median)
  cvr <- c(pregnant = par$req_cv[["preg"]], lactating = par$req_cv[["lact"]], other = par$req_cv[["other"]])
  st <- ifelse(is.na(people$status), "other", people$status)
  people$req_median <- med[st]; people$req_sdlog <- sdlog_from_cv(cvr[st])
  people$req <- ifelse(people$wra, rlnorm(nrow(people), log(people$req_median), people$req_sdlog), NA)
  people$inadequate <- ifelse(people$wra, people$fe_abs < people$req, NA)
  people$gap <- people$fe_abs - people$req # discrepancy: intake - requirement
  wra <- people[people$wra, ]
  truth <- list(
    prev_all = mean(wra$inadequate),
    prev_region = tapply(wra$inadequate, wra$region, mean),
    mean_fe_abs = mean(wra$fe_abs), frac_iron_other = frac_other,
    fe_true = fe_true, unit_true = unit_true, price_r = price_r, out_frac = out_frac, cluster = cluster
  )

  # 7. observation process (what the survey records)
  nv <- sample(par$n_visits[1]:par$n_visits[2], H, TRUE)
  obs <- vector("list", H); oth <- vector("list", H)
  omit_k <- pmin(par$omit_rate * ft$omit_mult, 0.9)
  for (h in seq_len(H)) {
    weeks <- sort(sample.int(W, nv[h])); r <- region[h]
    for (v in seq_along(weeks)) {
      q <- Qbuy[h, weeks[v], ]; keep <- which(q > 0 & runif(K) > omit_k)
      if (length(keep)) {
        as_unit <- runif(length(keep)) < par$p_unit_reported
        e <- rnorm(length(keep), 0, par$qty_err_sd)
        kg_rep <- q[keep] * exp(e)
        val <- ifelse(as_unit, pmax(0.5, round(2 * kg_rep / unit_true[r, keep]) / 2), round(kg_rep, 1))
        obs[[h]] <- rbind(obs[[h]], data.frame(hh = h, visit = v, week = weeks[v], region = r, food = keep,
          type = ifelse(as_unit, "unit", "kg"), value = val,
          spend = q[keep] * price_r[r, keep] * exp(rnorm(length(keep), 0, par$spend_err_sd))))
      }
      o_rep <- if (runif(1) < par$omit_rate / 2) 0 else Ospend_rec[h, weeks[v]] * exp(rnorm(1, 0, par$spend_err_sd))
      oth[[h]] <- rbind(oth[[h]], data.frame(hh = h, visit = v, week = weeks[v], region = r, other_spend = o_rep))
    }
  }
  hh <- data.frame(hh = seq_len(H), region = region, n_people = hh_n, need = hh_need,
                   log_inc = log_inc, z_inc = z_inc, cluster = cluster, n_visits = nv)
  list(par = par, hh = hh, people = people, wra = wra,
       obs = do.call(rbind, obs), other = do.call(rbind, oth), truth = truth)
}
