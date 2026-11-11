# Tunable settings and the food table
# NOTE: iron / energy / price values are PLAUSIBLE PLACEHOLDERS (order of magnitude
# of West African food composition), to be replaced by the West African FCT / GLSS data.

# Food table (12 food groups; K is derived from this table)
# kcal_100g: energy per 100 g purchased; fe_100g: iron mg per 100 g (mean);
# heme: fraction of iron that is heme; price: GHS/kg baseline; bulk_share: relative
# chance an item is bought in bulk (sharing-prone)
food_table <- function() {
  data.frame(
    food     = c("rice","maize_flour","millet_sorghum","roots_tubers","legumes_nuts",
                 "leafy_veg","other_veg","fruit","meat_poultry","fish","dairy_eggs","oil_sugar"),
    kcal_100g= c(360, 350, 350, 130, 380, 30, 30, 50, 200, 130, 130, 700),
    fe_100g  = c(1.0, 2.5, 4.5, 0.8, 6.0, 3.0, 1.0, 0.4, 2.0, 1.5, 1.2, 0.1),
    fe_cv    = c(.25, .25, .25, .30, .25, .35, .35, .35, .25, .30, .30, .50), # composition uncertainty
    heme     = c(0, 0, 0, 0, 0, 0, 0, 0, .40, .30, .10, 0),
    price    = c(8, 6, 6, 3, 9, 4, 5, 4, 28, 22, 18, 12),
    unit_kg  = c(1.2, 1.0, 1.0, 2.5, 0.9, 0.3, 0.4, 0.8, 1.0, 0.8, 0.5, 0.9), # kg per "bowl/bag/bunch"
    unit_sd  = c(.15, .15, .15, .25, .15, .30, .30, .25, .15, .20, .20, .15),  # between-region sd (log)
    # diet shape: baseline energy share, regional log-tilts, income slope (per sd of log income)
    base     = c(.12, .14, .06, .24, .07, .01, .01, .03, .03, .03, .02, .24),
    reg_North= c(-.3,  .2, 1.2,-.8,  .6, .4, .1,-.3,-.2,-.7,-.2, 0),
    reg_East = c(-.2,  .5, .2, .2,  .1, 0, 0, 0, 0,-.2, 0, 0),
    reg_South= c( .2, -.2,-.5, .4,  0,  0, 0, .3, 0, .7, 0, 0),
    reg_West = c( .5, -.3,-.5, .2,  0,  0, 0, .2, 0, .3, 0, 0),
    inc_slope= c(.3, -.1,-.3,-.3,   0,  0, .1, .4, .5, .3, .5, .1),
    bulk_prone=c(.8, .8, .8, .2, .6, 0, 0, 0, .1, .1, .1, .7),   # prob. multiplier for bulk-buying events
    omit_mult =c(.5, .5, .5, .6, .8,1.6,1.4,1.6, 1, 1, 1.2, .8), # relative omission rate per food
    stringsAsFactors = FALSE
  )
}

# Default simulation parameters
default_params <- function() {
  ft <- food_table()
  list(
    # design
    n_hh = 400, n_regions = 4, n_visits = c(2, 4), # visits per household drawn uniformly in this range
    n_weeks_truth = 12, # latent weeks simulated per household (usual intake = mean over them)
    region_names = c("North","East","South","West"),
    ft = ft, K = nrow(ft),
    # population
    hh_size_mean = 5.5,
    p_pregnant = 0.10, p_lactating = 0.15,
    # diet / purchase
    energy_ratio_sd = 0.20,    # sd (log) of purchased energy / household need (visit level)
    purchase_ratio_mean = 1.0, # energy available to residents / energy need (own production ignored here)
    other_markup = 2.5,        # prepared food costs this many times the raw-food price per kg
    other_conc = 15,           # Dirichlet concentration of the "other" spend composition around its centre
    other_center = c(0.22, 0.10, 0.04, 0.06, 0.07, 0.05, 0.08, 0.02, 0.14, 0.12, 0.03, 0.07), # spend shares (sum=1)
    diet_conc = 25,            # Dirichlet concentration of energy shares (lower = more heterogeneity)
    visit_diet_conc = 60,      # week-to-week diet variation around household mean diet
    p_other = 0.05,            # share of expenditure going to "general food/other"
    p_other_sd = 0.02,         # sd of that share across households/visits
    # sharing / networks
    p_cluster = 0.20,       # proportion of households belonging to a sharing cluster
    cluster_size = c(2, 4), # cluster size range
    share_out = 0.15,       # fraction of a cluster-member's purchase consumed within network (beyond own hh)
    nonres_out = 0.07,      # fraction leaving hh to non-resident relatives (all hh)
    bulk_buyer_prob = 0.5,  # prob. a cluster bulk purchase is recorded at ONE buyer (vs split)
    # within-household allocation
    women_bias = 0,  # log-multiplier on women's shares (0 = none, <0 = discriminate)
    alloc_conc = 60, # Dirichlet noise for allocation (lower = noisier)
    # physiology / iron
    min_purchase_kg = 0.05, # weekly amounts below this are treated as "not bought"
    absorption = 0.14,      # constant fraction of dietary iron absorbed (non-heme + heme mix average)
    absorption_sd = 0,      # sd (log) of individual absorption (0 = constant; step 2 of the project)
    req_pregnant_median = 3.2, req_lact_median = 1.0, req_other_median = 1.46, # absorbed iron mg/day
    req_cv = c(other = 0.30, preg = 0.30, lact = 0.30),
    # survey accuracy (observation model)
    omit_rate = 0.08,      # prob. an item purchased in a visit is forgotten
    qty_err_sd = 0.25,     # sd (log) multiplicative recall error on quantities
    p_unit_reported = 0.5, # prob. quantity reported in local unit (bowl, bag) rather than kg
    unit_prior_sd = 0.30,  # analyst's prior sd (log) on the true unit conversion
    spend_err_sd = 0.15,   # sd (log) on reported expenditure
    price_sd_region = 0.20 # sd (log) of region-specific prices
  )
}
