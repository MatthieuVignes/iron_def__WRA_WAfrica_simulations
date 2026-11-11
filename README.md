# Iron-intake inadequacy from household food-purchase surveys: a simulation study

This is a simulation framework (R + Stan) to quantify how accurately we can estimate 
the prevalence of **inadequate absorbed iron intake in women of reproductive age (WRA)** 
from a household food-purchase survey, and in particular how much accuracy is lost by:

1. ignoring the "general food / other" category (prepared meals, take-away) and
2. assuming a (random) composition for it,

under realistic survey error (omissions, recall error, bowl/bag units), food sharing between
neighbouring households, and non-resident feeding.

The generating model is inspired by a Ghana household survey (4 fictional regions, 
~400 households, 2-4 weekly visits). **All numerical inputs are placeholders** (see 
[Limitations](#limitations)) to be replaced them with real food-composition, price 
and unit-weight data by an expert.

## What the project does

* **Generates** households and a latent weekly food demand over 12 food groups 
 plus an "other" (prepared food) category whose composition is a random mix of the 12; 
 sharing between neighbouring households, bulk buying recorded under one buyer, and 
 food sent to non-resident relatives; allocation within the household; and 
 each woman's true absorbed iron intake and individual requirement.
* **Derives a noisy survey** from it: omitted items, recall error, quantities in 
 local units (bowl, bag) with unknown regional conversions, and noisy spending.
* **Compares estimators**, from naive to usual-intake to Bayesian (Stan), on bias, 
 RMSE and 95% credible-interval coverage of the true prevalence.
* **Varies** the share of expenditure on "other", survey accuracy, sharing structure and 
 sample size, to answer how much accuracy is lost by ignoring "other", by assuming a 
 composition for it, and by survey error and sharing.
* **Reports** the results in a Quarto document (to be turned into a tutorial?)

**Headline result of the shipped (small) run:** ignoring "other" costs little, when it 
 is about 5% of spending (about 2 percentage points) and matters from roughly 10% upwards 
 (about 8 points at 20%). Imputing "other" removes that dependence on its share, but the 
 overall level of the estimate is driven more by systematic survey errors (food leaving 
 the household, omissions) than by "other", and a larger sample does not shrink it.
 The shipped results use 12 replicates per scenario (6 for Bayesian ones), so differences 
 below about 3 points are within Monte-Carlo noise; treat them as illustrative until rerun 
 with more replicates. (and placeholder inputs are replaced)

> **Scope.** Food data inform *inadequate iron intake*, not iron *deficiency* (ferritin, anaemia),
> which also depends on infection, inflammation and blood loss. The estimand is the prevalence of
> absorbed iron intake below individual requirement, by the probability approach.

## Contents

```
R/00_params.R        parameters (default_params) and the 12-group food table (food_table)
R/01_simulate.R      generative model + survey observation model  -> simulate_survey()
R/02_estimators.R    naive and usual-intake estimators (A, A+, B0, B1, B1w)
R/03_bayes.R         modular Bayesian estimator calling Stan (C0, C1, C1w) -> est_bayes()
R/04_study.R         scenario library, replicate runner, parallel + resumable study driver
stan/intake_hier.stan  hierarchical model of latent usual household iron intake
install_deps.R       installs the required CRAN packages
run_study.R          command-line entry point (quick / full)
05_illustration.R    one worked example for the report figures
06_bias_budget.R     switch error sources off one at a time (paired seeds)
report.qmd           Quarto report (rationale, methods, results)
results/             saved outputs used by the report
docs/                TUTORIAL_GUIDE.md and a runnable Quarto-book template (docs/tutorial_template)
```

## Requirements

* R >= 4.2, a C++ toolchain (needed by `rstan`), and [Quarto](https://quarto.org) >= 1.4 to render the report
* R packages: `rstan`, `dplyr`, `tidyr`, `ggplot2`, `scales`, `future`, `furrr`, `knitr`, `rmarkdown`

Notes: the Stan model is compiled once (~1 minute) and cached next to the `.stan` file.
On Debian/Ubuntu, if compilation fails with "Boost not found", install `libboost-dev` and make sure the
`BH` R package provides its headers (a plain `install.packages("BH")` from CRAN does).

## Quick start

Run everything from the project root (paths are relative).

```bash
Rscript run_study.R quick        # ~3 min sanity check (2 scenarios, a few replicates)
Rscript run_study.R full         # complete study; resumable, see below
Rscript 05_illustration.R        # worked example for the report (~1 min)
Rscript 06_bias_budget.R         # error-source decomposition (~2 min)
quarto render report.qmd         # builds report.html from results/*.rds
```

A single simulated survey from the R console:

```r
for (f in c("R/00_params.R", "R/01_simulate.R", "R/02_estimators.R")) source(f)
s <- simulate_survey(default_params(), seed = 1)
s$truth$prev_all                 # true prevalence among simulated WRA
est_naive(s)$prev_all            # estimate from the survey alone, ignoring "other"
est_eb_other(s)$prev_all         # usual-intake estimate including "other"
```

### Running the full study within a time budget

Each replicate is cached in `results/cache_<mode>/`; a chunk stops starting new replicates after the
given number of minutes, and the next call picks up where it left off:

```bash
Rscript run_study.R full 8       # run for ~8 minutes, then stop
Rscript run_study.R full 8       # ... repeat until it prints "Done"
```

Replicates per scenario and the Stan sampler settings are set in `run_study.R` (`cfg`) and
`run_replicate()`. The numbers shipped in `results/` use 12 replicates (6 for Bayesian 
scenarios) on 2 cores; increase them for high (publication?) precision.

## Estimators

| Label | Description |
|---|---|
| A | Nominal unit conversions and iron table, household mean intake, ignore "other" |
| A+ | A plus iron imputed for "other" from its spend and a fixed assumed composition |
| B0 / B1 | Usual-intake (empirical-Bayes) estimator, without / with "other"; corrects assumed omission and outflow |
| B1w | B1 with a deliberately wrong composition for "other" |
| C0 / C1 / C1w | Bayesian (rstan): prior-based propagation of unit, composition, "other", omission and outflow uncertainty (stage 1), then a hierarchical Stan model of latent usual intake (stage 2) |

Metrics: bias, RMSE and 95% credible-interval coverage of the true prevalence, overall and by region.

## Main parameters to vary

| Parameter | Meaning |
|---|---|
| `n_hh`, `n_visits`, `n_regions` | design (visits drawn uniformly in the range) |
| `p_other`, `other_center`, `other_conc` | share of expenditure on "other" and its composition |
| `omit_rate`, `qty_err_sd`, `p_unit_reported`, `spend_err_sd` | survey accuracy |
| `p_cluster`, `share_out`, `bulk_buyer_prob`, `nonres_out` | sharing network and non-resident feeding |
| `alloc_conc`, `women_bias` | within-household allocation |
| `absorption`, `absorption_sd` | iron absorption (constant now; `absorption_sd > 0` gives woman-specific values) |
| `food_table()` | number of food groups K and their composition, price, unit weight and diet shape |

Scenarios are defined in `make_scenarios()` (`R/04_study.R`); add your own by extending it.

## Limitations

* **Placeholder inputs:** iron content, prices, bowl/bag weights, requirements and 
 diet shapes are approximate, not taken from the West African food composition table or 
 from survey data.
* Own production, gifts and meals eaten out are not modelled; all intake comes from purchases.
* Bioavailability is constant by default; diet- or person-dependent absorption is a 
 natural extension.
* All members of a household eat the same diet (allocation changes quantity, not food type).
* The Bayesian model is modular ("cut"): the data cannot correct the calibration priors, 
 so interval honesty depends on those priors.
* The shipped results come from few replicates; differences of a few percentage points 
 are within Monte-Carlo noise.

## Reproducibility

Seeds are fixed per replicate and shared across scenarios (common random numbers), and 
the "true" composition table, prices and unit conversions use a separate `table_seed`, 
so scenario comparisons are paired. Parallel runs use `future`/`furrr`.

## Tutorial companion

`report.qmd` is written as a results report. 
[`docs/TUTORIAL_GUIDE.md`](docs/TUTORIAL_GUIDE.md) explains how to turn it into a 
tutorial that walks through the simulations, and 
[`docs/tutorial_template/`](docs/tutorial_template/) is a runnable Quarto book with 
one complete chapter (`cd docs/tutorial_template && quarto render`).

## Citation and licence

Add a `LICENSE` (for example MIT for code and CC-BY for the report) and, if you wish, a `CITATION.cff` before making the repository public.
