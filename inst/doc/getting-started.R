## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(spliv)
set.seed(42)
n <- 240
z <- rnorm(n)
w <- rnorm(n)
exposure <- pnorm(w)
inactive <- seq_len(n) <= n / 2
x <- ifelse(inactive, 0, 1) * z + 0.4 * w + rnorm(n)
y <- 1.2 * x + 0.25 * w + 0.15 * exposure * z + rnorm(n)
d <- data.frame(y, x, z, w, exposure, inactive)
f <- y ~ x + w | z + w

## ----baseline-----------------------------------------------------------------
baseline <- spliv(f, d, vcov = "hc1")
baseline$estimates
uniform <- spliv(f, d, method = "uci", delta = 0.20, vcov = "hc1",
                 grid = list(steps = 9))
uniform$estimates

## ----pattern------------------------------------------------------------------
pattern <- spliv_pattern(
  name = "Exposure pattern", pattern = ~ exposure,
  rationale = "The alternative channel is stronger at higher exposure.",
  variables_used = "exposure", pattern_type = "theory_defined",
  normalize = "max_abs"
)
spliv_eval_pattern(pattern, d)[1:5]
patterned_uci <- spliv(f, d, method = "uci", delta = 0.20, vcov = "hc1",
                       violation_pattern = pattern, grid = list(steps = 9))
patterned_ltz <- spliv(f, d, method = "ltz", delta = 0.20, vcov = "hc1",
                       violation_pattern = pattern)
patterned_uci$estimates
patterned_ltz$estimates

## ----paths--------------------------------------------------------------------
path <- spliv_sensitivity_path(
  f, d, method = "uci", delta_grid = seq(0, 0.30, by = 0.05),
  vcov = "hc1", violation_pattern = pattern
)
head(path)
spliv_tipping_point(path)

## ----bpe----------------------------------------------------------------------
design <- bpe_design(
  name = "Theory-defined inactive subset", subset = ~ inactive,
  rationale = "The treatment channel is absent in the inactive subset.",
  variables_used = "inactive", subset_type = "theory_defined",
  pre_specified = TRUE,
  transportability_rationale = "The subset direct effect is informative for the target sample."
)
# Illustrative synthetic margin: 0.25 residual treatment SD per
# one-residual-SD instrument shift. In substantive work, pre-specify the
# margin rather than tuning it to obtain BPE eligibility.
bpe_margin <- 0.25
validation <- bpe_validate_design(
  f, d, design = design, vcov = "hc1",
  bpe_min_n_S = 40,
  bpe_equiv_margin = bpe_margin
)
validation[c("n_S", "equivalence_passed", "eligibility_passed")]

bpe_fit <- spliv(f, d, method = "bpe", bpe_design = design,
                 vcov = "hc1", bpe_min_n_S = 40,
                 bpe_equiv_margin = bpe_margin)
bpe_fit$estimates

