# ------------------------------------------------------------------------------
# Normality tests used in the study.
#
# We implement three procedures applied to a sample of K sample-means:
#   1. Kolmogorov-Smirnov against a *fully specified* normal (mu = mu_T,
#      sigma^2 = sigma2_T / n), so the p-value is exact (no Lilliefors fix).
#   2. Shapiro-Wilk via the built-in stats::shapiro.test.
#   3. D'Agostino-Pearson omnibus K^2 test, hand-coded from D'Agostino (1970)
#      and Anscombe & Glynn (1983). Returns the standardised skewness Z1 and
#      kurtosis Z2 alongside the p-value so we can interpret the test through
#      the lens of the empirical third / fourth moments.
# ------------------------------------------------------------------------------

# Unbiased moment-based skewness sqrt(b1) and kurtosis b2 as defined in
# D'Agostino (1970). Population moments (divide by K, not K-1).
sample_skewness <- function(x) {
  K  <- length(x)
  xc <- x - mean(x)
  m2 <- mean(xc^2)
  m3 <- mean(xc^3)
  m3 / m2^1.5
}

sample_kurtosis <- function(x) {
  K  <- length(x)
  xc <- x - mean(x)
  m2 <- mean(xc^2)
  m4 <- mean(xc^4)
  m4 / m2^2
}

# Kolmogorov-Smirnov against N(mu, sigma) with KNOWN parameters.
# `mu` and `sigma` here are mu_T and sqrt(sigma2_T / n) of the sample-mean
# under the parent distribution - i.e. the theoretical limit law.
ks_normal_test <- function(x, mu, sigma) {
  res <- suppressWarnings(stats::ks.test(x, "pnorm", mean = mu, sd = sigma))
  list(statistic = unname(res$statistic), p.value = unname(res$p.value))
}

# Shapiro-Wilk - estimates mu/sigma internally, valid for 3 <= K <= 5000.
sw_normal_test <- function(x) {
  res <- stats::shapiro.test(x)
  list(statistic = unname(res$statistic), p.value = unname(res$p.value))
}

# ---- D'Agostino-Pearson K^2 omnibus test -------------------------------------
# Reference formulas: D'Agostino (1970) for the skewness transform,
# Anscombe & Glynn (1983) for the kurtosis transform, combined into the
# omnibus K^2 = Z1^2 + Z2^2 ~ chi^2(2) under H0.
#
# Validity caveat: the skewness transform assumes K >= 8 and the kurtosis
# transform is recommended for K >= 20. Below 20 the p-value is approximate
# and we warn the caller; we still return a value (NA for the few degenerate
# small-K cases) so the power table can be filled.

# Standardised skewness Z1(sqrt(b1)) - D'Agostino (1970).
dagostino_z_skew <- function(g1, K) {
  if (K < 8) return(NA_real_)
  mu2     <- 6 * (K - 2) / ((K + 1) * (K + 3))
  gamma2  <- 36 * (K - 7) * (K^2 + 2*K - 5) /
             ((K - 2) * (K + 5) * (K + 7) * (K + 9))
  W2      <- -1 + sqrt(2 * gamma2 + 4)
  delta   <- 1 / sqrt(log(sqrt(W2)))
  alpha   <- sqrt(2 / (W2 - 1))
  y       <- g1 / sqrt(mu2)
  delta * log(y / alpha + sqrt((y / alpha)^2 + 1))
}

# Standardised kurtosis Z2(b2) - Anscombe & Glynn (1983).
anscombe_z_kurt <- function(b2, K) {
  if (K < 20) return(NA_real_)
  mu1   <- 3 * (K - 1) / (K + 1)
  mu2   <- 24 * K * (K - 2) * (K - 3) /
           ((K + 1)^2 * (K + 3) * (K + 5))
  gamma1 <- 6 * (K^2 - 5*K + 2) / ((K + 7) * (K + 9)) *
            sqrt(6 * (K + 3) * (K + 5) /
                 (K * (K - 2) * (K - 3)))
  A     <- 6 + 8 / gamma1 * (2 / gamma1 + sqrt(1 + 4 / gamma1^2))
  x     <- (b2 - mu1) / sqrt(mu2)
  num   <- 1 - 2 / (9 * A)
  den   <- 1 + x * sqrt(2 / (A - 4))
  # Cube root must preserve sign; only den == 0 is degenerate.
  ratio <- if (den == 0) NA_real_ else (1 - 2 / A) / den
  z     <- (num - sign(ratio) * abs(ratio)^(1/3)) /
           sqrt(2 / (9 * A))
  z
}

# Omnibus D'Agostino-Pearson K^2 test. Returns the standardised components
# so callers can report skewness / kurtosis Z-values too.
dagostino_test <- function(x) {
  K  <- length(x)
  g1 <- sample_skewness(x)
  b2 <- sample_kurtosis(x)
  Z1 <- dagostino_z_skew(g1, K)
  Z2 <- anscombe_z_kurt(b2, K)
  if (is.na(Z1) || is.na(Z2)) {
    return(list(statistic = NA_real_, p.value = NA_real_,
                Z1 = Z1, Z2 = Z2, g1 = g1, b2 = b2))
  }
  Ksq <- Z1^2 + Z2^2
  list(
    statistic = Ksq,
    p.value   = stats::pchisq(Ksq, df = 2, lower.tail = FALSE),
    Z1        = Z1,
    Z2        = Z2,
    g1        = g1,
    b2        = b2
  )
}
