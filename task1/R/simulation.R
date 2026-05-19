# ------------------------------------------------------------------------------
# Monte-Carlo routines for the normality-test power study.
#
# A single replication generates K sample means - each computed from n iid
# Exp(1) draws - and runs the three normality tests on this K-vector. Power
# at level alpha is the fraction of replications whose p-value falls below
# alpha. We also record the empirical skewness / kurtosis so we can comment
# on the components driving the D'Agostino omnibus statistic.
# ------------------------------------------------------------------------------

# Generate one batch of K sample means of size n drawn from `rng`.
# Drawing the n*K observations as a single matrix is faster than a Python-style
# inner loop.
draw_sample_means <- function(n, K, rng) {
  x   <- rng(n * K)
  mat <- matrix(x, nrow = n, ncol = K)
  colMeans(mat)
}

# One Monte-Carlo replication: run all three tests on a single K-sample,
# return p-values and the skewness / kurtosis used by D'Agostino.
single_replication <- function(n, K, rng, mu, sigma2) {
  means <- draw_sample_means(n, K, rng)
  sigma_mean <- sqrt(sigma2 / n)

  ks <- ks_normal_test(means, mu = mu, sigma = sigma_mean)
  sw <- sw_normal_test(means)
  da <- dagostino_test(means)

  list(
    ks_p  = ks$p.value,
    sw_p  = sw$p.value,
    da_p  = da$p.value,
    g1    = da$g1,
    b2    = da$b2,
    Z1    = da$Z1,
    Z2    = da$Z2
  )
}

# Run `n_sim` replications for fixed (n, K) and aggregate into a one-row tibble
# of empirical powers and average shape estimates.
run_power_cell <- function(n, K, n_sim, rng, mu, sigma2, alpha) {
  ks_p <- numeric(n_sim)
  sw_p <- numeric(n_sim)
  da_p <- numeric(n_sim)
  g1   <- numeric(n_sim)
  b2   <- numeric(n_sim)
  Z1   <- numeric(n_sim)
  Z2   <- numeric(n_sim)

  for (i in seq_len(n_sim)) {
    r       <- single_replication(n, K, rng, mu, sigma2)
    ks_p[i] <- r$ks_p
    sw_p[i] <- r$sw_p
    da_p[i] <- r$da_p
    g1[i]   <- r$g1
    b2[i]   <- r$b2
    Z1[i]   <- r$Z1
    Z2[i]   <- r$Z2
  }

  data.frame(
    n               = n,
    K               = K,
    n_sim           = n_sim,
    alpha           = alpha,
    power_ks        = mean(ks_p < alpha,  na.rm = TRUE),
    power_sw        = mean(sw_p < alpha,  na.rm = TRUE),
    power_da        = mean(da_p < alpha,  na.rm = TRUE),
    mean_skew       = mean(g1, na.rm = TRUE),
    mean_kurt       = mean(b2, na.rm = TRUE),
    mean_excess_kurt = mean(b2 - 3, na.rm = TRUE),
    mean_Z1         = mean(Z1, na.rm = TRUE),
    mean_Z2         = mean(Z2, na.rm = TRUE)
  )
}

# Sweep over a grid of (n, K) values - this is the headline experiment table.
run_power_grid <- function(sample_sizes, k_sizes, n_sim, rng, mu, sigma2, alpha) {
  rows <- list()
  idx  <- 1L
  for (K in k_sizes) {
    for (n in sample_sizes) {
      message(sprintf("  n = %3d, K = %3d ...", n, K))
      rows[[idx]] <- run_power_cell(n, K, n_sim, rng, mu, sigma2, alpha)
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

# Generate the raw p-value samples (one full Monte-Carlo run for a single
# (n, K) pair). Used to inspect the *distribution* of p-values, not just the
# rejection rate, on the diagnostic plots.
run_pvalue_dump <- function(n, K, n_sim, rng, mu, sigma2) {
  ks_p <- numeric(n_sim)
  sw_p <- numeric(n_sim)
  da_p <- numeric(n_sim)
  for (i in seq_len(n_sim)) {
    r       <- single_replication(n, K, rng, mu, sigma2)
    ks_p[i] <- r$ks_p
    sw_p[i] <- r$sw_p
    da_p[i] <- r$da_p
  }
  data.frame(
    n   = n,
    K   = K,
    rep = seq_len(n_sim),
    KS  = ks_p,
    SW  = sw_p,
    DA  = da_p
  )
}
