# ------------------------------------------------------------------------------
# Entry point for the normality-test power study (Project 2 / Task 1).
#
# Run from the task root with:
#     Rscript R/main.R
#
# Produces:
#   plots/*.png     - figures used in the report
#   output/*.csv    - power grid + p-value dump for the focus (n, K) cells
# ------------------------------------------------------------------------------

source("R/config.R")
source("R/tests.R")
source("R/simulation.R")
source("R/plots.R")

set.seed(RANDOM_SEED)
dir.create(PLOTS_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Parent distribution: Exp(1). Same family as Project 1/Task 2 so the
# theoretical moments carry over directly.
rng_exp <- function(k) rexp(k, rate = EXP_RATE)

# ---- Headline experiment: power grid over (n, K) ----------------------------
message("Running power grid (n x K) ...")
power_grid <- run_power_grid(
  sample_sizes = SAMPLE_SIZES,
  k_sizes      = K_SIZES,
  n_sim        = N_SIMULATIONS,
  rng          = rng_exp,
  mu           = EXP_MEAN,
  sigma2       = EXP_VAR,
  alpha        = ALPHA
)

write.csv(power_grid,
          file = file.path(OUTPUT_DIR, "power_grid.csv"),
          row.names = FALSE)

message("\n== Power grid (rejection rates at alpha = ", ALPHA, ") ==")
print(power_grid[, c("n", "K", "power_ks", "power_sw", "power_da",
                     "mean_skew", "mean_excess_kurt")],
      digits = 3, row.names = FALSE)

# ---- p-value dump for a small / medium / large n at K_DEFAULT ---------------
# These three replications give the reader a sense of the *distribution* of
# p-values - not just the rejection rate.
message("\nDumping p-values for selected n values at K = ", K_DEFAULT, " ...")
set.seed(RANDOM_SEED + 1L)
pdump_small  <- run_pvalue_dump(n = 5L,   K = K_DEFAULT,
                                n_sim = N_SIMULATIONS, rng = rng_exp,
                                mu = EXP_MEAN, sigma2 = EXP_VAR)
set.seed(RANDOM_SEED + 2L)
pdump_medium <- run_pvalue_dump(n = 30L,  K = K_DEFAULT,
                                n_sim = N_SIMULATIONS, rng = rng_exp,
                                mu = EXP_MEAN, sigma2 = EXP_VAR)
set.seed(RANDOM_SEED + 3L)
pdump_large  <- run_pvalue_dump(n = 200L, K = K_DEFAULT,
                                n_sim = N_SIMULATIONS, rng = rng_exp,
                                mu = EXP_MEAN, sigma2 = EXP_VAR)

write.csv(rbind(pdump_small, pdump_medium, pdump_large),
          file = file.path(OUTPUT_DIR, "pvalue_dump.csv"),
          row.names = FALSE)

# ---- One illustrative batch of K means at several n values ------------------
# Snapshot, not aggregated - to show what the tests actually "see".
set.seed(RANDOM_SEED + 100L)
illustrative_means <- list()
for (n in c(2L, 5L, 30L, 200L)) {
  illustrative_means[[as.character(n)]] <-
    draw_sample_means(n = n, K = K_DEFAULT, rng = rng_exp)
}

# ---- Figures -----------------------------------------------------------------
message("\nRendering figures ...")

save_plot(
  plot_means_density(
    illustrative_means, K = K_DEFAULT,
    title    = "Standaryzowane średnie X_bar_n - jedna realizacja"),
  "01_means_density.png"
)

save_plot(
  plot_power_vs_n(power_grid, K_focus = K_DEFAULT, alpha = ALPHA,
    title    = sprintf("Moce testów w funkcji n (K = %d, alpha = %.2f)",
                       K_DEFAULT, ALPHA),
    subtitle = "Linia przerywana: poziom istotności alpha (oczekiwana moc gdy H0 prawdziwa)"),
  "02_power_vs_n.png"
)

save_plot(
  plot_power_facets_K(power_grid, alpha = ALPHA,
    title    = "Wpływ wielkości próby testowej K na moc",
    subtitle = "Większe K = wiekszy zasób informacji dla testu = większa moc"),
  "03_power_facets_K.png"
)

save_plot(
  plot_shape_convergence(power_grid, K_focus = K_DEFAULT,
    title    = "Skośność i nadwyżka kurtozy X_bar_n - empiryczne vs. teoretyczne",
    subtitle = "Teoria: skew = 2/sqrt(n), excess kurt = 6/n (rozkład Exp(1))"),
  "04_shape_convergence.png"
)

save_plot(
  plot_dagostino_components(power_grid, K_focus = K_DEFAULT,
    title    = "Standaryzowane składowe testu D'Agostino-Pearson K^2",
    subtitle = "Linie przerywane: kwantyle +/- 1.96 (krytyczne wartości N(0,1) przy alpha = 0.05)"),
  "05_dagostino_components.png"
)

save_plot(
  plot_pvalue_distribution(pdump_small, alpha = ALPHA,
    title    = "Rozkład p-value przy n = 5",
    subtitle = sprintf("K = %d, alpha = %.2f - silna niezgodność z N: nagromadzenie p-value blisko 0",
                       K_DEFAULT, ALPHA)),
  "06_pvalues_n5.png"
)

save_plot(
  plot_pvalue_distribution(pdump_medium, alpha = ALPHA,
    title    = "Rozkład p-value przy n = 30",
    subtitle = sprintf("K = %d - regime pośredni: częściowy pile-up przy zera",
                       K_DEFAULT)),
  "07_pvalues_n30.png"
)

save_plot(
  plot_pvalue_distribution(pdump_large, alpha = ALPHA,
    title    = "Rozkład p-value przy n = 200",
    subtitle = sprintf("K = %d - CTG zadziałało: p-value przybliżają rozkład Uniform(0,1)",
                       K_DEFAULT)),
  "08_pvalues_n200.png"
)

message("Done. Outputs written to '", PLOTS_DIR, "/' and '", OUTPUT_DIR, "/'.")
