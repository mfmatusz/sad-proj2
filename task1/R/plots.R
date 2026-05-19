# ------------------------------------------------------------------------------
# Visualisation helpers for the normality-test power study.
# Same theme + palette as Project 1/Task 2 for visual continuity.
# ------------------------------------------------------------------------------

library(ggplot2)

# Shared theme (mirrors Project 1)
ntp_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "grey30"),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

# Fixed colour for each test - reused on every plot for legibility.
TEST_COLOURS <- c(KS = "#4C78A8", SW = "#F58518", DA = "#54A24B")
TEST_LABELS  <- c(KS = "Kołmogorow-Smirnow",
                  SW = "Shapiro-Wilk",
                  DA = "D'Agostino-Pearson K^2")

save_plot <- function(plot_obj, filename, width = PLOT_WIDTH_IN,
                      height = PLOT_HEIGHT_IN, dpi = PLOT_DPI) {
  path <- file.path(PLOTS_DIR, filename)
  ggsave(path, plot_obj, width = width, height = height, dpi = dpi)
  invisible(path)
}

# Long-format reshape of the power-grid table for ggplot.
power_long <- function(grid_df) {
  data.frame(
    n      = rep(grid_df$n, 3),
    K      = rep(grid_df$K, 3),
    test   = factor(rep(c("KS", "SW", "DA"), each = nrow(grid_df)),
                    levels = c("KS", "SW", "DA")),
    power  = c(grid_df$power_ks, grid_df$power_sw, grid_df$power_da)
  )
}

#' Power as a function of n for a fixed K (the headline figure).
plot_power_vs_n <- function(grid_df, K_focus, alpha,
                            title, subtitle = NULL) {
  df <- power_long(grid_df[grid_df$K == K_focus, ])
  ggplot(df, aes(x = n, y = power, colour = test, shape = test)) +
    geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.8) +
    scale_x_log10(breaks = sort(unique(df$n))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_colour_manual(values = TEST_COLOURS, labels = TEST_LABELS) +
    scale_shape_manual(values = c(KS = 16, SW = 17, DA = 15),
                       labels = TEST_LABELS) +
    labs(title = title, subtitle = subtitle,
         x = "Liczność próby n (skala log)",
         y = sprintf("Moc testu (alpha = %.2f)", alpha),
         colour = "Test", shape = "Test") +
    ntp_theme()
}

#' Power vs n faceted by K, so we can see how K influences detection.
plot_power_facets_K <- function(grid_df, alpha, title, subtitle = NULL) {
  df <- power_long(grid_df)
  df$K_label <- factor(paste0("K = ", df$K),
                       levels = paste0("K = ", sort(unique(df$K))))
  ggplot(df, aes(x = n, y = power, colour = test, shape = test)) +
    geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.4) +
    facet_wrap(~ K_label) +
    scale_x_log10(breaks = sort(unique(df$n))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_colour_manual(values = TEST_COLOURS, labels = TEST_LABELS) +
    scale_shape_manual(values = c(KS = 16, SW = 17, DA = 15),
                       labels = TEST_LABELS) +
    labs(title = title, subtitle = subtitle,
         x = "Liczność próby n (skala log)",
         y = sprintf("Moc testu (alpha = %.2f)", alpha),
         colour = "Test", shape = "Test") +
    ntp_theme()
}

#' Empirical skewness / kurtosis of X_bar_n vs. their theoretical curves.
#' For X ~ Exp(1): skew(X_bar_n) = 2/sqrt(n), excess-kurt(X_bar_n) = 6/n.
plot_shape_convergence <- function(grid_df, K_focus, title, subtitle = NULL) {
  df <- grid_df[grid_df$K == K_focus, ]
  long <- data.frame(
    n        = rep(df$n, 2),
    value    = c(df$mean_skew, df$mean_excess_kurt),
    quantity = factor(rep(c("skewness", "excess kurtosis"),
                          each = nrow(df)),
                      levels = c("skewness", "excess kurtosis"))
  )
  theory <- data.frame(
    n        = rep(df$n, 2),
    value    = c(EXP_SKEW / sqrt(df$n), EXP_KURT_E / df$n),
    quantity = factor(rep(c("skewness", "excess kurtosis"),
                          each = nrow(df)),
                      levels = c("skewness", "excess kurtosis"))
  )
  ggplot(long, aes(x = n, y = value)) +
    geom_line(data = theory, aes(linetype = "teoretyczna"),
              colour = "#E45756", linewidth = 0.9) +
    geom_line(aes(linetype = "empiryczna"),
              colour = "#4C78A8", linewidth = 0.9) +
    geom_point(colour = "#4C78A8", size = 2.2) +
    facet_wrap(~ quantity, scales = "free_y") +
    scale_x_log10(breaks = sort(unique(df$n))) +
    scale_linetype_manual(values = c("teoretyczna" = "dashed",
                                     "empiryczna"  = "solid")) +
    labs(title = title, subtitle = subtitle,
         x = "Liczność próby n (skala log)", y = NULL,
         linetype = "Krzywa") +
    ntp_theme()
}

#' Mean of the standardised D'Agostino components (Z1 = skew, Z2 = kurt).
#' Tells us *which* component is doing the work in the omnibus test.
plot_dagostino_components <- function(grid_df, K_focus, title, subtitle = NULL) {
  df <- grid_df[grid_df$K == K_focus, ]
  long <- data.frame(
    n         = rep(df$n, 2),
    Z         = c(df$mean_Z1, df$mean_Z2),
    component = factor(rep(c("Z1 (skośność)", "Z2 (kurtoza)"),
                           each = nrow(df)),
                       levels = c("Z1 (skośność)", "Z2 (kurtoza)"))
  )
  ggplot(long, aes(x = n, y = Z, colour = component, shape = component)) +
    geom_hline(yintercept = qnorm(0.975), linetype = "dashed",
               colour = "grey60") +
    geom_hline(yintercept = -qnorm(0.975), linetype = "dashed",
               colour = "grey60") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.8) +
    scale_x_log10(breaks = sort(unique(df$n))) +
    scale_colour_manual(values = c("Z1 (skośność)" = "#4C78A8",
                                   "Z2 (kurtoza)"  = "#F58518")) +
    scale_shape_manual(values = c("Z1 (skośność)" = 16,
                                  "Z2 (kurtoza)"  = 17)) +
    labs(title = title, subtitle = subtitle,
         x = "Liczność próby n (skala log)",
         y = "Średnia standaryzowanej składowej",
         colour = NULL, shape = NULL) +
    ntp_theme()
}

#' Density of one batch of K sample means (a single illustrative replication).
#' Useful as a "what is being tested?" picture in the report.
plot_means_density <- function(means_by_n, K, title, subtitle = NULL) {
  df <- do.call(rbind, lapply(names(means_by_n), function(key) {
    means <- means_by_n[[key]]
    n     <- as.integer(key)
    data.frame(
      n     = n,
      n_lbl = paste0("n = ", n),
      value = (means - 1) / sqrt(1 / n)  # mu = 1, sigma^2 = 1
    )
  }))
  df$n_lbl <- factor(df$n_lbl,
                     levels = paste0("n = ", sort(unique(df$n))))
  ggplot(df, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 18, fill = "#4C78A8", colour = "white",
                   alpha = 0.85) +
    stat_function(fun = dnorm, colour = "#E45756", linewidth = 0.9) +
    facet_wrap(~ n_lbl, scales = "free_y") +
    labs(title = title,
         subtitle = subtitle %||% sprintf(
           "Jedna realizacja: K = %d standaryzowanych średnich; czerwona linia: N(0,1).", K),
         x = "Standaryzowana średnia",
         y = "Gęstość") +
    ntp_theme()
}

# Tiny null-coalescing helper used above (avoids loading rlang just for `%||%`).
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Histograms of empirical p-values from one Monte-Carlo run.
#' Under H0 p-values are Uniform[0,1]; piling near zero ⇔ high power.
plot_pvalue_distribution <- function(pdump_df, alpha, title, subtitle = NULL) {
  long <- data.frame(
    test    = factor(rep(c("KS", "SW", "DA"), each = nrow(pdump_df)),
                     levels = c("KS", "SW", "DA")),
    p.value = c(pdump_df$KS, pdump_df$SW, pdump_df$DA)
  )
  ggplot(long, aes(x = p.value, fill = test)) +
    geom_histogram(bins = 30, colour = "white", alpha = 0.9) +
    geom_vline(xintercept = alpha, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ test, labeller = labeller(test = TEST_LABELS)) +
    scale_fill_manual(values = TEST_COLOURS, guide = "none") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(title = title, subtitle = subtitle,
         x = "p-value", y = "Liczność") +
    ntp_theme()
}
