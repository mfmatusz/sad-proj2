# wykresy - Problem 2

library(ggplot2)
library(dplyr)
library(tidyr)

PLOT_DIR <- "plots"

# ---- pomocnicze ustawienia stylu ----
motyw <- theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank()
  )

# ---- Wykres 1: rozklad statystyki T pod H0 i H1 ----
# ilustruje zasade dzialania detektora
wykres_rozklad_T <- function(N, A, sigma2, alpha) {
  eta  <- prog_krytyczny(N, sigma2, alpha)
  sd_T <- sqrt(sigma2 * N)

  x_min <- min(-4 * sd_T, A * N - 4 * sd_T)
  x_max <- max( 4 * sd_T, A * N + 4 * sd_T)
  t_seq <- seq(x_min, x_max, length.out = 600)

  df <- data.frame(
    t       = rep(t_seq, 2),
    gestosc = c(dnorm(t_seq, mean = 0,     sd = sd_T),
                dnorm(t_seq, mean = A * N, sd = sd_T)),
    hipoteza = rep(c("H\u2080: brak sygna\u0142u", "H\u2081: sygna\u0142 obecny"),
                   each = length(t_seq))
  )

  # obszar odrzucenia (pod H0, prawy ogon)
  df_h0  <- df[df$hipoteza == "H\u2080: brak sygna\u0142u", ]
  df_rej <- df_h0[df_h0$t >= eta, ]

  p <- ggplot(df, aes(x = t, y = gestosc, color = hipoteza, linetype = hipoteza)) +
    geom_ribbon(data = df_rej,
                aes(x = t, ymin = 0, ymax = gestosc),
                inherit.aes = FALSE,
                fill = "#d62728", alpha = 0.25) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = eta, linetype = "dashed", color = "black", linewidth = 0.7) +
    annotate("text", x = eta, y = max(df$gestosc) * 1.05,
             label = paste0("\u03b7 = ", round(eta, 1)),
             hjust = -0.1, size = 3.2) +
    annotate("text", x = eta + 0.5 * sd_T, y = max(df$gestosc) * 0.15,
             label = paste0("P(FA) = ", alpha),
             color = "#d62728", size = 3, hjust = 0) +
    scale_color_manual(values = c("H\u2080: brak sygna\u0142u" = "#1f77b4",
                                  "H\u2081: sygna\u0142 obecny" = "#ff7f0e")) +
    scale_linetype_manual(values = c("H\u2080: brak sygna\u0142u" = "solid",
                                     "H\u2081: sygna\u0142 obecny" = "solid")) +
    labs(title = "Rozk\u0142ad statystyki korelacyjnej T",
         subtitle = paste0("N = ", N, ", A = ", A, ", \u03c3\u00b2 = ", sigma2,
                           ", \u03b1 = ", alpha),
         x = "Warto\u015b\u0107 statystyki T = \u27e8x, s\u27e9",
         y = "G\u0119sto\u015b\u0107",
         color = NULL, linetype = NULL) +
    motyw
  ggsave(file.path(PLOT_DIR, "p2_01_rozklad_T.png"), p,
         width = 7, height = 4, dpi = 150)
}

# ---- Wykres 2: moc vs A dla roznych N (sigma2 stale) ----
wykres_moc_vs_A_N <- function(siatka, sigma2, alpha) {
  df <- siatka
  df$N <- factor(df$N, levels = sort(unique(df$N)),
                 labels = paste0("N = ", sort(unique(df$N))))

  # dlugi format: emp i teor jako osobne serie
  df_long <- df %>%
    pivot_longer(cols = c(moc_emp, moc_teor),
                 names_to  = "typ",
                 values_to = "moc") %>%
    mutate(typ = ifelse(typ == "moc_emp", "empiryczna", "teoretyczna"))

  p <- ggplot(df_long, aes(x = A, y = moc, color = N, linetype = typ)) +
    geom_hline(yintercept = alpha, linetype = "dotted", color = "grey50", linewidth = 0.6) +
    geom_line(linewidth = 0.8) +
    geom_point(data = df_long[df_long$typ == "empiryczna", ],
               size = 1.6) +
    scale_linetype_manual(values = c("empiryczna" = "solid", "teoretyczna" = "dashed")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Moc detektora w funkcji amplitudy A",
         subtitle = paste0("\u03c3\u00b2 = ", sigma2, ", \u03b1 = ", alpha,
                           " (linia przerywana: teoretyczna)"),
         x = "Amplituda sygna\u0142u A",
         y = "Moc (P wykrycia)",
         color = "D\u0142ugo\u015b\u0107 kodu N",
         linetype = NULL) +
    motyw
  ggsave(file.path(PLOT_DIR, "p2_02_moc_vs_A_N.png"), p,
         width = 7.5, height = 4.5, dpi = 150)
}

# ---- Wykres 3: moc vs A dla roznych sigma2 (N stale) ----
wykres_moc_vs_A_sigma2 <- function(siatka, N, alpha) {
  df <- siatka
  df$sigma2 <- factor(df$sigma2, levels = sort(unique(df$sigma2)),
                      labels = paste0("\u03c3\u00b2 = ", sort(unique(df$sigma2))))

  df_long <- df %>%
    pivot_longer(cols = c(moc_emp, moc_teor),
                 names_to  = "typ",
                 values_to = "moc") %>%
    mutate(typ = ifelse(typ == "moc_emp", "empiryczna", "teoretyczna"))

  p <- ggplot(df_long, aes(x = A, y = moc, color = sigma2, linetype = typ)) +
    geom_hline(yintercept = alpha, linetype = "dotted", color = "grey50", linewidth = 0.6) +
    geom_line(linewidth = 0.8) +
    geom_point(data = df_long[df_long$typ == "empiryczna", ],
               size = 1.6) +
    scale_linetype_manual(values = c("empiryczna" = "solid", "teoretyczna" = "dashed")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Moc detektora w funkcji amplitudy A",
         subtitle = paste0("N = ", N, ", \u03b1 = ", alpha,
                           " (linia przerywana: teoretyczna)"),
         x = "Amplituda sygna\u0142u A",
         y = "Moc (P wykrycia)",
         color = "Wariancja szumu",
         linetype = NULL) +
    motyw
  ggsave(file.path(PLOT_DIR, "p2_03_moc_vs_A_sigma2.png"), p,
         width = 7.5, height = 4.5, dpi = 150)
}

# ---- Wykres 4: moc vs N dla malej amplitudy (slaby sygnal) ----
wykres_moc_vs_N <- function(N_vec, A_vec_small, sigma2, alpha, M) {
  siatka <- expand.grid(N = N_vec, A = A_vec_small)
  siatka$moc_teor <- mapply(moc_teoretyczna,
                             N = siatka$N, A = siatka$A,
                             MoreArgs = list(sigma2 = sigma2, alpha = alpha))
  siatka$moc_emp  <- mapply(moc_empiryczna,
                             N = siatka$N, A = siatka$A,
                             MoreArgs = list(sigma2 = sigma2, alpha = alpha, M = M))

  siatka$A <- factor(siatka$A,
                     labels = paste0("A = ", sort(unique(siatka$A))))

  df_long <- siatka %>%
    pivot_longer(cols = c(moc_emp, moc_teor),
                 names_to  = "typ",
                 values_to = "moc") %>%
    mutate(typ = ifelse(typ == "moc_emp", "empiryczna", "teoretyczna"))

  p <- ggplot(df_long, aes(x = N, y = moc, color = A, linetype = typ)) +
    geom_hline(yintercept = alpha, linetype = "dotted", color = "grey50", linewidth = 0.6) +
    geom_line(linewidth = 0.8) +
    geom_point(data = df_long[df_long$typ == "empiryczna", ], size = 1.8) +
    scale_x_log10(breaks = N_vec) +
    scale_linetype_manual(values = c("empiryczna" = "solid", "teoretyczna" = "dashed")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Moc detektora w funkcji d\u0142ugo\u015bci kodu N",
         subtitle = paste0("\u03c3\u00b2 = ", sigma2, ", \u03b1 = ", alpha,
                           "  (przypadek s\u0142abego sygna\u0142u A \u226a \u03c3)"),
         x = "D\u0142ugo\u015b\u0107 kodu N (skala log)",
         y = "Moc (P wykrycia)",
         color = "Amplituda",
         linetype = NULL) +
    motyw
  ggsave(file.path(PLOT_DIR, "p2_04_moc_vs_N.png"), p,
         width = 7.5, height = 4.5, dpi = 150)
}

# ---- Wykres 5: rozklad p-value pod H0 (weryfikacja rozmiaru testu) ----
wykres_pvalue_H0 <- function(N, sigma2, alpha, M) {
  s <- generuj_kod(N)
  eta <- prog_krytyczny(N, sigma2, alpha)
  sd_T <- sqrt(sigma2 * N)

  T_vals <- replicate(M, {
    x <- generuj_sygnal(s, A = 0, sigma2, hipoteza = "H0")
    statystyka_T(x, s)
  })
  p_vals <- 1 - pnorm(T_vals, mean = 0, sd = sd_T)

  df <- data.frame(p_value = p_vals)
  alpha_emp <- mean(p_vals < alpha)

  p <- ggplot(df, aes(x = p_value)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, fill = "#1f77b4", color = "white", alpha = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    annotate("text", x = 0.5, y = 1.2,
             label = paste0("rozmiar empiryczny = ", round(alpha_emp, 4),
                            "\n(nominalne \u03b1 = ", alpha, ")"),
             size = 3.2) +
    labs(title = "Rozk\u0142ad p-value pod H\u2080",
         subtitle = paste0("N = ", N, ", \u03c3\u00b2 = ", sigma2,
                           ", M = ", M, " replikacji"),
         x = "p-value", y = "G\u0119sto\u015b\u0107") +
    motyw
  ggsave(file.path(PLOT_DIR, "p2_05_pvalue_H0.png"), p,
         width = 6, height = 3.8, dpi = 150)
}
