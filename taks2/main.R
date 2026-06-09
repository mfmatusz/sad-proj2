# main.R - uruchamia caly eksperyment dla Problemu 2
# wywolanie: Rscript R/main.R  (z katalogu glownego projektu)

source("R/config.R")
source("R/signal.R")
source("R/detector.R")
source("R/simulation.R")
source("R/plots.R")

cat("Problem 2: detekcja sygnalu pseudoszumowego\n")
cat("============================================\n")

# --- 1. weryfikacja rozmiaru testu ---
cat("Sprawdzam rozmiar testu pod H0...\n")
alpha_emp <- rozmiar_empiryczny(N = N_BASE, sigma2 = SIGMA2_BASE,
                                alpha = ALPHA, M = M)
cat(sprintf("  Rozmiar empiryczny: %.4f  (nominalne alfa = %.4f)\n",
            alpha_emp, ALPHA))

# --- 2. rozklad statystyki T (ilustracja zasady dzialania) ---
cat("Generuje wykres rozkladu T...\n")
wykres_rozklad_T(N = N_BASE, A = A_BASE, sigma2 = SIGMA2_BASE, alpha = ALPHA)

# --- 3. moc vs A: rozne N, stale sigma2 ---
cat("Obliczam moc vs A dla roznych N (moze chwile zajac)...\n")
siatka_N <- siatka_moc_vs_A(N_vec   = N_VEC,
                             A_vec   = A_VEC,
                             sigma2  = SIGMA2_BASE,
                             alpha   = ALPHA,
                             M       = M)
write.csv(siatka_N, "output/p2_moc_vs_A_N.csv", row.names = FALSE)
wykres_moc_vs_A_N(siatka_N, sigma2 = SIGMA2_BASE, alpha = ALPHA)

# --- 4. moc vs A: rozne sigma2, stale N ---
cat("Obliczam moc vs A dla roznych sigma2...\n")
siatka_s2 <- siatka_moc_vs_sigma2(N       = N_BASE,
                                   A_vec   = A_VEC,
                                   sigma2_vec = SIGMA2_VEC,
                                   alpha   = ALPHA,
                                   M       = M)
write.csv(siatka_s2, "output/p2_moc_vs_A_sigma2.csv", row.names = FALSE)
wykres_moc_vs_A_sigma2(siatka_s2, N = N_BASE, alpha = ALPHA)

# --- 5. moc vs N: slaby sygnal ---
cat("Obliczam moc vs N dla slabego sygnalu...\n")
A_SMALL <- c(0.25, 0.5, 1.0)
wykres_moc_vs_N(N_vec = c(10, 50, 200, 500, 1000),
                A_vec_small = A_SMALL,
                sigma2 = SIGMA2_BASE,
                alpha  = ALPHA,
                M      = M)

# --- 6. rozklad p-value pod H0 ---
cat("Generuje wykres p-value pod H0...\n")
wykres_pvalue_H0(N = N_BASE, sigma2 = SIGMA2_BASE, alpha = ALPHA, M = M)

cat("\nGotowe! Wykresy zapisane w: plots/\n")
cat("Dane zapisane w: output/\n")
