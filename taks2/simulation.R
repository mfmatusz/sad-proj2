# symulacje Monte Carlo - moc detektora

source("R/config.R")
source("R/signal.R")
source("R/detector.R")

# oblicza empiryczna moc dla jednej kombinacji (N, A, sigma2)
# moc = P(odrzucenie H0 | H1 prawdziwa)
moc_empiryczna <- function(N, A, sigma2, alpha, M, ziarno_kodu = 42) {
  s <- generuj_kod(N, ziarno = ziarno_kodu)
  trafienia <- replicate(M, {
    x <- generuj_sygnal(s, A, sigma2, hipoteza = "H1")
    decyzja(x, s, sigma2, alpha)
  })
  mean(trafienia)
}

# oblicza empiryczny rozmiar testu (alpha_emp) pod H0
rozmiar_empiryczny <- function(N, sigma2, alpha, M, ziarno_kodu = 42) {
  s <- generuj_kod(N, ziarno = ziarno_kodu)
  alarmy <- replicate(M, {
    x <- generuj_sygnal(s, A = 0, sigma2, hipoteza = "H0")
    decyzja(x, s, sigma2, alpha)
  })
  mean(alarmy)
}

# oblicza teoretyczna moc (analitycznie)
# T ~ N(A*N, sigma2*N) pod H1, prog = qnorm(1-alpha)*sqrt(sigma2*N)
moc_teoretyczna <- function(N, A, sigma2, alpha) {
  eta  <- prog_krytyczny(N, sigma2, alpha)
  mean_H1 <- A * N
  sd_T    <- sqrt(sigma2 * N)
  1 - pnorm(eta, mean = mean_H1, sd = sd_T)
}

# buduje siatke wynikow: moc vs A dla roznych N (bazowe sigma2)
siatka_moc_vs_A <- function(N_vec, A_vec, sigma2, alpha, M) {
  wyniki <- expand.grid(N = N_vec, A = A_vec)
  wyniki$moc_emp  <- mapply(moc_empiryczna,
                            N = wyniki$N, A = wyniki$A,
                            MoreArgs = list(sigma2 = sigma2, alpha = alpha, M = M))
  wyniki$moc_teor <- mapply(moc_teoretyczna,
                            N = wyniki$N, A = wyniki$A,
                            MoreArgs = list(sigma2 = sigma2, alpha = alpha))
  wyniki
}

# buduje siatke wynikow: moc vs A dla roznych sigma2 (bazowe N)
siatka_moc_vs_sigma2 <- function(N, A_vec, sigma2_vec, alpha, M) {
  wyniki <- expand.grid(sigma2 = sigma2_vec, A = A_vec)
  wyniki$moc_emp  <- mapply(moc_empiryczna,
                            A = wyniki$A, sigma2 = wyniki$sigma2,
                            MoreArgs = list(N = N, alpha = alpha, M = M))
  wyniki$moc_teor <- mapply(moc_teoretyczna,
                            A = wyniki$A, sigma2 = wyniki$sigma2,
                            MoreArgs = list(N = N, alpha = alpha))
  wyniki
}
