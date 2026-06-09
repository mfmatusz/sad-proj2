# generowanie sygnalu i obserwacji

# generuje deterministyczny kod pseudoszumowy z ziarna
# zwraca wektor N wartosci {-1, +1}
generuj_kod <- function(N, ziarno = 42) {
  set.seed(ziarno)
  sample(c(-1, 1), size = N, replace = TRUE)
}

# generuje probe losowa x(i) pod H0 lub H1
# H0: x = w,   H1: x = A*s + w
# w ~ N(0, sigma2)
generuj_sygnal <- function(s, A, sigma2, hipoteza = "H1") {
  N <- length(s)
  w <- rnorm(N, mean = 0, sd = sqrt(sigma2))
  if (hipoteza == "H0") {
    return(w)
  } else {
    return(A * s + w)
  }
}
