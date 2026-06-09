# detektor optymalny (test ilorazu wiarygodnosci / korelacyjny)
#
# Model:
#   H0: x(i) = w(i),          w(i) ~ N(0, sigma2)
#   H1: x(i) = A*s(i) + w(i), A nieznane, ale > 0 zakladamy detekcje jednostronna
#
# Funkcja wiarygodnosci (log) dla H1 wzgledem H0:
#   LLR = sum_i [ x(i)*A*s(i)/sigma2 - (A*s(i))^2 / (2*sigma2) ]
#
# Statystyka wystarczajaca (GLRT / korelacja dopasowana):
#   T = sum_i x(i)*s(i)
#
# Pod H0: T ~ N(0, sigma2 * sum(s^2)) = N(0, sigma2 * N)
# Pod H1: T ~ N(A * sum(s^2), sigma2 * sum(s^2)) = N(A*N, sigma2*N)
#
# Prog krytyczny dla poziomu istotnosci alpha (H1: A > 0):
#   eta = qnorm(1 - alpha) * sqrt(sigma2 * N)
#
# Odrzucamy H0 gdy T > eta

# oblicza statystyke korelacyjna T = <x, s>
statystyka_T <- function(x, s) {
  sum(x * s)
}

# wyznacza prog krytyczny dla zadanego alpha, sigma2, N
prog_krytyczny <- function(N, sigma2, alpha) {
  qnorm(1 - alpha) * sqrt(sigma2 * N)
}

# podejmuje decyzje: TRUE = wykryto sygnal (odrzucono H0)
decyzja <- function(x, s, sigma2, alpha) {
  T_obs <- statystyka_T(x, s)
  eta   <- prog_krytyczny(length(s), sigma2, alpha)
  T_obs > eta
}
