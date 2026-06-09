# parametry eksperymentu - Problem 2: detekcja sygnalu pseudoszumowego

set.seed(20260519)

# dlugosc kodu pseudoszumowego
N_VEC <- c(10, 50, 200, 1000)

# amplitudy sygnal/szum do analizy mocy
A_VEC <- seq(0, 3, by = 0.25)

# wariancja szumu gaussowskiego
SIGMA2_VEC <- c(0.5, 1, 4)

# poziom istotnosci = prawdopodobienstwo falszywego alarmu
ALPHA <- 0.001

# liczba replikacji Monte Carlo
M <- 10000

# bazowe parametry do ilustracji
N_BASE    <- 50
SIGMA2_BASE <- 1
A_BASE    <- 0.5
