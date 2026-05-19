# ------------------------------------------------------------------------------
# Global configuration and constants for the normality-test power study.
# Continuation of Project 1, Task 2: instead of *illustrating* the CLT we now
# *test* the convergence of the sample-mean distribution to N(mu, sigma^2/n).
# ------------------------------------------------------------------------------

# Reproducibility
RANDOM_SEED <- 20260519L

# Sample sizes "n" (number of variables being averaged).
# Span from severe non-normality (n=2) to a regime where CLT has kicked in.
SAMPLE_SIZES <- c(2L, 3L, 5L, 10L, 20L, 30L, 50L, 100L, 200L)

# "K" - the sample size handed to each normality test (number of sample means).
# Per the task spec K should be "kilkadziesiąt"; we sweep three representative
# values so we can comment on how K influences detection power.
K_SIZES <- c(20L, 50L, 100L)
K_DEFAULT <- 50L  # value highlighted in headline tables / power curves

# Monte-Carlo replications used to estimate test power.
N_SIMULATIONS <- 5000L

# Significance level for every test.
ALPHA <- 0.05

# Parent distribution: Exp(1), strongly right-skewed - the same family as in
# Project 1/Task 2, which lets us reuse the theoretical moments below.
EXP_RATE <- 1.0
EXP_MEAN <- 1.0   # E[X] = 1/rate
EXP_VAR  <- 1.0   # Var[X] = 1/rate^2

# Theoretical shape parameters of Exp(1):
#   skew(X)            = 2
#   excess kurt(X)     = 6
# Under iid sampling: skew(X_bar_n) = 2/sqrt(n), excess-kurt(X_bar_n) = 6/n.
EXP_SKEW   <- 2.0
EXP_KURT_E <- 6.0

# Output directories (relative to the task root)
PLOTS_DIR  <- "plots"
OUTPUT_DIR <- "output"

# Plot aesthetics (mirrors Project 1/Task 2 for visual consistency)
PLOT_WIDTH_IN  <- 9
PLOT_HEIGHT_IN <- 6
PLOT_DPI       <- 150
HIST_BINS      <- 40
