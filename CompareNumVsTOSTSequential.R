## ======================================================================
## Single-file script: compare sequential TOST-E vs multiplied-numeraires
## Used for figures in appendix E.2 of Koobs and Koning (2026)
## Edit ONLY the cfg block.
## ======================================================================

# Compare asymmetric sequential TOST-E vs multiplied-numeraire process
# in a Gaussian i.i.d. DGP.
#
# DGP: X_i ~ N(mu_true, sigma^2)
# Null: H0: mu <= Delta_minus OR mu >= Delta_plus
# Alternative used for LR construction: mu_alt (inside margins)
#
# Outputs:
#   1) CSV summary with expected e-values and rejection probabilities by n
#   2) PNG with two panels (expected e-value, rejection probability)

# For upper pane, use mu_true = 0 and for lower pane, set mu_true = 0.3

set.seed(123)

cfg <- list(
  reps = 200000L,
  n_min = 2L,
  n_max = 75L,
  alpha = 0.05,
  sigma = 1.0,
  mu_true = 0.3,       # DGP mean (inside margins for power checks)
  mu_alt = 0,        # reference alternative in LR numerator
  Delta_minus = -0.60, # asymmetric left margin
  Delta_plus = 0.40,   # asymmetric right margin
  out_csv = "asymmetric_seq_tost_vs_numeraire.csv",
  out_png = NULL       # set to filename (e.g. "plot.png") to save to file
)

if (!(cfg$Delta_minus < cfg$mu_alt && cfg$mu_alt < cfg$Delta_plus)) {
  stop("Require Delta_minus < mu_alt < Delta_plus.")
}
if (cfg$n_min < 1L || cfg$n_max < cfg$n_min) {
  stop("Require 1 <= n_min <= n_max.")
}

# Stable log(c*exp(a) + (1-c)*exp(b))
log_mix2 <- function(log_a, log_b, c_weight) {
  if (c_weight <= 0 || c_weight >= 1) stop("c_weight must be in (0,1).")
  u <- log(c_weight) + log_a
  v <- log(1 - c_weight) + log_b
  m <- pmax(u, v)
  out <- m + log(exp(u - m) + exp(v - m))
  both_neg_inf <- is.infinite(u) & (u < 0) & is.infinite(v) & (v < 0)
  out[both_neg_inf] <- -Inf
  out
}

# Numerically stable mean(exp(x)) from log-values
mean_exp_from_logs <- function(x) {
  m <- max(x)
  exp(m) * mean(exp(x - m))
}

# One-step numeraire calibration for asymmetric boundaries.
# We solve A(c) - B(c) = 0 where:
#   A(c) = E_{Delta_minus}[ q / (c p_- + (1-c) p_+) ]
#   B(c) = E_{Delta_plus }[ q / (c p_- + (1-c) p_+) ]
# and then scale by K = max(A(c*), B(c*)) to guarantee boundary control.
calibrate_numeraire <- function(mu_alt, sigma, Delta_minus, Delta_plus) {
  # Finite bounds avoid underflow pathologies in extreme tails.
  lo <- min(mu_alt, Delta_minus, Delta_plus) - 12 * sigma
  hi <- max(mu_alt, Delta_minus, Delta_plus) + 12 * sigma

  log_p_alt <- function(x) dnorm(x, mean = mu_alt, sd = sigma, log = TRUE)
  log_p_m <- function(x) dnorm(x, mean = Delta_minus, sd = sigma, log = TRUE)
  log_p_p <- function(x) dnorm(x, mean = Delta_plus, sd = sigma, log = TRUE)

  A_integrand <- function(x, c_weight) {
    l_alt <- log_p_alt(x)
    l_m <- log_p_m(x)
    l_mix <- log_mix2(l_m, log_p_p(x), c_weight)
    out <- exp(l_alt - l_mix + l_m)
    out[!is.finite(out)] <- 0
    out
  }
  B_integrand <- function(x, c_weight) {
    l_alt <- log_p_alt(x)
    l_p <- log_p_p(x)
    l_mix <- log_mix2(log_p_m(x), l_p, c_weight)
    out <- exp(l_alt - l_mix + l_p)
    out[!is.finite(out)] <- 0
    out
  }

  A_of_c <- function(c_weight) {
    integrate(
      f = function(x) A_integrand(x, c_weight),
      lower = lo, upper = hi, rel.tol = 1e-9, stop.on.error = FALSE
    )$value
  }
  B_of_c <- function(c_weight) {
    integrate(
      f = function(x) B_integrand(x, c_weight),
      lower = lo, upper = hi, rel.tol = 1e-9, stop.on.error = FALSE
    )$value
  }
  g <- function(c_weight) A_of_c(c_weight) - B_of_c(c_weight)

  eps <- 1e-6
  g_lo <- g(eps)
  g_hi <- g(1 - eps)
  if (g_lo * g_hi > 0) {
    stop("Could not bracket root for c* in (0,1).")
  }

  c_star <- uniroot(g, c(eps, 1 - eps), tol = 1e-12)$root
  A_star <- A_of_c(c_star)
  B_star <- B_of_c(c_star)
  K_scale <- max(A_star, B_star)

  list(c_star = c_star, K_scale = K_scale, A_star = A_star, B_star = B_star)
}

cal <- calibrate_numeraire(
  mu_alt = cfg$mu_alt,
  sigma = cfg$sigma,
  Delta_minus = cfg$Delta_minus,
  Delta_plus = cfg$Delta_plus
)

cat(sprintf("Calibrated c* = %.6f\n", cal$c_star))
cat(sprintf("Boundary checks (raw): A(c*) = %.6f, B(c*) = %.6f\n", cal$A_star, cal$B_star))
cat(sprintf("Scale factor K = %.6f\n", cal$K_scale))

# Simulate data matrix: rows = repetitions, cols = time
X <- matrix(
  rnorm(cfg$reps * cfg$n_max, mean = cfg$mu_true, sd = cfg$sigma),
  nrow = cfg$reps,
  ncol = cfg$n_max
)

# One-step log LR components
log_p_alt <- dnorm(X, mean = cfg$mu_alt, sd = cfg$sigma, log = TRUE)
log_p_L <- dnorm(X, mean = cfg$Delta_minus, sd = cfg$sigma, log = TRUE)
log_p_R <- dnorm(X, mean = cfg$Delta_plus, sd = cfg$sigma, log = TRUE)

log_e_L <- log_p_alt - log_p_L
log_e_R <- log_p_alt - log_p_R
log_e_mix_raw <- log_p_alt - log_mix2(log_p_L, log_p_R, cal$c_star)
log_e_mix <- log_e_mix_raw - log(cal$K_scale)

# Cumulative log-processes
log_E_L <- t(apply(log_e_L, 1, cumsum))
log_E_R <- t(apply(log_e_R, 1, cumsum))
log_E_TOST <- pmin(log_E_L, log_E_R)
log_E_MIX <- t(apply(log_e_mix, 1, cumsum))

n_grid <- cfg$n_min:cfg$n_max
log_E_TOST <- log_E_TOST[, n_grid, drop = FALSE]
log_E_MIX <- log_E_MIX[, n_grid, drop = FALSE]

# Expected e-value curves
avg_e_tost <- apply(log_E_TOST, 2, mean_exp_from_logs)
avg_e_mix <- apply(log_E_MIX, 2, mean_exp_from_logs)

# Sequential rejection probabilities by time n: P(max_{k<=n} E_k >= 1/alpha)
log_thr <- log(1 / cfg$alpha)
hit_tost <- log_E_TOST >= log_thr
hit_mix <- log_E_MIX >= log_thr

cum_hit_tost <- t(apply(hit_tost, 1, cummax))
cum_hit_mix <- t(apply(hit_mix, 1, cummax))

rej_tost <- colMeans(cum_hit_tost)
rej_mix <- colMeans(cum_hit_mix)

summary_df <- data.frame(
  n = n_grid,
  avg_e_tost = avg_e_tost,
  avg_e_numeraire = avg_e_mix,
  rej_prob_tost = rej_tost,
  rej_prob_numeraire = rej_mix
)

write.csv(summary_df, cfg$out_csv, row.names = FALSE)
cat(sprintf("Saved summary to %s\n", cfg$out_csv))

# Plot (ggplot2).
# If out_png is NULL, this plots to the current graphics device
# (e.g., RStudio's bottom-right Plots pane).
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Please install ggplot2 to generate the plots.")
}
library(ggplot2)

cols <- c("Symmetric LR" = "#E69F00", "TOST" = "#0072B2")

plot_df_rej <- data.frame(
  x = n_grid,
  TOST = rej_tost,
  Ratio = rej_mix
)

plot_df_e <- data.frame(
  x = n_grid,
  TOST = avg_e_tost,
  Ratio = avg_e_mix
)

p_rej <- ggplot(plot_df_rej, aes(x = x)) +
  geom_line(aes(y = Ratio, color = "Symmetric LR"), linewidth = 1.2) +
  geom_line(aes(y = TOST, color = "TOST"), linewidth = 1.2) +
  scale_color_manual(values = cols) +
  scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 1, 0.25)) +  # add this
  labs(x = "Sample size", y = "Rejection Probability") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

if (cfg$mu_true == 0) {
  p_e <- ggplot(plot_df_e, aes(x = x)) +
    geom_line(aes(y = Ratio, color = "Symmetric LR"), linewidth = 1.2) +
    geom_line(aes(y = TOST, color = "TOST"), linewidth = 1.2) +
    scale_color_manual(values = cols) +
    labs(x = "Sample size", y = "E-value") +
    scale_y_continuous(limits = c(0, 1250), breaks = seq(0, 1250, 250)) + 
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
} else if (cfg$mu_true == 0.3) {
  p_e <- ggplot(plot_df_e, aes(x = x)) +
    geom_line(aes(y = Ratio, color = "Symmetric LR"), linewidth = 1.2) +
    geom_line(aes(y = TOST, color = "TOST"), linewidth = 1.2) +
    scale_color_manual(values = cols) +
    labs(x = "Sample size", y = "E-value") +
    scale_y_continuous(limits = c(0, 70), breaks = seq(0, 60, 20)) + 
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
}


draw_two_panel <- function(left_plot, right_plot) {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(left_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(right_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
}

plot_to_file <- !is.null(cfg$out_png) && nzchar(cfg$out_png)
if (plot_to_file) {
  png(cfg$out_png, width = 1800, height = 800, res = 150)
}
draw_two_panel(p_e, p_rej)

if (plot_to_file) {
  dev.off()
  cat(sprintf("Saved figure to %s\n", cfg$out_png))
} else {
  cat("Plotted figure to current graphics device.\n")
}

# Console snapshot
cat("\nHead of summary table:\n")
print(head(summary_df, 8), row.names = FALSE)
cat("\nTail of summary table:\n")
print(tail(summary_df, 8), row.names = FALSE)

