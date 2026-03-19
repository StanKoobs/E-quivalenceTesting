## ======================================================================
## Single-file script: compare z-test (NP and e) and t-test (TOST and e)
## Used for Figure 2 of Koobs and Koning (2026)
## Edit ONLY the cfg block.
## ======================================================================

# Symmetric equivalence curves:
# - z-test e-curve at h = 0.5
# - z-test binary test (NP)
# - t-test TOST-E (noncentral-t based, shifted alternative on mean scale)
# - t-test binary TOST (central-t based)

cfg <- list(
  # Fixed realized statistic Xbar
  xbar = 0.10,

  # Data summary used for plotting fixed-realization curves
  n = 10,
  sigma = 1,
  s_obs = 1.0,  # observed sample SD used in t-statistics

  # Delta and alpha grids
  Delta_min = 0.000,
  Delta_max = 2.000,
  Delta_step = 0.001,
  alpha_min = 1e-4,
  alpha_max = 0.500,
  alpha_step = 0.0001,
  alpha_plot_min_left = 2e-3,

  # Curves
  h_vec = c(0.5),
  include_tost_h05 = TRUE,
  # TOST-E mode choices (coarse-tested):
  # - shifted          : shift = tost_shift / (1 - h)
  # - shifted_const    : shift = tost_shift
  # - shifted_prop     : shift = tost_prop * Delta
  # - shifted_prop_cap : shift = min(tost_prop * Delta, tost_cap)
  # - shifted_hybrid   : shift = tost_shift + tost_prop * Delta
  tost_alt_modes = c("shifted_prop"),
  tost_h = 0.5,
  tost_shift = 0.10,
  tost_prop = 0.50,
  tost_cap = 0.20,
  include_h1 = TRUE,
  include_tost_binary = TRUE,
  alpha0 = 0.05,

  # Quadrature for normalization
  quad_m = 4001,

  # Plot limits
  y_left_max = 1.70,
  y_right_max = 80,
  Delta_plot_max = 1.20,
  show_infinite_at_top = FALSE,
  log_x_left = TRUE,
  log_x_right = FALSE,

  # Draw in active device (RStudio pane) when interactive
  show_in_device = interactive(),

  # Save files
  save_png = "Figs/Section7_tmean_toste_n10_xbar010_wideDelta.png",
  save_pdf = "Figs/Section7_tmean_toste_n10_xbar010_wideDelta.pdf"
)

stopifnot(cfg$n >= 2, cfg$sigma > 0, cfg$alpha0 > 0, cfg$alpha0 < 1)
stopifnot(is.finite(cfg$s_obs), cfg$s_obs > 0)
stopifnot(cfg$Delta_max > cfg$Delta_min, cfg$Delta_step > 0)
stopifnot(cfg$Delta_plot_max > cfg$Delta_min, cfg$Delta_plot_max <= cfg$Delta_max)
stopifnot(cfg$alpha_max > cfg$alpha_min, cfg$alpha_step > 0, cfg$alpha_min > 0)
stopifnot(cfg$alpha_plot_min_left >= cfg$alpha_min, cfg$alpha_plot_min_left < cfg$alpha_max)
stopifnot(all(cfg$h_vec >= 0), all(cfg$h_vec < 1))
stopifnot(cfg$quad_m >= 201)
stopifnot(all(cfg$tost_alt_modes %in% c("shifted", "shifted_const", "shifted_prop", "shifted_prop_cap", "shifted_hybrid")))
stopifnot(is.finite(cfg$tost_shift), cfg$tost_shift >= 0)
stopifnot(is.finite(cfg$tost_prop), cfg$tost_prop >= 0)
stopifnot(is.finite(cfg$tost_cap), cfg$tost_cap >= 0)
stopifnot(is.finite(cfg$tost_h), cfg$tost_h > 0, cfg$tost_h < 1)

sd_x <- cfg$sigma / sqrt(cfg$n)
nu <- cfg$n - 1
Delta_grid <- seq(cfg$Delta_min, cfg$Delta_max, by = cfg$Delta_step)
alpha_grid <- seq(cfg$alpha_min, cfg$alpha_max, by = cfg$alpha_step)

# Densities for Xbar under mu
d0 <- function(x) dnorm(x, mean = 0, sd = sd_x)
dL <- function(x, Delta) dnorm(x, mean = -Delta, sd = sd_x)
dR <- function(x, Delta) dnorm(x, mean = Delta, sd = sd_x)

# h=0 direct symmetric LR
e0 <- function(x, Delta) d0(x) / (0.5 * dL(x, Delta) + 0.5 * dR(x, Delta))

# Quadrature nodes for E_{mu=Delta}
u_quad <- (seq_len(cfg$quad_m) - 0.5) / cfg$quad_m
z_quad <- qnorm(u_quad)

norm_const <- function(Delta, h) {
  k <- 1 / (1 - h)
  x_nodes <- Delta + sd_x * z_quad
  mean(e0(x_nodes, Delta)^k)
}

make_h_curve <- function(h) {
  k <- 1 / (1 - h)
  z_vec <- vapply(Delta_grid, function(D) norm_const(D, h), numeric(1))
  (e0(cfg$xbar, Delta_grid)^k) / z_vec
}

# t-statistics for testing equivalence on the mean scale with unknown sigma
t_left_obs <- function(Delta) sqrt(cfg$n) * (cfg$xbar + Delta) / cfg$s_obs
t_right_obs <- function(Delta) sqrt(cfg$n) * (Delta - cfg$xbar) / cfg$s_obs

safe_dt <- function(x, df, ncp = 0) pmax(dt(x, df = df, ncp = ncp), .Machine$double.xmin)

tost_lambda <- function(Delta, h, mode) {
  shift <- switch(
    mode,
    "shifted" = cfg$tost_shift / (1 - h),
    "shifted_const" = cfg$tost_shift,
    "shifted_prop" = cfg$tost_prop * Delta,
    "shifted_prop_cap" = pmin(cfg$tost_prop * Delta, cfg$tost_cap),
    "shifted_hybrid" = cfg$tost_shift + cfg$tost_prop * Delta,
    stop("Unsupported TOST alternative mode")
  )
  sqrt(cfg$n) * shift / cfg$s_obs
}

q_left_t <- function(t, Delta, h, mode) {
  safe_dt(t, df = nu, ncp = tost_lambda(Delta, h, mode))
}

q_right_t <- function(t, Delta, h, mode) {
  safe_dt(t, df = nu, ncp = tost_lambda(Delta, h, mode))
}

# TOST-E from one-sided t-based e-values:
# e_L = q_L(T_L)/f0(T_L), e_R = q_R(T_R)/f0(T_R), then min(e_L, e_R),
# with h-power normalization under the central t boundary model.
make_tost_curve <- function(h, mode) {
  k <- 1 / (1 - h)
  t_nodes <- qt(u_quad, df = nu)
  f0_nodes <- safe_dt(t_nodes, df = nu, ncp = 0)

  e0_left <- vapply(Delta_grid, function(D) {
    tl <- t_left_obs(D)
    q_left_t(tl, D, h, mode) / safe_dt(tl, df = nu, ncp = 0)
  }, numeric(1))

  e0_right <- vapply(Delta_grid, function(D) {
    tr <- t_right_obs(D)
    q_right_t(tr, D, h, mode) / safe_dt(tr, df = nu, ncp = 0)
  }, numeric(1))

  z_left <- vapply(Delta_grid, function(D) {
    mean((q_left_t(t_nodes, D, h, mode) / f0_nodes)^k)
  }, numeric(1))

  z_right <- vapply(Delta_grid, function(D) {
    mean((q_right_t(t_nodes, D, h, mode) / f0_nodes)^k)
  }, numeric(1))

  e_left <- (e0_left^k) / z_left
  e_right <- (e0_right^k) / z_right
  pmin(e_left, e_right)
}

# Binary TOST from one-sided t-tests at level alpha0,
# encoded as e-value in {0, 1/alpha0}.
make_tost_binary_curve <- function(alpha0) {
  p_left <- 1 - pt(t_left_obs(Delta_grid), df = nu)
  p_right <- 1 - pt(t_right_obs(Delta_grid), df = nu)
  reject <- pmax(p_left, p_right) <= alpha0
  ifelse(reject, 1 / alpha0, 0)
}

# Build e-curves
e_curves <- list()
for (h in cfg$h_vec) {
  nm <- "z_h05"
  cat(sprintf("Computing direct curve for h = %.3f ...\n", h))
  e_curves[[nm]] <- make_h_curve(h)
}
if (isTRUE(cfg$include_tost_h05)) {
  for (mode in cfg$tost_alt_modes) {
    nm <- paste0("tost_e_", mode)
    cat(sprintf("Computing TOST-E curve for h = %.3f (mode=%s) ...\n", cfg$tost_h, mode))
    e_curves[[nm]] <- make_tost_curve(cfg$tost_h, mode)
  }
}

# NP (h=1) encoded as {0,1/alpha0}
c_alpha <- function(Delta, alpha) {
  f <- function(c) {
    pnorm(c, mean = Delta, sd = sd_x) - pnorm(-c, mean = Delta, sd = sd_x) - alpha
  }
  uniroot(f, interval = c(0, Delta + 12 * sd_x), tol = 1e-10)$root
}
if (isTRUE(cfg$include_h1)) {
  c_vec <- vapply(Delta_grid, function(D) c_alpha(D, cfg$alpha0), numeric(1))
  e_curves[["z_bin"]] <- ifelse(abs(cfg$xbar) < c_vec, 1 / cfg$alpha0, 0)
}
if (isTRUE(cfg$include_tost_binary)) {
  cat(sprintf("Computing binary TOST curve (alpha=%.3f) ...\n", cfg$alpha0))
  e_curves[["tost_bin"]] <- make_tost_binary_curve(cfg$alpha0)
}

# Inversion alpha -> Delta_hat(alpha)
delta_hat_from_e <- function(e_vec) {
  raw <- vapply(alpha_grid, function(a) {
    idx <- which(e_vec >= 1 / a)[1]
    if (is.na(idx)) Inf else Delta_grid[idx]
  }, numeric(1))
  # Numerical inversion can introduce tiny non-monotone wiggles; enforce
  # the theoretical non-increasing shape in alpha.
  out <- raw
  running <- Inf
  for (i in seq_along(out)) {
    if (is.finite(out[i])) {
      running <- min(running, out[i])
      out[i] <- running
    }
  }
  out
}
dhat_curves <- lapply(e_curves, delta_hat_from_e)

clip_left <- function(x, ymax, show_inf_top = FALSE) {
  out <- rep(NA_real_, length(x))
  inside <- is.finite(x) & (x <= ymax)
  out[inside] <- x[inside]
  if (isTRUE(show_inf_top)) {
    out[!is.finite(x)] <- ymax
  }
  out
}

show_inf_for_curve <- function(curve_name) {
  if (curve_name %in% c("z_bin", "tost_bin")) {
    FALSE
  } else {
    isTRUE(cfg$show_infinite_at_top)
  }
}

h_levels <- c(
  "z_h05",
  if (isTRUE(cfg$include_h1)) "z_bin",
  if (isTRUE(cfg$include_tost_h05)) paste0("tost_e_", cfg$tost_alt_modes),
  if (isTRUE(cfg$include_tost_binary)) "tost_bin"
)

base_cols <- c(
  "z_h05" = "#E69F00",   # z-test class (orange)
  "z_bin" = "#D55E00",   # z-test class (vermillion)
  "tost_bin" = "#009E73" # TOST class (bluish green)
)
cols <- setNames(rep(NA_character_, length(h_levels)), h_levels)
for (nm in h_levels) {
  if (startsWith(nm, "tost_e_")) {
    cols[nm] <- "#56B4E9"
  } else {
    cols[nm] <- base_cols[nm]
  }
}

ltys <- setNames(rep("solid", length(h_levels)), h_levels)
# Match linetype by class:
# z-class = solid, TOST-class = custom dashed
if ("z_h05" %in% h_levels) ltys["z_h05"] <- "solid"
if ("z_bin" %in% h_levels) ltys["z_bin"] <- "solid"
for (nm in h_levels) {
  if (startsWith(nm, "tost_e_")) ltys[nm] <- "dashed"
}
if ("tost_bin" %in% h_levels) ltys["tost_bin"] <- "dashed"

lwds <- setNames(rep(2.0, length(h_levels)), h_levels)
for (nm in h_levels) {
  if (startsWith(nm, "tost_e_")) lwds[nm] <- 1.75
}
if ("tost_bin" %in% h_levels) lwds["tost_bin"] <- 1.75

labels <- setNames(rep("", length(h_levels)), h_levels)
for (nm in h_levels) {
  if (nm == "z_h05") labels[nm] <- "z-test (e)"
  if (nm == "z_bin") labels[nm] <- "z-test (NP)"
  if (nm == "tost_bin") labels[nm] <- "t-test (TOST)"
  if (startsWith(nm, "tost_e_")) labels[nm] <- "t-test (e)"
}

draw_panels <- function() {
  draw_manual_dashed <- function(x, y, col, lwd, dash_len = 0.07, gap_len = 0.05) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 2) return(invisible(NULL))
    rr <- rle(ok)
    ends <- cumsum(rr$lengths)
    starts <- ends - rr$lengths + 1
    for (k in seq_along(rr$values)) {
      if (!rr$values[k]) next
      i1 <- starts[k]; i2 <- ends[k]
      xs <- x[i1:i2]; ys <- y[i1:i2]
      if (length(xs) < 2) next
      seglen <- sqrt(diff(xs)^2 + diff(ys)^2)
      s <- c(0, cumsum(seglen))
      total <- tail(s, 1)
      if (!is.finite(total) || total <= 0) next
      period <- dash_len + gap_len
      dash_starts <- seq(0, total, by = period)
      for (ds in dash_starts) {
        de <- min(ds + dash_len, total)
        if (de <= ds) next
        mids <- s[s > ds & s < de]
        sq <- c(ds, mids, de)
        xq <- approx(s, xs, xout = sq, ties = "ordered")$y
        yq <- approx(s, ys, xout = sq, ties = "ordered")$y
        lines(xq, yq, col = col, lwd = lwd)
      }
    }
    invisible(NULL)
  }

  # Tight layout: closer x-axis titles and smaller gap between the two panels.
  layout(matrix(c(1, 2), nrow = 1, byrow = TRUE), widths = c(1, 1))
  par(mgp = c(1.15, 0.35, 0), tcl = -0.20, xaxs = "i", yaxs = "i", lend = 1)
  par(cex.axis = 1.18, cex.lab = 1.30)

  # Left: alpha -> Delta_hat_alpha
  par(mar = c(2.9, 3.8, 0.9, 0.25))
  first <- h_levels[1]
  alpha_plot <- pmax(alpha_grid, cfg$alpha_min)
  left_xmin <- cfg$alpha_plot_min_left
  x_left <- if (isTRUE(cfg$log_x_left)) log10(alpha_plot) else alpha_plot
  xlim_left <- if (isTRUE(cfg$log_x_left)) c(log10(left_xmin), log10(cfg$alpha_max)) else c(left_xmin, cfg$alpha_max)
  plot(NA, xlab = expression(alpha), ylab = "", xaxt = "n",
       xlim = xlim_left, ylim = c(0, cfg$y_left_max), type = "n")
  if (isTRUE(cfg$log_x_left)) {
    left_ticks <- c(0.005, 0.02, 0.05, 0.1, 0.2, 0.5)
    axis(
      1,
      at = log10(left_ticks),
      labels = c("0.005", "0.02", "0.05", "0.1", "0.2", "0.5"),
      gap.axis = -1
    )
  } else {
    axis(
      1,
      at = c(seq(0, 0.4, by = 0.1), 0.49),
      labels = c(sprintf("%.1f", seq(0, 0.4, by = 0.1)), "")
    )
  }
  mtext(expression(Delta), side = 2, line = 2.1)
  for (hs in h_levels) {
    yv <- clip_left(dhat_curves[[hs]], cfg$y_left_max, show_inf_for_curve(hs))
    if (ltys[hs] == "solid") {
      lines(x_left, yv, lwd = lwds[hs], col = cols[hs], lty = "solid")
    } else {
      draw_manual_dashed(x_left, yv, col = cols[hs], lwd = lwds[hs])
    }
  }
  # For binary curves, show the jump to infinity as an in-frame vertical rise
  # up to the plot boundary (without arrows outside the frame).
  for (hs in intersect(c("z_bin", "tost_bin"), h_levels)) {
    yv <- dhat_curves[[hs]]
    j <- which(is.finite(yv))[1]
    if (!is.na(j) && j > 1) {
      xj <- if (isTRUE(cfg$log_x_left)) log10(alpha_plot[j]) else alpha_plot[j]
      yj <- min(yv[j], cfg$y_left_max)
      if (ltys[hs] == "solid") {
        segments(xj, yj, xj, cfg$y_left_max, col = cols[hs], lty = "solid", lwd = lwds[hs])
      } else {
        draw_manual_dashed(c(xj, xj), c(yj, cfg$y_left_max), col = cols[hs], lwd = lwds[hs],
                           dash_len = 0.06, gap_len = 0.045)
      }
    }
  }
  # Right: Delta -> e_Delta
  par(mar = c(2.9, 3.2, 0.9, 0.55))
  Delta_plot <- Delta_grid
  make_sym_curve <- function(e_vec) {
    x_sym <- c(-rev(Delta_plot), Delta_plot)
    y_sym <- c(rev(e_vec), e_vec)
    list(x = x_sym, y = y_sym)
  }
  sym_first <- make_sym_curve(e_curves[[first]])
  plot(
    sym_first$x, sym_first$y,
    type = "l", lwd = lwds[first], col = cols[first], lty = ltys[first],
    xlab = expression(Delta), ylab = "", xaxt = "n",
    xlim = c(-cfg$Delta_plot_max, cfg$Delta_plot_max), ylim = c(0, cfg$y_right_max),
    log = if (isTRUE(cfg$log_x_right)) "x" else ""
  )
  if (isTRUE(cfg$log_x_right)) {
    axis(1, at = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5),
         labels = c("0.001", "0.005", "0.01", "0.02", "0.05", "0.1", "0.2", "0.5"))
  } else {
    axis(1, at = seq(-1.0, 1.0, by = 0.2))
  }
  mtext(expression(epsilon[Delta]), side = 2, line = 1.8, cex = 1.45)
  if (length(h_levels) > 1) {
    for (hs in h_levels[-1]) {
      sym_hs <- make_sym_curve(e_curves[[hs]])
      lines(sym_hs$x, sym_hs$y, lwd = lwds[hs], col = cols[hs], lty = ltys[hs])
    }
  }

  usr <- par("usr")
  legend(
    x = usr[1] + 0.24 * (usr[2] - usr[1]),
    y = usr[4],
    legend = unname(labels),
    col = cols,
    lty = ltys,
    lwd = lwds,
    cex = 1.24,
    bty = "n",
    xjust = 0,
    yjust = 1
  )
  layout(1)
}

if (isTRUE(cfg$show_in_device)) {
  draw_panels()
}

if (!is.null(cfg$save_png)) {
  dir.create(dirname(cfg$save_png), recursive = TRUE, showWarnings = FALSE)
  png(cfg$save_png, width = 1550, height = 620, res = 180)
  draw_panels()
  dev.off()
}

if (!is.null(cfg$save_pdf)) {
  dir.create(dirname(cfg$save_pdf), recursive = TRUE, showWarnings = FALSE)
  pdf(cfg$save_pdf, width = 8.4, height = 3.8)
  draw_panels()
  dev.off()
}

cat("Saved plot files:\n")
cat(" -", cfg$save_png, "\n")
cat(" -", cfg$save_pdf, "\n")
cat(sprintf("DGP: xbar=%.3f, n=%d, sigma=%.3f, s_obs=%.3f (sd_x=%.3f)\n", cfg$xbar, cfg$n, cfg$sigma, cfg$s_obs, sd_x))
cat(sprintf("h values: %s\n", paste(cfg$h_vec, collapse = ", ")))
cat(sprintf("TOST-E h=0.5 included: %s\n", ifelse(isTRUE(cfg$include_tost_h05), "yes", "no")))
cat(sprintf("TOST-E alternative modes: %s\n", paste(cfg$tost_alt_modes, collapse = ", ")))
cat(sprintf("TOST-E h parameter: %.3f\n", cfg$tost_h))
cat(sprintf("TOST shift parameter: %.3f\n", cfg$tost_shift))
cat(sprintf("TOST proportional shift parameter: %.3f\n", cfg$tost_prop))
cat(sprintf("TOST shift cap parameter: %.3f\n", cfg$tost_cap))
if ("z_bin" %in% names(e_curves)) {
  cat(sprintf("delta_star(z-binary): %.6f\n", Delta_grid[which(e_curves[["z_bin"]] > 0)[1]]))
}
