## ======================================================================
## Single-file script: B/W plot of equivalence curve and e-value curve
## Used for figure in introduction of Koobs and Koning (2026)
## Edit ONLY the cfg block.
## ======================================================================

suppressPackageStartupMessages(library(ggplot2))

# For the left pane in Figure 1, use the following settings below:
# plot_php = TRUE, swap_axes = TRUE,value_limits  = c(0, 0.5),    

# For the right pane in Figure 1, use the following settings below:
# plot_php = FALSE, swap_axes = FALSE,value_limits  = c(0, 50),  

# ---- 0) Easy parameters to tune ----
cfg <- list(
  # Model / experiment
  x      = .2,
  sigma  = .1,
  delta  = -.1,
  alpha  = 0.05,
  
  # Grid over mu
  mu_min  = 0,
  mu_max  = .7,
  mu_step = 0.001,
  
  # h curves (must be < 1; optionally add the h=1 limit curve)
  h_vec            = c(0.5),
  include_h1_limit = TRUE,
  
  # Capping (defaults to 1/alpha if NULL)
  caps   = c(1000),  # NULL, length 1, or length(h_vec)
  h1_cap = NULL,  # NULL -> 1/alpha
  
  # What to plot
  panel         = "capped",   # "root", "capped", "both"
  plot_php      = TRUE,     # TRUE: plot 1/e_mu ; FALSE: plot e_mu
  swap_axes     = TRUE,     # TRUE: x=value, y=mu ; FALSE: x=mu, y=value
  use_log_scale = FALSE,    # log10 on value axis
  value_limits  = c(0, 0.5),     # limits on value axis (x if swap_axes else y)
  log_floor     = 1e-6,     # display clamp (php/log)
  
  # Legend
  show_legend     = FALSE,   # TRUE / FALSE / "auto" (auto hides if only 1 h)
  legend_position = "bottom", # only used if show_legend resolves to TRUE
  legend_nrow     = 1,
  
  # Styling (paper-like B/W)
  base_size    = 22,
  base_family  = "serif",
  line_width   = 0.7,
  facet_labels = c(root = "Uncapped", capped = "Capped"),
  linetype_map = NULL,      # optional named vector, e.g. c("0"="solid","0.5"="dashed","1"="dotted")
  
  # Saving (optional)
  save_path = NULL,         # e.g. "fig_root.pdf"
  width     = 3.35,
  height    = 2.4,
  dpi       = 300
)

# ---- 1) Numerically-stable helpers ----
LOG_MAX <- log(.Machine$double.xmax)
LOG_MIN <- log(.Machine$double.xmin)

clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

logspace_add <- function(a, b) {
  # returns log(exp(a)+exp(b)) stably, allowing -Inf
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

# ---- 2) Capped LR boost (stable: root finding in log space) ----
# Returns log(b_C) rather than b_C to avoid overflow.
compute_log_boost <- function(mu, sigma, k, C) {
  stopifnot(length(mu) == 1, is.finite(mu), mu >= 0)
  stopifnot(length(sigma) == 1, is.finite(sigma), sigma > 0)
  stopifnot(length(k) == 1, is.finite(k), k > 0)
  stopifnot(length(C) == 1, is.finite(C), C > 0)
  
  delta <- (k * mu) / sigma
  logC  <- log(C)
  
  # Solve lhs(z) = 1/C  <=>  log(lhs(z)) + log(C) = 0
  g <- function(z) {
    log_t1 <- (-delta * (z - delta / 2)) + pnorm(z - delta, log.p = TRUE)
    log_t2 <- pnorm(z, lower.tail = FALSE, log.p = TRUE)  # log(1 - Phi(z))
    log_lhs <- logspace_add(log_t1, log_t2)
    log_lhs + logC
  }
  
  # Bracket root
  lower <- -10
  upper <-  10
  it <- 0
  while (g(lower) * g(upper) > 0 && it < 30) {
    lower <- lower * 1.5
    upper <- upper * 1.5
    it <- it + 1
  }
  if (g(lower) * g(upper) > 0) stop("compute_log_boost: failed to bracket root (try different parameters).")
  
  z_sol <- uniroot(g, lower = lower, upper = upper)$root
  
  # log(b_C) = log(C) - delta*(z - delta/2)
  log_b <- logC - delta * (z_sol - delta / 2)
  log_b
}

# ---- 3) Root family in log space ----
log_e_root_mu_h <- function(x, mu_vec, h, delta, sigma) {
  stopifnot(h < 1)
  dnorm(x, mean = mu_vec + delta / (1 - h), sd = sigma, log = TRUE) -
    dnorm(x, mean = mu_vec, sd = sigma, log = TRUE)
}

# ---- 4) Data builder (short, stable) ----
make_data <- function(cfg) {
  stopifnot(is.finite(cfg$x), length(cfg$x) == 1)
  stopifnot(is.finite(cfg$sigma), cfg$sigma > 0)
  stopifnot(is.finite(cfg$alpha), cfg$alpha > 0, cfg$alpha < 1)
  stopifnot(is.finite(cfg$delta), cfg$delta != 0)
  stopifnot(is.numeric(cfg$h_vec), length(cfg$h_vec) >= 1, all(is.finite(cfg$h_vec)), all(cfg$h_vec < 1))
  stopifnot(cfg$mu_min <= cfg$mu_max, cfg$mu_step > 0)
  
  mu_grid <- seq(cfg$mu_min, cfg$mu_max, by = cfg$mu_step)
  if (length(mu_grid) == 0) stop("Empty mu grid. Check mu_min/mu_max/mu_step.")
  
  default_cap <- 1 / cfg$alpha
  caps <- cfg$caps
  if (is.null(caps)) {
    caps <- rep(default_cap, length(cfg$h_vec))
  } else if (length(caps) == 1) {
    caps <- rep(as.numeric(caps), length(cfg$h_vec))
  } else {
    stopifnot(length(caps) == length(cfg$h_vec))
    caps <- as.numeric(caps)
  }
  h1_cap <- if (is.null(cfg$h1_cap)) default_cap else as.numeric(cfg$h1_cap)
  
  # Precompute log-boost per h (mu for boost is abs(delta), matching your earlier setup)
  log_boosts <- mapply(
    function(h, cap) compute_log_boost(mu = abs(cfg$delta), sigma = cfg$sigma, k = 1 / (1 - h), C = cap),
    h = cfg$h_vec, cap = caps
  )
  stopifnot(all(is.finite(log_boosts)))
  
  # Build root + capped
  root_list <- vector("list", length(cfg$h_vec))
  cap_list  <- vector("list", length(cfg$h_vec))
  
  for (i in seq_along(cfg$h_vec)) {
    h   <- cfg$h_vec[i]
    cap <- caps[i]
    logb <- log_boosts[i]
    
    log_e_root <- log_e_root_mu_h(cfg$x, mu_grid, h, cfg$delta, cfg$sigma)
    
    # Root (finite via clamping)
    e_root <- exp(clamp(log_e_root, LOG_MIN, LOG_MAX))
    root_list[[i]] <- data.frame(mu = mu_grid, h = h, family = "root", e_value = e_root)
    
    # Capped: log(e) = min(logb + log_e_root, log(cap)), then exp safely
    log_e_cap <- pmin(logb + log_e_root, log(cap))
    e_cap <- exp(clamp(log_e_cap, LOG_MIN, LOG_MAX))
    cap_list[[i]] <- data.frame(mu = mu_grid, h = h, family = "capped", e_value = e_cap)
  }
  
  out <- rbind(do.call(rbind, root_list), do.call(rbind, cap_list))
  
  # Optional h=1 limit curve (step function)
  if (isTRUE(cfg$include_h1_limit)) {
    indicator <- if (cfg$delta > 0) {
      as.numeric(cfg$x > qnorm(1 - cfg$alpha, mean = mu_grid, sd = cfg$sigma))
    } else {
      as.numeric(cfg$x < qnorm(cfg$alpha, mean = mu_grid, sd = cfg$sigma))
    }
    out <- rbind(out, data.frame(mu = mu_grid, h = 1, family = "capped", e_value = h1_cap * indicator))
  }
  
  stopifnot(all(is.finite(out$e_value)))
  stopifnot(all(out$e_value >= 0))
  out
}

# ---- 5) Plot builder (paper B/W) ----
paper_theme <- function(cfg) {
  theme_classic(base_size = cfg$base_size, base_family = cfg$base_family) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
      axis.ticks.length = grid::unit(2.2, "pt"),
      axis.title = element_text(margin = margin(t = 6, r = 6)),
      legend.position = cfg$legend_position,
      legend.title = element_text(face = "plain"),
      legend.key.width = grid::unit(1.6, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "plain")
    )
}

bw_linetypes <- function(n) {
  types <- c("solid", "dashed", "dotdash", "dotted", "twodash", "longdash")
  rep(types, length.out = n)
}

make_plot <- function(df, cfg) {
  # panel selection
  if (cfg$panel != "both") df <- df[df$family == cfg$panel, , drop = FALSE]
  
  df$family <- factor(df$family, levels = c("root", "capped"))
  
  # Stable ordering for h (numeric sort if possible)
  h_chr <- as.character(df$h)
  h_num <- suppressWarnings(as.numeric(unique(h_chr)))
  if (!any(is.na(h_num))) {
    df$h <- factor(h_chr, levels = as.character(sort(unique(as.numeric(h_chr)))))
  } else {
    df$h <- factor(h_chr, levels = sort(unique(h_chr)))
  }
  
  show_leg <- cfg$show_legend
  if (identical(show_leg, "auto")) show_leg <- (nlevels(df$h) > 1)
  
  # Value transform for plotting (display clamp)
  e_plot <- pmax(df$e_value, cfg$log_floor)
  value  <- if (isTRUE(cfg$plot_php)) 1 / e_plot else if (isTRUE(cfg$use_log_scale)) e_plot else df$e_value
  df$x_plot <- if (isTRUE(cfg$swap_axes)) value else df$mu
  df$y_plot <- if (isTRUE(cfg$swap_axes)) df$mu else value
  
  value_lab <- if (isTRUE(cfg$plot_php)) expression(alpha) else expression(epsilon[Delta])
  x_lab <- if (isTRUE(cfg$swap_axes)) value_lab else expression(Delta)
  y_lab <- if (isTRUE(cfg$swap_axes)) expression(Delta) else value_lab
  
  hs <- levels(df$h)
  if (!is.null(cfg$linetype_map)) {
    lt <- cfg$linetype_map
    miss <- setdiff(hs, names(lt))
    if (length(miss) > 0) {
      fill <- bw_linetypes(length(miss)); names(fill) <- miss
      lt <- c(lt, fill)
    }
    lt <- lt[hs]
  } else {
    lt <- setNames(bw_linetypes(length(hs)), hs)
  }
  
  p <- ggplot(df, aes(x = x_plot, y = y_plot, linetype = h, group = h)) +
    geom_path(color = "black", linewidth = cfg$line_width, lineend = "round") +
    scale_linetype_manual(values = lt) +
    labs(x = x_lab, y = y_lab, linetype = expression(h), title = NULL, subtitle = NULL, caption = NULL) +
    guides(linetype = if (isTRUE(show_leg)) guide_legend(nrow = cfg$legend_nrow, byrow = TRUE) else "none") +
    paper_theme(cfg)
  
  if (!isTRUE(show_leg)) p <- p + theme(legend.position = "none")
  
  if (cfg$panel == "both") {
    p <- p + facet_grid(family ~ ., switch = "y", labeller = as_labeller(cfg$facet_labels))
  }
  
  # reference line at value=1
  # if (isTRUE(cfg$swap_axes)) {
  #   p <- p + geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.35, color = "black")
  # } else {
  #   p <- p + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35, color = "black")
  # }
  
  # log scale on value axis
  if (isTRUE(cfg$use_log_scale)) {
    if (isTRUE(cfg$swap_axes)) p <- p + scale_x_log10() else p <- p + scale_y_log10()
  }
  
  # limits on value axis
  if (!is.null(cfg$value_limits)) {
    if (isTRUE(cfg$swap_axes)) p <- p + coord_cartesian(xlim = cfg$value_limits)
    else p <- p + coord_cartesian(ylim = cfg$value_limits)
  }
  
  p
}

# ---- 6) Run ----
df <- make_data(cfg)
p  <- make_plot(df, cfg)
print(p)

if (!is.null(cfg$save_path)) {
  ggsave(
    filename = cfg$save_path,
    plot = p,
    width = cfg$width,
    height = cfg$height,
    dpi = cfg$dpi,
    device = if (grepl("\\.pdf$", cfg$save_path, ignore.case = TRUE)) grDevices::cairo_pdf else NULL
  )
}

