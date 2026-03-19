## ======================================================================
## Single-file script: STP2/STP3 regions on the 2-simplex (mu plot only)
## Used for geometric illustration in appendix of Koobs and Koning (2026)
## Edit ONLY the cfg block.
## ======================================================================

suppressPackageStartupMessages(library(plotly))

# ---- 0) Easy parameters to tune ----
cfg <- list(
  # Fixed rows (mu_1 and mu_2)
  p_mu1 = c(16, 7, 1) / 24,
  p_mu2 = c(12, 6, 6) / 24,

  # Extra mu points to display (optional)
  p_mu3 = c(6, 6, 12) / 24,
  p_mu4 = c(2, 4, 18) / 24,

  # Grid resolution over simplex (higher = smoother, slower)
  simplex_resolution = 350,

  # Plot style
  stp2_only_color = "#E69F00",
  stp3_color = "#009E73",
  simplex_color = "#d0d8e8",
  mu_color = "black",
  marker_size_region = 3.5,
  marker_size_mu = 7,

  # Labels
  x_lab = "P(X=1)",
  y_lab = "P(X=2)",
  z_lab = "P(X=3)",
  legend_stp2 = "STP2 only",
  legend_stp3 = "STP3 (subset STP2)",

  # Output
  show_in_viewer = TRUE,
  save_html = "Figs/STP2_STP3_mu_simplex.html"
)

stopifnot(length(cfg$p_mu1) == 3, length(cfg$p_mu2) == 3)
stopifnot(abs(sum(cfg$p_mu1) - 1) < 1e-12, abs(sum(cfg$p_mu2) - 1) < 1e-12)
stopifnot(all(cfg$p_mu1 >= 0), all(cfg$p_mu2 >= 0))
stopifnot(cfg$simplex_resolution >= 50)

# ---- 1) TP checks for candidate third row p ----
is_stp2 <- function(p, p1, p2) {
  M <- rbind(p1, p2, p)

  row_pairs <- combn(1:3, 2, simplify = FALSE)
  col_pairs <- combn(1:3, 2, simplify = FALSE)

  for (rp in row_pairs) {
    for (cp in col_pairs) {
      if (det(M[rp, cp, drop = FALSE]) < 0) return(FALSE)
    }
  }
  TRUE
}

is_stp3 <- function(p, p1, p2) {
  M <- rbind(p1, p2, p)
  det(M) >= 0 && is_stp2(p, p1, p2)
}

# ---- 2) Build simplex grid ----
r <- cfg$simplex_resolution
grid <- expand.grid(
  x = seq(0, 1, length.out = r),
  y = seq(0, 1, length.out = r)
)
grid <- subset(grid, x + y <= 1)
grid$z <- 1 - grid$x - grid$y

tp_mat <- t(apply(grid[, c("x", "y", "z")], 1, function(v) {
  p <- as.numeric(v)
  c(stp2 = is_stp2(p, cfg$p_mu1, cfg$p_mu2), stp3 = is_stp3(p, cfg$p_mu1, cfg$p_mu2))
}))

grid$stp2 <- tp_mat[, "stp2"]
grid$stp3 <- tp_mat[, "stp3"]
grid$region <- ifelse(grid$stp3, "stp3", ifelse(grid$stp2, "stp2_only", "none"))

stp2_only <- subset(grid, region == "stp2_only")
stp3_pts <- subset(grid, region == "stp3")

# ---- 3) Mu points shown on top of regions ----
mu_mat <- rbind(cfg$p_mu1, cfg$p_mu2, cfg$p_mu3, cfg$p_mu4)
mu_df <- data.frame(
  x = mu_mat[, 1],
  y = mu_mat[, 2],
  z = mu_mat[, 3],
  mu = paste0("\u03BC", 1:4)
)

# ---- 4) Simplex face ----
simplex <- data.frame(
  x = c(1, 0, 0),
  y = c(0, 1, 0),
  z = c(0, 0, 1)
)


# ---- 5) Plot ----
plot_ly() %>%
  # 1. Simplex triangle surface
  add_trace(type = 'mesh3d',
            x = simplex$x, y = simplex$y, z = simplex$z,
            i = c(0), j = c(1), k = c(2),
            opacity = 0.15, color = I("#d0d8e8"),
            showlegend = FALSE) %>%
  
  # 2. TP₂ only region (orange points)
  add_trace(type = 'scatter3d', mode = 'markers',
            x = tp2_only$x, y = tp2_only$y, z = tp2_only$z,
            marker = list(size = 4, color = '#E69F00',
                          symbol = 'circle'),
            name = "TP₂ only") %>%
  
  # 3. TP₃ region (green points)
  add_trace(type = 'scatter3d', mode = 'markers',
            x = tp3_points$x, y = tp3_points$y, z = tp3_points$z,
            marker = list(size = 4, color = '#009E73',
                          symbol = 'circle'),
            name = "TP₃ (⊆ TP₂)") %>%
  
  # 4. θ points (added last to sit on top)
  add_trace(type = 'scatter3d', mode = 'markers+text',
            x = A[,1], y = A[,2], z = A[,3],
            marker = list(size = 8, color = 'black',
                          symbol = 'circle',
                          line = list(color = 'white', width = 2)),
            text = theta, textposition = 'top center',
            showlegend = FALSE) %>%
  
  # 5. Layout (same as before)
  layout(
    scene = list(
      xaxis = list(title = "P(X=1)", range = c(0, 1),
                   titlefont = list(size = 16), tickfont = list(size = 14)),
      yaxis = list(title = "P(X=2)", range = c(0, 1),
                   titlefont = list(size = 16), tickfont = list(size = 14)),
      zaxis = list(title = "P(X=3)", range = c(0, 1),
                   titlefont = list(size = 16), tickfont = list(size = 14)),
      aspectmode = "cube"
    ),
    legend = list(
      x = 1.02,         # Just outside the plot area on the right
      xanchor = "left",
      y = 0.9,
      font = list(size = 14),
      itemsizing = "constant",
      bordercolor = "black",
      borderwidth = 0
    ),
    margin = list(l = 0, r = 0, t = 0, b = 0),  # Remove outer whitespace
    title = NULL
  )







