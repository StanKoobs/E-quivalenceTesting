## ======================================================================
## Single-file script: compare TOST e-values vs symmetric LR e-values 
## in sequential t-test setting
## Used for Figure 4 in appendix of Koobs and Koning (2026)
## ======================================================================


tfunc <- function(x) {
  n <- length(x)
  if (n < 2) {
    stop("n < 2")
  }
  
  s <- sd(x)
  return(sqrt(n) * mean(x) / s)
}


set.seed(123)

M <- 200000
LRvec <- rep(0, M)
beta2vec <- rep(0, M)
LRvec2 <- rep(0, M)

for (m in 1:M) {
  # Parameters
  n <- 10
  mu <- 0
  sigma <- 1
  delta0 <- 0    # under H1
  delta1 <- -1   # null part 1
  delta3 <- 1    # null part 2
  
  # Generate data under delta = 1 (mu = 1, sigma = 1)
  x <- rnorm(n, mean = mu, sd = sigma)
  
  # Sample statistics
  xbar <- mean(x)
  s <- sd(x)
  t_stat <- xbar / (s / sqrt(n))
  df <- n - 1
  
  u <- x / sqrt(sum(x^2))
  
  beta2_stat <- t_stat^2 / (t_stat^2 + n - 1)
  beta2vec[m] <- beta2_stat
  
  # Noncentrality parameters
  lambda0 <- sqrt(n) * delta0
  lambda1 <- sqrt(n) * delta1
  lambda3 <- sqrt(n) * delta3
  
  # Compute densities of t-statistic under each lambda
  f0 <- dt(t_stat, df = df, ncp = lambda0)
  f1 <- dt(t_stat, df = df, ncp = lambda1)
  f3 <- dt(t_stat, df = df, ncp = lambda3)
  
  b0 <- dbeta(beta2_stat, 1/2, df / 2)
  b1 <- dbeta(beta2_stat, 1/2, df / 2, ncp = n * delta3^2)
  
  F0 <- df(t_stat^2, 1, df)
  F1 <- df(t_stat^2, 1, df, ncp = n * delta3^2)
  
  # Likelihood ratio (point vs mixture)
  p_null <- 0.5  # weight in the mixture
  f_null <- p_null * f1 + (1 - p_null) * f3
  #LR <- f0 / f_null
  LR <- min(f0 / f1, f0 / f3)
  
  LRvec[m] <- LR
  
  LRvec2[m] <- F0 / F1
}

mean(LRvec)
mean(LRvec2)


#### Code for LR process
M <- 50000
n <- 50
LRvec <- matrix(1, M, n)
beta2vec <- rep(0, M)
LRvec2 <- matrix(1, M, n)

for (m in 1:M) {
  # Parameters
  mu <- 0
  sigma <- 1
  delta0 <- 0    # under H1
  delta1 <- -0.5   # null part 1
  delta3 <- 0.5    # null part 2
  xvec <- rnorm(1, mu, sigma)
  for (i in 2:n) {
    xvec <- c(xvec, rnorm(1, mean = mu, sd = sigma))
    
    # Sample statistics
    xbar <- mean(xvec)
    s <- sd(xvec)
    t_stat <- xbar / (s / sqrt(i))
    df <- i - 1
    
    u <- xvec / sqrt(sum(xvec^2))
    
    beta2_stat <- t_stat^2 / (t_stat^2 + i - 1)
    beta2vec[m] <- beta2_stat
    
    # Noncentrality parameters
    lambda0 <- sqrt(i) * delta0
    lambda1 <- sqrt(i) * delta1
    lambda3 <- sqrt(i) * delta3
    
    df <- i - 1
    
    # Compute densities of t-statistic under each lambda
    f0 <- dt(t_stat, df = df, ncp = lambda0)
    f1 <- dt(t_stat, df = df, ncp = lambda1)
    f3 <- dt(t_stat, df = df, ncp = lambda3)
    
    b0 <- dbeta(beta2_stat, 1/2, df / 2)
    b1 <- dbeta(beta2_stat, 1/2, df / 2, ncp = i * delta3^2)
    
    F0 <- df(t_stat^2, 1, df)
    F1 <- df(t_stat^2, 1, df, ncp = i * delta3^2)
    
    # Likelihood ratio (point vs mixture)
    p_null <- 0.5  # weight in the mixture
    f_null <- p_null * f1 + (1 - p_null) * f3
    #LR <- f0 / f_null
    LR <- min(f0 / f1, f0 / f3)
    
    LRvec[m, i] <- LR
    
    LRvec2[m, i] <- F0 / F1
  }
  # Generate data under delta = 1 (mu = 1, sigma = 1)
}

colMeans(LRvec)
colMeans(LRvec2)

rejvec <- rep(0, n)
rejvec2 <- rep(0, n)
rejmat <- matrix(0, M, n)
rejmat2 <- matrix(0, M, n)
for (i in 1:n) {
  for (m in 1:M) {
    rejmat[m, i] <- any(LRvec[m, 1:i] > 20)
    rejmat2[m, i] <- any(LRvec2[m, 1:i] > 20)
  }
  rejvec[i] <- sum(rejmat[,i]) / nrow(rejmat)
  rejvec2[i] <-  sum(rejmat2[,i]) / nrow(rejmat2)
}




library(ggplot2)

#n <- 50  # define sample size range
df <- data.frame(
  x = 1:n,
  TOST = rejvec,
  Ratio = rejvec2
)

df2 <- data.frame(
  x = 1:50,
  TOST = colMeans(LRvec),
  Ratio = colMeans(LRvec2)
)

cols <- c("Symmetric LR" = "#E69F00",
          "TOST"         = "#0072B2")

library(ggplot2)

ggplot(df, aes(x = x)) +
  geom_line(aes(y = Ratio, color = "Symmetric LR"), linewidth = 1.2) +
  geom_line(aes(y = TOST,  color = "TOST"), linewidth = 1.2) +
  scale_color_manual(values = cols) +   # Assign custom colors
  labs(
    x = "Sample size",
    y = "Rejection Probability"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",           # <- Remove legend completely
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )


ggplot(df2, aes(x = x)) +
  geom_line(aes(y = Ratio, color = "Symmetric LR"), linewidth = 1.2) +
  geom_line(aes(y = TOST,  color = "TOST"), linewidth = 1.2) +
  scale_color_manual(values = cols) +   # Assign custom colors
  labs(
    x = "Sample size",
    y = "E-value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",           # <- Remove legend completely
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )



# Output
cat("t-statistic:", round(t_stat, 3), "\n")
cat("Density under H1 (delta = 0):", round(f0, 5), "\n")
cat("Mixture density under H0 (delta = ±1):", round(f_null, 5), "\n")

