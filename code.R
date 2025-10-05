library(ggplot2)
library(dplyr)

#true regression function
true_m <- function(x) {
  sin((x^3 + 0.1)^(-1))
}

# function to estimate h_AMISE 
calculate_h_amise <- function(X, Y, N) {
  n <- length(X)
  
  #split data into N blocks
  if (N == 1) {
    block_indices <- rep(1, n)
  } else {
    block_indices <- cut(seq_along(X), breaks = N, labels = FALSE)
  }
  
  theta_sum <- 0
  sigma_sum <- 0
  
  #loop over blocks
  for (j in 1:N) {
    block_mask <- block_indices == j
    X_block <- X[block_mask]
    Y_block <- Y[block_mask]
    n_block <- length(X_block)
    
    if (n_block >= 6) {
      #fit quartic polynomial
      fit <- lm(Y_block ~ X_block + I(X_block^2) + I(X_block^3) + I(X_block^4))
      beta <- coef(fit)
      
      #calculate second derivative
      for (i in 1:n_block) {
        x_val <- X_block[i]
        second_deriv <- 2 * beta[3] + 6 * beta[4] * x_val + 12 * beta[5] * x_val^2
        theta_sum <- theta_sum + second_deriv^2
      }
      
      #calculate residuals for sigma^2
      predictions <- predict(fit)
      sigma_sum <- sigma_sum + sum((Y_block - predictions)^2)
    }
  }
  
  #compute theta_22 and sigma^2 estimates
  theta_22_hat <- theta_sum / n
  sigma2_hat <- sigma_sum / (n - 5 * N)
  supp_X <- diff(range(X))
  
  #h_AMISE formula with quartic kernel constant 35
  h_amise <- n^(-1/5) * ((35 * sigma2_hat * supp_X) / theta_22_hat)^(1/5)
  
  return(h_amise)
}

#function to find optimal N using Cp criterion
find_optimal_N <- function(X, Y) {
  n <- length(X)
  N_max <- max(min(floor(n / 20), 5), 1)
  
  if (N_max == 1) return(1)
  
  rss_values <- numeric(N_max)
  
  for (N in 1:N_max) {
    if (N == 1) {
      block_indices <- rep(1, n)
    } else {
      block_indices <- cut(seq_along(X), breaks = N, labels = FALSE)
    }
    
    rss_N <- 0
    for (j in 1:N) {
      block_mask <- block_indices == j
      X_block <- X[block_mask]
      Y_block <- Y[block_mask]
      
      if (length(X_block) >= 6) {
        fit <- lm(Y_block ~ X_block + I(X_block^2) + I(X_block^3) + I(X_block^4))
        predictions <- predict(fit)
        rss_N <- rss_N + sum((Y_block - predictions)^2)
      }
    }
    rss_values[N] <- rss_N
  }
  
  #Cp criterion computation
  cp_values <- rep(NA, N_max)
  for (N in 1:(N_max - 1)) {
    cp_values[N] <- (rss_values[N] / (rss_values[N_max] / (n - 5 * N_max))) - (n - 10 * N)
  }
  
  #return index of minimum Cp
  which.min(cp_values[1:(N_max - 1)])
}

#simulation
set.seed(123)
n_values <- c(100, 500, 1000, 2000)
alpha_beta_values <- list(c(1, 1), c(2, 2), c(2, 5), c(6,1), c(1,6))

results <- data.frame()
optimal_N_results <- data.frame()

for (n in n_values) {
  for (ab in alpha_beta_values) {
    alpha <- ab[1]
    beta <- ab[2]
    shape_label <- paste0("Beta(", alpha, ",", beta, ")")
    
    #generate data
    X <- rbeta(n, alpha, beta)
    Y <- true_m(X) + rnorm(n, 0, 1)
    
    #find optimal N
    optimal_N <- find_optimal_N(X, Y)
    
    #store optimal N
    optimal_N_results <- rbind(optimal_N_results, data.frame(
      n = n,
      alpha = alpha,
      beta = beta,
      shape_label = shape_label,
      optimal_N = optimal_N
    ))
    
    #calculate h_AMISE for N = 1,...,5
    for (N in 1:5) {
      h_amise <- calculate_h_amise(X, Y, N)
      
      results <- rbind(results, data.frame(
        n = n,
        N = N,
        alpha = alpha,
        beta = beta,
        shape_label = shape_label,
        h_amise = h_amise,
        is_optimal = (N == optimal_N)
      ))
    }
  }
}

#plot h_AMISE vs N
p1 <- results %>%
  filter(n == 1000) %>%
  group_by(N, shape_label) %>%
  summarise(mean_h = mean(h_amise, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = N, y = mean_h, color = shape_label)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "h_AMISE vs Number of Blocks (N)",
       x = "Number of Blocks (N)",
       y = "h_AMISE",
       color = "Distribution") +
  theme_minimal()

ggsave("h_amise_vs_N.png", p1, width = 8, height = 6, dpi = 300)

#plot optimal N vs n
p2 <- optimal_N_results %>%
  ggplot(aes(x = n, y = optimal_N, color = shape_label)) +
  geom_point(size = 3) +
  labs(title = "Optimal N vs Sample Size",
       x = "Sample Size (n)",
       y = "Optimal Number of Blocks (N)",
       color = "Distribution") +
  theme_minimal()

ggsave("optimal_N_vs_n.png", p2, width = 8, height = 6, dpi = 300)

#plot optimal bandwidth by distribution
p3 <- results %>%
  filter(is_optimal) %>%
  group_by(shape_label) %>%
  summarise(mean_h = mean(h_amise, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = shape_label, y = mean_h, fill = shape_label)) +
  geom_col(alpha = 0.8) +
  labs(title = "Optimal Bandwidth by Distribution",
       x = "Beta Distribution",
       y = "h_AMISE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("bandwidth_by_distribution.png", p3, width = 8, height = 6, dpi = 300)

p4 <- results %>%
  filter(is_optimal == TRUE) %>% # only plot optimal bandwidths for clarity
  group_by(n, shape_label) %>%
  summarise(mean_h = mean(h_amise, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = n, y = mean_h, color = shape_label)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = unique(results$n)) +
  labs(title = "h_AMISE vs Sample Size (n)",
       x = "Sample Size (n)",
       y = "Optimal h_AMISE",
       color = "Distribution") +
  theme_minimal()

ggsave("h_amise_vs_nsample.png", p4, width = 8, height = 6, dpi = 300)
