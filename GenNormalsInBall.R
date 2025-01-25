library(ggplot2)

ntrials <- 2000 # from 500 - go to 5000
dim.min <- 3
dim.max <- 200 # from 100
dim.range <- dim.max - dim.min + 1

######################################
#### Functions to test randomness ####
######################################

# Function to calculate vector norms
vector_norm <- function(v) {
  sqrt(rowSums(v^2))
}

# Function to calculate pairwise distances
pairwise_distances <- function(vectors) {
  n <- nrow(vectors)
  distances <- matrix(NA, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      distances[i, j] <- sqrt(sum((vectors[i, ] - vectors[j, ])^2))
    }
  }
  distances[lower.tri(distances)] <- distances[upper.tri(distances)]
  diag(distances) <- 0
  return(distances[upper.tri(distances)])
}

# Main function to test randomness of unit vectors
test_unit_vectors <- function(vectors) {
  n <- nrow(vectors)
  d <- ncol(vectors)
  
  # 1. Norm Test
  norms <- vector_norm(vectors)
  if (any(abs(norms - 1) > 1e-6)) {
    cat("Norm Test: Failed. Some vectors are not unit vectors.\n")
  } else {
    cat("Norm Test: Passed. All vectors are unit vectors.\n")
  }
  
  # 2. Pairwise Distance Test
  cat("Pairwise Distance Test:\n")
  distances <- pairwise_distances(vectors)
  mean_dist <- mean(distances)
  cat(sprintf("Mean Pairwise Distance: %.4f\n", mean_dist))
  
  # Compare with expected mean distance for uniform distribution (Monte Carlo approximation)
  expected_mean <- sqrt(2) # true as d approaches infinity
  cat(sprintf("Expected Mean Distance: %.4f\n", expected_mean))
  
  if (abs(mean_dist - expected_mean) > 0.1) {
    cat("Pairwise Distance Test: Failed. Distribution deviates from uniformity.\n")
  } else {
    cat("Pairwise Distance Test: Passed.\n")
  }
}

######################################
#### Main program to generate RVs ####
######################################

# Generating random unit vectors using the scaling of Normal vectors
times_normal <- rep (0, dim.max)
for (j in 1:dim.max) {
  # Time each iteration for each dimension
  t <- system.time({
  Z <- matrix(rnorm(j * ntrials), nrow = j, ncol = ntrials)
  Z <- apply(Z, 2, function(col) col / sqrt(sum(col^2)))
  })
  # Save the time taken for each dimension
  times_normal[j] <- t[["user.self"]]/ntrials
}
print("Normal generation time to generate dim.max")
print(times_normal[dim.max])
print("Testing randomness:")
test_unit_vectors(t(Z))

# Generating random unit vectors using spherical coordinates
times_spherical <- rep (0, dim.range)
for (j in dim.min:dim.max) {
  # Time each iteration for each dimension
  t <- system.time({
    # Initialize the matrix to store the generated random vectors
    Z <- matrix(0, j, ntrials)
    # Generate the thetas
    theta <- matrix(0, j - 1, ntrials)
    accepted_samples <- numeric(ntrials)
    for (i in 1:(j-2)){
      k <- 0
      while (k < ntrials) {
        U <- runif(ntrials - k, 0, 2 * pi)
        temp <- runif(ntrials - k, 0, 1)
        accepted <- temp <= sin(U)^(j - i - 1)
        num_accepted <- sum(accepted)
        if (num_accepted != 0){
          accepted_samples[(k + 1):(k + num_accepted)] <- U[accepted]
        }
        k <- k + num_accepted
      }
      theta[i, ] <- accepted_samples
    }
    theta[j - 1, ] <- runif(ntrials, 0, 2 * pi)
    # Calculate Z's from thetas
    Z[j, ] <- cos(theta[1, ])
    prod_sin <- sin(theta[1, ])
    if (j > 2) {
      for (i in (j - 1):2) {
        Z[i, ] <- prod_sin * cos(theta[j - i + 1, ])
        prod_sin <- prod_sin * sin(theta[j - i + 1, ])
      }
    }
    Z[1, ] <- prod_sin
  })
  # Save the time taken for each dimension
  times_spherical[j] <- t[["user.self"]]/ntrials
}
print("Spherical generation time to generate dim.max")
print(times_spherical[dim.max])
print("Testing randomness:")
test_unit_vectors(t(Z))

# Generating random unit vectors using the Beta distribution from Proposition 6
times_beta <- rep (0, dim.max)
for (j in dim.min:dim.max) {
  # Time each iteration for each dimensions
  t <- system.time({
    # Initialize the matrix to store the generated random vectors
    Z <- matrix(0, j, ntrials)
    # For d = 1, we have the trivial case
    Z[1, ] <- sample(c(-1, 1), ntrials, replace = TRUE)
    # Proposition 4 Step 1
    r <- rep(1, ntrials)
    # Proposition 4 Step 2
    for (i in 2:(j - 2)) {
      X <- rbeta(ntrials, (i - 1) / 2, (i - 1) / 2)
      Z[i, ] <- r * (2 * X - 1)
      r <- sqrt(r^2 - Z[i, ]^2)
    }
    Z[1, ] <- r * Z[1, ]
    r <- sqrt(r^2 - Z[1, ]^2)
    # Proposition 4 Step 3
    theta <- runif(ntrials, 0, 2 * pi)
    Z[j - 1, ] <- r * cos(theta)
    Z[j, ] <- r * sin(theta)
  })
  # Save the time taken for each dimension
  times_beta[j] <- t[["user.self"]]/ntrials
}
print("Beta generation time to generate dim.max")
print(times_beta[dim.max])
print("Testing randomness:")
test_unit_vectors(t(Z))

# Generating random unit vectors using the Gamma distribution from Proposition 7
times_gamma <- rep (0, dim.max)
for (j in dim.min:dim.max) {
  # Time each iteration for each dimension
  t <- system.time({
    # Initialize the matrix to store the generated random vectors
    Z <- matrix(0, j, ntrials)
    # For d = 1, we have the trivial case
    Z[1, ] <- sample(c(-1, 1), ntrials, replace = TRUE)
    # Proposition 4 Step 1
    r <- rep(1, ntrials)
    # Proposition 4 Step 2
    for (i in 2:(j - 2)) {
      X <- rgamma(ntrials, (i - 1) / 2, (i - 1) / 2)
      Y <- rgamma(ntrials, (i - 1) / 2, (i - 1) / 2)
      Z[i, ] <- r * ((X - Y) / (X + Y))
      r <- sqrt(r^2 - Z[i, ]^2)
    }
    Z[1, ] <- r * Z[1, ]
    r <- sqrt(r^2 - Z[1, ]^2)
    # Proposition 4 Step 3
    theta <- runif(ntrials, 0, 2 * pi)
    Z[j - 1, ] <- r * cos(theta)
    Z[j, ] <- r * sin(theta)
  })
  # Save the time taken for each dimension
  times_gamma[j] <- t[["user.self"]]/ntrials
}
print("Gamma generation time to generate dim.max")
print(times_gamma[dim.max])
print("Testing randomness:")
test_unit_vectors(t(Z))

# Generating random unit vectors using the ratio-of-uniforms method 1
times_ratio1 <- rep (0, dim.range)
for (j in 4:dim.max) {
  # Initialize the matrix to store the generated random vectors
  Z <- matrix(0, j, ntrials)
  # For d = 1:3, we have the trivial case
  Z[1, ] <- sample(c(-1, 1), ntrials, replace = TRUE)
  Z[2, ] <- runif(ntrials, 0, 1)
  Z[2, ] <- 1/2 + asin(Z[2, ]) / pi
  Z[3, ] <- 1/2 + 2 * runif(ntrials, 0, 1)
  # Time each iteration for each dimension
  t <- system.time({
    # Proposition 4 Step 1
    r0 <- rep(1, ntrials)
    # Proposition 4 Step 2
    if (j - 2 >= 4) {
      for (i in 4:(j - 2)) {
      r <- (i - 3) / 2
      mbound <- sqrt(r^r / (r + 1)^(r + 1))
      accepted_samples <- numeric(ntrials)
      k <- 0
      while (k < ntrials) {
        U <- runif(ntrials - k, 0, 1)
        V <- runif(ntrials - k, -mbound, mbound)
        accepted <- V^2 <= U^2 * (1 - U^(2 / r))
        num_accepted <- sum(accepted)
        if (num_accepted != 0){
          accepted_samples[(k + 1):(k + num_accepted)] <- V[accepted] / U[accepted]
        }
        k <- k + num_accepted
      }
      Z[i, ] <- r0 * accepted_samples
      r0 <- sqrt(r0^2 - Z[i, ]^2)
      }
      Z[1, ] <- r0 * Z[1, ]
      r0 <- sqrt(r0^2 - Z[1, ]^2)
      Z[2, ] <- r0 * Z[2, ]
      r0 <- sqrt(r0^2 - Z[2, ]^2)
      Z[3, ] <- r0 * Z[3, ]
      r0 <- sqrt(r0^2 - Z[3, ]^2)
    }
    # Proposition 4 Step 3
    theta <- runif(ntrials, 0, 2 * pi)
    Z[j - 1, ] <- r0 * cos(theta)
    Z[j, ] <- r0 * sin(theta)
  })
  # Save the time taken for each dimension
  times_ratio1[j] <- t[["user.self"]]/ntrials
}
print("Ratio of Uniforms Method 1 time to generate dim.max")
print(times_ratio1[dim.max])
print("Testing randomness:")
test_unit_vectors(t(Z))

# Generating random unit vectors using the ratio-of-uniforms method 2
times_ratio2 <- rep (0, dim.range)
for (j in 4:dim.max) {
  # Initialize the matrix to store the generated random vectors
  Z <- matrix(0, j, ntrials)
  # For d = 1, we have the trivial case
  Z[1, ] <- sample(c(-1, 1), ntrials, replace = TRUE)
  Z[2, ] <- runif(ntrials, 0, 1)
  Z[2, ] <- 1/2 + asin(Z[2, ]) / pi
  Z[3, ] <- 1/2 + 2 * runif(ntrials, 0, 1)
  # Time each iteration for each dimension
  t <- system.time({
    # Proposition 4 Step 1
    r0 <- rep(1, ntrials)
    # Proposition 4 Step 2
    if ((j - 2) >= 4) {
      for (i in 4:(j - 2)) {
        r <- (i - 3) / 2
        Omega <- r * gamma( r / 2) * sqrt(pi) / (2 * gamma((r + 3) / 2))
        M <- (r / (r + 1))^((r + 1) / 2)
        constant <- 2 / (M * Omega)
        accepted_samples <- numeric(ntrials)
        k <- 0
        while (k < ntrials) {
          U <- runif(ntrials - k, 0, 1)
          temp <- runif(ntrials - k, 0, 1)
          accepted <- temp <= constant * U * sqrt(1 - U^(2 / r))
          num_accepted <- sum(accepted)
          if (num_accepted != 0){ 
            accepted_samples[(k + 1):(k + num_accepted)] <- U[accepted]
          }
          k <- k + num_accepted
        }
        Z[i, ] <- accepted_samples
        vbound <- Z[i, ] * sqrt((1 - Z[i, ]^(2 / r)))
        V <- runif(ntrials, -vbound, vbound)
        Z[i, ] <- r0 * V / Z[i, ]
        r0 <- sqrt(r0^2 - Z[i, ]^2)
     }
     Z[1, ] <- r0 * Z[1, ]
      r0 <- sqrt(r0^2 - Z[1, ]^2)
      Z[2, ] <- r0 * Z[2, ]
      r0 <- sqrt(r0^2 - Z[2, ]^2)
      Z[3, ] <- r0 * Z[3, ]
      r0 <- sqrt(r0^2 - Z[3, ]^2)
    }
    # Proposition 4 Step 3
    theta <- runif(ntrials, 0, 2 * pi)
    Z[j - 1, ] <- r0 * cos(theta)
    Z[j, ] <- r0 * sin(theta)
  })
# Save the time taken for each dimension
  times_ratio2[j] <- t[["user.self"]]/ntrials
}
print("Ratio of Uniforms Method 2 time to generate dim.max")
print(times_ratio2[dim.max])
print("Testing randomness:")
test_unit_vectors(t(Z))

#Plot the time series of computation times excluding Spherical
# Create a data frame from the times vector
times_df <- data.frame(
  Dimension = 1:dim.max,
  Time = c(times_beta,
           times_gamma,
           times_normal,
           #times_spherical,
           times_ratio1,
           times_ratio2),
  Type = rep(c("Beta",
               "Gamma",
               "Normal",
               #"Spherical",
               "Ratio M1",
               "Ratio M2"), each = dim.max)
)

plot <- ggplot(times_df, aes(x = Dimension, y = Time, color = Type)) +  #linetype = Type
  geom_line() +
  geom_point(aes(shape = Type)) +
  labs(title = "Time Series of Computation Times",
       x = "Dimension",
       y = "Elapsed Time (seconds)") +
  theme_minimal() #+
  #scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
  #scale_shape_manual(values = c(16, 17, 18, 19, 20, 21))

print(plot)

#Plot the time series of computation times including Spherical
# Create a data frame from the times vector
times_df <- data.frame(
  Dimension = 1:dim.max,
  Time = c(times_beta,
           times_gamma,
           times_normal,
           times_spherical,
           times_ratio1,
           times_ratio2),
  Type = rep(c("Beta",
               "Gamma",
               "Normal",
               "Spherical",
               "Ratio M1",
               "Ratio M2"), each = dim.max)
)

plot <- ggplot(times_df, aes(x = Dimension, y = Time, color = Type)) + 
  geom_line() +
  geom_point(aes(shape = Type)) +
  labs(title = "Time Series of Computation Times",
       x = "Dimension",
       y = "Elapsed Time (seconds)") +
  theme_minimal() #+
  #scale_color_manual(values = c("red", "blue", "green", "purple", "orange", "brown")) +
  #scale_shape_manual(values = c(16, 17, 18, 19, 20, 21))

print(plot)

#Plot the time series of computation times including Spherical and zooming in for d <= 50
# Create a data frame from the times vector
times_df <- data.frame(
  Dimension = 1:50,
  Time = c(times_beta[1:50],
           times_gamma[1:50],
           times_normal[1:50],
           times_spherical[1:50],
           times_ratio1[1:50],
           times_ratio2[1:50]),
  Type = rep(c("Beta",
               "Gamma",
               "Normal",
               "Spherical",
               "Ratio M1",
               "Ratio M2"), each = 50)
)

plot <- ggplot(times_df, aes(x = Dimension, y = Time, color = Type)) +  #linetype = Type
  geom_line() +
  geom_point(aes(shape = Type)) +
  labs(title = "Time Series of Computation Times",
       x = "Dimension",
       y = "Elapsed Time (seconds)") +
  theme_minimal() #+
  #scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
  #scale_shape_manual(values = c(16, 17, 18, 19, 20, 21))

print(plot)

#Plot the time series of computation times including Spherical and zooming in for d <= 50
# Create a data frame from the times vector
times_df <- data.frame(
  Dimension = 1:50,
  Time = c(times_beta[1:50],
           times_gamma[1:50],
           times_normal[1:50],
           #times_spherical[1:50],
           times_ratio1[1:50],
           times_ratio2[1:50]),
  Type = rep(c("Beta",
               "Gamma",
               "Normal",
               #"Spherical",
               "Ratio M1",
               "Ratio M2"), each = 50)
)

plot <- ggplot(times_df, aes(x = Dimension, y = Time, color = Type)) +  #linetype = Type
  geom_line() +
  geom_point(aes(shape = Type)) +
  labs(title = "Time Series of Computation Times",
       x = "Dimension",
       y = "Elapsed Time (seconds)") +
  theme_minimal() #+
  #scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
  #scale_shape_manual(values = c(16, 17, 18, 19, 20, 21))

print(plot)

print("Difference between Normal and Beta:")
print(times_normal[dim.max]-times_beta[dim.max])