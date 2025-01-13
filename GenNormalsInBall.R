library(ggplot2)

ntrials <- 5000 # from 500
dim.min <- 3
dim.max <- dim.min + 97 # from 100
dim.range <- dim.max - dim.min + 1

# for dim = 1, we have the trivial case
Z <- matrix(0, dim.max, ntrials)
Z[1, ] <- sample(c(-1, 1), ntrials, replace = TRUE)
Z[dim.max - 1, ] <- rep(0, ntrials)
Z[dim.max, ] <- rep(0, ntrials)

# Gamma generation
times_gamma <- rep (0, dim.range)
for (j in dim.min:dim.max) {
  # initialize the r
  r <- rep(1, ntrials)

  # time each iteration for each dimensions
  t <- system.time({
    for (i in 2:(j - 2)) {
      X <- rgamma(ntrials, (i - 1) / 2, (i - 1) / 2)
      Y <- rgamma(ntrials, (i - 1) / 2, (i - 1) / 2)
      Z[i, ] <- r * ((X - Y) / (X + Y))
      r <- sqrt(r^2 - Z[i, ]^2)
    }

    ## Generate a random angle uniformly distributed between 0 and 2π
    theta <- runif(ntrials, 0, 2 * pi)

    # Compute the coordinates (x, y) on the unit circle
    Z[j - 1, ] <- r * cos(theta)
    Z[j, ] <- r * sin(theta)
  })
  # preseve the time taken for each dimension
  times_gamma[j] <- t[["user.self"]]/ntrials
}
p_values <- apply(Z, 2, function(x) shapiro.test(x)$p.value)
print("Gamma generation p-values summary")
print(summary(p_values))

# Beta generation
times_beta <- rep (0, dim.range)
for (j in dim.min:dim.max) {
  # initialize the r
  r <- rep(1, ntrials)

  # time each iteration for each dimensions
  t <- system.time({
    for (i in 2:(j - 2)) {
      X <- rbeta(ntrials, (i - 1) / 2, (i - 1) / 2)
      Z[i, ] <- r * (2 * X - 1)
      r <- sqrt(r^2 - Z[i, ]^2)
    }

    ## Generate a random angle uniformly distributed between 0 and 2π
    theta <- runif(ntrials, 0, 2 * pi)

    # Compute the coordinates (x, y) on the unit circle
    Z[j - 1, ] <- r * cos(theta)
    Z[j, ] <- r * sin(theta)
  })

  # preseve the time taken for each dimension
  times_beta[j] <- t[["user.self"]]/ntrials
}

p_values <- apply(Z, 2, function(x) shapiro.test(x)$p.value)
print("Beta generation p-values summary")
print(summary(p_values))

# Normal random vector generation
times_normal <- rep (0, dim.range)
for (j in dim.min:dim.max) {

  # time each iteration for each dimensions
  t <- system.time({
    Z[j, ] <- rnorm(ntrials)
    Z[j, ] <- Z[j, ] / sqrt(sum(Z[j, ]^2))
  })

  # preseve the time taken for each dimension
  times_normal[j] <- t[["user.self"]]/ntrials
}

p_values <- apply(Z, 2, function(x) shapiro.test(x)$p.value)
print("Normal generation p-values summary")
print(summary(p_values))

# Create a data frame from the times vector
times_df <- data.frame(
  Dimension = 1:dim.max,
  Time = c(times_beta,
           times_gamma,
           times_normal),
  Type = rep(c("Beta",
               "Gamma",
               "Normal"), each = dim.max)
)

plot <- ggplot(times_df, aes(x = Dimension, y = Time, color = Type)) + 
  geom_line() +
  geom_point() +
  labs(title = "Time Series of Computation Times",
       x = "Dimension",
       y = "Elapsed Time (seconds)") +
  theme_minimal()

print(plot)