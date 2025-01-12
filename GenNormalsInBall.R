library(ggplot2)

ntrials <- 50000 # from 500
dim.min <- 3
dim.max <- dim.min + 97 # from 100
dim.range <- dim.max - dim.min + 1

# for dim = 1, we have the trivial case
Z <- matrix(0, dim.max, ntrials)
Z[1, ] <- sample(c(-1, 1), ntrials, replace = TRUE)
Z[dim.max - 1, ] <- rep(0, ntrials)
Z[dim.max, ] <- rep(0, ntrials)


# Beta generation
times_beta <- rep (0, dim.range)
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
  times_beta[j] <- t[["user.self"]]/ntrials
}




# Create a data frame from the times vector
times_df <- data.frame(
  Dimension = 1:dim.max,
  Time = times_beta
)

# Create the time series plot using ggplot2
plot <- ggplot(times_df, aes(x = Dimension, y = Time)) +
  geom_line() +
  geom_point() +
  labs(title = "Time Series of Computation Times",
       x = "Dimension",
       y = "Elapsed Time (seconds)") +
  theme_minimal()
print(plot)