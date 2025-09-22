# Packages ----------------------------------------------------------------
packages <- c("tidyverse", "matrixStats", "MASS", "VGAM", "data.table", "ggpubr")

invisible(
  lapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      message(sprintf("Package '%s' is not installed. Please install
          the package to run the code.", pkg))
    }
  })
)

# Theme for plotting ------------------------------------------------------
theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.key.size = unit(1.5, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(2.5, "cm"),
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    strip.text.x = element_text(size = 20, colour = "black", angle = 0, face = "bold"),
    strip.text.y = element_text(size = 15, colour = "black", angle = 270, face = "bold")
  )

# Helper functions --------------------------------------------------------


# Function simulating nsim samples of either Beta-distributed data or, if cut = TRUE, ordinal categorical data
# Input:
## - ncat:      # Number of categories (used only if cut = TRUE)
## - nx:        # Sample size for group X
## - ny:        # Sample size for group Y
## - alphax:    # Alpha parameter for the beta distribution (group X)
## - alphay:    # Alpha parameter for the beta distribution (group Y)
## - betax:     # Beta parameter for the beta distribution (group X)
## - betay:     # Beta parameter for the beta distribution (group Y)
## - nsim:      # Number of simulation replicates (default = 1e5)
## - cut:       # If TRUE, transform beta values into ordinal Likert-type categories
## - seed:      # Random seed for reproducibility
# Output:
## A list with two elements:
## - sX: A matrix of shape [nx × nsim], either beta or ordinal categorical data for group X
## - sY: A matrix of shape [ny × nsim], either beta or ordinal categorical data for group Y
simulateTwoSamplesOrdinal <- function(ncat, nx, ny, alphax, alphay, betax, betay,
                                      nsim = 1e5, cut = F, seed = 123) {
  set.seed(seed)

  result <- list()

  nx <- round(nx, 0)
  ny <- round(ny, 0)

  sX <- rbeta(nx * nsim, shape1 = alphax, shape2 = betax)
  sY <- rbeta(ny * nsim, shape1 = alphay, shape2 = betay)

  if (cut) {
    breaks <- seq(0, 1, length.out = ncat + 1)
    result$sX <- matrix(as.integer(cut(sX, breaks = breaks, labels = 1:ncat)),
      ncol = nsim, nrow = nx
    )
    result$sY <- matrix(as.integer(cut(sY, breaks = breaks, labels = 1:ncat)),
      ncol = nsim, nrow = ny
    )
  } else {
    result$sX <- matrix(sX, ncol = nsim, nrow = nx)
    result$sY <- matrix(sY, ncol = nsim, nrow = ny)
  }
  return(result)
}

# Function simulating nsim samples of normal, exponential, Poisson, or Laplace data
# with corresponding parameter settings for two independent groups
# Input:
## - distx:     # Distribution type for group X ("normal", "exponential", "poisson", "laplace")
## - disty:     # Distribution type for group Y (same options as distx)
## - nx:        # Sample size for group X
## - ny:        # Sample size for group Y
## - mux:       # Mean for group X (used if distx is "normal" or "laplace")
## - muy:       # Mean for group Y (used if disty is "normal" or "laplace")
## - sdx:       # Standard deviation or scale parameter for group X (used if distx is "normal" or "laplace")
## - sdy:       # Standard deviation or scale parameter for group Y (used if distx is "normal" or "laplace")
## - lambdax:   # Rate for exponential or Poisson distribution (group X)
## - lambday:   # Rate for exponential or Poisson distribution (group Y)
## - nsim:      # Number of simulation replicates (default = 1e5)
## - seed:      # Random seed for reproducibility
# Output:
## A list with two elements:
## - sX: A matrix of shape [nx × nsim], simulated data for group X
## - sY: A matrix of shape [ny × nsim], simulated data for group Y
simulateTwoSamples <- function(distx = c("normal", "exponential", "poisson", "laplace"),
                               disty = c("normal", "exponential", "poisson", "laplace"),
                               nx, ny, mux, muy, sdx, sdy, lambdax = NA, lambday = NA,
                               nsim = 1e5, seed = NA) {
  set.seed(seed)
  result <- list()

  # To prevent sample sizes to become float numbers
  nx <- round(nx, 0)
  ny <- round(ny, 0)

  sX <- switch(distx,
    "normal" = rnorm(nx * nsim, mean = mux, sd = sdx),
    "laplace" = rlaplace(nx * nsim, location = mux, scale = sdx),
    "exponential" = rexp(nx * nsim, rate = lambdax),
    "poisson" = rpois(nx * nsim, lambda = lambdax)
  )

  sY <- switch(disty,
    "normal" = rnorm(ny * nsim, mean = muy, sd = sdy),
    "laplace" = rlaplace(ny * nsim, location = muy, scale = sdy),
    "exponential" = rexp(ny * nsim, rate = lambday),
    "poisson" = rpois(ny * nsim, lambda = lambday)
  )

  result$sX <- matrix(sX, ncol = nsim, nrow = nx)
  result$sY <- matrix(sY, ncol = nsim, nrow = ny)

  return(result)
}

# Compute true theta in case of ordered categorical data (Appendix A.3)
# Input:
## - c: Vector containing cutpoints for discretization
## - alphax: Alpha parameter for beta distribution (group X)
## - alphay: Alpha parameter for beta distribution (group Y)
## - betax:  Beta parameter for beta  distribution (group X)
## - betay:  Beta parameter for beta distribution (group Y)
# Output:
## - theta: corresponding true Mann-Whitney Effect
calculateTheta <- function(c = seq(0, 1, 0.2), alphay, betay, alphax, betax) {
  K <- length(c) - 1

  s <- numeric(K)
  b <- numeric(K)

  for (i in 2:(K + 1)) {
    term1 <- pbeta(c[i], alphay, betay) - pbeta(c[i - 1], alphay, betay)
    inner_sum <- sum(diff(pbeta(c[1:i - 1], alphax, betax)))
    s[i - 1] <- term1 * inner_sum
  }

  for (i in 1:K) {
    term2 <- pbeta(c[i + 1], alphax, betax) - pbeta(c[i], alphax, betax)
    term3 <- pbeta(c[i + 1], alphay, betay) - pbeta(c[i], alphay, betay)
    b[i] <- 0.5 * term2 * term3
  }

  p <- sum(s) + sum(b)
  return(p)
}


# Simulation functions -----------------------------------------------------

# Type-1 Error Rate; This function iterates through a tibble of all simulation settings in the script
# simulation.R
# Input:
## - Nxy:        # Total sample size
## - nx:         # Sample size for group X
## - ny:         # Sample size for group Y
## - sigmaX:     # Standard deviation for group X
## - sigmaY:     # Standard deviation for group Y
## - niveau:     # Significance level (e.g., 0.05)
## - seed:       # Random seed for reproducibility
## - t1e_BM:     # t1e for BM-Test (kept empty, argument just for technical reasons)
## - t1E_C2:     # t1e for C2-Test (kept empty, argument just for technical reasons)
## - t1e_permu:  # t1e for permutation test (kept empty, argument just for technical reasons)
## - dist:       # Distribution type
## - alphax:     # Alpha parameter for beta distribution (group X)
## - alphay:     # Alpha parameter for beta distribution (group Y)
## - betax:      # Beta parameter for beta  distribution (group X)
## - betay:      # Beta parameter for beta distribution (group Y)
## - likert:     # If TRUE; Beta-Distribution is dichotomized
## - lambda:     # Poisson or exponential parameter (if applicable)
## - iterationTime: # Runtime of the simulation
## - nsim:       # Number of simulation runs (default 1e5)
# Output:
## Data.frame containing the input parameters as well as the simulated
# type-1 error rate values for all three methods
simulateT1E <- function(Nxy = NA, nx = NA, ny = NA, sigmaX = NA, sigmaY = NA,
                        niveau = NA, seed = NA, t1e_BM = NA, t1E_C2 = NA, t1e_permu = NA,
                        dist = NA, alphax = NA, alphay = NA, betax = NA, betay = NA,
                        likert = NA, lambda = NA, iterationTime = NA, nsim = 1e5) {
  # Start measuring execution time
  start <- Sys.time()

  # Simulate data depending on the specified distribution
  if (dist == "normal" | dist == "poisson" | dist == "exponential" | dist == "laplace") {
    # For continuous or count data, simulate from standard parametric distributions
    data <- simulateTwoSamples(
      distx = dist,
      disty = dist,
      nx = nx,
      ny = ny,
      mux = 0,
      muy = 0,
      sdx = sigmaX,
      sdy = sigmaY,
      lambdax = lambda,
      lambday = lambda,
      nsim = nsim,
      seed = seed
    )
  } else if (dist == "beta") {
    # For ordinal data generated from a beta distribution, use a separate simulation function
    data <- simulateTwoSamplesOrdinal(
      ncat = 5, # Number of Likert-scale categories
      nx = nx,
      ny = ny,
      alphax = alphax,
      alphay = alphay,
      betax = betax,
      betay = betay,
      nsim = nsim,
      cut = likert, # If TRUE, beta-distribution is discretized according to ncat
      seed = seed
    )
  }

  # Apply statistical tests and get rejection results
  results <- getTestResults(
    X = data$sX, # Sample X
    Y = data$sY, # Sample Y
    permutation = TRUE, # Include permutation test
    alpha = niveau, # Significance level
    seed = seed,
    nperm = 1e4 # Number of permutations
  )

  # Aggregate Type I error rate (mean rejection under H0) by test type
  errorRates <- aggregate(rejH0 ~ test, FUN = "mean", data = results)

  # Extract error rates for the three tests
  t1e_BM <- errorRates[1, 2] # Type I error for baseline method
  t1e_C2 <- errorRates[2, 2] # Type I error for second method (e.g., C2)
  t1e_permu <- errorRates[3, 2] # Type I error for permutation test

  # Stop timing
  end <- Sys.time()

  # Compile all simulation parameters and results into a data frame
  df <- data.frame(
    Nxy = Nxy,
    nx = nx,
    ny = ny,
    sigmaX = sigmaX,
    sigmaY = sigmaY,
    niveau = niveau,
    seed = seed,
    t1e_BM = t1e_BM,
    t1e_C2 = t1e_C2,
    t1e_permu = t1e_permu,
    dist = dist,
    alphax = alphax,
    alphay = alphay,
    betax = betax,
    betay = betay,
    likert = likert,
    lambda = lambda,
    iterationTime = as.numeric(end - start) # Total runtime
  )

  # Return the results
  return(df)
}


# Power and coverage simulation; This function iterates through a tibble of all simulation settings in the script
# simulation.R
# Input:
## - Nxy:         # Total sample size
## - nx:          # Sample size for group X
## - ny:          # Sample size for group Y
## - mux:         # Mean for group X (used in normal distributions)
## - muy:         # Mean for group Y (used in normal distributions)
## - sigmaX:      # Standard deviation for group X
## - sigmaY:      # Standard deviation for group Y
## - niveau:      # Significance level (e.g., 0.05)
## - seed:        # Random seed for reproducibility
## - theta:       # True effect size (used to evaluate CI coverage)
## - pwr_BM:      # Power for BM-Test (kept empty, filled by function)
## - pwr_C2:      # Power for C2-Test (kept empty, filled by function)
## - pwr_permu:   # Power for permutation test (kept empty, filled by function)
## - cov_BM:      # Coverage for BM-Test (kept empty, filled by function)
## - cov_C2:      # Coverage for C2-Test (kept empty, filled by function)
## - cov_permu:   # Coverage for permutation test (kept empty, filled by function)
## - dist:        # Distribution type ("normal", "exponential", "beta")
## - alphax:      # Alpha parameter for beta distribution (group X)
## - alphay:      # Alpha parameter for beta distribution (group Y)
## - betax:       # Beta parameter for beta distribution (group X)
## - betay:       # Beta parameter for beta distribution (group Y)
## - likert:      # If TRUE; beta distribution is categorized into Likert scale
## - lambdax:     # Rate parameter for group X (used in exponential distribution)
## - lambday:     # Rate parameter for group Y (used in exponential distribution)
## - iterationTime: # Runtime of the simulation (kept empty, filled by function)
## - separation:  # Estimated proportion of separations between the two groups (kept empty, filled by function)
## - nsim:        # Number of simulation runs (default 1e5)
# Output:
## Data.frame containing the input parameters along with:
## - Estimated power for BM-Test, C2-Test, and Permutation Test
## - Estimated confidence interval coverage rates for all three methods
## - Estimated proportion of separations between the two groups
simulatePWRCOV <- function(Nxy = NA, nx = NA, ny = NA, mux = NA, muy = NA,
                           sigmaX = NA, sigmaY = NA, niveau = NA, seed = NA, theta = NA,
                           pwr_BM = NA, pwr_C2 = NA, pwr_permu = NA,
                           cov_BM = NA, cov_C2 = NA, cov_permu = NA, dist = NA,
                           alphax = NA, alphay = NA, betax = NA, betay = NA,
                           likert = NA, lambdax = NA, lambday = NA, iterationTime = NA,
                           separation = NA, nsim = 1e5) {
  # Start timing the simulation run
  start <- Sys.time()

  # Generate data depending on the distribution specified
  if (dist == "normal" | dist == "exponential") {
    # For continuous data, use simulateTwoSamples with given means and SDs/lambdas
    data <- simulateTwoSamples(
      distx = dist,
      disty = dist,
      nx = nx,
      ny = ny,
      mux = mux,
      muy = muy,
      sdx = sigmaX,
      sdy = sigmaY,
      lambdax = lambdax,
      lambday = lambday,
      nsim = nsim,
      seed = seed
    )
  } else if (dist == "beta") {
    # For ordinal data from beta distributions, use specialized simulation
    data <- simulateTwoSamplesOrdinal(
      ncat = 5, # Number of Likert-scale categories
      nx = nx,
      ny = ny,
      alphax = alphax,
      alphay = alphay,
      betax = betax,
      betay = betay,
      nsim = nsim,
      cut = likert, # If TRUE, beta-distribution is discretized according to ncat
      seed = seed
    )
  }

  # Run statistical tests on the simulated data
  results <- getTestResults(
    X = data$sX,
    Y = data$sY,
    permutation = TRUE, # Include permutation test
    alpha = niveau,
    seed = seed,
    nperm = 1e4 # Number of permutations for the test
  )

  # Estimate power: proportion of rejections under H1
  rejRates <- aggregate(rejH0 ~ test, FUN = "mean", data = results)
  pwr_BM <- rejRates[1, 2] # Power for BM-Test
  pwr_C2 <- rejRates[2, 2] # Power for C2-Test
  pwr_permu <- rejRates[3, 2] # Power for Permutation test

  # Filter result data for each test method
  BM <- filter(results, test == "BM-Test")
  C2 <- filter(results, test == "C2-Test")
  Permu <- filter(results, test == "Studentized Permutation test")

  # Estimate coverage probability: proportion of times true parameter is inside CI
  cov_BM <- mean(BM$lower < theta & theta < BM$upper)
  cov_C2 <- mean(C2$lower < theta & theta < C2$upper)
  cov_permu <- mean(Permu$lower < theta & theta < Permu$upper)

  # Stop timing
  end <- Sys.time()

  # Compile all parameters and results into a data frame
  df <- data.frame(
    Nxy = Nxy,
    nx = nx,
    ny = ny,
    theta = theta, # True parameter value
    dist = dist,
    sigmaX = sigmaX,
    sigmaY = sigmaY,
    mux = mux,
    muy = muy,
    niveau = niveau,
    seed = seed,
    pwr_BM = pwr_BM,
    pwr_C2 = pwr_C2,
    pwr_permu = pwr_permu,
    cov_BM = cov_BM,
    cov_C2 = cov_C2,
    cov_permu = cov_permu,
    alphax = alphax,
    alphay = alphay,
    betax = betax,
    betay = betay,
    likert = likert,
    lambdax = lambdax,
    lambday = lambday,
    iterationTime = as.numeric(end - start), # Runtime
    separation = mean(results$separation) # Estimated fraction of separations
  )

  # Return simulation results
  return(df)
}


# Plot functions ----------------------------------------------------------

# Wrapper function for all plots in the type-1 error rate simulation
# Input: data
# Output: Type-1 error rate plot
plotT1E <- function(data) {
  # Extract unique 'niveau' value from data (e.g., 0.05, 0.01, etc.)
  niveau_val <- unique(data$niveau)

  # Extract unique distribution type from data (e.g., "normal", "beta", etc.)
  dist <- unique(data$dist)

  # Define y-axis settings depending on 'niveau' (significance level)
  y_settings <- switch(as.character(niveau_val),
                       "0.05" = list(limits = c(0.035, 0.06), breaks = seq(0.035, 0.065, 0.005), hline = 0.05),
                       "0.01" = list(limits = c(0.000, 0.02), breaks = seq(0.000, 0.020, 0.005), hline = 0.01),
                       "0.005" = list(limits = c(0.000, 0.01), breaks = seq(0.000, 0.010, 0.002), hline = 0.005),
                       "0.001" = list(limits = c(0.000, 0.0033), breaks = seq(0.000, 0.003, 0.001), hline = 0.001)
  )

  # Choose the facet variable name dynamically depending on distribution
  # This variable will be used for rows in facet_grid()
  facet_var <- switch(dist,
                      "normal" = "facetNL",
                      "laplace" = "facetNL",
                      "beta" = "facetB",
                      "ordinal" = "facetB",
                      "exponential" = "facetPE",
                      "poisson" = "facetPE"
  )

  # Start building the ggplot
  p <- ggplot(data, aes(x = Nxy, y = t1e, group = test, color = test, shape = test)) +
    facet_grid(rows = vars(!!rlang::sym(facet_var)), cols = vars(r_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    xlab(expression(N)) +
    ylab("Simulated Type I Error Rate") +
    scale_color_manual(
      values = c("C2-Test" = "seagreen4", "BM-Test" = "orange2", "Permutation Test" = "blue2"),
      labels = c("BM-Test" = expression(T[N]^BM), "C2-Test" = expression(C^2), "Permutation Test" = expression(T[perm]^BM))
    ) +
    scale_shape_manual(
      values = c("BM-Test" = 4, "C2-Test" = 15, "Permutation Test" = 6),
      labels = c("BM-Test" = expression(T[N]^BM), "C2-Test" = expression(C^2), "Permutation Test" = expression(T[perm]^BM))
    ) +
    scale_x_continuous(
      limits = c(20, 150),
      breaks = unique(data$Nxy)
    ) +
    scale_y_continuous(
      limits = y_settings$limits,
      breaks = y_settings$breaks
    ) +
    geom_hline(yintercept = y_settings$hline, color = "grey", linetype = "dashed", linewidth = 1.2) +
    theme

  # Return the plot
  return(p)
}

# Wrapper function for power and coverage plots
# Input:
#   - data: Data frame with simulation results
#   - type: One of "power" or "coverage" to control y-axis and reference lines
# Output: power or coverage plot
plotPOWCOV <- function(data, type = c("power", "coverage")) {
  type <- match.arg(type)

  # Choose y variable and axis label depending on type
  y_var <- if (type == "power") "pwr" else "cov"
  y_label <- if (type == "power") "Simulated Power" else "Simulated Coverage"

  # Start building the base ggplot
  base_plot <- ggplot(data, aes(
    x = theta, y = .data[[y_var]],
    group = test, color = test, linetype = test
  )) +
    geom_line(linewidth = 1.2) +
    xlab(expression(theta)) +
    ylab(y_label) + # Y-axis label (power or coverage)
    scale_color_manual(
      values = c("C2-Test" = "seagreen4", "BM-Test" = "orange2", "Permutation Test" = "blue2"),
      labels = c(
        "BM-Test" = expression(T[N]^BM),
        "C2-Test" = expression(C^2),
        "Permutation Test" = expression(T[perm]^BM)
      )
    ) +
    scale_linetype_manual(
      values = c("solid", "dashed", "dotted"),
      labels = c(
        "BM-Test" = expression(T[N]^BM),
        "C2-Test" = expression(C^2),
        "Permutation Test" = expression(T[perm]^BM)
      )
    ) +
    theme

  # Add type-specific y-axis settings and horizontal reference lines
  if (type == "power") {
    # Power plot settings
    base_plot +
      scale_y_continuous(
        limits = c(0.0, 1),
        breaks = seq(0.0, 1, 0.2)
      ) +
      geom_hline(
        yintercept = 0.8,
        color = "grey", linetype = "dashed", linewidth = 1.2
      ) +
      scale_x_continuous(limits = c(0.5, 0.9), breaks = seq(0.4, 0.9, 0.1))
  } else {
    # Coverage plot settings
    base_plot +
      scale_y_continuous(
        limits = c(0.65, 1),
        breaks = c(seq(0.6, 1, 0.1), 0.95)
      ) +
      geom_hline(
        yintercept = 0.95,
        color = "grey", linetype = "dashed", linewidth = 1.2
      ) +
      geom_hline(
        yintercept = 0.925,
        color = "grey", linetype = "dotted", linewidth = 1
      ) +
      geom_hline(
        yintercept = 0.975,
        color = "grey", linetype = "dotted", linewidth = 1
      ) +
      scale_x_continuous(limits = c(0.5, 0.95), breaks = seq(0.4, 0.95, 0.1))
  }
}
