# Function applying the Brunner-Munzel test, the C2 test, and (optionally) a Studentized permutation test
# to (simulated) data from two groups. Returns confidence intervals, test statistics, p-values, and effect
# estimates for each test across nsim simulations.
# Input:
## - X:           # Matrix [nx × nsim]; simulated data for group X (each column = one dataset)
## - Y:           # Matrix [ny × nsim]; simulated data for group Y (each column = one dataset)
## - alpha:       # Significance level (default = 0.05)
## - permutation: # If TRUE, also perform the Studentized permutation test (default = TRUE)
## - nperm:       # Number of permutations used in the permutation test (default = 10000)
## - seed:        # Seed for reproducibility (default = random draw from 1:100000) if permutation=TRUE
# Output:
## A data.table containing for each test and simulation:
## - test:        Name of the test ("BM-Test", "C2-Test", or "Studentized Permutation test")
## - lower:       Lower bound of the confidence interval for θ
## - theta:       Point estimate of the effect (θ)
## - upper:       Upper bound of the confidence interval for θ
## - variance:    Variance estimate used for test statistic
## - ties:        Estimated proportion of ties (used in C2-test)
## - statistic:   Test statistic value (BM Tn, C2, or permutation Tn)
## - pval:        P-value of the test
## - rejH0:       Logical; TRUE if H0 rejected at level α
## - separation:  TRUE if θ equals 0 or 1 (boundary case)
## - sample:      Label of the simulation sample (e.g., "Sample 1", "Sample 2", ...)
getTestResults <- function(X, Y, alpha = 0.05, permutation = TRUE,
                           nperm = 10000, seed = sample(1:100000, 1)) {
  # Convert data to matrix
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  # Since each col represents a sample, the number of cols should be equal
  if (ncol(X) != ncol(Y)) {
    stop("number of samples should be equal for both groups")
  }

  # Extract number of cols
  d <- ncol(X)

  # Template for results
  results <- data.table(
    test = character(),
    lower = numeric(),
    theta = numeric(),
    upper = numeric(),
    variance = numeric(),
    ties = numeric(),
    statistic = numeric(),
    pval = numeric(),
    rejH0 = logical(),
    separation = logical()
  )

  # Loop iterating through all samples
  for (i in 1:d) {
    # Extracting sample i
    x <- X[, i]
    y <- Y[, i]

    # Compute lengths and combined length
    nx <- length(x)
    ny <- length(y)
    Nxy <- nx + ny
    m <- min(nx, ny)

    # Combine and rank the data
    xy <- c(x, y)
    rxy <- rank(xy)
    rx <- rank(x)
    ry <- rank(y)
    rxry <- c(rx, ry)
    rxmax <- rank(x, ties.method = "max")
    rxmin <- rank(x, ties.method = "min")
    rxymax <- rank(xy, ties.method = "max")[1:nx]
    rxymin <- rank(xy, ties.method = "min")[1:nx]

    # Calculate placements and their variance
    pl <- rxy - rxry
    mplx <- mean(pl[1:nx])
    mply <- mean(pl[(nx + 1):Nxy])
    plx <- (1 / ny) * (rxy[1:nx] - rx)
    ply <- (1 / nx) * (rxy[(nx + 1):Nxy] - ry)
    vx <- var(plx)
    vy <- var(ply)

    # Effect estimate in Equation (2)
    theta <- mean(ply)

    # Brunner-Munzel Test -----------------------------------------------------

    # DeLong variance estimator in Equation (4)
    sigmaBF <- Nxy * (vx / nx + vy / ny)

    # Degrees of freedom in Equation (10)
    df <- (vx / nx + vy / ny)^2 / (vx^2 / (nx^2 * (nx - 1)) + vy^2 / (ny^2 * (ny - 1)))

    # Exception handling
    df[is.na(df)] <- 1000
    thetaBM <- theta
    theta0 <- theta == 0 # Indicator; True if theta == 0
    theta1 <- theta == 1 # Indicator; True if theta == 1

    # Effect according to appendix A.4
    thetaBM[theta0] <- if (any(x %in% y)) 1 / (2 * nx * ny) else 1 / (nx * ny)
    thetaBM[theta1] <- if (any(x %in% y)) (2 * nx * ny - 1) / (2 * nx * ny) else (nx * ny - 1) / (nx * ny)
    sigmaBF[sigmaBF == 0] <- if (any(x %in% y)) Nxy / (2 * nx^2 * ny^2) else 2 * Nxy / (nx^2 * ny^2)
    se <- sqrt(sigmaBF / Nxy)

    # Test statistic in Equation (9)
    Tn_BM <- sqrt(Nxy) * (thetaBM - 1 / 2) / sqrt(sigmaBF)
    pval_BM <- 2 * min(pt(Tn_BM, df = df), 1 - pt(Tn_BM, df = df))

    # Confidence interval in Equation (12)
    lower_BM <- theta - se * qt(1 - alpha / 2, df = df)
    upper_BM <- theta + se * qt(1 - alpha / 2, df = df)

    # C2-Test -----------------------------------------------------------------

    # Estimator for the probability of ties in Equation (6)
    tau <- 1 / ny * (mean(rxymax) - mean(rxymin) - (mean(rxmax) - mean(rxmin)))

    # Unbiased variance estimator in Equation (7)
    sigmaN <- 1 / (nx * (nx - 1) * ny * (ny - 1)) *
      (sum((pl - c(rep(mplx, nx), rep(mply, ny)))^2) -
        nx * ny * (theta * (1 - theta) - 0.25 * tau))

    # Estimator of q in Equation (18)
    q <- sigmaN / (theta * (1 - theta))

    # Test statistic of the new test in Equation (20)
    C2 <- (4 / q) * (theta - 1 / 2)^2

    # Exception handling; Test statistic in Equation (22)
    C2[sigmaN == 0] <- 4 * m * (theta - 1 / 2)^2
    pval_C2 <- 1 - pchisq(C2, df = 1)

    # 1-alpha quantile of a Chi(1)-distribution
    cAlpha <- qchisq(1 - alpha, df = 1)

    # Confidence interval formula in Equation (21)
    ctr <- (2 * theta + q * cAlpha) / (2 * (1 + q * cAlpha))
    margin <- (sqrt(cAlpha) * sqrt(4 * q * theta * (1 - theta) + q^2 * cAlpha)) /
      (2 * (1 + q * cAlpha))
    lower_r <- ctr - margin
    upper_r <- ctr + margin

    # Confidence interval formula in Equation (22)
    lower_BK <- (2 * m * theta + cAlpha) / (2 * (m + cAlpha)) - (sqrt(cAlpha) *
      sqrt(4 * m * theta * (1 - theta) + cAlpha)) / (2 * (m + cAlpha))
    upper_BK <- (2 * m * theta + cAlpha) / (2 * (m + cAlpha)) + (sqrt(cAlpha) *
      sqrt(4 * m * theta * (1 - theta) + cAlpha)) / (2 * (m + cAlpha))

    # Summary -----------------------------------------------------------------

    # Save results from Brunner-Munzel test
    results_BM <- data.table(
      test = "BM-Test",
      lower = lower_BM,
      theta = thetaBM,
      upper = upper_BM,
      variance = sigmaBF,
      ties = tau,
      statistic = Tn_BM,
      pval = pval_BM,
      rejH0 = pval_BM < alpha,
      separation = theta0 | theta1
    )

    # Save results for C2-Test
    results_C2 <- data.table(
      test = "C2-Test",
      lower = ifelse(theta0 | theta1, lower_BK, lower_r),
      theta = theta,
      upper = ifelse(theta0 | theta1, upper_BK, upper_r),
      variance = sigmaN * Nxy,
      ties = tau,
      statistic = C2,
      pval = pval_C2,
      rejH0 = pval_C2 < alpha,
      separation = theta0 | theta1
    )

    # Add test results to the output table
    results <- rbind(results, results_BM, results_C2)

    # If permutation==TRUE, the studentized permutation test by Neubert-Brunner
    # (2007) is run and added
    if (permutation) {
      results_perm <- getPermutationResults(
        x = x, y = y, Tn_BM = Tn_BM,
        theta = thetaBM, se = se, nperm = nperm, alpha = alpha,
        seed = seed, separation = theta0 | theta1
      )
      results <- rbind(results, results_perm)
    }
  }

  # Results are ordered by test and samples are labeled
  results <- results[order(results$test), ]
  results$sample <- if (permutation) {
    rep(paste("Sample", 1:d), 3)
  } else {
    rep(paste("Sample", 1:d), 2)
  }

  # Return all test results
  return(results)
}

