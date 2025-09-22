# Function computing results of the Studentized permutation test as described by Neubert & Brunner (2007)
# Input:
## - x:           # Numeric vector of observations for group X
## - y:           # Numeric vector of observations for group Y
## - Tn_BM:       # Test statistic from the Brunner-Munzel test (used as reference)
## - theta:       # Effect estimate (mean placement difference)
## - se:          # Standard error used for confidence interval calculation
## - nperm:       # Number of permutations
## - alpha:       # Significance level (e.g., 0.05)
## - seed:        # Random seed for reproducibility
## - separation:  # Logical; TRUE if effect size is 0 or 1 (boundary case)
# Output:
## A data.table containing:
## - test:        Name of the test ("Studentized Permutation test")
## - lower:       Lower bound of the confidence interval
## - theta:       Effect estimate
## - upper:       Upper bound of the confidence interval
## - variance:    Always NA (not estimated in permutation test, inherited from BM-Test)
## - ties:        Always NA (not used in this test)
## - statistic:   Brunner-Munzel statistic used as reference (to compute p-value)
## - pval:        P-value from the permutation distribution
## - rejH0:       Logical; TRUE if H0 rejected at level alpha
## - separation:  Logical flag passed from input
getPermutationResults <- function(x, y, Tn_BM, theta, se, nperm, alpha, seed,
                            separation) {
  set.seed(seed)

  # Prepare variables and data
  nx <- length(x)
  ny <- length(y)
  Nxy <- nx + ny
  xy <- c(x, y)
  P <- replicate(nperm, sample(1:Nxy))
  xyP <- matrix(xy[P], ncol = nperm)

  # Separate x and y from permutations
  xP <- xyP[1:nx, ]
  yP <- xyP[(nx + 1):Nxy, ]

  # Compute ranks
  rxP <- t(colRanks(xP, ties.method = "average"))
  ryP <- t(colRanks(yP, ties.method = "average"))
  rxyP <- t(colRanks(xyP, ties.method = "average"))

  # Compute rank differences
  rxryP <- rbind(rxP, ryP)
  plP <- rxyP - rxryP
  plxP <- (1 / ny) * (rxyP[1:nx, ] - rxP)
  plyP <- (1 / nx) * (rxyP[(nx + 1):Nxy, ] - ryP)

  # Variance calculations
  vxP <- colVars(plxP)
  vyP <- colVars(plyP)
  sigmaBFP <- Nxy * (vxP / nx + vyP / ny)
  sigmaBFP[sigmaBFP == 0] <- Nxy / (2 * nx * ny)

  # Compute test statistics
  thetaP <- colMeans(plyP)
  T_perm <- sqrt(Nxy) * (thetaP - 1 / 2) / sqrt(sigmaBFP)

  # p-value and rejection decision according to Equation (13)
  pval <- 2 * min(mean(T_perm <= Tn_BM), mean(T_perm >= Tn_BM))
  pval <- ifelse(pval==0,1/(nperm+1),pval)
  rejH0 <- pval < alpha

  # Confidence interval according to Equation (14)
  c1 <- quantile(T_perm, 1 - alpha / 2)
  c2 <- quantile(T_perm, alpha / 2)
  lower <- theta - c1 * se
  upper <- theta - c2 * se

  # Create results as a tibble
  results_perm <- data.table(
    test = "Studentized Permutation test",
    lower = lower,
    theta = theta,
    upper = upper,
    variance = NA,
    ties = NA,
    statistic = Tn_BM,
    pval = pval,
    rejH0 = rejH0,
    separation = separation
  )

  return(results_perm)
}
