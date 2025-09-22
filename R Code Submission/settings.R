# Load required libraries and import theme
source("functions/util.R")

# Type I Error Rate -------------------------------------------------------

# Seed: set in order to ensure that the individual seeds per simulation run remain
# the same
set.seed(20112024)

nx <- c(15, 15, 30, 30, 20, 40, 45, 30, 60, 60, 40, 80, 75, 50, 100)
ny <- c(15, 30, 15, 30, 40, 20, 45, 60, 30, 60, 80, 40, 75, 100, 50)
nLevels <- c(0.001, 0.005, 0.01, 0.05) # Significance levels
Nxy <- nx + ny

## Normal ------------------------------------------------------------------

# Settings for the normal distribution
sigmaX <- 1
sigmaY <- c(1, 3)
nCases <- length(nx) * length(sigmaY) * length(nLevels)

# Create the tibble to iterate through in the simulation
df_normal <- tibble(
  Nxy = rep(Nxy, length(sigmaY) * length(nLevels)),
  nx = rep(nx, length(sigmaY) * length(nLevels)),
  ny = rep(ny, length(sigmaY) * length(nLevels)),
  sigmaX = rep(sigmaX, nCases),
  sigmaY = rep(sigmaY, each = length({{ nx }}) * length(nLevels)),
  niveau = rep(rep(nLevels, each = length({{ nx }})), length({{ sigmaY }})),
  seed = ceiling(runif(nCases, min = 1, max = 1e5)),
  t1e_BM = numeric(nCases),
  t1e_C2 = numeric(nCases),
  t1e_permu = numeric(nCases),
  dist = "normal"
)

df_normal <- df_normal[order(df_normal$niveau, df_normal$sigmaY, df_normal$Nxy), ]

## Laplace ------------------------------------------------------------------

# Laplace settings are the same as normal, hence copy the tibble above and rename "dist"
df_laplace <- df_normal
df_laplace$dist <- "laplace"
df_laplace <- df_laplace[order(df_laplace$niveau, df_laplace$sigmaY, df_laplace$Nxy), ]

## Ordinal -----------------------------------------------------------------

# Settings for the orderd categorical data
alphax <- c(1, 5, 1, 2)
alphay <- c(1, 2, 5, 2)
betax <- c(1, 5, 1, 5)
betay <- c(1, 2, 5, 5)
nCases <- length(nx) * length(alphax) * length(nLevels)

# Create the tibble to iterate through in the simulation
df_ordinal <- tibble(
  Nxy = rep(Nxy, length(alphax) * length(nLevels)),
  nx = rep(nx, length(alphax) * length(nLevels)),
  ny = rep(ny, length(alphax) * length(nLevels)),
  alphax = rep(alphax, each = length({{ nx }}) * length(nLevels)),
  alphay = rep(alphay, each = length({{ nx }}) * length(nLevels)),
  betax = rep(betax, each = length({{ nx }}) * length(nLevels)),
  betay = rep(betay, each = length({{ nx }}) * length(nLevels)),
  niveau = rep(rep(nLevels, each = length({{ nx }})), times = length({{ alphax }})),
  seed = ceiling(runif(nCases, min = 1, max = 1e5)),
  t1e_BM = numeric(nCases),
  t1e_C2 = numeric(nCases),
  t1e_permu = numeric(nCases),
  dist = "beta",
  likert = TRUE # Indicator; TRUE if data should be discretized
)


df_ordinal <- df_ordinal[order(df_ordinal$niveau, df_ordinal$alphax, df_ordinal$betay, df_ordinal$Nxy), ]


## Beta --------------------------------------------------------------------

# Beta settings are the same as ordinal, hence copy the tibble above and rename "dist"
df_beta <- df_ordinal
df_beta$likert <- FALSE # Set likert==FALSE to avoid discretization
df_beta <- df_beta[order(df_beta$niveau, df_beta$alphax, df_beta$betay, df_beta$Nxy), ]


## Exponential -----------------------------------------------------------------

# Settings for the exponential distribution
lambda <- 1
nCases <- length(nx) * length(lambda) * length(nLevels)

# Create the tibble to iterate through in the simulation
df_exp <- tibble(
  Nxy = rep(Nxy, length(lambda) * length(nLevels)),
  nx = rep(nx, length(lambda) * length(nLevels)),
  ny = rep(ny, length(lambda) * length(nLevels)),
  lambda = rep(lambda, each = length({{ nx }}) * length(nLevels)),
  niveau = rep(nLevels, times = length({{ nx }}) * length({{ lambda }})),
  seed = ceiling(runif(nCases, min = 1, max = 1e5)),
  t1e_BM = numeric(nCases),
  t1e_C2 = numeric(nCases),
  t1e_permu = numeric(nCases),
  dist = "exponential"
)

## Poisson -----------------------------------------------------------------

# Poisson settings are the same as exponential, hence copy the tibble above and rename "dist"
df_pois <- df_exp
df_pois$dist <- "poisson"

## Merging everything ------------------------------------------------------

# Take all the above tibbles and join them
df_t1e <- bind_rows(df_normal, df_laplace, df_beta, df_ordinal, df_exp, df_pois)

# Add the iteration time to monitor runtime length as additional column
df_t1e$iterationTime <- vector("numeric", nrow(df_t1e))


# Power and Coverage --------------------------------------------------------------

# Seed: set in order to ensure that the individual seeds per simulation run remain
# the same
set.seed(03122024)

# Parameter settings for the power simulation
nx <- c(15, 15, 30)
ny <- c(15, 30, 15)
niveau <- c(0.05)
theta <- c(seq(0.5, 0.9, 0.05),0.91,0.92,0.93,0.94,0.95)

## Normal -----------------------------------------------------------------

# Parameters for the normal
sigmaX <- 1
sigmaY <- c(1, 3)
mux <- 0 # mux is fixed to 0 and muy computed according to the target theta

nCases <- length(nx) * length(niveau) * length(sigmaY) * length(theta)
df_normal <- expand.grid(
  nx = nx,
  ny = ny,
  sigmaX = sigmaX,
  sigmaY = sigmaY,
  niveau = niveau,
  theta = theta
) %>%
  unique() %>%
  filter((nx == 15 & ny == 15) | (nx == 15 & ny == 30) | (nx == 30 & ny == 15)) %>%
  mutate(
    Nxy = nx + ny,
    mux = mux,
    muy = qnorm(theta) * sqrt(sigmaX^2 + sigmaY^2), # muy is computed to match theta
    seed = ceiling(runif(nCases, min = 1, max = 1e5)),
    pwr_BM = numeric(nCases),
    pwr_C2 = numeric(nCases),
    pwr_permu = numeric(nCases),
    cov_BM = numeric(nCases),
    cov_C2 = numeric(nCases),
    cov_permu = numeric(nCases),
    dist = "normal"
  ) %>%
  dplyr::select(
    Nxy, nx, ny, theta, dist, sigmaX, sigmaY, mux, muy, niveau, seed,
    pwr_BM, pwr_C2, pwr_permu, cov_BM, cov_C2, cov_permu
  ) %>%
  arrange(theta)


## Ordinal --------------------------------------------------------------

# Max possible value is 0.9
theta <- seq(0.5, 0.9, 0.05)


# Parameters for the ordered categorical data
alphax <- 1
betax <- 1
betay <- 1
alphay <- numeric(length(theta))

# alphay is computed to match a target theta using the calculateTheta function
# from the util.R script; See also Appendix A.3
for (i in 1:length(theta)) {
  alphay[i] <- uniroot(f = function(x) {
    g <- calculateTheta(
      c = seq(0, 1, 0.2), alphay = x, betay = betay,
      alphax = alphax, betax = betax
    ) - theta[i]
    return(g)
  }, interval = c(0, 1e5), extendInt = "yes")$root
}

nCases <- length(nx) * length(niveau) * length(alphay)
df_ordinal <- expand.grid(
  nx = nx,
  ny = ny,
  theta = theta,
  niveau = niveau
) %>%
  filter((nx == 15 & ny == 15) | (nx == 15 & ny == 30) | (nx == 30 & ny == 15)) %>%
  unique() %>%
  mutate(
    Nxy = nx + ny,
    seed = ceiling(runif(nCases, min = 1, max = 1e5)),
    alphax = alphax,
    betax = betax,
    betay = betay,
    alphay = rep(alphay, each = 3),
    pwr_BM = numeric(nCases),
    pwr_C2 = numeric(nCases),
    pwr_permu = numeric(nCases),
    cov_BM = numeric(nCases),
    cov_C2 = numeric(nCases),
    cov_permu = numeric(nCases),
    dist = "beta",
    likert = TRUE
  ) %>%
  dplyr::select(
    Nxy, nx, ny, theta, dist, alphax, alphay, betax, betay, likert,
    niveau, seed, pwr_BM, pwr_C2, pwr_permu, cov_BM, cov_C2, cov_permu
  ) %>%
  arrange(theta)

## Exponential ----------------------------------------------------------------

theta <- c(seq(0.5, 0.9, 0.05),0.91,0.92,0.93,0.94,0.95)

# Parameters for the ordered categorical data
lambdax <- 1

nCases <- length(nx) * length(niveau) * length(theta)
df_exp <- expand.grid(
  nx = nx,
  ny = ny,
  theta = theta,
  niveau = niveau
) %>%
  unique() %>%
  filter((nx == 15 & ny == 15) | (nx == 15 & ny == 30) | (nx == 30 & ny == 15)) %>%
  mutate(
    Nxy = nx + ny,
    seed = ceiling(runif(nCases, min = 1, max = 1e5)),
    lambdax = lambdax,
    lambday = lambdax / (1 / (1 - theta) - 1), # lambday computed to match a target theta
    pwr_BM = numeric(nCases),
    pwr_C2 = numeric(nCases),
    pwr_permu = numeric(nCases),
    cov_BM = numeric(nCases),
    cov_C2 = numeric(nCases),
    cov_permu = numeric(nCases),
    dist = "exponential"
  ) %>%
  dplyr::select(
    Nxy, nx, ny, theta, dist, lambdax, lambday, niveau, seed, pwr_BM,
    pwr_C2, pwr_permu, cov_BM, cov_C2, cov_permu
  ) %>%
  arrange(theta)


## Merging everything ------------------------------------------------------

# All power simulation settings are merged
df_pwr_cov <- bind_rows(df_normal, df_ordinal, df_exp)

# Add the iteration time to monitor runtime length as additional column
df_pwr_cov$iterationTime <- vector("numeric", nrow(df_pwr_cov))

# Add seperation as column to monitor the fraction of completely separated samples
df_pwr_cov$separation <- vector("numeric", nrow(df_pwr_cov))


