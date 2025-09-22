# Load required libraries, functions and settings
source("functions/util.R")
source("functions/getTestResults.R")
source("functions/getStudentizedPermutationTest.R")
source("settings.R")

# Import the simulation results
t1e_results <- read.csv("./results/df_t1e.csv")
pwr_cov_results <- read.csv("./results/df_pwr_cov.csv")

# Type I Error Rate -------------------------------------------------------

## Setting 187 -------------------------------------------------

df_t1e187 <- df_t1e[187,]

res187 <- simulateT1E(
  Nxy = df_t1e187$Nxy[1],
  nx = df_t1e187$nx[1],
  ny = df_t1e187$ny[1],
  sigmaX = df_t1e187$sigmaX[1],
  sigmaY = df_t1e187$sigmaY[1],
  niveau = df_t1e187$niveau[1],
  seed = df_t1e187$seed[1],
  dist = df_t1e187$dist[1],
  alphax = df_t1e187$alphax[1],
  alphay = df_t1e187$alphay[1],
  betax = df_t1e187$betax[1],
  betay = df_t1e187$betay[1],
  likert = df_t1e187$likert[1],
  lambda = df_t1e187$lambda[1],
  iterationTime = NA,
  nsim = 1e5
)

# Should be equal
t1e_results[187, c("t1e_BM", "t1e_C2", "t1e_permu")]
res187[1, c("t1e_BM", "t1e_C2", "t1e_permu")]

## Setting 286 -------------------------------------------------

df_t1e286 <- df_t1e[286,]

res286 <- simulateT1E(
  Nxy = df_t1e286$Nxy[1],
  nx = df_t1e286$nx[1],
  ny = df_t1e286$ny[1],
  sigmaX = df_t1e286$sigmaX[1],
  sigmaY = df_t1e286$sigmaY[1],
  niveau = df_t1e286$niveau[1],
  seed = df_t1e286$seed[1],
  dist = df_t1e286$dist[1],
  alphax = df_t1e286$alphax[1],
  alphay = df_t1e286$alphay[1],
  betax = df_t1e286$betax[1],
  betay = df_t1e286$betay[1],
  likert = df_t1e286$likert[1],
  lambda = df_t1e286$lambda[1],
  iterationTime = NA,
  nsim = 1e5
)

# Should be equal
t1e_results[286, c("t1e_BM", "t1e_C2", "t1e_permu")]
res286[1, c("t1e_BM", "t1e_C2", "t1e_permu")]


## Setting 601 -------------------------------------------------

df_t1e601 <- df_t1e[601,]

res601 <- simulateT1E(
    Nxy = df_t1e601$Nxy[1],
    nx = df_t1e601$nx[1],
    ny = df_t1e601$ny[1],
    sigmaX = df_t1e601$sigmaX[1],
    sigmaY = df_t1e601$sigmaY[1],
    niveau = df_t1e601$niveau[1],
    seed = df_t1e601$seed[1],
    dist = df_t1e601$dist[1],
    alphax = df_t1e601$alphax[1],
    alphay = df_t1e601$alphay[1],
    betax = df_t1e601$betax[1],
    betay = df_t1e601$betay[1],
    likert = df_t1e601$likert[1],
    lambda = df_t1e601$lambda[1],
    iterationTime = NA,
    nsim = 1e5
)

# Should be equal
t1e_results[601, c("t1e_BM", "t1e_C2", "t1e_permu")]
res601[1, c("t1e_BM", "t1e_C2", "t1e_permu")]

# Power and Coverage ------------------------------------------------------

# Setting 31 --------------------------------------------------------------

df_pwr_cov31 <- df_pwr_cov[31,]

res31 <- simulatePWRCOV(
  Nxy = df_pwr_cov31$Nxy[1],
  nx = df_pwr_cov31$nx[1],
  ny = df_pwr_cov31$ny[1],
  theta = df_pwr_cov31$theta[1],
  mux = df_pwr_cov31$mux[1],
  muy = df_pwr_cov31$muy[1],
  sigmaX = df_pwr_cov31$sigmaX[1],
  sigmaY = df_pwr_cov31$sigmaY[1],
  niveau = df_pwr_cov31$niveau[1],
  seed = df_pwr_cov31$seed[1],
  dist = df_pwr_cov31$dist[1],
  alphax = df_pwr_cov31$alphax[1],
  alphay = df_pwr_cov31$alphay[1],
  betax = df_pwr_cov31$betax[1],
  betay = df_pwr_cov31$betay[1],
  likert = df_pwr_cov31$likert[1],
  lambdax = df_pwr_cov31$lambdax[1],
  lambday = df_pwr_cov31$lambday[1],
  iterationTime = NA,
  separation = NA,
  nsim = 1e5
)

pwr_cov_results[31, c("pwr_BM", "pwr_C2", "pwr_permu", "cov_BM", "cov_C2", "cov_permu")]
res31[1, c("pwr_BM", "pwr_C2", "pwr_permu", "cov_BM", "cov_C2", "cov_permu")]

