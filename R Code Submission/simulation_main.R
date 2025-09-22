# Load required libraries, functions and settings
source("functions/util.R")
source("functions/getTestResults.R")
source("functions/getStudentizedPermutationTest.R")
source("settings.R")

# Type I Error Rate -------------------------------------------------------
# Running this loop produces the type-1 error rate results from the paper
for (i in 1:nrow(df_t1e)) {
  res <- simulateT1E(
    Nxy = df_t1e$Nxy[i],
    nx = df_t1e$nx[i],
    ny = df_t1e$ny[i],
    sigmaX = df_t1e$sigmaX[i],
    sigmaY = df_t1e$sigmaY[i],
    niveau = df_t1e$niveau[i],
    seed = df_t1e$seed[i],
    dist = df_t1e$dist[i],
    alphax = df_t1e$alphax[i],
    alphay = df_t1e$alphay[i],
    betax = df_t1e$betax[i],
    betay = df_t1e$betay[i],
    likert = df_t1e$likert[i],
    lambda = df_t1e$lambda[i],
    iterationTime = NA,
    nsim = 1e5
  )

  df_t1e$t1e_BM[i] <- res$t1e_BM
  df_t1e$t1e_C2[i] <- res$t1e_C2
  df_t1e$t1e_permu[i] <- res$t1e_permu
  df_t1e$iterationTime[i] <- res$iterationTime

  print(paste(i, "/", nrow(df_t1e)))
}

# Save results
write.csv(df_t1e, file = "results/df_t1e.csv")

# Power and coverage -------------------------------------------------------------------
# Running this loop produces the power and coverage results from the paper
for (i in 1:nrow(df_pwr_cov)) {
  res <- simulatePWRCOV(
    Nxy = df_pwr_cov$Nxy[i],
    nx = df_pwr_cov$nx[i],
    ny = df_pwr_cov$ny[i],
    theta = df_pwr_cov$theta[i],
    mux = df_pwr_cov$mux[i],
    muy = df_pwr_cov$muy[i],
    sigmaX = df_pwr_cov$sigmaX[i],
    sigmaY = df_pwr_cov$sigmaY[i],
    niveau = df_pwr_cov$niveau[i],
    seed = df_pwr_cov$seed[i],
    dist = df_pwr_cov$dist[i],
    alphax = df_pwr_cov$alphax[i],
    alphay = df_pwr_cov$alphay[i],
    betax = df_pwr_cov$betax[i],
    betay = df_pwr_cov$betay[i],
    likert = df_pwr_cov$likert[i],
    lambdax = df_pwr_cov$lambdax[i],
    lambday = df_pwr_cov$lambday[i],
    iterationTime = NA,
    separation = NA,
    nsim = 1e5
  )

  df_pwr_cov$pwr_BM[i] <- res$pwr_BM
  df_pwr_cov$pwr_C2[i] <- res$pwr_C2
  df_pwr_cov$pwr_permu[i] <- res$pwr_permu
  df_pwr_cov$cov_BM[i] <- res$cov_BM
  df_pwr_cov$cov_C2[i] <- res$cov_C2
  df_pwr_cov$cov_permu[i] <- res$cov_permu
  df_pwr_cov$iterationTime[i] <- res$iterationTime
  df_pwr_cov$separation[i] <- res$separation

  print(paste(i, "/", nrow(df_pwr_cov)))
}

write.csv(df_pwr_cov, file = "results/df_pwr_cov.csv")
