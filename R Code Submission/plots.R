# Load required libraries and import theme
source("functions/util.R")

# Settings of the Simulation Study ----------------------------------------

## H0 ----------------------------------------------------------------------

x_vals <- seq(-15, 15, length.out = 1000)
x_vals_pois <- seq(0,999,1)

data_H0 <- data.frame(
  x = rep(x_vals, 2),
  x_pois = rep(x_vals_pois, 2),
  density_norm = c(dnorm(x_vals, mean = 0, sd = 1), dnorm(x_vals, mean = 0, sd = 3)),
  density_laplace = c(dlaplace(x_vals, 0, 1), dlaplace(x_vals, 0, 3)),
  density_beta5511 = c(dbeta(x_vals, 5, 5), dbeta(x_vals, 1, 1)),
  density_beta5522 = c(dbeta(x_vals, 5, 5), dbeta(x_vals, 2, 2)),
  density_beta2525 = c(dbeta(x_vals, 2, 5), dbeta(x_vals, 2, 5)),
  density_poisson = c(dpois(x_vals_pois, 1), dpois(x_vals_pois, 1)),
  density_exponential = c(dexp(x_vals,1), dexp(x_vals, 1)),
  distribution_norm = factor(rep(c("N(0,1)", "N(0,9)"), each = length(x_vals))),
  distribution_laplace = factor(rep(c("L(0,1)", "L(0,3)"), each = length(x_vals))),
  distribution_beta5511 = factor(rep(c("B(5,5)", "B(1,1)"), each = length(x_vals))),
  distribution_beta5522 = factor(rep(c("B(5,5)", "B(2,2)"), each = length(x_vals))),
  distribution_beta2525 = factor(rep(c("B(2,5)", "B(2,5)"), each = length(x_vals))),
  distribution_poisson = factor(rep(c("Pois(1)", "Pois(1)"), each = length(x_vals))),
  distribution_exponential = factor(rep(c("Exp(1)", "Exp(1)"), each = length(x_vals)))
)


# Plots - Homoscedasticity

beta25 <- ggplot(data_H0, aes(x = x, y = density_beta2525, color = distribution_beta2525)) +
  geom_line(size = 1) +
  labs(title = "Density plots of B(2,5)",
       x = "x",
       y = "Density",
       color = "Distribution") +
  xlim(c(0,1)) +
  theme +
  theme(legend.position = "none")

exponential <- ggplot(data_H0, aes(x = x, y = density_exponential, color = distribution_exponential)) +
  geom_line(size = 1) +
  labs(title = "Density plots of Exp(1)",
       x = "x",
       y = "Density",
       color = "Distribution") +
  xlim(c(0,15)) +
  theme +
  theme(legend.position = "none")

poisson <- ggplot(data_H0, aes(x = x_pois, y = density_poisson, color = distribution_poisson)) +
  geom_line(size = 1) +
  labs(title = "Density plots of Pois(1)",
       x = "x",
       y = "Density",
       color = "Distribution") +
  xlim(c(0,15)) +
  theme +
  theme(legend.position = "none")

ggarrange(beta25, exponential, poisson, ncol=3,nrow=1) %>%
  ggsave(
    filename = sprintf("./plots/Supplement/section1.1_table3-hom.eps"),
    device = "eps", width = 18, height = 6
  )


# Plots - Heteroscedasticity
norm <- ggplot(data_H0, aes(x = x, y = density_norm, color = distribution_norm)) +
  geom_line(size = 1) +
  labs(title = "Density plots of N(0,1) and N(0,9)",
       x = "x",
       y = "Density",
       color = "Distribution") +
  theme

beta11 <- ggplot(data_H0, aes(x = x, y = density_beta5511, color = distribution_beta5511)) +
  geom_line(size = 1) +
  labs(title = "Density plots of B(5,5) and B(1,1)",
       x = "x",
       y = "Density",
       color = "Distribution") +
  xlim(c(0,1)) +
  theme

beta22 <- ggplot(data_H0, aes(x = x, y = density_beta5522, color = distribution_beta5522)) +
  geom_line(size = 1) +
  labs(title = "Density plots of B(5,5) and B(2,2)",
       x = "x",
       y = "Density",
       color = "Distribution") +
  xlim(c(0,1)) +
  theme

laplace <- ggplot(data_H0, aes(x = x, y = density_laplace, color = distribution_laplace)) +
  geom_line(size = 1) +
  labs(title = "Density plots of L(0,1) and L(0,3)",
       x = "x",
       y = "Density",
       color = "Distribution") +
  theme

ggarrange(norm, beta11, beta22, laplace, ncol=2,nrow=2) %>%
  ggsave(
    filename = sprintf("./plots/Supplement/section1.1_table3-het.pdf"),
    device = "pdf", width = 16.3, height = 14
  )

## H1 ----------------------------------------------------------------------

x_vals <- seq(-15, 15, length.out = 1000)

data_H1 <- data.frame(
  x = rep(x_vals, 2),
  density_norm_hom = c(dnorm(x_vals, mean = 0, sd = 1), dnorm(x_vals, mean = 2, sd = 1)),
  distribution_norm_hom = factor(rep(c("N(0,1)", "N(2,1)"), each = length(x_vals))),
  density_norm_het = c(dnorm(x_vals, mean = 0, sd = 3), dnorm(x_vals, mean = 2, sd = 1)),
  distribution_norm_het = factor(rep(c("N(0,9)", "N(2,9)"), each = length(x_vals))),
  density_beta11 = c(dbeta(x_vals, 1, 1), dbeta(x_vals, 2, 1)),
  distribution_beta11 = factor(rep(c("B(1,1)", "B(2,1)"), each = length(x_vals))),
  density_exp = c(dexp(x_vals, 1), dexp(x_vals, 0.5)),
  distribution_exp = factor(rep(c("Exp(1)", "Exp(0.5)"), each = length(x_vals)))
)

norm_hom <- ggplot(data_H1, aes(x = x, y = density_norm_hom, color = distribution_norm_hom)) +
  geom_line(size = 1) +
  labs(title = expression(paste("Density plots of N(0,1) and N(2,1) - ", theta, " = 0.92")),
       x = "x",
       y = "Density",
       color = "Distribution") +
  theme


norm_het <- ggplot(data_H1, aes(x = x, y = density_norm_het, color = distribution_norm_het)) +
  geom_line(size = 1) +
  labs(title = expression(paste("Density plots of N(0,1) and N(2,9) - ", theta, " = 0.74")),
       x = "x",
       y = "Density",
       color = "Distribution") +
  theme

beta11 <- ggplot(data_H1, aes(x = x, y = density_beta11, color = distribution_beta11)) +
  geom_line(size = 1) +
  labs(title = expression(paste("Density plots of B(1,1) and B(2,1) - ", theta, " = 0.66")),
       x = "x",
       y = "Density",
       color = "Distribution") +
  xlim(c(0,1)) +
  theme

exp <- ggplot(data_H1, aes(x = x, y = density_exp, color = distribution_exp)) +
  geom_line(size = 1) +
  labs(title = expression(paste("Density plots of Exp(1) and Exp(0.5) - ", theta, " = 2/3")),
       x = "x",
       y = "Density",
       color = "Distribution") +
  scale_color_manual(
    values = c("#00BFC4", "#F8766D")
  ) +
  xlim(c(0,10)) +
  theme

ggarrange(norm_hom, norm_het, beta11, exp, ncol=2,nrow=2) %>%
  ggsave(
    filename = sprintf("./plots/Supplement/section1.2_table4.eps"),
    device = "eps", width = 16.3, height = 14
  )


# Type-I Error Rate -------------------------------------------------------

t1e_results <- read.csv("./results/df_t1e.csv") %>%
  pivot_longer(cols = c(t1e_BM, t1e_C2, t1e_permu), names_to = "test", values_to = "t1e") %>%
  mutate(
    temp = nx,
    nx = if_else(dist %in% c("ordinal", "beta") & alphax == 1 &
      betax == 1 & alphay == 5, ny, nx),
    ny = if_else(dist %in% c("ordinal", "beta") & alphax == 1 &
      betax == 1 & alphay == 5, temp, ny),
    alphax = if_else(dist %in% c("ordinal", "beta") & alphax == 1 &
      betax == 1 & alphay == 5, 5, alphax),
    betax = if_else(dist %in% c("ordinal", "beta") & alphax == 5 &
      betax == 1 & alphay == 5, 5, betax),
    alphay = if_else(dist %in% c("ordinal", "beta") & alphax == 5 &
      betax == 5 & alphay == 5, 1, alphay),
    betay = if_else(dist %in% c("ordinal", "beta") & alphax == 5 &
      betax == 5 & alphay == 1, 1, betay)
  ) %>%
  dplyr::select(-temp) %>%
  mutate(
    dist = case_when(
      is.na(likert) ~ dist,
      likert == TRUE ~ "ordinal",
      TRUE ~ dist
    ),
    r = round(nx / ny, 2),
    r_label = paste("n1/n2 = ", r),
    test = case_when(
      test == "t1e_BM" ~ "BM-Test",
      test == "t1e_C2" ~ "C2-Test",
      test == "t1e_permu" ~ "Permutation Test"
    ),
    facetNL = paste0("(sd1, sd2) = ", "(", sigmaX, ",", sigmaY, ")"),
    facetB = paste0("(", alphax, ",", betax, ") , (", alphay, ",", betay, ")"),
    facetPE = paste0("Rate = ", lambda)
  )


## Figure 1 in Section 4.1 -------------------------------------

Figure1_left <- t1e_results %>%
  filter(
    nx > ny, sigmaX == 1, sigmaY == 3, test == "BM-Test",
    niveau == 0.05, dist == "normal"
  ) %>%
  ggplot(aes(x = Nxy, y = t1e, group = test, color = test, shape = test)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  xlab(expression(N)) +
  ylab("Simulated Type I Error Rate") +
  scale_color_manual(
    values = c("BM-Test" = "orange2"),
    labels = c("BM-Test" = expression(T[N]^BM))
  ) +
  scale_shape_manual(
    values = c("BM-Test" = 4),
    labels = c("BM-Test" = expression(T[N]^BM))
  ) +
  scale_x_continuous(
    limits = c(20, 150),
    breaks = unique(t1e_results$Nxy)
  ) +
  scale_y_continuous(
    limits = c(0.035, 0.06),
    breaks = seq(0.035, 0.065, 0.005)
  ) +
  geom_hline(yintercept = 0.05, color = "grey", linetype = "dashed", linewidth = 1.2) +
  theme +
  labs(title = expression(paste(alpha, " = 0.05")))

Figure1_right <- t1e_results %>%
  filter(
    nx > ny, sigmaX == 1, sigmaY == 3, test == "BM-Test",
    niveau == 0.005, dist == "normal"
  ) %>%
  ggplot(aes(x = Nxy, y = t1e, group = test, color = test, shape = test)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  xlab(expression(N)) +
  ylab("Simulated Type I Error Rate") +
  scale_color_manual(
    values = c("BM-Test" = "orange2"),
    labels = c("BM-Test" = expression(T[N]^BM))
  ) +
  scale_shape_manual(
    values = c("BM-Test" = 4),
    labels = c("BM-Test" = expression(T[N]^BM))
  ) +
  scale_x_continuous(
    limits = c(20, 150),
    breaks = unique(t1e_results$Nxy)
  ) +
  scale_y_continuous(
    limits = c(0.000, 0.01),
    breaks = seq(0.000, 0.010, 0.002)
  ) +
  geom_hline(yintercept = 0.005, color = "grey", linetype = "dashed", linewidth = 1.2) +
  theme +
  labs(title = expression(paste(alpha, " = 0.005")))

Figure1 <- ggarrange(Figure1_left, Figure1_right, ncol = 2)

ggsave(filename = "./plots/Main Paper/Figure1.eps", device = "eps", width = 14, height = 6)


## Figures in Section 5 -------------------------------------------------------------------------

### Normal ------------------------------------------------------------------

niveaus <- c(0.05, 0.01, 0.005, 0.001)

for (i in seq_along(niveaus)) {
  p_hom <- t1e_results %>%
    filter(dist == "normal", niveau == niveaus[i], sigmaX == sigmaY, nx >= ny) %>%
    plotT1E()

  p_het <- t1e_results %>%
    filter(dist == "normal", niveau == niveaus[i], sigmaX != sigmaY) %>%
    plotT1E()

  # Combine plots and save
  ggpubr::ggarrange(p_hom, p_het, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
    ggsave(
      filename = sprintf("./plots/Supplement/section2.1.%s_normal_t1e.eps", format(i, scientific = FALSE)),
      device = "eps", width = 16.3, height = 11
    )

  # This is a main paper plot
  if(niveaus[i] == 0.005) {
    ggpubr::ggarrange(p_hom, p_het, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
      ggsave(
        filename = "./plots/Main Paper/Figure2.eps",
        device = "eps", width = 16.3, height = 11
      )
  }

}

### Beta --------------------------------------------------------------------

for (i in seq_along(niveaus)) {
  a1b1 <- t1e_results %>%
    filter(
      dist == "beta", niveau == niveaus[i], alphax == 1, betax == 1, alphay == 1,
      betay == 1, nx >= ny
    ) %>%
    plotT1E()

  a2b5 <- t1e_results %>%
    filter(
      dist == "beta", niveau == niveaus[i], alphax == 2, betax == 5, alphay == 2,
      betay == 5, nx >= ny
    ) %>%
    plotT1E()

  a15b15a21b21 <- t1e_results %>%
    filter(
      dist == "beta", niveau == niveaus[i], alphax == 5, betax == 5,
      alphay == 1, betay == 1
    ) %>%
    plotT1E()

  a15b15b22b22 <- t1e_results %>%
    filter(
      dist == "beta", niveau == niveaus[i], alphax == 5, betax == 5,
      alphay == 2, betay == 2
    ) %>%
    plotT1E()

  ggpubr::ggarrange(a1b1, a2b5, a15b15a21b21, a15b15b22b22,
                    nrow = 4,
                    common.legend = TRUE, legend = "bottom"
  ) %>%
    ggsave(
      filename = sprintf("./plots/Supplement/section2.2.%s_beta_t1e.eps", format(i, scientific = FALSE)),
      device = "eps", width = 16.3, height = 18
    )
}

### Ordinal -----------------------------------------------------------------

for (i in seq_along(niveaus)) {
  a1b1 <- t1e_results %>%
    filter(
      dist == "ordinal", niveau == niveaus[i], alphax == 1, betax == 1,
      alphay == 1, betay == 1, nx >= ny
    ) %>%
    plotT1E()

  a2b5 <- t1e_results %>%
    filter(
      dist == "ordinal", niveau == niveaus[i], alphax == 2, betax == 5,
      alphay == 2, betay == 5, nx >= ny
    ) %>%
    plotT1E()

  a15b15a21b21 <- t1e_results %>%
    filter(
      dist == "ordinal", niveau == niveaus[i], alphax == 5, betax == 5,
      alphay == 1, betay == 1
    ) %>%
    plotT1E()

  a15b15b22b22 <- t1e_results %>%
    filter(
      dist == "ordinal", niveau == niveaus[i], alphax == 5, betax == 5,
      alphay == 2, betay == 2
    ) %>%
    plotT1E()

  # Combine plots and save
  ggpubr::ggarrange(a1b1, a2b5, a15b15a21b21, a15b15b22b22,
                    nrow = 4,
                    common.legend = TRUE, legend = "bottom"
  ) %>%
    ggsave(
      filename = sprintf("./plots/Supplement/section2.3.%s_ordinal_t1e.eps", format(i, scientific = FALSE)),
      device = "eps", width = 16.3, height = 18
    )

  # This is a main paper plot
  if(niveaus[i] == 0.005) {
    ggpubr::ggarrange(a1b1, a2b5, a15b15a21b21, a15b15b22b22,
                      nrow = 4,
                      common.legend = TRUE, legend = "bottom"
    ) %>%
      ggsave(
        filename = "./plots/Main Paper/Figure3.eps",
        device = "eps", width = 16.3, height = 18
      )

  }

}

### Poisson -----------------------------------------------------------------

for (i in seq_along(niveaus)) {
  t1e_results %>%
    filter(dist == "poisson", niveau == niveaus[i], nx >= ny) %>%
    plotT1E() %>%
    ggsave(
      filename = sprintf("./plots/Supplement/section2.4.%s_poisson_t1e.eps", format(i, scientific = FALSE)),
      device = "eps", width = 16.3, height = 7
    )
}

### Exponential -------------------------------------------------------------

for (i in seq_along(niveaus)) {
  t1e_results %>%
    filter(dist == "exponential", niveau == niveaus[i], nx >= ny) %>%
    plotT1E() %>%
    ggsave(
      filename = sprintf("./plots/Supplement/section2.5.%s_exponential_t1e.eps", format(i, scientific = FALSE)),
      device = "eps", width = 16.3, height = 7
    )
}

### Laplace -----------------------------------------------------------------

for (i in seq_along(niveaus)) {
  p_hom <- t1e_results %>%
    filter(dist == "laplace", niveau == niveaus[i], sigmaX == sigmaY, nx >= ny) %>%
    plotT1E()

  p_het <- t1e_results %>%
    filter(dist == "laplace", niveau == niveaus[i], sigmaX != sigmaY) %>%
    plotT1E()

  # Combine plots and save
  ggpubr::ggarrange(p_hom, p_het, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
    ggsave(
      filename = sprintf("./plots/Supplement/section2.6.%s_laplace_t1e.eps", format(i, scientific = FALSE)),
      device = "eps", width = 16.3, height = 11
    )
}


# Power -------------------------------------------------------------------

pwr_cov_results <- read.csv("./results/df_pwr_cov.csv") %>%
  pivot_longer(
    cols = c(pwr_BM, pwr_C2, pwr_permu, cov_BM, cov_C2, cov_permu),
    names_to = c(".value", "test"),
    names_sep = "_"
  ) %>%
  mutate(
    r = round(nx / ny, 2),
    r_label = paste("n1/n2 = ", r),
    test = case_when(
      test == "BM" ~ "BM-Test",
      test == "C2" ~ "C2-Test",
      test == "permu" ~ "Permutation Test"
    ),
    facetN = paste0("(sd1, sd2) = ", "(", sigmaX, ",", sigmaY, ")")
  )

## Normal ------------------------------------------------------------------

plot_normal_equal <- pwr_cov_results %>%
  filter(dist == "normal", nx >= ny, sigmaX == sigmaY) %>%
  plotPOWCOV(type = "power") +
  facet_grid(rows = vars(facetN), cols = vars(r_label))

plot_normal_unequal <- pwr_cov_results %>%
  filter(dist == "normal", sigmaX != sigmaY) %>%
  plotPOWCOV(type = "power") +
  facet_grid(rows = vars(facetN), cols = vars(r_label))

ggpubr::ggarrange(plot_normal_equal, plot_normal_unequal, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
  ggsave(
    filename = "./plots/Supplement/section3.1_normal_pwr.eps",
    device = "eps", width = 16.3, height = 11
  )

# Save also in the plots belonging to the main paper
ggpubr::ggarrange(plot_normal_equal, plot_normal_unequal, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
  ggsave(
    filename = "./plots/Main Paper/Figure4.eps",
    device = "eps", width = 16.3, height = 11
  )

## Ordinal -----------------------------------------------------------------

pwr_cov_results %>%
  filter(dist == "beta") %>%
  plotPOWCOV(type = "power") +
  facet_grid(cols = vars(r_label))

ggsave(
  filename = "./plots/Supplement/section3.2_ordinal_pwr.eps",
    device = "eps", width = 16.3, height = 6
)


## Exponential -------------------------------------------------------------

pwr_cov_results %>%
  filter(dist == "exponential") %>%
  plotPOWCOV(type = "power") +
  facet_grid(cols = vars(r_label))

ggsave(
  filename = "./plots/Supplement/section3.3_exponential_pwr.eps",
    device = "eps", width = 16.3, height = 6
)

# Coverage -------------------------------------------------------------------

## Normal ------------------------------------------------------------------

plot_normal_equal <- pwr_cov_results %>%
  filter(dist == "normal", nx >= ny, sigmaX == sigmaY) %>%
  plotPOWCOV(type = "coverage") +
  facet_grid(rows = vars(facetN), cols = vars(r_label))

plot_normal_unequal <- pwr_cov_results %>%
  filter(dist == "normal", sigmaX != sigmaY) %>%
  plotPOWCOV(type = "coverage") +
  facet_grid(rows = vars(facetN), cols = vars(r_label))

ggpubr::ggarrange(plot_normal_equal, plot_normal_unequal, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
  ggsave(
    filename = "./plots/Supplement/section4.1_normal_coverage.eps",
    device = "eps", width = 16.3, height = 10
  )

# Save also in the plots belonging to the main paper
ggpubr::ggarrange(plot_normal_equal, plot_normal_unequal, nrow = 2, common.legend = TRUE, legend = "bottom") %>%
  ggsave(
    filename = "./plots/Main Paper/Figure5.eps",
    device = "eps", width = 16.3, height = 10
  )


