# Import functions
source("functions/util.R")
source("functions/getTestResults.R")
source("functions/getStudentizedPermutationTest.R")

# Data example in Table 1 -------------------------------------------------

group1 <- c(1956, 3828, 2051, 3721, 3233, 2000, 4000, 4428, 2603, 2370)
group2 <- c(820, 3364, 1957, 1851, 2984, 744, 2044)

getTestResults(group2, group1, seed=6082025)

# Pain Scores example in Section 2 and 6 ----------------------------------

# Pain Scores
intervention <- rep(1:5, c(16, 5, 0, 1, 0))
control <- rep(1:5, c(4, 1, 5, 7, 2))

# Combine data into a data frame
data <- data.frame(
  PainScore = c(intervention, control),
  Group = rep(
    c("Specific Suction (n=22)", "Control (n=19)"),
    c(length(intervention), length(control))
  )
)

# Equal distance between the points in the boxplot
data$position <- ave(data$PainScore, data$Group, data$PainScore,
  FUN = function(x) seq(-0.3, 0.3, length.out = length(x))
)
data$position[c(22, 27)] <- 0

# Create boxplot
ggplot(data, aes(x = Group, y = PainScore)) +
  stat_boxplot(geom = "errorbar", width = 0.6) +
  geom_boxplot(outlier.colour = "white", fill = "white", colour = "black") +
  geom_point(aes(x = as.numeric(factor(Group)) + position, y = PainScore),
    size = 4, shape = 16
  ) +
  ylim(c(0.5, 5.5)) +
  labs(x = "Treatment Group", y = "Pain Score") +
  labs(
    x = "Treatment Group",
    y = "Pain Score"
  ) +
  theme +
  theme(legend.position = "none")

# Save boxplot
ggsave(filename = "./plots/Main Paper/Data_Table1.eps", device = "eps", width = 8, height = 7)

# Conduct test
getTestResults(
  X = intervention, Y = control, alpha = 0.05, permutation = T,
  nperm = 1e4, seed = 05022025
)

