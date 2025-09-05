library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

# Create data frames from the tables
create_data <- function(rho) {
  if(rho == 0.1) {
    data <- data.frame(
      pi1 = rep(c(rep(0.25, 9), 2)),
                Delta = rep(c(rep(0.10, 3), rep(0.15, 3), 4), 2),
                ni = rep(c(25, 50, 100), 12),
                Wald_ECP = c(0.940, 0.942, 0.948, 0.939, 0.935, 0.943, 0.932, 0.939, 0.945,
                             0.935, 0.945, 0.941, 0.943, 0.945, 0.942, 0.922, 0.948, 0.944),
                Wald_MIW = c(0.365, 0.261, 0.186, 0.371, 0.266, 0.189, 0.374, 0.267, 0.191,
                             0.397, 0.284, 0.202, 0.393, 0.281, 0.200, 0.385, 0.276, 0.196),
                Likelihood_ECP = c(0.913, 0.898, 0.889, 0.902, 0.885, 0.889, 0.889, 0.887, 0.665,
                                   0.901, 0.905, 0.907, 0.905, 0.913, 0.949, 0.907, 0.824, 0.442),
                Likelihood_MIW = c(0.255, 0.192, 0.139, 0.271, 0.195, 0.126, 0.273, 0.182, 0.106,
                                   0.250, 0.185, 0.136, 0.259, 0.189, 0.128, 0.262, 0.180, 0.114),
                Score_ECP = c(0.830, 0.813, 0.796, 0.830, 0.782, 0.736, 0.812, 0.751, 0.657,
                              0.777, 0.766, 0.744, 0.764, 0.738, 0.949, 0.726, 0.698, 0.781),
                Score_MIW = c(0.282, 0.221, 0.170, 0.332, 0.242, 0.180, 0.332, 0.250, 0.180,
                              0.314, 0.251, 0.195, 0.374, 0.278, 0.209, 0.363, 0.288, 0.213),
                MOVER1_ECP = c(0.947, 0.946, 0.951, 0.946, 0.945, 0.945, 0.942, 0.946, 0.945,
                               0.940, 0.949, 0.943, 0.947, 0.951, 0.946, 0.940, 0.951, 0.946),
                MOVER1_MIW = c(0.359, 0.259, 0.185, 0.363, 0.263, 0.188, 0.365, 0.264, 0.189,
                               0.384, 0.279, 0.201, 0.380, 0.276, 0.198, 0.374, 0.272, 0.195),
                MOVER2_ECP = c(0.948, 0.947, 0.953, 0.948, 0.946, 0.945, 0.944, 0.947, 0.946,
                               0.941, 0.950, 0.943, 0.948, 0.952, 0.946, 0.941, 0.951, 0.946),
                MOVER2_MIW = c(0.361, 0.260, 0.185, 0.365, 0.263, 0.188, 0.367, 0.265, 0.190,
                               0.384, 0.279, 0.201, 0.380, 0.276, 0.198, 0.375, 0.272, 0.195),
                SA_ECP = c(0.949, 0.945, 0.950, 0.948, 0.953, 0.954, 0.949, 0.953, 0.945,
                           0.945, 0.951, 0.954, 0.950, 0.948, 0.951, 0.938, 0.944, 0.949),
                SA_MIW = c(0.359, 0.259, 0.185, 0.364, 0.263, 0.188, 0.368, 0.265, 0.189,
                           0.390, 0.281, 0.201, 0.385, 0.278, 0.198, 0.378, 0.272, 0.195)
      )
  } else if(rho == 0.5) {
    data <- data.frame(
      pi1 = rep(c(rep(0.25, 9), 2)),
                Delta = rep(c(rep(0.10, 3), rep(0.15, 3), rep(0.20, 3)), 2),
                            ni = rep(c(25, 50, 100), 12),
                            Wald_ECP = c(0.941, 0.956, 0.944, 0.937, 0.948, 0.952, 0.934, 0.939, 0.940,
                                         0.939, 0.941, 0.946, 0.940, 0.948, 0.941, 0.938, 0.939, 0.947),
                            Wald_MIW = c(0.423, 0.303, 0.216, 0.431, 0.309, 0.219, 0.435, 0.312, 0.222,
                                         0.465, 0.332, 0.236, 0.458, 0.328, 0.233, 0.450, 0.321, 0.228),
                            Likelihood_ECP = c(0.919, 0.929, 0.904, 0.908, 0.907, 0.878, 0.895, 0.883, 0.665,
                                               0.894, 0.893, 0.903, 0.905, 0.907, 0.899, 0.906, 0.897, 0.741),
                            Likelihood_MIW = c(0.295, 0.224, 0.168, 0.318, 0.239, 0.166, 0.329, 0.236, 0.145,
                                               0.292, 0.216, 0.160, 0.301, 0.221, 0.151, 0.303, 0.211, 0.140),
                            Score_ECP = c(0.846, 0.840, 0.807, 0.833, 0.827, 0.785, 0.825, 0.803, 0.657,
                                          0.774, 0.766, 0.749, 0.762, 0.737, 0.698, 0.745, 0.715, 0.782),
                            Score_MIW = c(0.310, 0.242, 0.188, 0.341, 0.268, 0.203, 0.364, 0.280, 0.205,
                                          0.347, 0.274, 0.215, 0.380, 0.303, 0.233, 0.406, 0.316, 0.237),
                            MOVER1_ECP = c(0.947, 0.956, 0.944, 0.947, 0.952, 0.953, 0.943, 0.942, 0.943,
                                           0.947, 0.947, 0.948, 0.945, 0.949, 0.944, 0.944, 0.945, 0.949),
                            MOVER1_MIW = c(0.418, 0.303, 0.216, 0.424, 0.307, 0.219, 0.426, 0.309, 0.221,
                                           0.450, 0.326, 0.234, 0.444, 0.323, 0.231, 0.438, 0.317, 0.227),
                            MOVER2_ECP = c(0.949, 0.957, 0.945, 0.949, 0.953, 0.953, 0.945, 0.943, 0.943,
                                           0.949, 0.947, 0.948, 0.946, 0.950, 0.944, 0.945, 0.946, 0.950),
                            MOVER2_MIW = c(0.421, 0.304, 0.217, 0.426, 0.308, 0.220, 0.429, 0.310, 0.221,
                                           0.450, 0.326, 0.234, 0.445, 0.323, 0.231, 0.439, 0.317, 0.227),
                            SA_ECP = c(0.944, 0.945, 0.952, 0.944, 0.941, 0.944, 0.937, 0.942, 0.947,
                                       0.936, 0.951, 0.944, 0.941, 0.944, 0.948, 0.950, 0.941, 0.945),
                            SA_MIW = c(0.414, 0.301, 0.215, 0.419, 0.304, 0.218, 0.432, 0.307, 0.220,
                                       0.448, 0.326, 0.234, 0.443, 0.322, 0.231, 0.435, 0.316, 0.226)
                )
  } else if(rho == 0.9) {
    data <- data.frame(
      pi1 = rep(c(rep(0.25, 9), 2)),
                Delta = rep(c(rep(0.10, 3), rep(0.15, 3), rep(0.20, 3)), 2),
                ni = rep(c(25, 50, 100), 12),
                Wald_ECP = c(0.937, 0.945, 0.949, 0.933, 0.939, 0.951, 0.943, 0.939, 0.951,
                             0.943, 0.937, 0.947, 0.938, 0.945, 0.949, 0.943, 0.942, 0.948),
                Wald_MIW = c(0.479, 0.343, 0.244, 0.487, 0.348, 0.248, 0.492, 0.352, 0.250,
                             0.523, 0.374, 0.266, 0.516, 0.369, 0.263, 0.506, 0.362, 0.257),
                Likelihood_ECP = c(0.927, 0.931, 0.927, 0.917, 0.916, 0.905, 0.917, 0.898, 0.874,
                                   0.898, 0.880, 0.879, 0.884, 0.889, 0.901, 0.882, 0.896, 0.890),
                Likelihood_MIW = c(0.330, 0.252, 0.193, 0.360, 0.276, 0.205, 0.377, 0.284, 0.193,
                                   0.327, 0.242, 0.179, 0.333, 0.247, 0.176, 0.332, 0.239, 0.160),
                Score_ECP = c(0.861, 0.858, 0.852, 0.848, 0.840, 0.818, 0.846, 0.820, 0.777,
                              0.760, 0.739, 0.720, 0.734, 0.719, 0.697, 0.733, 0.704, 0.739),
                Score_MIW = c(0.337, 0.264, 0.207, 0.373, 0.293, 0.228, 0.396, 0.312, 0.233,
                              0.370, 0.284, 0.216, 0.400, 0.310, 0.203, 0.425, 0.292, 0.124),
                MOVER1_ECP = c(0.944, 0.947, 0.951, 0.946, 0.945, 0.951, 0.952, 0.942, 0.952,
                               0.953, 0.940, 0.948, 0.948, 0.950, 0.950, 0.953, 0.947, 0.951),
                MOVER1_MIW = c(0.472, 0.341, 0.243, 0.477, 0.345, 0.247, 0.481, 0.348, 0.249,
                               0.506, 0.368, 0.264, 0.500, 0.363, 0.260, 0.492, 0.357, 0.256),
                MOVER2_ECP = c(0.946, 0.948, 0.952, 0.948, 0.946, 0.952, 0.954, 0.944, 0.953,
                               0.954, 0.942, 0.948, 0.950, 0.950, 0.960, 0.955, 0.948, 0.952),
                MOVER2_MIW = c(0.475, 0.342, 0.244, 0.480, 0.346, 0.247, 0.483, 0.349, 0.249,
                               0.506, 0.368, 0.264, 0.501, 0.363, 0.260, 0.494, 0.358, 0.256),
                SA_ECP = c(0.940, 0.942, 0.952, 0.941, 0.951, 0.943, 0.942, 0.949, 0.951,
                           0.945, 0.942, 0.947, 0.942, 0.938, 0.944, 0.933, 0.942, 0.955),
                SA_MIW = c(0.458, 0.335, 0.241, 0.464, 0.340, 0.245, 0.469, 0.343, 0.247,
                           0.497, 0.364, 0.262, 0.491, 0.359, 0.259, 0.485, 0.353, 0.254)
      )
  }
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(gridExtra)
  
  # Create data frame for ρ = 0.9
  data_rho9 <- data.frame(
    pi1 = rep(c(rep(0.25, 9), rep(0.50, 9))),
    Delta = rep(c(rep(0.10, 3), rep(0.15, 3), rep(0.20, 3)), 2),
    ni = rep(c(25, 50, 100), 6),
    Wald_ECP = c(0.937, 0.945, 0.949, 0.933, 0.939, 0.951, 0.943, 0.939, 0.951,
                 0.943, 0.937, 0.947, 0.938, 0.945, 0.949, 0.943, 0.942, 0.948),
    Wald_MIW = c(0.479, 0.343, 0.244, 0.487, 0.348, 0.248, 0.492, 0.352, 0.250,
                 0.523, 0.374, 0.266, 0.516, 0.369, 0.263, 0.506, 0.362, 0.257),
    Likelihood_ECP = c(0.927, 0.931, 0.927, 0.917, 0.916, 0.905, 0.917, 0.898, 0.874,
                       0.898, 0.880, 0.879, 0.884, 0.889, 0.901, 0.882, 0.896, 0.890),
    Likelihood_MIW = c(0.330, 0.252, 0.193, 0.360, 0.276, 0.205, 0.377, 0.284, 0.193,
                       0.327, 0.242, 0.179, 0.333, 0.247, 0.176, 0.332, 0.239, 0.160),
    Score_ECP = c(0.861, 0.858, 0.852, 0.848, 0.840, 0.818, 0.846, 0.820, 0.777,
                  0.760, 0.739, 0.720, 0.734, 0.719, 0.697, 0.733, 0.704, 0.739),
    Score_MIW = c(0.337, 0.264, 0.207, 0.373, 0.293, 0.228, 0.396, 0.312, 0.233,
                  0.370, 0.284, 0.216, 0.400, 0.310, 0.203, 0.425, 0.292, 0.124),
    MOVER1_ECP = c(0.944, 0.947, 0.951, 0.946, 0.945, 0.951, 0.952, 0.942, 0.952,
                   0.953, 0.940, 0.948, 0.948, 0.950, 0.950, 0.953, 0.947, 0.951),
    MOVER1_MIW = c(0.472, 0.341, 0.243, 0.477, 0.345, 0.247, 0.481, 0.348, 0.249,
                   0.506, 0.368, 0.264, 0.500, 0.363, 0.260, 0.492, 0.357, 0.256),
    MOVER2_ECP = c(0.946, 0.948, 0.952, 0.948, 0.946, 0.952, 0.954, 0.944, 0.953,
                   0.954, 0.942, 0.948, 0.950, 0.950, 0.960, 0.955, 0.948, 0.952),
    MOVER2_MIW = c(0.475, 0.342, 0.244, 0.480, 0.346, 0.247, 0.483, 0.349, 0.249,
                   0.506, 0.368, 0.264, 0.501, 0.363, 0.260, 0.494, 0.358, 0.256),
    SA_ECP = c(0.940, 0.942, 0.952, 0.941, 0.951, 0.943, 0.942, 0.949, 0.951,
               0.945, 0.942, 0.947, 0.942, 0.938, 0.944, 0.933, 0.942, 0.955),
    SA_MIW = c(0.458, 0.335, 0.241, 0.464, 0.340, 0.245, 0.469, 0.343, 0.247,
               0.497, 0.364, 0.262, 0.491, 0.359, 0.259, 0.485, 0.353, 0.254)
  )
  
  # Reshape data to long format
  data_long <- data_rho9 %>%
    pivot_longer(cols = -c(pi1, Delta, ni),
                 names_to = c("Method", "Metric"),
                 names_sep = "_",
                 values_to = "Value") %>%
    pivot_wider(names_from = "Metric", values_from = "Value")
  
  # Convert to factors for better plotting
  data_long <- data_long %>%
    mutate(
      pi1 = factor(pi1),
      Delta = factor(Delta),
      ni = factor(ni, levels = c(25, 50, 100), ordered = TRUE),
      Method = factor(Method, levels = c("Wald", "Likelihood", "Score", "MOVER1", "MOVER2", "SA"))
    )
  
  ## Create Plots for ρ = 0.9
  
  # Plot 1: ECP by Method and Sample Size (faceted by pi1 and Delta)
  p1 <- ggplot(data_long, aes(x = ni, y = ECP, color = Method, group = Method)) +
    geom_point(size = 3) +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    facet_grid(pi1 ~ Delta, labeller = label_both) +
    labs(title = "Empirical Coverage Probability (ρ = 0.9)",
         subtitle = "Red dashed line indicates ideal 95% coverage",
         x = "Sample Size", y = "ECP") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  # Plot 2: MIW by Method and Sample Size (faceted by pi1 and Delta)
  p2 <- ggplot(data_long, aes(x = ni, y = MIW, color = Method, group = Method)) +
    geom_point(size = 3) +
    geom_line() +
    facet_grid(pi1 ~ Delta, labeller = label_both) +
    labs(title = "Mean Interval Width (ρ = 0.9)",
         x = "Sample Size", y = "MIW") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  # Plot 3: ECP vs MIW scatter plot (colored by Method, shaped by sample size)
  p3 <- ggplot(data_long, aes(x = MIW, y = ECP, color = Method, shape = ni)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    facet_grid(pi1 ~ Delta, labeller = label_both) +
    labs(title = "ECP vs MIW Trade-off (ρ = 0.9)",
         subtitle = "Ideal methods would be in top-left corner (high coverage, low width)",
         x = "Mean Interval Width", y = "Empirical Coverage Probability") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Plot 4: Boxplot of ECP by Method
  p4 <- ggplot(data_long, aes(x = Method, y = ECP, fill = Method)) +
    geom_boxplot() +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    labs(title = "Distribution of ECP by Method (ρ = 0.9)",
         x = "Method", y = "ECP") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 5: Heatmap of average ECP by Method and Sample Size
  p5_data <- data_long %>%
    group_by(Method, ni) %>%
    summarize(avg_ECP = mean(ECP), .groups = "drop")
  
  p5 <- ggplot(p5_data, aes(x = ni, y = Method, fill = avg_ECP)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", 
                         midpoint = 0.95, limits = c(0.8, 1)) +
    geom_text(aes(label = round(avg_ECP, 3)), color = "black", size = 3) +
    labs(title = "Average ECP by Method and Sample Size (ρ = 0.9)",
         x = "Sample Size", y = "Method", fill = "Average ECP") +
    theme_bw()
  
  # Plot 6: Bar plot comparing average ECP and MIW
  p6_data <- data_long %>%
    group_by(Method) %>%
    summarize(avg_ECP = mean(ECP), avg_MIW = mean(MIW), .groups = "drop") %>%
    pivot_longer(cols = c(avg_ECP, avg_MIW), names_to = "Metric", values_to = "Value")
  
  p6 <- ggplot(p6_data, aes(x = Method, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(data = data.frame(Metric = "avg_ECP", yint = 0.95), 
               aes(yintercept = yint), linetype = "dashed", color = "red") +
    labs(title = "Average Performance Metrics by Method (ρ = 0.9)",
         x = "Method", y = "Value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Arrange plots in a grid
  grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)