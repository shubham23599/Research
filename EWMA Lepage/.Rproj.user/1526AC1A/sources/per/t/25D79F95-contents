library(ggplot2)
library(dplyr)
library(tidyr)

# Data
x  <- seq(-5, 5, length.out = 2000)
x_pos <- seq(0.001, 5, length.out = 1500)

df_common <- data.frame(
  x = x,
  Normal   = dnorm(x, 0, 1),
  Logistic = dlogis(x, 0, 1),
  Laplace  = 0.5 * exp(-abs(x))
)

df_lognorm <- data.frame(
  x = x_pos,
  Lognormal = dlnorm(x_pos, 0, 1)
)

df_all <- bind_rows(
  pivot_longer(df_common, -x, names_to="Distribution", values_to="Density"),
  pivot_longer(df_lognorm, -x, names_to="Distribution", values_to="Density")
)

# Colors
cols <- c(
  "Laplace"   = "#B22222",
  "Logistic"  = "#1E90FF",
  "Lognormal" = "#228B22",
  "Normal"    = "#4B0082"
)

# Labels with slight gap
labs <- c(
  "Normal"    = expression("  " * N(theta, delta)),
  "Lognormal" = expression("  " * LN(theta, delta)),
  "Laplace"   = expression("  " * L(theta, delta)),
  "Logistic"  = expression("  " * Log(theta, delta))
)

# ⬇️ CLEAN MINIMAL GRAPH WITH AXES + NO BACKGROUND
ggplot(df_all, aes(x, Density, color = Distribution)) +
  geom_line(size = 1) +
  scale_color_manual(
    values = cols,
    labels = labs,
    limits = c("Normal", "Lognormal", "Laplace", "Logistic")  # ORDER SET HERE
  ) +
  labs(x = "v", y = "Probability density function", color = "") +
  theme_minimal(base_size = 16) +
  theme(
    # --- Legend inside plot ---
    legend.position = c(0.95, 0.90),
    legend.justification = c(1, 1),
    legend.text = element_text(size = 16),
    legend.key.width = unit(1.2, "cm"),
    
    # --- Axes lines ---
    axis.line = element_line(color = "black", linewidth = 0.7),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    # --- Remove background ---
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
