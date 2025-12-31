library(ggplot2)

# Sample data
data <- data.frame(
  Sample = 1:10,
  Statistic = c(5, 6, 5.5, 7, 6.8, 6.5, 6.1, 5.7, 7.2, 6)
)

# Control limits
CL <- mean(data$Statistic)  # Center Line
UCL <- CL + 1  # Upper Control Limit
LCL <- CL - 1  # Lower Control Limit

# Plot
ggplot(data, aes(x = Sample, y = Statistic)) +
  geom_line(aes(x = Sample, y = Statistic), color = "black", size = 1) +
  geom_point(color = "black", size = 2) +
  geom_segment(aes(x = min(Sample), xend = max(Sample), y = CL, yend = CL), color = "black", size = 1) +
  geom_segment(aes(x = min(Sample), xend = max(Sample), y = UCL, yend = UCL), color = "black", size = 1) +
  geom_segment(aes(x = min(Sample), xend = max(Sample), y = LCL, yend = LCL), color = "black", size = 1) +
  annotate("text", x = max(data$Sample) + 1.5, y = UCL, label = "UCL", hjust = 1, vjust = -0.5, color = "black", fontface = "bold") +
  annotate("text", x = max(data$Sample) + 1.5, y = CL, label = "CL", hjust = 1, vjust = -0.5, color = "black", fontface = "bold") +
  annotate("text", x = max(data$Sample) + 1.5, y = LCL, label = "LCL", hjust = 1, vjust = 1.5, color = "black", fontface = "bold") +
  labs(x = "Sample number or time", y = "Control Statistic") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
