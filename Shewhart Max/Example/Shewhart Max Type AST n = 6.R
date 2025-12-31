library(DescTools)
library(ggplot2)
source("Distance_Statistics.R")

# Read the dataset
df=read.csv("Example//Example Data AST.csv")
head(df)

# Split Phase I and Phase II data
Phase_1 <- df[df$Phase == "I", paste0("X", 1:5)]
Phase_2 <- df[df$Phase == "II", paste0("X", 1:5)]

k <- 6  # subgroup size

## Convert Phase_2 to row-wise vector (same as before)
vec <- as.vector(t(as.matrix(Phase_2)))

## Number of samples
nS <- length(vec) - k + 1

## Preallocate matrix (fast)
Phase2_samples <- matrix(NA_real_, nrow = nS, ncol = k)

## Fill using vectorized indexing
for (j in 1:k) {
  Phase2_samples[, j] <- vec[j:(j + nS - 1)]
}

## Convert to data frame
Phase2_samples <- as.data.frame(Phase2_samples)
colnames(Phase2_samples) <- paste0("X", 1:k)

## Optional: add sample index
Phase2_samples$Sample <- seq_len(nS)

Phase2_samples

# Flatten Phase I data into a single vector 'u'
u <- unlist(lapply(paste0("X", 1:5), function(x) Phase_1[[x]]))
print(paste0("m = ", length(u)))

# Wilcoxon–Rank–Sum and Ansari–Bradley (WRSAB) statistic
WRSAB_MAX <- function(u, v) {
  m <- length(u); n <- length(v); N <- m + n
  m1 <- 0.5 * (n * (N + 1))
  s1 <- (1 / 12) * (m * n * (N + 1))
  
  if (N %% 2 == 0) {
    m2 <- (n * N) / 4
    s2 <- (1 / 48) * (m * n * (N^2 - 4) / (N - 1))
  } else {
    m2 <- (n * (N^2 - 1)) / (4 * N)
    s2 <- (1 / 48) * (m * n * (N + 1) * (N^2 + 3)) / N^2
  }
  
  S1 <- abs((sum(rep(c(0, 1), c(m, n)) * rank(c(u, v))) - m1) / sqrt(s1))
  S2 <- abs((sum(abs((rank(c(u, v))) - 0.5 * (N + 1)) * rep(c(0, 1), c(m, n))) - m2) / sqrt(s2))
  S <- max(S1, S2)
  
  return(c(S1, S2, S))
}

# Van der Waerden and Mood Test (VWM) statistic
VWM_MAX <- function(u, v) {
  m <- length(u); n <- length(v); N <- m + n
  Ind <- rep(c(0, 1), c(m, n))
  
  m2 <- (n * (N^2 - 1)) / 12
  s2 <- (1 / 180) * (m * n * (N + 1) * (N^2 - 4))
  
  norm_scores <- qnorm(rank(c(u, v)) / (N + 1))
  XN <- sum(Ind * norm_scores)
  
  s1 <- (m * n) / (N * (N - 1)) * sum(norm_scores^2)
  
  S1 <- abs(XN / sqrt(s1))
  S2 <- abs((sum((rank(c(u, v)) - 0.5 * (N + 1))^2 * Ind) - m2) / sqrt(s2))
  S <- max(S1, S2)
  
  return(c(S1, S2, S))
}

# Generate statistics for Phase II using Phase I as baseline
generate_plot_df <- function(phase_2, u) {
  plot_df <- data.frame(
    Sample = integer(),
    WRS_Statistic = numeric(),
    AB_Statistic = numeric(),
    WRSAB_Max = numeric(),
    VW_Statistic = numeric(),
    Mood_Statistic = numeric(),
    VWM_Max = numeric(),
    Lepage = numeric(),
    Cucconi = numeric(),
    SLVM = numeric(),
    Wilcoxon_p = numeric(),
    VanWaerden_p = numeric(),
    AnsariBradley_p = numeric(),
    Mood_p = numeric()
  )
  
  for (t in 1:nrow(phase_2)) {
    v <- as.vector(unlist(phase_2[t, ]))
    test1 <- WRSAB_MAX(u, v)
    test2 <- VWM_MAX(u, v)
    lepage <- WRSAB(u, v)
    cucconi <- SC(u, v)
    SLVM <- VWM(u,v)
    plot_df[t, ] <- c(
      t,
      test1[1],
      test1[2],
      test1[3],
      test2[1],
      test2[2],
      test2[3],
      lepage,cucconi,SLVM,
      wilcox.test(u, v)$p.value,
      VanWaerdenTest(list(v, u))$p.value,
      ansari.test(u, v)$p.value,
      mood.test(u, v)$p.value
    )
  }
  
  return(plot_df)
}

# Run the function
plot_df <- generate_plot_df(Phase2_samples, u)
print(plot_df)



# Define the plot function (UDF)
create_plot <- function(data, chart_name,lower_limit, upper_limit,stat_name,adj=0.2) {
  dataset_name <- ""
  shifts <- 1:nrow(data)
  ARL <- data[[chart_name]]
  chart <- rep("Statistic", length(shifts))
  plot_data <- data.frame(chart, shifts, ARL)
  colnames(plot_data) <- c("Charts", "Shifts", "ARL")
  
  # Generate the plot
  plot <- ggplot(plot_data, aes(x = Shifts, y = ARL)) +
    geom_line(color = "black", size = 1.1) +
    geom_point(color = "black", size = 2.2) +
    theme_bw() +
    geom_hline(yintercept = upper_limit, size = 0.9) +
    geom_hline(yintercept = lower_limit) +
    theme(
      legend.position = c(0.2, 0.9),
      legend.title = element_blank(),
      axis.title.y = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      text = element_text(size = 20),
      axis.text = element_text(size = 10),
      axis.ticks.length = unit(0.45, "cm")
    ) +
    xlab("Sample Number") +
    ylab(stat_name) +
    scale_x_continuous(breaks = seq(0, nrow(data), by = 2)) +
    #scale_y_continuous(limits = c(floor(min(ARL)), ceiling(max(ARL))),
    #                   breaks = seq(floor(min(ARL)), ceiling(max(ARL)), by = 0.5)) +
    geom_text(aes(x = nrow(data) + 1, label = "UCL", y = upper_limit - adj), size = 5)
  return(plot)
}

plot_df = plot_df[1:25,]
# Generate and display the plot
ARL_graph <- create_plot(plot_df,"WRSAB_Max",0,2.975,bquote(M["W,A"]),0.2)
print(ARL_graph)

ARL_graph <- create_plot(plot_df,"VWM_Max",0,3.070,bquote(M["V,M"]),0.2)
print(ARL_graph)

# Generate and display the plot
ARL_graph <- create_plot(plot_df,"Lepage",0,11.25,bquote(Lepage),0.5)
print(ARL_graph)

ARL_graph <- create_plot(plot_df,"Cucconi",0,6.05,bquote(Cucconi),0.35)
print(ARL_graph)

ARL_graph <- create_plot(plot_df,"SLVM",0,13.200,bquote(L["V,M"]),0.75)
print(ARL_graph)


# Save the final dataset
output_file <- "Shewhart Max Example AST Statistics Updated n=6 .csv"
write.csv(plot_df, output_file)
