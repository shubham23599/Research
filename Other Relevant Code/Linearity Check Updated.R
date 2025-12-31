# Define the test statistic functions
WRS = function(u, v) {
  m = length(u); n = length(v); N = m + n
  S = sum(rep(c(0, 1), c(m, n)) * rank(c(u, v)))
  S
}

MW = function(u, v) {
  m = length(u); n = length(v); N = m + n
  S = sum(sapply(1:m, function(i) u[i] > v))
  S
}

VW = function(u, v) {
  m = length(u); n = length(v); N = m + n
  Ind = rep(c(0, 1), c(m, n))
  XN = sum(Ind * qnorm(rank(c(u, v)) / (N + 1)))
  XN
}

Mood = function(u, v) {
  m = length(u); n = length(v); N = m + n
  S = sum(((rank(c(u, v))) - 1/2 * (N + 1))^2 * rep(c(0, 1), c(m, n)))
  S
}

AB = function(u, v) {
  m = length(u); n = length(v); N = m + n
  S = sum(abs((rank(c(u, v))) - 1/2 * (N + 1)) * rep(c(0, 1), c(m, n)))
  S
}


corr_matrix_scatter <- function(simulation = 1000, location = 0, scale = 1) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(officer)
  library(flextable)
  library(GGally)
  library(progress)
  
  m_values <- c(50, 75, 100, 125, 150, 175, 200, 300)
  n_values <- c(5, 10, 15, 20, 25, 30, 50)
  
  # Initialize Word document
  doc <- read_docx()
  
  doc <- doc %>%
    body_add_par("Introduction", style = "heading 1") %>%
    body_add_par(paste(
      "This document shows the correlation matrix, p-value matrix, and scatter plot matrix for different combinations of m and n based on",
      simulation,
      "simulations, assuming the process is in control (IC), i.e., under the alternative hypothesis (Ha), with location (",
      location,
      ") and scale (",
      scale,
      ") while IC location and scale parameters are (0,1). Here, m denotes the reference sample size, and n denotes the testing or monitoring sample size. The statistics calculated are:",
      sep = " "
    ), style = "Normal") %>%
    body_add_par("1. WRS - Wilcoxon Rank Sum Test", style = "Normal") %>%
    body_add_par("2. MW - Mann-Whitney Test", style = "Normal") %>%
    body_add_par("3. VW - Vander Waerden Test", style = "Normal") %>%
    body_add_par("4. AB - Ansari-Bradley Test", style = "Normal") %>%
    body_add_par("5. Mood - Mood's Median Test", style = "Normal") %>%
    body_add_par("For simulation, I used the normal distribution due to the nonparametric setup. In Lepage (1971), it was proven that under the null hypothesis (H0), the WRS and AB statistics are uncorrelated for all m and n.", style = "Normal") %>%
    body_add_par("Reference: Lepage, Y. (1971). A combination of Wilcoxon's and Ansari-Bradley's statistics. Biometrika, 58(1), 213-217.", style = "Normal") %>%
    body_add_par("This document also provides correlation matrices and scatter plot matrices to analyze the relationship between different statistical tests.", style = "Normal") %>%
    body_add_par("The following sections detail the results for each combination of m and n, as well as the location and scale parameters.", style = "Normal")
  
  # Add Index Page
  doc <- doc %>%
    body_add_par("Index", style = "heading 1") %>%
    body_add_par("This document contains results for the following combinations of m, n, location, and scale:", style = "Normal")
  
  # Create an index list for combinations of m, n, location, and scale
  index_list <- list()
  
  for (m in m_values) {
    for (n in n_values) {
      index_list <- c(index_list, paste("Combination of m =", m, ", n =", n, ", location =", location, ", scale =", scale))
    }
  }
  
  # Add the index list to the document
  for (entry in index_list) {
    doc <- doc %>%
      body_add_par(entry, style = "Normal")
  }
  
  # Initialize progress bar
  total_combinations <- length(m_values) * length(n_values)
  pb <- progress_bar$new(
    format = "  Generating [:bar] :percent (:current/:total) ETA: :eta",
    total = total_combinations,
    clear = FALSE,
    width = 60
  )
  
  # Add page break for results
  doc <- doc %>%
    body_add_par("Results", style = "heading 1") %>%
    body_add_par("The following sections present the correlation matrices, p-value matrices, and scatter plot matrices for each combination of m, n, location, and scale.", style = "Normal")
  
  for (m in m_values) {
    for (n in n_values) {
      pb$tick()  # Update progress bar
      
      # Add section header
      doc <- doc %>% body_add_par(paste("Results for m =", m, ", n =", n, ", location =", location, ", scale =", scale), style = "heading 2")
      
      # Generate data
      data <- data.frame(WRS = numeric(simulation), 
                         MW = numeric(simulation), 
                         VW = numeric(simulation), 
                         AB = numeric(simulation), 
                         Mood = numeric(simulation))
      
      for (i in 1:simulation) {
        u <- rnorm(m, mean = 0, sd = 1)
        v <- rnorm(n, mean = location, sd = scale)
        data$WRS[i] <- WRS(u, v)
        data$MW[i] <- MW(u, v)
        data$VW[i] <- VW(u, v)
        data$AB[i] <- AB(u, v)
        data$Mood[i] <- Mood(u, v)
      }
      
      # Calculate correlation and p-value matrices
      correlation_pvalues <- function(data) {
        n <- ncol(data)
        p_mat <- matrix(NA, n, n)
        colnames(p_mat) <- colnames(data)
        rownames(p_mat) <- colnames(data)
        for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
            test <- cor.test(data[[i]], data[[j]])
            p_mat[i, j] <- test$p.value
            p_mat[j, i] <- test$p.value
          }
        }
        diag(p_mat) <- NA
        return(p_mat)
      }
      
      corr_matrix <- cor(data)
      p_matrix <- correlation_pvalues(data)
      
      # Add correlation matrix
      corr_matrix_df <- as.data.frame(round(corr_matrix, 4))
      corr_matrix_df <- cbind(Statistic = rownames(corr_matrix_df), corr_matrix_df)
      ft_corr <- as_flextable(corr_matrix_df) %>%
        colformat_double(j = 2:ncol(corr_matrix_df), digits = 4)
      doc <- doc %>% body_add_par("Correlation Matrix", style = "heading 3") %>%
        body_add_flextable(ft_corr)
      
      # Add p-value matrix
      p_matrix_df <- as.data.frame(round(p_matrix, 4))
      p_matrix_df <- cbind(Statistic = rownames(p_matrix_df), p_matrix_df)
      ft_p <- as_flextable(p_matrix_df) %>%
        colformat_double(j = 2:ncol(p_matrix_df), digits = 4)
      doc <- doc %>% body_add_par("P-Value Matrix", style = "heading 3") %>%
        body_add_flextable(ft_p)
      
      # Scatterplot matrix
      scatter_matrix <- ggpairs(data, title = paste("Scatterplot Matrix for m =", m, ", n =", n, ", location =", location, ", scale =", scale))
      plot_file <- tempfile(fileext = ".png")
      ggsave(plot_file, scatter_matrix, width = 12, height = 12)
      doc <- doc %>% body_add_par("Scatter Plot Matrix:", style = "heading 3") %>%
        body_add_img(src = plot_file, width = 6.69, height = 8)
      unlink(plot_file) # Remove temp file
    }
  }
  
  print(doc, target = "Test_Statistics_Summary_IC.docx")
  message("Word document created: Test_Statistics_Summary.docx")
}

# Example call
corr_matrix_scatter(simulation = 50000, location = 0, scale = 1)





