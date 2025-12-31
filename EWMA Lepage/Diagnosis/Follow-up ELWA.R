library(openxlsx)
library(dplyr)

run_lepage_simulation_once <- function(
    mu_shift, sigma_shift,
    Simulations = 10000,
    lambda = 0.1, mu_lep = 2, sd_lep = 2,
    m = 100, n = 5, tau = 100,
    mu = 0, sigma = 1,
    K = 2.920,
    alpha = 1/500
) {
  set.seed(123)
  
  N <- m + n
  Indicator <- rep(c(0,1), c(m,n))
  
  # Control limit
  H <- round(mu_lep + K * sd_lep * sqrt(lambda/(2 - lambda)), 3)
  alpha_prime <- 1 - sqrt(1 - alpha)
  
  # --- Mean/variance formulae ---
  m1 <- n * (N + 1) / 2
  s1 <- m * n * (N + 1) / 12
  
  if (N %% 2 == 0) {
    m2 <- (n * N) / 4
    s2 <- (1 / 48) * (m * n * (N^2 - 4) / (N - 1))
  } else {
    m2 <- (n * (N^2 - 1)) / (4 * N)
    s2 <- (1 / 48) * (m * n * (N + 1) * (N^2 + 3)) / N^2
  }
  
  sqrt_s1 <- sqrt(s1)
  sqrt_s2 <- sqrt(s2)
  
  location <- scale <- both <- false <- 0
  results <- vector("list", Simulations)
  
  # =====================================================
  #                  MAIN LOOP
  # =====================================================
  for (j in seq_len(Simulations)) {
    u <- rnorm(m, mu, sigma)
    Z <- mu_lep
    i <- 0
    test_sample <- list()
    
    # -------- warm-up --------
    while (i < tau) {
      v <- rnorm(n, mu, sigma)
      test_sample[[i + 1]] <- v
      
      R <- rank(c(u, v))
      WRS <- ((sum(R * Indicator) - m1) / sqrt_s1)^2
      AB <- ((sum(abs(R - 0.5 * (N + 1)) * Indicator) - m2) / sqrt_s2)^2
      S <- WRS + AB
      
      Z <- lambda * S + (1 - lambda) * Z
      Z <- max(mu_lep, Z)
      
      if (Z > H) {
        i <- 0
        Z <- mu_lep
      } else {
        i <- i + 1
      }
    }
    
    # -------- OOC --------
    RL <- 0
    while (Z <= H) {
      v <- rnorm(n, mu_shift, sigma_shift)
      test_sample[[length(test_sample) + 1]] <- v
      
      R <- rank(c(u, v))
      WRS <- ((sum(R * Indicator) - m1) / sqrt_s1)^2
      AB <- ((sum(abs(R - 0.5 * (N + 1)) * Indicator) - m2) / sqrt_s2)^2
      S <- WRS + AB
      
      Z <- lambda * S + (1 - lambda) * Z
      Z <- max(mu_lep, Z)
      
      RL <- RL + 1
    }
    
    # -------- Classification --------
    last_v <- test_sample[[length(test_sample)]]
    
    p1 <- wilcox.test(u, last_v)$p.value
    p2 <- ansari.test(u, last_v)$p.value
    
    classify <- function(p1, p2) {
      if (p1 <= alpha_prime & p2 > alpha_prime) return("location")
      if (p1 > alpha_prime & p2 <= alpha_prime) return("scale")
      if (p1 <= alpha_prime & p2 <= alpha_prime) return("both")
      return("none")
    }
    
    class <- classify(p1, p2)
    
    if (class == "none") {
      if (length(test_sample) >= 2) {
        combined <- c(test_sample[[length(test_sample)]],
                      test_sample[[length(test_sample)-1]])
        p1b <- wilcox.test(u, combined)$p.value
        p2b <- ansari.test(u, combined)$p.value
        class <- classify(p1b, p2b)
        p1 <- p1b; p2 <- p2b
      }
    }
    
    if (class == "location") location <- location + 1 else
      if (class == "scale") scale <- scale + 1 else
        if (class == "both") both <- both + 1 else
          false <- false + 1
    
    # Save replication results
    results[[j]] <- data.frame(
      Replication = j,
      OOC_DETECTED_POINT = RL,
      p_wilcox = p1,
      p_ansari = p2,
      Classification = class
    )
  }
  
  # ---------- Summary ----------
  results_df <- bind_rows(results)
  
  ARL <- mean(results_df$OOC_DETECTED_POINT)
  ARL_median <- median(results_df$OOC_DETECTED_POINT)
  ARL_sd <- sd(results_df$OOC_DETECTED_POINT)
  
  actual_shift <- if (mu_shift != mu & sigma_shift == sigma) {
    "Location"
  } else if (mu_shift == mu & sigma_shift != sigma) {
    "Scale"
  } else if (mu_shift != mu & sigma_shift != sigma) {
    "Both"
  } else {
    "None"
  }
  
  final_df <- data.frame(
    IC_LOCATION = mu,
    IC_SCALE = sigma,
    SHIFT_LOCATION = mu_shift,
    SHIFT_SCALE = sigma_shift,
    Actual_Shift = actual_shift,
    ARL_mean = ARL,
    ARL_median = ARL_median,
    ARL_sd = ARL_sd,
    location_count = location,
    scale_count = scale,
    simultaneous_count = both,
    false_count = false,
    location_rate = location / Simulations,
    scale_rate = scale / Simulations,
    simultaneous_rate = both / Simulations,
    false_rate = false / Simulations
  )
  
  return(list(results_df = results_df, final_df = final_df))
}


# ================================================================
#                BATCH FUNCTION FOR MULTIPLE COMBINATIONS
# ================================================================

run_lepage_simulation_batch <- function(
    theta_values, delta_values,
    results_excel = "All_Results.xlsx",
    summary_excel = "Summary_Final.xlsx",
    IC_mu = 0,        # IC location mean
    IC_sigma = 1      # IC scale (std dev)
) {
  
  wb_results <- createWorkbook()
  wb_summary <- createWorkbook()
  
  summary_all <- data.frame()
  
  for (θ in theta_values) {
    for (δ in delta_values) {
      
      # ==========================================================
      #                SKIP IC CASE (θ = µ AND δ = σ)
      # ==========================================================
      if (θ == IC_mu && δ == IC_sigma) {
        cat("Skipping IC case θ =", θ, "δ =", δ, "\n")
        next   # Skip this combination
      }
      
      cat("Running θ =", θ, ", δ =", δ, "...\n")
      
      out <- run_lepage_simulation_once(
        mu_shift = θ,
        sigma_shift = δ
      )
      
      results_df <- out$results_df
      final_df <- out$final_df
      
      sheet_name <- paste0("theta_", θ, "_delta_", δ)
      
      addWorksheet(wb_results, sheet_name)
      writeData(wb_results, sheet_name, results_df)
      
      summary_all <- dplyr::bind_rows(summary_all, final_df)
    }
  }
  
  # Save final summary
  addWorksheet(wb_summary, "Summary")
  writeData(wb_summary, "Summary", summary_all)
  
  saveWorkbook(wb_results, results_excel, overwrite = TRUE)
  saveWorkbook(wb_summary, summary_excel, overwrite = TRUE)
  
  cat("Batch complete. Files saved.\n")
}

theta_vals <- c(0,0.5,1,1.5,2)     # location shifts
delta_vals <- c(1,1.25,1.5,2)     # scale shifts

run_lepage_simulation_batch(
  theta_values = theta_vals,
  delta_values = delta_vals,
  results_excel = "ELWA_Results_All.xlsx",
  summary_excel = "ELWA_Summary_All.xlsx"
)

