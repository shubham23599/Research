###############################################################################
# EWMA Cucconi Control Chart Simulation
#
# Purpose:
#   Compute the Zero-State Out-of-Control Average Run Length (OOC ZSARL)
#   of the EWMA Cucconi EC(p,q) control chart under location/scale shifts.
#
# Notes:
#   To simulate for other distributions:
#       - Normal:       rnorm()
#       - Log-normal:   rlnorm()
#       - Laplace:      rdoublex()  (from smoothmest)
#       - Logistic:     rlogis()
#
# Algorithm 2 implemented below (Steps 1–9).
###############################################################################

# ------------------------------ INITIAL SETUP --------------------------------
rm(list = ls())
library(smoothmest)      # for Laplace (double-exponential) random numbers

###############################################################################
# Step 1: Set parameters m, n, λ, θ ≠ 0, δ ≠ 1 and warm-up length τ
###############################################################################
Simulations <- 50000          # Number of Monte Carlo replications (B)
lambda <- 0.05              # EWMA smoothing parameter
mu_cuc <- 1                 # IC mean of Cucconi statistic
sd_cuc <- 1                 # IC SD of Cucconi statistic
m <- 100                    # Reference sample size
n <- 5                      # Monitoring sample size
tau <- 0                    # Steady-state warm-up length (τ)

# IC distribution parameters (θ = 0, δ = 1)
mu <- 0
sigma <- 1

N <- m + n                  # Combined sample size
Indicator <- rep(c(0, 1), c(m, n))  # Group indicator for ranks

###############################################################################
# Step 2: Determine K_(p,q) using Algorithm 1 and compute H_(p,q)
###############################################################################
K <- 2.417                   # Charting constant obtained from Algorithm 1
H <- round(mu_cuc + K * sd_cuc * sqrt(lambda / (2 - lambda)), 3)

###############################################################################
# Test statistic components: Cucconi U and V
###############################################################################
m2 <- (n * (N + 1) * (2 * N + 1)) / 6
s2 <- (m * n * (N + 1) * (2 * N + 1) * (8 * N + 11)) / 180
rho <- ((2 * (N^2 - 4)) / ((2 * N + 1) * (8 * N + 11))) - 1

###############################################################################
# Step 3 & Step 5 shift settings: location shift δ and scale shift σ_g
###############################################################################
delta <- c(0.0, 0.25, 0.5, 1.0, 2.0)     # Location shifts
sg0   <- c(1.0, 1.25, 1.5, 2.0)          # Scale shifts

###############################################################################
# Output storage
###############################################################################
ARL <- numeric()

data_arl <- data.frame(
  "Delta" = numeric(),
  "Sg"    = numeric(),
  "ARL"   = numeric()
  )

count <- 1


###############################################################################
# ------------------------------ MAIN SIMULATION ------------------------------
###############################################################################
for (sg in sg0) {                  # Loop over scale shifts
  for (d in 1:length(delta)) {     # Loop over location shifts
    
    RL <- numeric(Simulations)     # Store run length for each replication
    
    ############################################################################
    # Step 9: Repeat Steps 3–8 for B replications
    ############################################################################
    for (j in 1:Simulations) {
      
      ##########################################################################
      # Step 3: Generate an IC reference sample U (θ = 0, δ = 1)
      ##########################################################################
      u <- rnorm(m, mu, sigma)
      
      Z <- mu_cuc         # Initialize EWMA statistic
      steady <- 1         # Warm-up counter
      
      ##########################################################################
      # Step 4: Warm-up (steady-state) phase
      #         Generate IC samples V, update EC(p,q),
      #         If EC > H, restart warm-up.
      ##########################################################################
      while (steady < (tau+1)) {
        
        v <- rnorm(n, mu, sigma)   # IC monitoring sample
        
        # Compute combined ranks
        Ranks <- rank(c(u, v))
        
        # Cucconi U and V components
        U <- (sum(Indicator * Ranks^2) - m2) / sqrt(s2)
        V <- (sum(Indicator * (N + 1 - Ranks)^2) - m2) / sqrt(s2)
        
        # Cucconi statistic (standardized quadratic form)
        S <- (U^2 + V^2 - 2 * rho * U * V) / (2 * (1 - rho^2))
        
        # EWMA update
        Z <- lambda * S + (1 - lambda) * Z
        Z <- max(mu_cuc, Z)        # Avoid downward drift
        
        # Restart if control limit exceeded
        if (Z > H) {
          steady <- 1
          Z <- mu_cuc
        } else {
          steady <- steady + 1
        }
      }
      
      ##########################################################################
      # Step 5–8: Monitoring phase under OOC distribution
      # Step 7: Continue until EC > H
      # Step 8: Record run length
      ##########################################################################
      while (Z <= H & RL[j] <= 5000) {
        
        # Step 5: Generate OOC sample V
        v <- rnorm(n, mu + delta[d], sg)
        
        # Compute ranks
        Ranks <- rank(c(u, v))
        
        # Cucconi components
        U <- (sum(Indicator * Ranks^2) - m2) / sqrt(s2)
        V <- (sum(Indicator * (N + 1 - Ranks)^2) - m2) / sqrt(s2)
        S <- (U^2 + V^2 - 2 * rho * U * V) / (2 * (1 - rho^2))
        
        # Step 6: Update EC(p,q)
        Z <- lambda * S + (1 - lambda) * Z
        Z <- max(mu_cuc, Z)
        
        RL[j] <- RL[j] + 1   # Step 8: Increase run-length count
      }
    } # End simulation loop
    
    ############################################################################
    # Compute performance measures for each shift combination
    ############################################################################
    ARL[d]  <- mean(RL)

  }
  
  # Store results for this scale shift
  d1 <- data.frame(
    "Delta" = delta,
    "Sg"    = sg,
    "ARL"   = ARL
  )
  
  print(d1)
  
  data_arl[count:(count + length(delta) - 1), ] <- d1
  count <- count + length(delta)
}

###############################################################################
# Save all results
###############################################################################
print(data_arl)
filename <- paste0(lambda, "_EWMA_Cucconi_ARL.csv")
write.csv(data_arl, filename, row.names = FALSE)
paste0("\n✅ Simulation complete: Results saved as ", filename)
###############################################################################
