###############################################################################
# EWMA Lepage (VW–Mood) Control Chart based Simulation
#
# Purpose:
#   Compute the Zero-State Average Run Length (ZSARL) and Steady-State 
#   Average Run Length (SSARL) for the EWMA Lepage EL(V,M) 
#   control chart under joint location and scale shifts.
#
# Description:
#   This code implements the simulation procedure described in Algorithm 2 
#   (Steps 1–9) of the  manuscript. The procedure evaluates the 
#   performance of the EWMA Lepage chart by generating monitoring samples from 
#   a specified distribution and computing the EL(V,M) statistic sequentially 
#   until a signal occurs.
#
# Distribution Options:
#       - Normal:        rnorm()
#       - Log-normal:    rlnorm()
#       - Laplace:       rdoublex()   (from smoothmest package)
#       - Logistic:      rlogis()
#
# Usage Notes:
#   - To compute ZSARL, set tau = 0
#   - To compute SSARL, set tau = 100 
#
#   The simulation is repeated across multiple replications to obtain stable 
#   estimates of both ZSARL and SSARL for the given shift scenario.
###############################################################################

rm(list = ls())
library(smoothmest)

###############################################################################
# Step 1: Set parameters m, n, λ, θ ≠ 0, δ ≠ 1 and warm-up length τ
###############################################################################
Simulations <- 50000          # Number of Monte Carlo replications (B)
lambda <- 0.05                # EWMA smoothing parameter
mu_lep <- 2                   # IC mean of Lepage statistic
sd_lep <- 2                   # IC SD of Lepage statistic
m <- 100                      # Reference sample size
n <- 5                        # Monitoring sample size
tau <- 100                    # Steady-state warm-up length (τ)

# IC distribution parameters (θ = 0, δ = 1)
mu <- 0
sigma <- 1

N <- m + n                    # Combined sample size
Indicator <- rep(c(0, 1), c(m, n))  # Group indicator for VW and Mood statistics

###############################################################################
# Step 2: Determine K_(p,q) using Algorithm 1 and compute H_(p,q)
###############################################################################
K <- 2.461
H <- round(mu_lep + K * sd_lep * sqrt(lambda / (2 - lambda)), 3)

###############################################################################
# Test statistic components: mean and variance of Mood test
###############################################################################
m2 <- (n * (N^2 - 1)) / 12
s2 <- (1 / 180) * (m * n * (N + 1) * (N^2 - 4))

###############################################################################
# Step 3 & 5: Shift combinations 
###############################################################################
delta <- c(0.00, 0.25, 0.50, 1.00, 2.00)
sg0   <- c(1.00, 1.25, 1.50, 2.00)

###############################################################################
# Output storage
###############################################################################
ARL  <- numeric()

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
      
      Z <- mu_lep         # Initialize EWMA statistic
      steady <- 1         # Warm-up counter
      
      ##########################################################################
      # Step 4: Warm-up (steady-state) phase
      #         Generate IC samples V, update EL(p,q),
      #         If EL > H, restart warm-up.
      ##########################################################################
      while (steady < (tau+1)){
        
        v <- rnorm(n, mu, sigma)   # IC monitoring sample
        
        # Compute combined ranks
        Ranks <- rank(c(u, v))
        
        VW <- ((sum(Indicator * qnorm(Ranks/(N + 1)))) /
                 sqrt((m * n)/(N * (N - 1)) * sum(qnorm(Ranks/(N + 1))^2)))^2
        
        Mood <- ((sum((Ranks - 0.5*(N+1))^2 * Indicator) - m2) / sqrt(s2))^2
        
        S <- VW + Mood
        
        Z <- lambda * S + (1 - lambda) * Z
        Z <- max(mu_lep, Z)
        
        # Reset warm-up if control limit crossed
        if (Z > H) {
          steady <- 1
          Z <- mu_lep
        } else {
          steady <- steady + 1
        }
      }
      
      ##########################################################################
      # Step 5–8: Monitoring phase under OOC distribution
      # Step 7: Continue until EL > H
      # Step 8: Record run length
      ##########################################################################
      while (Z <= H & RL[j] <= 5000) {
        # Step 5: Generate OOC sample V
        v <- rnorm(n, mu + delta[d], sg)
        Ranks <- rank(c(u, v))
        
        VW <- ((sum(Indicator * qnorm(Ranks/(N + 1)))) /
                 sqrt((m * n)/(N * (N - 1)) * sum(qnorm(Ranks/(N + 1))^2)))^2
        
        Mood <- ((sum((Ranks - 0.5*(N+1))^2 * Indicator) - m2) / sqrt(s2))^2
        
        S <- VW + Mood
        
        # Step 6: Update EL(p,q)
        Z <- lambda * S + (1 - lambda) * Z
        Z <- max(mu_lep, Z)
        
        RL[j] <- RL[j] + 1   # Step 8: Increase run-length count
      }
    }
    
    ############################################################################
    # Performance statistics
    ############################################################################
    ARL[d]  <- mean(RL)

  }
  
  # Store block results for current sg
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
# Save results
###############################################################################
print(data_arl)
filename <- paste0(lambda, "_EWMA_ARL_VW_Mood.csv")
write.csv(data_arl, filename, row.names = FALSE)
paste0("\n✅ Simulation complete: Results saved as ", filename)
###############################################################################
