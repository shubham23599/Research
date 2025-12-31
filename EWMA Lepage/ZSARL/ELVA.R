# EWMA Lepage based on V-AB Chart Simulation
# To simulate for other distributions, replace rnorm() with the appropriate function.
# For a log-normal distribution, use rlnorm().
# For a Laplace distribution, use rdoublex(). 
# For a logistic distribution, use rlogis().

# Clear workspace
rm(list=ls())
library('smoothmest')
# Set the number of simulations
Simulations <- 50000

# Mean and SD of Lepage or Lepage-type statistic
mu_lep <- 2
sd_lep <- 2

# Set the value of lambda 
lambda <- 0.05

# Charting constant
K <- 2.479

# Estimate the control limits for ARL calculations
H <- round(mu_lep + K*sd_lep*sqrt(lambda/(2-lambda)),3)

# Set parameters for IC distribution
mu <- 0 
sigma <- 1

m <- 100  # Reference sample size
n <- 5    # Monitoring Sample size
N <- m + n

# Indicator function
Indicator = rep(c(0, 1), c(m, n))  

# Calculate mean and variance for AB test
if (N %% 2 == 0) {
  m2 <- (n * N) / 4
  s2 <- 1/48 * (m * n * (N^2 - 4) / (N - 1))
} else {
  m2 <- (n * (N^2 - 1)) / (4 * N)
  s2 <- 1/48 * (m * n * (N + 1) * (N^2 + 3)) / N^2
}

# Define delta and sg values
delta <- c(seq(0, 1.5, .25), 2, 3) # For Location
sg0 <- seq(1, 2, 0.25)  # For Scale

# Initialize ARL vector
ARL <- numeric()

# Create an empty data frame for storing results
data_arl <- data.frame("Delta" = numeric(),
                       "Sg" = numeric(),
                       "ARL" = numeric())
count<-1
# Loop over sg values
for (sg in sg0) {
  # Loop over delta values
  for(d in 1:length(delta)){
    RL <- numeric(Simulations) # Initialize RL vector for each simulation
    # Loop over simulations
    for (j in 1:Simulations) {
      Z <- mu_lep
      u <- rnorm(m, mu, sigma) 
      
      # Run simulation until control limit is reached
      while (Z <= H & RL[j] <= 5000) {
        v <- rnorm(n, mu + delta[d], sg)
        Ranks <- rank(c(u, v))
        
        # Computing test statistic XN for VW and variance for VW
        XN <- sum(Indicator * (qnorm(Ranks/(N + 1))))
        s_1 <- (m * n) / (N * (N - 1)) * sum(qnorm(Ranks/(N + 1))^2)
        S <- ((XN - 0) / sqrt(s_1))^2  +
          ((sum(abs((Ranks - 1/2 * (N + 1))) * Indicator) - m2) / sqrt(s2))^2
        Z <- lambda*S + (1-lambda)*Z
        Z <- max(mu_lep,Z)
        RL[j] <- RL[j] + 1
      }
    }
    # Calculate ARL and store in ARL vector
    ARL[d] <- mean(RL)
  }
  d1 <- data.frame("delta"=delta,'Sg'=sg,"ARL"=ARL)
  print(d1)
  data_arl[count:(count+(length(delta)-1)),] <- d1
  count<-count+length(delta)
  
}

# Print the final data frame
print(data_arl)

# Save results to a CSV file
write.csv(data_arl, "L_EWMA_ARL_V_AB.csv", row.names = FALSE)
