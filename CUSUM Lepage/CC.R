# CUSUM Cucconi Chart Simulation
# To simulate for other distributions, replace rnorm() with the appropriate function.
# For a log-normal distribution, use rlnorm().
# For a Laplace distribution, use rdoublex(). 
# For a logistic distribution, use rlogis().

# Clear workspace
rm(list=ls())
library('smoothmest')

# Set the number of simulations
Simulations <- 50000

# Set the control limits for ARL calculations
H <- 12.471

# Mean of Lepage or Lepage-type statistic
mu_lep <- 1

# Set the value of k 
k <- 0

# Set parameters for IC distribution
mu <- 0 
sigma <- 1

m <- 100  # Reference sample size
n <- 5    # Monitoring Sample size
N <- m + n

# Indicator function
Indicator = rep(c(0, 1), c(m, n))  

# Calculate Required Constants for test statistic
m2 <-(n*(N+1)*(2*N+1))/6
s2 <- (m*n)/180*((N+1)*(2*N+1)*(8*N+11))
rho <-((2*(N^2-4))/((2*N+1)*(8*N+11)))-1

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
  for(d in 1:length(delta)) {
    RL <- numeric(Simulations) # Initialize RL vector for each simulation
    # Loop over simulations
    for (j in 1:Simulations) {
      C <- 0
      u <- rnorm(m, mu, sigma)
      
      # Run simulation until control limit is reached
      while (C <= H & RL[j] <= 5000) {
        v <- rnorm(n, mu + delta[d], sg)
        Ranks <- rank(c(u, v))
        U <- ((sum(Indicator*Ranks^2)-m2)/sqrt(s2))
        V <- ((sum(Indicator*(N+1-Ranks)^2)-m2)/sqrt(s2))
        S <-(U^2+V^2-2*rho*U*V)/(2*(1-rho^2))
        C <- max(0,C+(S-mu_lep)-k)
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
write.csv(data_arl, "L_CUSUM_Cucconi.csv", row.names = FALSE)
