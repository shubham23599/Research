# Shewhart Lepage based on V-M Chart Simulation
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
H <- 12.7

# Set parameters for IC distribution
mu <- 0 
sigma <- 1

m <- 100  # Reference sample size
n <- 5    # Monitoring Sample size
N <- m + n

# Indicator function
Indicator = rep(c(0, 1), c(m, n))  

# Calculate mean and variance for Mood test
m2 <- (n*(N^2-1))/(12)
s2 <- 1/180*(m*n*(N+1)*(N^2-4))

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
      S <- 0
      u <- rnorm(m, mu, sigma)
      
      # Run simulation until control limit is reached
      while (S <= H & RL[j] <= 5000) {
        v <- rnorm(n, mu + delta[d], sg)
        Ranks <- rank(c(u, v))
        
        # Computing test statistic XN for VW and variance for VW
        XN <- sum(Indicator * (qnorm(Ranks/(N + 1))))
        s_1 <- (m * n) / (N * (N - 1)) * sum(qnorm(Ranks/(N + 1))^2)
        S <- ((XN - 0) / sqrt(s_1))^2  +
          ((sum((Ranks-1/2*(N+1))^2*Indicator)-m2)/sqrt(s2))^2
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
write.csv(data_arl, "L_ShewARL_V_M.csv", row.names = FALSE)
