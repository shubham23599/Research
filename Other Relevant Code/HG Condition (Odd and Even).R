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
m <- 100
n <- 5
wrs <- 0
vw <- 0
ab <- 0
mood <- 0
N <- 10000
for(i in 1:N){
  u <- rlnorm(m)
  v <- rlnorm(n)
  wrs_ps <- WRS(u,v)-(n*(m+n+1))/2
  wrs_ns <- WRS(-u,-v)-(n*(m+n+1))/2
  if(round(wrs_ps,8) == -round(wrs_ns,8)) wrs <- wrs+1
  vw_ps <- VW(u,v)
  VW_ns <- VW(-u,-v)
  if(round(vw_ps,8) == round(-VW_ns,8)) vw <- vw+1
  
  ab_ps <- AB(u,v)
  ab_ns <- AB(u,v)
  if(round(ab_ps,8) == round(ab_ns,8)) ab <- ab+1
  mood_ps <- Mood(u,v)
  mood_ns <- Mood(-u,-v)
  if(round(mood_ps,8) == round(mood_ns,8)) mood <- mood+1
  
  

  }
wrs/N
vw/N
ab/N
mood/N