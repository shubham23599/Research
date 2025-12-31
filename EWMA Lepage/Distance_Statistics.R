ZWRS = function(u,v){
  m=length(u); n=length(v);N=m+n
  m1=1/2*(n*(N+1));s1=1/12*(m*n*(N+1))
  return((sum(rep(c(0,1),c(m,n))*rank(c(u,v)))-m1)/sqrt(s1))
}
ZAB = function(u,v){
  m=length(u); n=length(v);N=m+n
  if(N%%2==0 ){m2=(n*N)/4;s2=1/48*(m*n*(N^2-4)/(N-1)) }else{m2=(n*(N^2-1))/(4*N);s2=1/48*(m*n*(N+1)*(N^2+3))/N^2}
  return((sum(abs((rank(c(u,v)))-1/2*(N+1))*rep(c(0,1),c(m,n)))-m2)/sqrt(s2))
}
ZVW = function(u,v){
  m=length(u); n=length(v);N=m+n
  Ind = rep(c(0, 1), c(m, n))  
  XN = sum(Ind * (qnorm(rank(c(u, v)) / (N + 1))))
  s_1 = (m * n) / (N * (N - 1)) * sum(qnorm(rank(c(u, v)) / (N + 1))^2)
   return((XN - 0) / sqrt(s_1))
}
ZM = function(u,v){
  m=length(u); n=length(v);N=m+n
  m2=(n*(N^2-1))/(12);s2=1/180*(m*n*(N+1)*(N^2-4))
  return((sum(((rank(c(u,v)))-1/2*(N+1))^2*rep(c(0,1),c(m,n)))-m2)/sqrt(s2))
}

WRSAB =function(u,v){
  m=length(u); n=length(v);N=m+n
  m1=1/2*(n*(N+1));s1=1/12*(m*n*(N+1))
  if(N%%2==0 ){m2=(n*N)/4;s2=1/48*(m*n*(N^2-4)/(N-1)) }else{m2=(n*(N^2-1))/(4*N);s2=1/48*(m*n*(N+1)*(N^2+3))/N^2}
  S=((sum(rep(c(0,1),c(m,n))*rank(c(u,v)))-m1)/sqrt(s1))^2+((sum(abs((rank(c(u,v)))-1/2*(N+1))*rep(c(0,1),c(m,n)))-m2)/sqrt(s2))^2
  return(S)
}
WRSM =function(u,v){
  m=length(u); n=length(v);N=m+n
  m1=m*n/2;s1=1/12*(m*n*(N+1))
  m2=(n*(N^2-1))/(12);s2=1/180*(m*n*(N+1)*(N^2-4))
  S=((sum(rep(c(0,1),c(m,n))*rank(c(u,v)))-m1)/sqrt(s1))^2+((sum(((rank(c(u,v)))-1/2*(N+1))^2*rep(c(0,1),c(m,n)))-m2)/sqrt(s2))^2
}
MWAB = function(u,v){
  m=length(u); n=length(v);N=m+n
  m1=m*n/2
  s1=1/12*(m*n*(N+1))
  # Calculate mean and variance for AB test
  if(N%%2==0 ){ m2=(n*N)/4;s2=1/48*(m*n*(N^2-4)/(N-1))
  } else { m2=(n*(N^2-1))/(4*N); s2=1/48*(m*n*(N+1)*(N^2+3))/N^2}
  S=((sum(sapply(1:length(u),function(i) u[i]>v))-m1)/sqrt(s1))^2+
    ((sum(abs((rank(c(u,v)))-1/2*(N+1))*rep(c(0,1),c(m,n)))-m2)/sqrt(s2))^2
  return(S)
}
MWM = function(u,v){
  m=length(u); n=length(v);N=m+n
  m1=m*n/2;s1=1/12*(m*n*(N+1))
  m2=(n*(N^2-1))/(12);s2=1/180*(m*n*(N+1)*(N^2-4))
  S=((sum(sapply(1:length(u),function(i) u[i]>v))-m1)/sqrt(s1))^2+((sum(((rank(c(u,v)))-1/2*(N+1))^2*rep(c(0,1),c(m,n)))-m2)/sqrt(s2))^2
  return(S)
}
SC =function(u,v){
  m=length(u); n=length(v);N=m+n
  m2 = (n * (N + 1) * (2 * N + 1)) / 6  # Constant used in calculation
  s2 = (m * n) / 180 * ((N + 1) * (2 * N + 1) * (8 * N + 11))  # Constant used in calculation
  pho = ((2 * (N^2 - 4)) / ((2 * N + 1) * (8 * N + 11))) - 1  # Constant used in calculation
  U = ((sum(rep(c(0, 1), c(m, n)) * rank(c(u, v))^2) - m2) / sqrt(s2))  # Calculate U statistic
  V = ((sum(rep(c(0, 1), c(m, n)) * (N + 1 - rank(c(u, v)))^2) - m2) / sqrt(s2))  # Calculate V statistic
  C = (U^2 + V^2 - 2 * pho * U * V) / (2 * (1 - pho^2))  # Calculate S statistic
  return(C)
}
VWAB = function(u,v){
  m=length(u); n=length(v);N=m+n
  Ind = rep(c(0, 1), c(m, n))  
  if(N%%2==0 ){ m2=(n*N)/4;s2=1/48*(m*n*(N^2-4)/(N-1))
  } else { m2=(n*(N^2-1))/(4*N); s2=1/48*(m*n*(N+1)*(N^2+3))/N^2}
  XN = sum(Ind * (qnorm(rank(c(u, v)) / (N + 1))))
  s_1 = (m * n) / (N * (N - 1)) * sum(qnorm(rank(c(u, v)) / (N + 1))^2)
  # Computing the test statistic S
  S = ((XN - 0) / sqrt(s_1))^2 + ((sum(abs((rank(c(u, v))) - 1/2 * (N + 1)) * Ind) - m2) / sqrt(s2))^2
  return(S)
}
VWM = function(u,v){
  m=length(u); n=length(v);N=m+n
  Ind = rep(c(0, 1), c(m, n))  
  m2 <- (n * (N^2 - 1)) / (12)
  s2 <- 1 / 180 * (m * n * (N + 1) * (N^2 - 4))
  XN <- sum(Ind * (qnorm(rank(c(u, v)) / (N + 1))))
  s_1 <- (m * n) / (N * (N - 1)) * sum(qnorm(rank(c(u, v)) / (N + 1))^2)
  S <- ((XN - 0) / sqrt(s_1))^2 + ((sum(((rank(c(u, v))) - 1/2 * (N + 1))^2 * Ind) - m2) / sqrt(s2))^2
  return(S)
}