############################################################
## COMPLETE PHASE-I ANALYSIS FOR AST DATA
############################################################

## ---------------------------------------------------------
## 1. INSTALL & LOAD REQUIRED PACKAGES
## ---------------------------------------------------------
install.packages(c("dfphase1", "tseries", "moments", "ggplot2"))
library(dfphase1)
library(tseries)
library(moments)
library(ggplot2)

############################################################
## 2. READ DATA
############################################################
ast <- read.csv("HWMA Example Data AST.csv")

## Phase-I observations only
phase1 <- subset(ast, Phase == "I")

## Subgroup matrix (n = 5)
X <- as.matrix(phase1[, c("X1", "X2", "X3", "X4", "X5")])
m <- nrow(X)    # number of subgroups
n <- ncol(X)    # subgroup size

## Vector format
x_vec <- as.vector(t(X))

############################################################
## 3. DESCRIPTIVE STATISTICS
############################################################
mean_ast  <- mean(x_vec)
sd_ast    <- sd(x_vec)
se_ast    <- sd_ast / sqrt(length(x_vec))
skew_ast  <- skewness(x_vec)
kurt_ast  <- kurtosis(x_vec)

cat("\n--- Descriptive Statistics ---\n")
cat("Mean        :", mean_ast, "\n")
cat("Std Dev     :", sd_ast, "\n")
cat("Std Error   :", se_ast, "\n")
cat("Skewness    :", skew_ast, "\n")
cat("Kurtosis    :", kurt_ast, "\n")

############################################################
## 4. INDEPENDENCE CHECK
############################################################
cat("\n--- Ljung-Box Independence Test ---\n")
for (lag in c(1, 3, 5, 10)) {
  test <- Box.test(x_vec, lag = lag, type = "Ljung-Box")
  cat("Lag", lag, ": p-value =", round(test$p.value, 4), "\n")
}

############################################################
## 5. NORMALITY TESTS
############################################################
jb_test <- jarque.bera.test(x_vec)
sw_test <- shapiro.test(x_vec)

cat("\n--- Normality Tests ---\n")
cat("Jarque-Bera p-value :", round(jb_test$p.value, 4), "\n")
cat("Shapiro-Wilk p-value:", round(sw_test$p.value, 4), "\n")

############################################################
## 6. KERNEL DENSITY PLOT
############################################################
ggplot(data.frame(AST = x_vec), aes(x = AST)) +
  geom_density(fill = "lightblue", alpha = 0.6) +
  labs(
    title = "Kernel Density Estimate of Phase-I AST Data",
    x = "Average Service Time",
    y = "Density"
  ) +
  theme_minimal()

############################################################
## 7. RS/P PHASE-I ANALYSIS
############################################################
rsp_res <- rsp(
  y        = x_vec,
  plot     = TRUE,
  L        = 10000,
  seed     = 11642257,
  alpha    = 0.10
)

cat("\n--- RS/P Phase-I Results ---\n")
cat("Location p-value:", rsp_res$p.loc, "\n")
cat("Scale p-value   :", rsp_res$p.scale, "\n")

############################################################
## 8. SHEWHART-TYPE PHASE-I ANALYSIS
############################################################
## XbarS chart (parametric)
shewhart(
  x           = X,
  stat        = "XbarS",
  aggregation = "mean",
  plot        = TRUE,
  FAP         = 0.10,
  seed        = 11642257,
  L           = 1000
)

## Nonparametric
shewhart(X, stat="Rank", aggregation="mean", plot=TRUE, FAP=0.10,seed = 11642257,L= 1000)
shewhart(X, stat="lRank", aggregation="mean", plot=TRUE, FAP=0.10,seed = 11642257,L= 1000)
shewhart(X, stat="sRank", aggregation="mean", plot=TRUE, FAP=0.10,seed = 11642257,L= 1000)

############################################################
## END OF COMPLETE PHASE-I ANALYSIS
############################################################
