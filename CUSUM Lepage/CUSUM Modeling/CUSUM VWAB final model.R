# Load required libraries
library(lmtest)
library(nortest)
library(ggplot2)
library(gridExtra)

# Load the data
data <- read.csv("CUSUM Modeling//CUSUM VWAB ARL.csv")

# Define the proportion for the train-test split
set.seed(123)  # For reproducibility
train_index <- sample(1:nrow(data), round(nrow(data) * 0.8), replace = FALSE)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]
write.csv(train_data,"CUSUM Modeling//TRAIN_VWAB.csv")
write.csv(test_data,"CUSUM Modeling//TEST_VWAB.csv")
# Define a small constant to avoid log(0)
small_constant <- 0.00001

# Define the model formula with log transformations on the independent variables
formula <- H ~ log(k + small_constant) + log(m) + log(n) + log(ARL)

# Fit the linear model
model <- lm(formula, data = train_data)

# Summary of the model
model_summary <- summary(model)
print(model_summary)
# Assumption tests
# 1. Normality of residuals (Shapiro-Wilk test)
shapiro_p <- shapiro.test(residuals(model))$p.value

# 2. Homoscedasticity (Breusch-Pagan test)
bp_p <- bptest(model)$p.value

# 3. Autocorrelation (Durbin-Watson test)
dw_p <- dwtest(model)$p.value

# Calculate RMSE and MAPE on test data
# Predictions on test set
test_predictions <- predict(model, newdata = test_data)

# Actual values
test_actuals <- test_data$H



# RMSE (Root Mean Squared Error)
test_rmse <- sqrt(mean((test_actuals - test_predictions)^2))

# MAPE (Mean Absolute Percentage Error)
test_mape <- mean(abs((test_actuals - test_predictions) / test_actuals)) * 100

# Display results
cat("Model Summary:\n")
print(model_summary)

cat("\nShapiro-Wilk p-value for residuals:", round(shapiro_p, 4), "\n")
cat("Breusch-Pagan p-value for heteroscedasticity:", round(bp_p, 4), "\n")
cat("Durbin-Watson p-value for autocorrelation:", round(dw_p, 4), "\n")
cat("Test RMSE:", round(test_rmse, 4), "\n")
cat("Test MAPE:", round(test_mape, 2), "%\n")

# Residuals vs Fitted
p1 <- ggplot(data = as.data.frame(model$fitted.values), aes(x = model$fitted.values, y = model$residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals")

# Q-Q Plot
p2 <- ggplot(data = as.data.frame(model$residuals), aes(sample = model$residuals)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Normal Q-Q", x = "Theoretical Quantiles", y = "Standardized Residuals")

# Scale-Location Plot
p3 <- ggplot(data = as.data.frame(model$fitted.values), aes(x = model$fitted.values, y = sqrt(abs(model$residuals)))) +
  geom_point() +
  geom_smooth(se = FALSE, color = "blue") +
  labs(title = "Scale-Location", x = "Fitted values", y = "Sqrt(Standardized Residuals)")
# Extract leverage values and residuals
leverage_values <- hatvalues(model)
residuals <- model$residuals

# Create a data frame with these values
diagnostics_df <- data.frame(Leverage = leverage_values, Residuals = residuals)
# Leverage Plot
p4 <- ggplot(diagnostics_df, aes(x = Leverage, y = Residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Leverage vs Residuals", x = "Leverage", y = "Residuals")

# Arrange plots in a grid
grid.arrange(p1, p2, p3, p4, nrow = 2)

library(car)  # For vif()

# Calculate VIF for the model
vif_values <- vif(model)

# Define a threshold for VIF
vif_threshold <- 5

# Check if all VIF values are below the threshold
vif_flag <- ifelse(all(vif_values < vif_threshold), "Yes", "No")

# Create the final results dataframe
final_results <- data.frame(
  Test = c("Shapiro-Wilk (Normality of Residuals)",
           "Breusch-Pagan (Heteroscedasticity)",
           "Durbin-Watson (Autocorrelation)",
           "Test RMSE",
           "VIF (Variance Inflation Factor)"),
  P_Value = c(round(shapiro_p, 4),
              round(bp_p, 4),
              round(dw_p, 4),
              NA,  # RMSE does not have a p-value
              NA),  # VIF does not have a p-value
  Value = c(NA,  # p-values for assumption tests
            NA,
            NA,
            round(test_rmse, 4),  # Test RMSE
            paste0("k: ", round(vif_values["log(k + small_constant)"], 2),
                   ", m: ", round(vif_values["log(m)"], 2),
                   ", n: ", round(vif_values["log(n)"], 2),
                   ", ARL: ", round(vif_values["log(ARL)"], 2))
  ),
  Flag = c(ifelse(shapiro_p > 0.01, "Yes", "No"),
           ifelse(bp_p > 0.01, "Yes", "No"),
           ifelse(dw_p > 0.01, "Yes", "No"),
           NA,  # RMSE does not have a flag
           vif_flag)  # VIF flag
)

# Display the final results
View(final_results)
