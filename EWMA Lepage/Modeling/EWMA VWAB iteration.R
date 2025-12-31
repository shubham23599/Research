# Load required libraries
library(car)
library(MASS)
library(lmtest)
library(nortest)
library(caret)
data <- read.csv("Modeling//EWMA VWAB ARL.csv")
# Define the proportion for the train-test split
set.seed(123)  # For reproducibility
train_index <- sample(1:nrow(data), round(nrow(data) * 0.8), replace = FALSE)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]
write.csv(train_data,"E_TRAIN_VWAB.csv")
write.csv(test_data,"E_TEST_VWAB.csv")
# Define the main terms
terms <- c("k", "m", "n", "ARL")

# Generate only unique two-way interaction combinations (order 2)
two_way_interactions <- combn(terms, 2, simplify = FALSE)
two_way_interactions <- sapply(two_way_interactions, function(x) paste(x, collapse = ":"))

# Initialize a list to store the model types
model_types <- list(
  "With Intercept" = "H ~ k + m + n + ARL",
  "Without Intercept" = "H ~ k + m + n + ARL - 1"
)

# Add two-way interaction combinations to model types
for(interaction in two_way_interactions) {
  model_types[[paste("Interaction", interaction, sep = ": ")]] <- paste("H ~ k + m + n + ARL +", interaction)
}

# Initialize a dataframe to store results
results <- data.frame()

# Function to reverse transformations
reverse_transform <- function(x, transform_name) {
  if (transform_name == "log") {
    return(exp(x))
  } else if (transform_name == "sqrt") {
    return(x^2)
  } else if (transform_name == "square") {
    return(sqrt(x))
  } else if (transform_name == "normalize") {
    return(x * (max(train_data$H) - min(train_data$H)) + min(train_data$H))
  } else if (transform_name == "standardize") {
    return(x * sd(train_data$H) + mean(train_data$H))
  } else if (transform_name == "inverse") {
    return(1 / x)
  } else {
    return(x)
  }
}

# Define transformations
transformations <- list(
  none = function(x) x,
  log = function(x) log(x),    # Adding 1 to avoid log(0)
  sqrt = function(x) sqrt(x),
  square = function(x) x^2,
  normalize = function(x) (x - min(x)) / (max(x) - min(x)),
  standardize = function(x) (x - mean(x)) / sd(x),
  inverse = function(x) 1 / (x)
)

# Loop over transformations and model types
for(dep_transform_name in names(transformations)) {
  dep_transform <- transformations[[dep_transform_name]]
  
  for(ind_transform_name in names(transformations)) {
    ind_transform <- transformations[[ind_transform_name]]
    
    for(model_type_name in names(model_types)) {
      formula <- as.formula(model_types[[model_type_name]])
      
      # Apply transformations on train data
      transformed_train_data <- data.frame(
        k = ind_transform(train_data$k),
        m = ind_transform(train_data$m),
        n = ind_transform(train_data$n),
        ARL = ind_transform(train_data$ARL),
        H = dep_transform(train_data$H)
      )
      
      # Apply transformations on test data
      transformed_test_data <- data.frame(
        k = ind_transform(test_data$k),
        m = ind_transform(test_data$m),
        n = ind_transform(test_data$n),
        ARL = ind_transform(test_data$ARL),
        H = dep_transform(test_data$H)
      )
      
      # Fit the linear model on the training set
      model <- lm(formula, data = transformed_train_data)
      
      # Predict on the train and test sets
      train_predictions <- predict(model, newdata = transformed_train_data)
      test_predictions <- predict(model, newdata = transformed_test_data)
      
      # Retransform predictions back to original scale
      original_train_predictions <- reverse_transform(train_predictions, dep_transform_name)
      original_test_predictions <- reverse_transform(test_predictions, dep_transform_name)
      
      # Retransform the actual values back to original scale
      original_train_actuals <- reverse_transform(transformed_train_data$H, dep_transform_name)
      original_test_actuals <- reverse_transform(transformed_test_data$H, dep_transform_name)
      
      # Calculate performance metrics for both train and test sets
      train_rmse <- sqrt(mean((original_train_actuals - original_train_predictions)^2))
      train_mape <- mean(abs((original_train_actuals - original_train_predictions) / original_train_actuals)) * 100
      test_rmse <- sqrt(mean((original_test_actuals - original_test_predictions)^2))
      test_mape <- mean(abs((original_test_actuals - original_test_predictions) / original_test_actuals)) * 100
      
      # Calculate additional diagnostics for the model
      adj_r2 <- summary(model)$adj.r.squared
      shapiro_p <- shapiro.test(model$residuals)$p.value
      dw_p <- dwtest(model)$p.value
      bp_p <- bptest(model)$p.value
      vif_values <- vif(model)
      
      # Extract model coefficients and p-values
      coeffs <- summary(model)$coefficients
      coeff_p_values <- coeffs[, "Pr(>|t|)"]
      
      # Determine flags for significance and diagnostic tests
      all_coeff_significant <- all(coeff_p_values < 0.01)
      normality_pass <- shapiro_p > 0.01
      autocorrelation_pass <- dw_p > 0.01
      vif_pass <- all(vif_values < 5)
      heteroscedasticity_pass <- bp_p > 0.01
      
      # Final flag: All conditions must pass
      final_flag <- all(all_coeff_significant, normality_pass, autocorrelation_pass, vif_pass, heteroscedasticity_pass)
      
      # Construct the equation for the model
      equation <- paste(
        if (model_type_name == "With Intercept") {
          paste("H =", round(coeffs["(Intercept)", "Estimate"], 4))
        } else {
          "H ="
        },
        paste0(
          round(coeffs[!rownames(coeffs) %in% "(Intercept)", "Estimate"], 4),
          "*",
          rownames(coeffs)[!rownames(coeffs) %in% "(Intercept)"],
          collapse = " + "
        )
      )
      
      # Store results
      results <- rbind(results, data.frame(
        Model_Type = model_type_name,
        Dependent_Transformation = dep_transform_name,
        Independent_Transformation = ind_transform_name,
        Equation = equation,
        Adjusted_R2 = adj_r2,
        Train_RMSE = train_rmse,
        Train_MAPE = train_mape,
        Test_RMSE = test_rmse,
        Test_MAPE = test_mape,
        Normality_p_value = shapiro_p,
        DW_p_value = dw_p,
        BP_p_value = bp_p,
        VIF_k = vif_values['k'],
        VIF_m = vif_values['m'],
        VIF_n = vif_values['n'],
        VIF_ARL = vif_values['ARL'],
        Coeff_k = coeffs['k', 'Estimate'],
        Coeff_m = coeffs['m', 'Estimate'],
        Coeff_n = coeffs['n', 'Estimate'],
        Coeff_ARL = coeffs['ARL', 'Estimate'],
        P_value_k = coeff_p_values['k'],
        P_value_m = coeff_p_values['m'],
        P_value_n = coeff_p_values['n'],
        P_value_ARL = coeff_p_values['ARL'],
        All_Coefficients_Significant = all_coeff_significant,
        Normality_Pass = normality_pass,
        Autocorrelation_Pass = autocorrelation_pass,
        VIF_Pass = vif_pass,
        Heteroscedasticity_Pass = heteroscedasticity_pass,
        Final_Flag = final_flag
      ))
    }
  }
}

# Display the results
View(results)
