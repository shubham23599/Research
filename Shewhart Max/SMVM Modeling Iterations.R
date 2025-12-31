# Load required libraries
library(car)
library(MASS)
library(lmtest)
library(nortest)
library(caret)

# Load data
df <- read.csv("SMVM Modelling Data.csv")
data <- df # [df$a>0,]
# Define the proportion for the train-test split
set.seed(123)  # For reproducibility
train_index <- sample(1:nrow(data), round(nrow(data) * 0.9), replace = FALSE)
train_data <- data[train_index, ]
test_data <- data[-train_index, ]


# Define transformations
transformations <- list(
  none = function(x) x,
  log = function(x) log(x),    # Adding 1 to avoid log(0)
  sqrt = function(x) sqrt(x),
  square = function(x) x^2,
  normalize = function(x) (x - min(x)) / (max(x) - min(x)),
  standardize = function(x) (x - mean(x)) / sd(x),
  inverse = function(x) 1 / (x)  # Avoid division by zero
)

# Function to reverse transformations
reverse_transform <- function(x, transform_name) {
  if (transform_name == "log") {
    return(exp(x))
  } else if (transform_name == "sqrt") {
    return(x^2)
  } else if (transform_name == "square") {
    return(sqrt(x))
  } else if (transform_name == "normalize") {
    return(x * (max(train_data$a) - min(train_data$a)) + min(train_data$a))
  } else if (transform_name == "standardize") {
    return(x * sd(train_data$a) + mean(train_data$a))
  } else if (transform_name == "inverse") {
    return(1 / x)
  } else {
    return(x)
  }
}

# Automatically handle independent and dependent variables
ind_var <- c('m',"n","ARL" )
dep_var <- "H"

# Generate unique two-way interaction combinations
two_way_interactions <- combn(ind_var, 2, simplify = FALSE)
two_way_interactions <- sapply(two_way_interactions, function(x) paste(x, collapse = ":"))

# Define model types
model_types <- list(
  "With Intercept" = paste(dep_var, "~", paste(ind_var, collapse = " + ")),
  "Without Intercept" = paste(dep_var, "~", paste(ind_var, collapse = " + "), "- 1")
)

# Add interaction terms
for (interaction in two_way_interactions) {
  model_types[[paste("Interaction", interaction, sep = ": ")]] <- paste(dep_var, "~", paste(ind_var, collapse = " + "), "+", interaction)
}

# Initialize results dataframe
results <- data.frame()

# Loop through transformations and model types
for (dep_transform_name in names(transformations)) {
  dep_transform <- transformations[[dep_transform_name]]
  
  for (ind_transform_name in names(transformations)) {
    ind_transform <- transformations[[ind_transform_name]]
    
    for (model_type_name in names(model_types)) {
      formula <- as.formula(model_types[[model_type_name]])
      
      # Transform train and test data
      transformed_train_data <- train_data
      transformed_test_data <- test_data
      
      for (var in ind_var) {
        transformed_train_data[[var]] <- ind_transform(train_data[[var]])
        transformed_test_data[[var]] <- ind_transform(test_data[[var]])
      }
      transformed_train_data[[dep_var]] <- dep_transform(train_data[[dep_var]])
      transformed_test_data[[dep_var]] <- dep_transform(test_data[[dep_var]])
      
      # Fit the model
      model <- lm(formula, data = transformed_train_data)
      
      # Predictions
      train_predictions <- predict(model, newdata = transformed_train_data)
      test_predictions <- predict(model, newdata = transformed_test_data)
      
      # Reverse transformations
      original_train_predictions <- reverse_transform(train_predictions, dep_transform_name)
      original_test_predictions <- reverse_transform(test_predictions, dep_transform_name)
      original_train_actuals <- reverse_transform(transformed_train_data[[dep_var]], dep_transform_name)
      original_test_actuals <- reverse_transform(transformed_test_data[[dep_var]], dep_transform_name)
      
      # Metrics
      train_rmse <- sqrt(mean((original_train_actuals - original_train_predictions)^2))
      train_mape <- mean(abs((original_train_actuals - original_train_predictions) / original_train_actuals)) * 100
      test_rmse <- sqrt(mean((original_test_actuals - original_test_predictions)^2))
      test_mape <- mean(abs((original_test_actuals - original_test_predictions) / original_test_actuals)) * 100
      
      # Model diagnostics
      adj_r2 <- summary(model)$adj.r.squared
      shapiro_p <- shapiro.test(model$residuals)$p.value
      dw_p <- dwtest(model)$p.value
      bp_p <- bptest(model)$p.value
      vif_values <- vif(model)
      
      # Coefficients and significance
      coeffs <- summary(model)$coefficients
      coeff_p_values <- coeffs[, "Pr(>|t|)"]
      
      all_coeff_significant <- all(coeff_p_values < 0.01)
      normality_pass <- shapiro_p > 0.01
      autocorrelation_pass <- dw_p > 0.01
      vif_pass <- all(vif_values < 5)
      heteroscedasticity_pass <- bp_p > 0.01
      
      final_flag <- all(all_coeff_significant, normality_pass, autocorrelation_pass, vif_pass, heteroscedasticity_pass)
      
      # Equation
      equation <- paste(
        ifelse("(Intercept)" %in% rownames(coeffs), paste(dep_var," =", round(coeffs["(Intercept)", "Estimate"], 4)),paste0(dep_var," =")),
        paste(
          round(coeffs[rownames(coeffs) != "(Intercept)", "Estimate"], 4),
          "*",
          rownames(coeffs)[rownames(coeffs) != "(Intercept)"],
          collapse = " + "
        )
      )
      
      # Append results
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

# View results
View(results)
# Install the openxlsx package if not already installed
#install.packages("openxlsx")

# Load the package
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# Add sheets and write data
addWorksheet(wb, "Train Data")
writeData(wb, "Train Data", train_data)

addWorksheet(wb, "Test Data")
writeData(wb, "Test Data", test_data)

addWorksheet(wb, "Model")
writeData(wb, "Model", results)

# Save the Excel file
saveWorkbook(wb, "SMVM_Modelling_iterations.xlsx", overwrite = TRUE)
