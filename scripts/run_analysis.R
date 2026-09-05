arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop(
    "Usage: Rscript scripts/run_analysis.R path/to/anonymized_clinical_data.csv",
    call. = FALSE
  )
}

source("R/data_contract.R")
source("R/metrics.R")
source("R/models.R")
source("R/resampling.R")
source("R/threshold_model.R")
source("R/scorecard.R")

data <- utils::read.csv(arguments[[1L]], stringsAsFactors = FALSE)
data <- prepare_analysis_data(data)

mortality_predictors <- c(
  "extent_low", "albumin", "creatinine", "platelet_count"
)
missing_predictors <- setdiff(mortality_predictors, names(data))
if (length(missing_predictors) > 0L) {
  stop(
    "The mortality analysis requires: ",
    paste(missing_predictors, collapse = ", "),
    call. = FALSE
  )
}

mortality_formula <- stats::reformulate(
  mortality_predictors,
  response = "mortality"
)

cat("\nRepeated stratified validation for logistic regression\n")
logistic_validation <- repeated_stratified_cv(
  data,
  mortality_formula,
  v = 5,
  repeats = 50,
  seed = 2026,
  method = "logistic",
  threshold_method = "balanced"
)
print(logistic_validation$summary, row.names = FALSE)

cat("\nRepeated stratified validation for threshold model\n")
threshold_validation <- evaluate_threshold_model_cv(
  data,
  outcome = "mortality",
  predictors = mortality_predictors,
  v = 5,
  repeats = 50,
  seed = 2026
)
print(threshold_validation$summary, row.names = FALSE)

cat("\nFinal logistic model\n")
logistic_model <- fit_logistic_model(
  data,
  outcome = "mortality",
  predictors = mortality_predictors,
  selection = "none"
)
print(summary(logistic_model$model))
print(diagnose_logistic_fit(logistic_model))

cat("\nThreshold model and integer scorecard\n")
threshold_model <- fit_threshold_model(
  data,
  outcome = "mortality",
  predictors = mortality_predictors
)
training_risk <- predict(threshold_model, data)
risk_cutoff <- select_balanced_threshold(data$mortality, training_risk)$threshold
scorecard <- derive_integer_scorecard(threshold_model, risk_cutoff)
print(threshold_model$breakpoints)
print(scorecard$criteria)
cat("Integer score cutoff:", scorecard$integer_score_cutoff, "\n")
