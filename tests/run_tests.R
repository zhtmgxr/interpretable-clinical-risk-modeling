source("R/data_contract.R")
source("R/metrics.R")
source("R/models.R")
source("R/resampling.R")
source("R/threshold_model.R")
source("R/scorecard.R")

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

expect_equal <- function(actual, expected, tolerance = 1e-8, message = NULL) {
  equal <- isTRUE(all.equal(actual, expected, tolerance = tolerance))
  if (!equal) {
    if (is.null(message)) {
      message <- paste("Expected", expected, "but received", actual)
    }
    stop(message, call. = FALSE)
  }
}

make_synthetic_data <- function(n = 120L, seed = 11L) {
  set.seed(seed)
  extent <- sample(
    c("below_wrist", "wrist_to_elbow", "above_elbow"),
    n,
    replace = TRUE,
    prob = c(0.35, 0.4, 0.25)
  )
  albumin <- stats::rnorm(n, 37, 6)
  creatinine <- pmax(30, stats::rnorm(n, 135, 65))
  platelet_count <- pmax(25, stats::rnorm(n, 180, 70))
  linear_predictor <-
    -2.1 +
    0.9 * (extent != "below_wrist") +
    0.8 * (albumin < 35) +
    1.0 * (creatinine > 150) +
    0.6 * (platelet_count < 150)
  mortality <- stats::rbinom(n, 1, stats::plogis(linear_predictor))
  if (sum(mortality == 1) < 10L) {
    mortality[seq_len(10L)] <- 1L
  }
  if (sum(mortality == 0) < 10L) {
    mortality[seq_len(10L)] <- 0L
  }

  data.frame(
    patient_id = sprintf("synthetic_%03d", seq_len(n)),
    mortality = mortality,
    extent = extent,
    albumin = albumin,
    creatinine = creatinine,
    platelet_count = platelet_count,
    stringsAsFactors = FALSE
  )
}

cat("Testing data validation\n")
synthetic <- make_synthetic_data()
validate_clinical_data(synthetic)
prepared <- prepare_analysis_data(synthetic)
expect_true("extent_low" %in% names(prepared), "extent_low was not created.")
expect_true(
  all(prepared$extent_low %in% c(0, 1)),
  "extent_low must be binary."
)

cat("Testing binary metrics\n")
metrics <- binary_metrics(
  truth = c(1, 1, 0, 0),
  prediction = c(1, 0, 1, 0)
)
expect_equal(unname(metrics["accuracy"]), 0.5)
expect_equal(unname(metrics["sensitivity"]), 0.5)
expect_equal(unname(metrics["specificity"]), 0.5)

cat("Testing balanced cutoff selection\n")
cutoff <- select_balanced_threshold(
  truth = c(0, 0, 1, 1),
  probability = c(0.1, 0.3, 0.7, 0.9)
)
expect_true(cutoff$threshold > 0.3, "Cutoff should separate the two classes.")
expect_true(cutoff$threshold < 0.7, "Cutoff should separate the two classes.")

cat("Testing logistic model and diagnostics\n")
logistic <- fit_logistic_model(
  prepared,
  outcome = "mortality",
  predictors = c("extent_low", "albumin", "creatinine", "platelet_count")
)
logistic_probability <- predict(logistic, prepared)
expect_true(length(logistic_probability) == nrow(prepared), "Prediction length mismatch.")
diagnostics <- diagnose_logistic_fit(logistic)
expect_true(is.logical(diagnostics$converged), "Convergence flag must be logical.")

cat("Testing repeated stratified validation\n")
validation <- repeated_stratified_cv(
  prepared,
  mortality ~ extent_low + albumin + creatinine + platelet_count,
  v = 5,
  repeats = 3,
  seed = 19,
  method = "logistic"
)
expect_true(nrow(validation$per_repeat) == 3L, "Expected three validation repeats.")
expect_true(
  all(validation$summary$mean[1:4] >= 0 & validation$summary$mean[1:4] <= 1),
  "Validation metrics must be probabilities."
)

cat("Testing threshold model\n")
threshold_model <- fit_threshold_model(
  prepared,
  outcome = "mortality",
  predictors = c("albumin", "creatinine", "platelet_count"),
  candidate_probabilities = seq(0.2, 0.8, by = 0.1)
)
threshold_probability <- predict(threshold_model, prepared)
expect_true(length(threshold_probability) == nrow(prepared), "Threshold prediction mismatch.")
expect_true(length(threshold_model$breakpoints) == 3L, "Expected three breakpoints.")

cat("Testing scorecard derivation\n")
risk_cutoff <- select_balanced_threshold(
  prepared$mortality,
  threshold_probability
)$threshold
scorecard <- derive_integer_scorecard(threshold_model, risk_cutoff)
scored <- apply_integer_scorecard(scorecard, prepared)
expect_true(nrow(scored) == nrow(prepared), "Scorecard row count mismatch.")
expect_true(all(scored$total_score >= 0), "Score totals must be nonnegative.")

cat("Testing four variable study scorecard\n")
examples <- data.frame(
  extent = c("above_elbow", "below_wrist"),
  platelet_count = c(100, 300),
  albumin = c(40, 30),
  creatinine = c(200, 100),
  stringsAsFactors = FALSE
)
study_scores <- score_upper_limb_nf(examples)
expect_equal(study_scores$total_score, c(11, 0))
expect_equal(study_scores$predicted_high_risk, c(1, 0))

cat("All tests passed.\n")
