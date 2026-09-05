fit_logistic_model <- function(
  data,
  outcome,
  predictors,
  selection = c("none", "backward_aic")
) {
  selection <- match.arg(selection)
  required <- c(outcome, predictors)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing model columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  complete <- stats::complete.cases(data[, required, drop = FALSE])
  model_data <- data[complete, required, drop = FALSE]
  if (nrow(model_data) == 0L) {
    stop("No complete observations are available for the model.", call. = FALSE)
  }
  .assert_binary(model_data[[outcome]], outcome)
  if (length(unique(model_data[[outcome]])) != 2L) {
    stop("The outcome must contain both classes.", call. = FALSE)
  }

  formula <- stats::reformulate(predictors, response = outcome)
  initial <- stats::glm(
    formula,
    data = model_data,
    family = stats::binomial()
  )
  fitted <- if (selection == "backward_aic") {
    suppressWarnings(stats::step(initial, direction = "backward", trace = 0))
  } else {
    initial
  }

  structure(
    list(
      model = fitted,
      initial_model = initial,
      outcome = outcome,
      predictors = predictors,
      selection = selection,
      rows_used = which(complete)
    ),
    class = "clinical_logistic_model"
  )
}

predict.clinical_logistic_model <- function(object, newdata, ...) {
  stats::predict(object$model, newdata = newdata, type = "response", ...)
}

diagnose_logistic_fit <- function(model, coefficient_limit = 15, se_limit = 10) {
  if (!inherits(model, "clinical_logistic_model")) {
    stop("model must be a clinical_logistic_model.", call. = FALSE)
  }
  summary_table <- summary(model$model)$coefficients
  coefficients <- summary_table[, "Estimate"]
  standard_errors <- summary_table[, "Std. Error"]

  list(
    converged = isTRUE(model$model$converged),
    extreme_coefficients = names(coefficients)[
      is.finite(coefficients) & abs(coefficients) > coefficient_limit
    ],
    large_standard_errors = names(standard_errors)[
      is.finite(standard_errors) & standard_errors > se_limit
    ],
    rank_deficient = model$model$rank < length(stats::coef(model$model)),
    observations = stats::nobs(model$model),
    events = sum(model$model$y == 1)
  )
}

fit_multinomial_model <- function(data, outcome, predictors) {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("Install the nnet package to fit a multinomial model.", call. = FALSE)
  }
  required <- c(outcome, predictors)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing model columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  model_data <- data[
    stats::complete.cases(data[, required, drop = FALSE]),
    required,
    drop = FALSE
  ]
  model_data[[outcome]] <- factor(model_data[[outcome]])
  if (nlevels(model_data[[outcome]]) < 3L) {
    stop("The multinomial outcome must contain at least three classes.", call. = FALSE)
  }

  formula <- stats::reformulate(predictors, response = outcome)
  nnet::multinom(formula, data = model_data, trace = FALSE)
}
