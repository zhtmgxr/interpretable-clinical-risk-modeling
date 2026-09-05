make_stratified_folds <- function(truth, v = 5L, seed = NULL) {
  .validate_binary_vectors(truth)
  if (v < 2L) {
    stop("v must be at least 2.", call. = FALSE)
  }
  class_counts <- table(truth)
  if (any(class_counts < v)) {
    stop("Each outcome class must contain at least v observations.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  folds <- integer(length(truth))
  for (class_value in c(0, 1)) {
    indices <- sample(which(truth == class_value))
    folds[indices] <- rep(seq_len(v), length.out = length(indices))
  }
  folds
}

.fit_resampling_model <- function(data, formula, method) {
  if (method == "logistic") {
    return(stats::glm(formula, data = data, family = stats::binomial()))
  }
  stats::lm(formula, data = data)
}

.predict_resampling_model <- function(model, data, method) {
  if (method == "logistic") {
    return(as.numeric(stats::predict(model, newdata = data, type = "response")))
  }
  as.numeric(stats::predict(model, newdata = data))
}

.summarize_repeated_metrics <- function(per_repeat) {
  metric_names <- c(
    "accuracy", "sensitivity", "specificity", "balanced_accuracy",
    "precision", "npv", "fpr", "fdr", "fnr", "f1"
  )
  means <- vapply(metric_names, function(metric) {
    values <- per_repeat[[metric]]
    if (all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
  }, numeric(1))
  standard_deviations <- vapply(metric_names, function(metric) {
    values <- per_repeat[[metric]]
    if (sum(!is.na(values)) < 2L) NA_real_ else stats::sd(values, na.rm = TRUE)
  }, numeric(1))
  data.frame(
    metric = metric_names,
    mean = unname(means),
    standard_deviation = unname(standard_deviations),
    row.names = NULL
  )
}

repeated_stratified_cv <- function(
  data,
  formula,
  v = 5L,
  repeats = 50L,
  seed = 2026L,
  method = c("logistic", "linear_probability"),
  threshold_method = c("balanced", "fixed"),
  fixed_threshold = 0.5
) {
  method <- match.arg(method)
  threshold_method <- match.arg(threshold_method)
  outcome <- as.character(formula[[2L]])
  required <- all.vars(formula)
  model_data <- data[
    stats::complete.cases(data[, required, drop = FALSE]),
    required,
    drop = FALSE
  ]
  .assert_binary(model_data[[outcome]], outcome)

  repeat_metrics <- vector("list", repeats)
  prediction_records <- vector("list", repeats)

  for (repeat_id in seq_len(repeats)) {
    fold_id <- make_stratified_folds(
      model_data[[outcome]],
      v = v,
      seed = seed + repeat_id - 1L
    )
    truth <- probability <- prediction <- rep(NA_real_, nrow(model_data))
    thresholds <- numeric(v)

    for (fold in seq_len(v)) {
      training <- model_data[fold_id != fold, , drop = FALSE]
      testing <- model_data[fold_id == fold, , drop = FALSE]
      fitted <- .fit_resampling_model(training, formula, method)
      training_probability <- .predict_resampling_model(
        fitted, training, method
      )
      threshold <- if (threshold_method == "balanced") {
        select_balanced_threshold(
          training[[outcome]],
          training_probability
        )$threshold
      } else {
        fixed_threshold
      }
      testing_probability <- .predict_resampling_model(fitted, testing, method)
      positions <- which(fold_id == fold)
      truth[positions] <- testing[[outcome]]
      probability[positions] <- testing_probability
      prediction[positions] <- as.integer(testing_probability >= threshold)
      thresholds[fold] <- threshold
    }

    repeat_metrics[[repeat_id]] <- binary_metrics(
      truth = truth,
      prediction = prediction
    )
    prediction_records[[repeat_id]] <- data.frame(
      repeat_id = repeat_id,
      row_id = seq_len(nrow(model_data)),
      truth = truth,
      probability = probability,
      prediction = prediction,
      fold_id = fold_id,
      fold_threshold = thresholds[fold_id]
    )
  }

  per_repeat <- as.data.frame(do.call(rbind, repeat_metrics))
  rownames(per_repeat) <- NULL
  list(
    summary = .summarize_repeated_metrics(per_repeat),
    per_repeat = per_repeat,
    predictions = do.call(rbind, prediction_records),
    configuration = list(
      folds = v,
      repeats = repeats,
      seed = seed,
      method = method,
      threshold_method = threshold_method
    )
  )
}
