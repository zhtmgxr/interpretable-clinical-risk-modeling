.candidate_breakpoints <- function(x, probabilities) {
  if (!is.numeric(x)) {
    stop("Threshold predictors must be numeric.", call. = FALSE)
  }
  values <- unique(as.numeric(stats::quantile(
    x,
    probs = probabilities,
    na.rm = TRUE,
    names = FALSE,
    type = 2
  )))
  values[values > min(x, na.rm = TRUE) & values < max(x, na.rm = TRUE)]
}

.threshold_design <- function(data, predictors, breakpoints, outcome = NULL) {
  design <- as.data.frame(lapply(predictors, function(variable) {
    as.integer(data[[variable]] > breakpoints[[variable]])
  }))
  names(design) <- paste0("step_", make.names(predictors, unique = TRUE))
  if (!is.null(outcome)) {
    design[[outcome]] <- data[[outcome]]
  }
  design
}

.fit_threshold_at <- function(data, outcome, predictors, breakpoints) {
  design <- .threshold_design(data, predictors, breakpoints, outcome)
  formula <- stats::reformulate(
    paste0("step_", make.names(predictors, unique = TRUE)),
    response = outcome
  )
  stats::lm(formula, data = design)
}

fit_threshold_model <- function(
  data,
  outcome,
  predictors,
  candidate_probabilities = seq(0.15, 0.85, by = 0.05),
  max_iterations = 20L
) {
  required <- c(outcome, predictors)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing threshold model columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  model_data <- data[
    stats::complete.cases(data[, required, drop = FALSE]),
    required,
    drop = FALSE
  ]
  .assert_binary(model_data[[outcome]], outcome)
  if (length(unique(model_data[[outcome]])) != 2L) {
    stop("The outcome must contain both classes.", call. = FALSE)
  }
  if (any(!vapply(model_data[predictors], is.numeric, logical(1)))) {
    stop("Every threshold predictor must be numeric.", call. = FALSE)
  }

  candidates <- lapply(
    model_data[predictors],
    .candidate_breakpoints,
    probabilities = candidate_probabilities
  )
  empty_candidates <- names(candidates)[lengths(candidates) == 0L]
  if (length(empty_candidates) > 0L) {
    stop(
      "No interior breakpoint candidates for: ",
      paste(empty_candidates, collapse = ", "),
      call. = FALSE
    )
  }

  breakpoints <- vapply(candidates, stats::median, numeric(1))
  names(breakpoints) <- predictors
  history <- list()
  converged <- FALSE

  for (iteration in seq_len(max_iterations)) {
    previous <- breakpoints
    for (variable in predictors) {
      candidate_aic <- vapply(candidates[[variable]], function(value) {
        trial <- breakpoints
        trial[[variable]] <- value
        stats::AIC(.fit_threshold_at(model_data, outcome, predictors, trial))
      }, numeric(1))
      breakpoints[[variable]] <- candidates[[variable]][which.min(candidate_aic)]
    }
    fitted <- .fit_threshold_at(model_data, outcome, predictors, breakpoints)
    history[[iteration]] <- data.frame(
      iteration = iteration,
      aic = stats::AIC(fitted),
      t(breakpoints),
      check.names = FALSE
    )
    if (isTRUE(all.equal(previous, breakpoints, tolerance = 0))) {
      converged <- TRUE
      break
    }
  }

  fitted <- .fit_threshold_at(model_data, outcome, predictors, breakpoints)
  structure(
    list(
      model = fitted,
      outcome = outcome,
      predictors = predictors,
      breakpoints = breakpoints,
      candidates = candidates,
      history = do.call(rbind, history),
      converged = converged,
      rows_used = which(stats::complete.cases(data[, required, drop = FALSE]))
    ),
    class = "clinical_threshold_model"
  )
}

predict.clinical_threshold_model <- function(object, newdata, ...) {
  design <- .threshold_design(
    newdata,
    object$predictors,
    object$breakpoints
  )
  as.numeric(stats::predict(object$model, newdata = design, ...))
}

evaluate_threshold_model_cv <- function(
  data,
  outcome,
  predictors,
  v = 5L,
  repeats = 50L,
  seed = 2026L,
  candidate_probabilities = seq(0.15, 0.85, by = 0.1)
) {
  required <- c(outcome, predictors)
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
      fitted <- fit_threshold_model(
        training,
        outcome = outcome,
        predictors = predictors,
        candidate_probabilities = candidate_probabilities
      )
      training_probability <- predict(fitted, training)
      threshold <- select_balanced_threshold(
        training[[outcome]],
        training_probability
      )$threshold
      testing_probability <- predict(fitted, testing)
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
      candidate_probabilities = candidate_probabilities
    )
  )
}
