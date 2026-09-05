.safe_ratio <- function(numerator, denominator) {
  if (denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}

.validate_binary_vectors <- function(truth, probability = NULL) {
  if (!is.numeric(truth) || !all(truth %in% c(0, 1))) {
    stop("truth must be a numeric vector containing only 0 and 1.", call. = FALSE)
  }
  if (!is.null(probability)) {
    if (!is.numeric(probability) || length(probability) != length(truth)) {
      stop("probability must be numeric and match truth in length.", call. = FALSE)
    }
    if (any(!is.finite(probability))) {
      stop("probability must contain only finite values.", call. = FALSE)
    }
  }
}

binary_metrics <- function(
  truth,
  probability = NULL,
  prediction = NULL,
  threshold = 0.5
) {
  .validate_binary_vectors(truth, probability)

  if (is.null(prediction)) {
    if (is.null(probability)) {
      stop("Provide either probability or prediction.", call. = FALSE)
    }
    prediction <- as.integer(probability >= threshold)
  }
  if (length(prediction) != length(truth) || !all(prediction %in% c(0, 1))) {
    stop("prediction must match truth and contain only 0 and 1.", call. = FALSE)
  }

  tp <- sum(prediction == 1 & truth == 1)
  tn <- sum(prediction == 0 & truth == 0)
  fp <- sum(prediction == 1 & truth == 0)
  fn <- sum(prediction == 0 & truth == 1)

  sensitivity <- .safe_ratio(tp, tp + fn)
  specificity <- .safe_ratio(tn, tn + fp)
  precision <- .safe_ratio(tp, tp + fp)
  npv <- .safe_ratio(tn, tn + fn)

  c(
    accuracy = .safe_ratio(tp + tn, length(truth)),
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    precision = precision,
    npv = npv,
    fpr = .safe_ratio(fp, fp + tn),
    fdr = .safe_ratio(fp, fp + tp),
    fnr = .safe_ratio(fn, fn + tp),
    f1 = .safe_ratio(2 * tp, 2 * tp + fp + fn),
    tp = tp,
    tn = tn,
    fp = fp,
    fn = fn
  )
}

select_balanced_threshold <- function(truth, probability) {
  .validate_binary_vectors(truth, probability)
  unique_values <- sort(unique(probability))
  if (length(unique_values) == 1L) {
    candidates <- unique_values
  } else {
    candidates <- c(
      unique_values[1L],
      (unique_values[-1L] + unique_values[-length(unique_values)]) / 2,
      unique_values[length(unique_values)]
    )
  }

  evaluations <- lapply(candidates, function(candidate) {
    metrics <- binary_metrics(
      truth = truth,
      probability = probability,
      threshold = candidate
    )
    data.frame(
      threshold = candidate,
      sensitivity = unname(metrics["sensitivity"]),
      specificity = unname(metrics["specificity"]),
      balanced_accuracy = unname(metrics["balanced_accuracy"])
    )
  })
  table <- do.call(rbind, evaluations)
  table$difference <- abs(table$sensitivity - table$specificity)
  table$difference[!is.finite(table$difference)] <- Inf
  table$balanced_accuracy[!is.finite(table$balanced_accuracy)] <- -Inf
  ordering <- order(
    table$difference,
    -table$balanced_accuracy,
    table$threshold
  )
  best <- table[ordering[1L], , drop = FALSE]

  list(
    threshold = best$threshold[[1L]],
    metrics = best,
    candidates = table
  )
}
