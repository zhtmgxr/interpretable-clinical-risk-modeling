derive_integer_scorecard <- function(model, risk_cutoff, scale = 10) {
  if (!inherits(model, "clinical_threshold_model")) {
    stop("model must be a clinical_threshold_model.", call. = FALSE)
  }
  if (!is.numeric(risk_cutoff) || length(risk_cutoff) != 1L) {
    stop("risk_cutoff must be one numeric value.", call. = FALSE)
  }
  coefficient_names <- paste0(
    "step_",
    make.names(model$predictors, unique = TRUE)
  )
  coefficients <- stats::coef(model$model)[coefficient_names]
  names(coefficients) <- model$predictors
  if (any(!is.finite(coefficients))) {
    stop("Cannot derive a scorecard from nonfinite coefficients.", call. = FALSE)
  }

  adjusted_intercept <- unname(stats::coef(model$model)["(Intercept)"])
  negative <- coefficients < 0
  adjusted_intercept <- adjusted_intercept + sum(coefficients[negative])

  criteria <- data.frame(
    variable = model$predictors,
    comparison = ifelse(negative, "less_or_equal", "greater_than"),
    breakpoint = unname(model$breakpoints[model$predictors]),
    raw_weight = abs(unname(coefficients)),
    points = round(scale * abs(unname(coefficients))),
    stringsAsFactors = FALSE
  )
  continuous_cutoff <- scale * (risk_cutoff - adjusted_intercept)

  structure(
    list(
      criteria = criteria,
      risk_cutoff = risk_cutoff,
      scale = scale,
      adjusted_intercept = adjusted_intercept,
      continuous_score_cutoff = continuous_cutoff,
      integer_score_cutoff = ceiling(continuous_cutoff)
    ),
    class = "clinical_scorecard"
  )
}

apply_integer_scorecard <- function(scorecard, newdata) {
  if (!inherits(scorecard, "clinical_scorecard")) {
    stop("scorecard must be a clinical_scorecard.", call. = FALSE)
  }
  criteria <- scorecard$criteria
  missing_columns <- setdiff(criteria$variable, names(newdata))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing scorecard columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  contributions <- vapply(seq_len(nrow(criteria)), function(index) {
    values <- newdata[[criteria$variable[index]]]
    active <- if (criteria$comparison[index] == "greater_than") {
      values > criteria$breakpoint[index]
    } else {
      values <= criteria$breakpoint[index]
    }
    as.integer(active) * criteria$points[index]
  }, numeric(nrow(newdata)))
  if (is.null(dim(contributions))) {
    contributions <- matrix(contributions, nrow = nrow(newdata))
  }
  colnames(contributions) <- paste0(criteria$variable, "_points")
  total <- rowSums(contributions)

  data.frame(
    contributions,
    total_score = total,
    predicted_high_risk = as.integer(
      total >= scorecard$integer_score_cutoff
    ),
    row.names = NULL
  )
}

upper_limb_nf_scorecard <- function() {
  structure(
    list(
      criteria = data.frame(
        variable = c(
          "extent", "platelet_count", "albumin", "creatinine"
        ),
        comparison = c(
          "higher_extent", "less_or_equal", "greater_than", "greater_than"
        ),
        breakpoint = c(NA, 215, 34, 148),
        points = c(3, 1, 3, 4),
        stringsAsFactors = FALSE
      ),
      integer_score_cutoff = 8
    ),
    class = "upper_limb_nf_scorecard"
  )
}

score_upper_limb_nf <- function(newdata, scorecard = upper_limb_nf_scorecard()) {
  required <- c("extent", "platelet_count", "albumin", "creatinine")
  missing_columns <- setdiff(required, names(newdata))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing scorecard columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  extent_points <- ifelse(
    newdata$extent %in% c("wrist_to_elbow", "above_elbow"),
    3,
    0
  )
  platelet_points <- ifelse(newdata$platelet_count <= 215, 1, 0)
  albumin_points <- ifelse(newdata$albumin > 34, 3, 0)
  creatinine_points <- ifelse(newdata$creatinine > 148, 4, 0)
  total <- extent_points + platelet_points + albumin_points + creatinine_points

  data.frame(
    extent_points = extent_points,
    platelet_points = platelet_points,
    albumin_points = albumin_points,
    creatinine_points = creatinine_points,
    total_score = total,
    predicted_high_risk = as.integer(total >= scorecard$integer_score_cutoff),
    row.names = NULL
  )
}
