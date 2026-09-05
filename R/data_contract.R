clinical_schema <- function() {
  data.frame(
    variable = c(
      "patient_id", "mortality", "amputation_status",
      "amputation_type", "age", "extent", "organism",
      "blood_culture", "diabetes", "ckd", "cad", "smoker",
      "special_medication", "shock", "fever_38", "wcc",
      "hemoglobin", "platelet_count", "inr", "crp",
      "abnormal_liver_function", "albumin", "creatinine",
      "glucose", "ph"
    ),
    type = c(
      "character", "binary", "binary", "amputation_type",
      "numeric", "extent", "character", rep("binary", 5),
      "character", rep("binary", 2), rep("numeric", 5),
      "binary", rep("numeric", 4)
    ),
    required = c(TRUE, TRUE, rep(FALSE, 23)),
    stringsAsFactors = FALSE
  )
}

.assert_binary <- function(values, variable) {
  observed <- unique(values[!is.na(values)])
  if (!all(observed %in% c(0, 1))) {
    stop(variable, " must contain only 0, 1, or NA.", call. = FALSE)
  }
}

validate_clinical_data <- function(data, strict = FALSE) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame.", call. = FALSE)
  }

  schema <- clinical_schema()
  required <- schema$variable[schema$required]
  missing_required <- setdiff(required, names(data))
  if (length(missing_required) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyDuplicated(data$patient_id)) {
    stop("patient_id must be unique.", call. = FALSE)
  }
  if (any(is.na(data$patient_id)) || any(data$patient_id == "")) {
    stop("patient_id cannot be missing or empty.", call. = FALSE)
  }

  binary_variables <- intersect(
    schema$variable[schema$type == "binary"],
    names(data)
  )
  for (variable in binary_variables) {
    .assert_binary(data[[variable]], variable)
  }

  numeric_variables <- intersect(
    schema$variable[schema$type == "numeric"],
    names(data)
  )
  for (variable in numeric_variables) {
    if (!is.numeric(data[[variable]])) {
      stop(variable, " must be numeric.", call. = FALSE)
    }
    if (any(!is.finite(data[[variable]][!is.na(data[[variable]])]))) {
      stop(variable, " contains a nonfinite value.", call. = FALSE)
    }
  }

  if ("extent" %in% names(data)) {
    allowed_extent <- c("below_wrist", "wrist_to_elbow", "above_elbow")
    observed_extent <- unique(as.character(data$extent[!is.na(data$extent)]))
    invalid_extent <- setdiff(observed_extent, allowed_extent)
    if (length(invalid_extent) > 0L) {
      stop(
        "Invalid extent values: ",
        paste(invalid_extent, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if ("amputation_type" %in% names(data)) {
    allowed_type <- c("no_amputation", "digit_ray", "higher_level")
    observed_type <- unique(
      as.character(data$amputation_type[!is.na(data$amputation_type)])
    )
    invalid_type <- setdiff(observed_type, allowed_type)
    if (length(invalid_type) > 0L) {
      stop(
        "Invalid amputation_type values: ",
        paste(invalid_type, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (strict) {
    unknown <- setdiff(names(data), schema$variable)
    if (length(unknown) > 0L) {
      stop(
        "Columns not present in the data contract: ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
  }

  invisible(data)
}

prepare_analysis_data <- function(data) {
  validate_clinical_data(data)
  prepared <- data
  if ("extent" %in% names(prepared)) {
    prepared$extent_low <- as.integer(prepared$extent == "below_wrist")
  }
  prepared
}
