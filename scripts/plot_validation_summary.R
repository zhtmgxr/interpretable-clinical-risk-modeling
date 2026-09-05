arguments <- commandArgs(trailingOnly = TRUE)
input <- if (length(arguments) >= 1L) {
  arguments[[1L]]
} else {
  "results/four_variable_model_comparison.csv"
}
output <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  "figures/model_validation_summary.png"
}

results <- utils::read.csv(input, stringsAsFactors = FALSE)
labels <- gsub("_", " ", results$model)
metrics <- rbind(
  results$cv_sensitivity,
  results$cv_specificity,
  results$cv_accuracy
)
colnames(metrics) <- labels
rownames(metrics) <- c("Sensitivity", "Specificity", "Accuracy")

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
grDevices::png(output, width = 1800, height = 1050, res = 180)
graphics::par(mar = c(10, 5, 3, 1))
graphics::barplot(
  metrics,
  beside = TRUE,
  col = c("#B94A48", "#3E6B89", "#687864"),
  border = NA,
  ylim = c(0, 100),
  las = 2,
  ylab = "Repeated cross validation performance (%)",
  main = "Overall accuracy can conceal failure on mortality cases"
)
graphics::abline(h = seq(0, 100, by = 20), col = "#DDDDDD", lty = 3)
graphics::legend(
  "topright",
  legend = rownames(metrics),
  fill = c("#B94A48", "#3E6B89", "#687864"),
  border = NA,
  bty = "n"
)
grDevices::dev.off()
cat("Wrote", output, "\n")
