# Statistical design

## Outcome scarcity changes the objective

The cohort contains forty two patients and eight deaths. In this regime, a classifier can obtain about eighty one percent accuracy by predicting survival for everyone. Model comparison must therefore include mortality sensitivity and specificity rather than accuracy alone.

## Repeated stratified validation

Each validation repeat assigns mortality cases and survival cases to five folds separately. This preserves both classes as far as the sample permits. The split is repeated fifty times with deterministic seeds, and each metric is summarized across repeats.

The classification threshold is estimated from the training partition within each fold. It is never selected from the held out outcomes. This ordering is essential because choosing a threshold before cross validation would leak outcome information into evaluation.

Repeated cross validation reduces dependence on one fortunate partition. It does not make the observations independent, create new information, or remove the need for external validation.

## Model selection and inference

Backward AIC selection is available because it was part of the analytical workflow. The resulting coefficient table is useful for prediction and exploratory association analysis, but ordinary confidence intervals and P values do not account for the selection event. They should not be read as confirmatory evidence.

Rare categories can also cause complete or quasi complete separation. The diagnostic function reports failed convergence, extreme coefficients, large standard errors, and rank deficiency so that an apparently decisive estimate is not accepted silently.

## Thresholded linear model

For numeric predictors, the model is

```text
Y = beta 0 + sum of beta j × indicator(X j > psi j) + error
```

The implementation searches candidate breakpoints by coordinate descent. Candidate values come from interior empirical quantiles, which prevents thresholds at isolated minima or maxima. Given a set of breakpoints, coefficient estimation is ordinary least squares on binary threshold indicators.

This search is deliberately explicit. It makes the nonconvex breakpoint problem visible, records the optimization path, and allows the entire search to be repeated inside each training fold.

## From coefficients to a scorecard

A negative threshold coefficient is rewritten as a positive weight on the complementary condition. The adjusted coefficients are scaled and rounded to integer points. The transformed intercept is used to update the cutoff.

Both the continuous cutoff and final integer cutoff are returned. This matters because coefficient rounding changes the decision boundary and should not be hidden.

## Missing data

The provided analysis uses complete cases for the variables in each model. For a new study, missingness should be described explicitly. Any imputation must be estimated within the training fold to avoid leakage. Multiple imputation would also need to propagate uncertainty across the breakpoint and score derivation process.

## Intended interpretation

The implementation is suitable for methodological review and independent validation. It is not calibrated for clinical deployment. Any future clinical use would require external data, prospective validation, prespecified endpoints, decision curve analysis, calibration assessment, and governance appropriate to a medical decision support tool.
