# Interpretable Clinical Risk Modeling

[![R tests](https://github.com/zhtmgxr/interpretable-clinical-risk-modeling/actions/workflows/tests.yml/badge.svg)](https://github.com/zhtmgxr/interpretable-clinical-risk-modeling/actions/workflows/tests.yml)

An R implementation of interpretable risk modeling for small clinical cohorts, developed around a twenty year retrospective study of upper limb necrotizing fasciitis.

## Research question

Rare disease modeling creates a difficult statistical regime. The cohort is too small for flexible models to learn stable decision boundaries, the positive outcome is uncommon, and clinicians need a decision rule they can inspect under time pressure.

This project asks a practical question: how can a model preserve enough structure to detect mortality risk while remaining simple enough to audit and calculate?

The analysis combines classical regression, threshold based modeling, repeated stratified validation, and an additive scorecard. Its main lesson is that overall accuracy is not sufficient when the outcome is imbalanced. A model can appear accurate while failing almost every patient in the clinically important minority class.

## Study setting

1. Forty two consecutive patients collected over twenty years
2. Thirty four survivors and eight deaths
3. Mortality, amputation status, and amputation type as outcomes
4. Demographic, comorbidity, infection, presentation, and laboratory covariates
5. R 4.4.1 for the statistical analysis

Patient level records are not included. The repository contains a data contract, analysis functions, aggregate findings, and synthetic fixtures for testing.

## Evaluation insight

The most informative comparison uses four predictors and stratified five fold cross validation repeated fifty times.

| Model | Sensitivity | Specificity | Accuracy |
| --- | ---: | ---: | ---: |
| Linear probability model | 81.0% | 73.3% | 74.4% |
| Logistic regression | 70.0% | 77.6% | 76.3% |
| Thresholded linear model | 60.6% | 72.1% | 69.8% |
| Classification tree | 2.0% | 97.2% | 79.2% |
| Random forest | 11.4% | 92.5% | 77.0% |

The tree and random forest achieved competitive overall accuracy by predicting the majority survival class. Their mortality sensitivity shows that they were not clinically reliable. Logistic regression produced the strongest balance among the directly interpretable models. The thresholded model gave up some predictive performance in exchange for a score that can be evaluated without software.

The complete aggregate comparison is stored in [`results/four_variable_model_comparison.csv`](results/four_variable_model_comparison.csv).

## Statistical workflow

1. Validate the clinical schema and outcome encoding.
2. Summarize continuous and categorical variables by outcome.
3. Fit mortality logistic regression with optional backward AIC selection.
4. Diagnose nonconvergence and coefficient inflation caused by separation.
5. Fit multinomial logistic regression for amputation type when the optional dependency is available.
6. Estimate threshold effects with a coordinate search over candidate breakpoints.
7. Choose a classification cutoff by balancing sensitivity and specificity on training data only.
8. Evaluate every decision rule with repeated stratified cross validation.
9. Translate threshold coefficients into positive integer points and an explicit cutoff.

The validation code keeps cutoff selection inside each training fold. This prevents the held out fold from influencing the decision threshold.

## Repository structure

```text
R/
  data_contract.R             Clinical schema and validation
  metrics.R                   Confusion matrix metrics and cutoff selection
  models.R                    Logistic and multinomial regression
  resampling.R                Repeated stratified cross validation
  threshold_model.R           Breakpoint search and threshold regression
  scorecard.R                 Integer score derivation and application
scripts/
  run_analysis.R              End to end analysis entry point
  plot_validation_summary.R   Aggregate model comparison figure
data/
  data_dictionary.csv         Expected variables and encodings
results/
  four_variable_model_comparison.csv
docs/
  statistical_design.md       Assumptions, validation, and limitations
tests/
  run_tests.R                 Dependency free unit and integration tests
```

## Run the analysis

R 4.4.1 or a later compatible release is recommended.

```bash
Rscript tests/run_tests.R
Rscript scripts/run_analysis.R path/to/anonymized_clinical_data.csv
```

The analysis script expects the columns defined in [`data/data_dictionary.csv`](data/data_dictionary.csv). It prints repeated validation summaries, fits the final mortality model, derives an integer scorecard, and writes no patient level output by default.

## Scorecard design

For numeric covariates, the threshold model has the form

```text
predicted risk = intercept + sum of coefficient j × indicator(x j > breakpoint j)
```

Each negative coefficient is rewritten as a positive contribution on the complementary condition. The transformed coefficients are scaled and rounded into integer points. The classification cutoff is adjusted for the transformed intercept and rounded conservatively.

The code preserves both the continuous cutoff and the integer implementation so that rounding error remains visible.

## Statistical limitations

The cohort contains only eight deaths. Every estimate therefore has substantial uncertainty, and model selection can amplify that uncertainty. P values and confidence intervals obtained after AIC selection should not be interpreted as if the model had been fixed in advance. Repeated cross validation measures split sensitivity, but it does not replace external validation on an independent hospital cohort.

This repository is a research implementation. It is not a validated medical device and must not be used to make clinical decisions.
