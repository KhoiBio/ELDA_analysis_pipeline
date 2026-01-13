library(statmod)

# =============================================================================
# EXTREME LIMITING DILUTION ANALYSIS (ELDA) - COMPLETE WORKFLOW
# =============================================================================
# This script performs a comprehensive ELDA analysis including:
# - Stem cell frequency estimation with confidence intervals
# - Statistical comparison between treatment groups
# - Model validation tests (goodness of fit)
# - Pairwise group comparisons
# - Visualization of results
# =============================================================================

cat("\n")
cat("===============================================================================\n")
cat("           EXTREME LIMITING DILUTION ANALYSIS (ELDA) - START\n")
cat("===============================================================================\n\n")

# -----------------------------------------------------------------------------
# STEP 1: DATA ENTRY
# -----------------------------------------------------------------------------
# Purpose: Create a dataframe containing your experimental data
# Required columns:
#   - cells: number of cells plated per well
#   - positive: number of wells that showed colony growth
#   - tested: total number of wells tested at each dose
#   - group: treatment/condition name for each observation
# -----------------------------------------------------------------------------

cat("LOADING DATA\n")
cat("-----------------------------------------------------------------------------\n")

elda_data <- data.frame(
  cells = c(50, 40, 30, 20, 10, 50, 40, 30, 20, 10, 50, 40, 30, 20, 10),
  positive = c(10, 10, 10, 10, 8, 4, 3, 1, 0, 0, 6, 7, 0, 4, 1),
  tested = rep(10, 15),
  group = c(
    rep("GSC8-11_mock_DMSOctrl", 5),
    rep("GSC8-11_2Gyx3_ICM2.5uM", 5),
    rep("GSC8-11_2Gyx3_DMSOctrl", 5)
  )
)

# Convert group to factor for proper statistical analysis
elda_data$group <- factor(elda_data$group)

cat("Data loaded successfully!\n")
cat("  Number of observations:", nrow(elda_data), "\n")
cat("  Treatment groups:", nlevels(elda_data$group), "\n\n")

# -----------------------------------------------------------------------------
# RUN ELDA ANALYSIS
# -----------------------------------------------------------------------------
# Purpose: Calculate stem cell frequencies using maximum likelihood estimation
# Parameters:
#   - response: number of positive wells
#   - dose: number of cells plated
#   - tested: total wells tested
#   - group: treatment groups
#   - observed: FALSE uses bias-reduced estimates (more accurate)
# -----------------------------------------------------------------------------

cat("RUNNING ELDA ANALYSIS\n")
cat("-----------------------------------------------------------------------------\n")

elda_result <- elda(
  response = elda_data$positive,
  dose = elda_data$cells,
  tested = elda_data$tested,
  group = elda_data$group,
  observed = FALSE  # Use bias-reduced estimates
)

cat("ELDA analysis completed!\n\n")

# -----------------------------------------------------------------------------
# STEM CELL FREQUENCY ESTIMATES
# -----------------------------------------------------------------------------
# Purpose: Display the estimated stem cell frequencies with confidence intervals
# Interpretation: The "Estimate" tells you how many cells are needed to 
#                 contain one stem cell (lower = higher stem cell frequency)
# -----------------------------------------------------------------------------

cat("STEM CELL FREQUENCY ESTIMATES\n")
cat("-----------------------------------------------------------------------------\n")
cat("Interpretation: 'Estimate' = number of cells needed per stem cell\n")
cat("                Lower values indicate HIGHER stem cell content\n\n")

# Create formatted table
freq_table <- data.frame(
  Group = rownames(elda_result$CI),
  Estimate = round(elda_result$CI[, 2], 2),
  Lower_CI = round(elda_result$CI[, 3], 2),
  Upper_CI = round(elda_result$CI[, 1], 2)
)

print(freq_table, row.names = FALSE)
cat("\n\n")

# -----------------------------------------------------------------------------
# OVERALL GROUP COMPARISON
# -----------------------------------------------------------------------------
# Purpose: Test if there are statistically significant differences in stem 
#          cell frequency between ANY of the treatment groups
# Null Hypothesis: All groups have the same stem cell frequency
# -----------------------------------------------------------------------------

cat("OVERALL TEST FOR GROUP DIFFERENCES\n")
cat("-----------------------------------------------------------------------------\n")
cat("Tests: Are there significant differences between ANY groups?\n\n")

# Create formatted table for overall test
overall_table <- data.frame(
  Test = "Overall",
  Chisq = round(as.numeric(elda_result$test.difference[1]), 2),
  DF = as.numeric(elda_result$test.difference[3]),
  P_Value = formatC(as.numeric(elda_result$test.difference[2]), format = "e", digits = 2)
)

print(overall_table, row.names = FALSE)

# Interpretation
overall_p <- as.numeric(elda_result$test.difference[2])
cat("\nInterpretation: ")
if (overall_p < 0.001) {
  cat("*** HIGHLY SIGNIFICANT (p < 0.001) - Groups differ substantially\n\n\n")
} else if (overall_p < 0.05) {
  cat("** SIGNIFICANT (p < 0.05) - Groups show differences\n\n\n")
} else {
  cat("NOT SIGNIFICANT (p >= 0.05) - No evidence of group differences\n\n\n")
}

# -----------------------------------------------------------------------------
# GOODNESS OF FIT TESTS
# -----------------------------------------------------------------------------
# Purpose: Validate that the ELDA model assumptions are met and the fit is good
# Three tests:
#   1. Likelihood Ratio: Tests if groups have parallel slopes (single-hit model)
#   2. Goodness of Fit: Tests if the model fits the observed data adequately
#   3. Overdispersion: Tests for extra variability beyond binomial expectation
# -----------------------------------------------------------------------------

cat("MODEL VALIDATION - GOODNESS OF FIT TESTS\n")
cat("-----------------------------------------------------------------------------\n")
cat("These tests check if the ELDA model assumptions are valid\n\n")

# Prepare data for GLM analysis
# negative = number of wells WITHOUT colonies
negative <- elda_data$tested - elda_data$positive
y <- cbind(elda_data$positive, negative)

# Fit three nested models with increasing complexity:

# Model 1 (Full): Each group gets its own slope AND intercept
#         Formula: y ~ log(cells) * group  (expands to: log(cells) + group + log(cells):group)
#         Interpretation: Groups can have completely different dose-response curves
model_full <- glm(y ~ log(elda_data$cells) * elda_data$group, 
                  family = binomial(link = "cloglog"))

# Model 2 (Single-hit): All groups have SAME slope but different intercepts
#         Formula: y ~ log(cells) + group  (no interaction term)
#         Interpretation: Groups differ in frequency but follow same single-hit kinetics
#         This is the STANDARD ELDA assumption
model_single_hit <- glm(y ~ log(elda_data$cells) + elda_data$group, 
                        family = binomial(link = "cloglog"))

# Model 3 (Null): All groups are identical
#         Formula: y ~ log(cells)  (no group term at all)
#         Interpretation: No differences between groups
model_null <- glm(y ~ log(elda_data$cells), 
                  family = binomial(link = "cloglog"))

# Test 1: Likelihood Ratio Test
# Compares: Single-hit model vs Full model
# Question: Do we need different slopes, or are parallel slopes sufficient?
lr_test <- anova(model_single_hit, model_full, test = "Chisq")

# Test 2: Goodness of Fit (Deviance Test)
# Compares: Model predictions vs observed data
# Question: Does the single-hit model fit the data well?
deviance_test <- model_single_hit$deviance
df_resid <- model_single_hit$df.residual
p_gof <- pchisq(deviance_test, df_resid, lower.tail = FALSE)

# Test 3: Overdispersion Test (Pearson Chi-square)
# Compares: Observed variance vs expected binomial variance
# Question: Is there extra variability beyond what's expected?
fitted_probs <- fitted(model_single_hit)
pearson_resid <- (elda_data$positive - elda_data$tested * fitted_probs) / 
  sqrt(elda_data$tested * fitted_probs * (1 - fitted_probs))
pearson_chisq <- sum(pearson_resid^2)
p_pearson <- pchisq(pearson_chisq, df_resid, lower.tail = FALSE)

# Create formatted results table
gof_table <- data.frame(
  Test = c("Likelihood_Ratio", "Goodness_of_Fit", "Overdispersion"),
  Chisq = round(c(lr_test$Deviance[2], deviance_test, pearson_chisq), 2),
  DF = c(lr_test$Df[2], df_resid, df_resid),
  P_Value = formatC(c(lr_test$`Pr(>Chi)`[2], p_gof, p_pearson), 
                    format = "e", digits = 2)
)

print(gof_table, row.names = FALSE)

# Interpretation of each test
cat("\nInterpretations:\n")
cat("  Likelihood_Ratio: ", 
    ifelse(lr_test$`Pr(>Chi)`[2] > 0.05, 
           "PASS - Parallel slopes OK (single-hit model valid)", 
           "FAIL - Groups have different slopes"), "\n")
cat("  Goodness_of_Fit:  ", 
    ifelse(p_gof > 0.05, 
           "PASS - Model fits data well", 
           "FAIL - Poor model fit"), "\n")
cat("  Overdispersion:   ", 
    ifelse(p_pearson > 0.05, 
           "PASS - No extra variability detected", 
           "FAIL - Extra variability present"), "\n\n\n")

# -----------------------------------------------------------------------------
# PAIRWISE COMPARISONS
# -----------------------------------------------------------------------------
# Purpose: Compare each pair of groups directly to identify which specific
#          groups differ from each other
# Useful when overall test is significant to pinpoint differences
# -----------------------------------------------------------------------------

cat("PAIRWISE COMPARISONS BETWEEN GROUPS\n")
cat("-----------------------------------------------------------------------------\n")
cat("Direct comparison of each group pair\n\n")

# Get list of all groups
groups <- levels(elda_data$group)

# Initialize empty dataframe to store results
pairwise_results <- data.frame()

# Loop through all pairs of groups
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    
    # Extract data for just these two groups
    sub_data <- subset(elda_data, group %in% c(groups[i], groups[j]))
    
    # Run ELDA on this pair
    pw_result <- limdil(sub_data$positive, sub_data$cells, 
                        sub_data$tested, sub_data$group)
    
    # Store results in a row
    row <- data.frame(
      Group1 = groups[i],
      Group2 = groups[j],
      Chisq = round(as.numeric(pw_result$test.difference[1]), 2),
      DF = as.numeric(pw_result$test.difference[3]),
      P_Value = formatC(as.numeric(pw_result$test.difference[2]), 
                        format = "e", digits = 2)
    )
    
    # Add row to results table
    pairwise_results <- rbind(pairwise_results, row)
  }
}

# Display results
print(pairwise_results, row.names = FALSE)
cat("\n\n")

# -----------------------------------------------------------------------------
# VISUALIZATION
# -----------------------------------------------------------------------------
# Purpose: Generate a log-fraction plot showing dose-response curves
# The plot shows:
#   - Each group as a separate line
#   - Dots are observed data points
#   - Lines are model fits
#   - Parallel lines indicate valid single-hit model
# -----------------------------------------------------------------------------

cat("GENERATING VISUALIZATION\n")
cat("-----------------------------------------------------------------------------\n")
cat("Creating log-fraction plot...\n\n")

plot(elda_result, main = "Limiting Dilution Analysis (ELDA)\nLog-Fraction Plot")

cat("Plot generated!\n")
cat("  - Each line represents one treatment group\n")
cat("  - Dots show observed data points\n")
cat("  - Lines should be parallel if single-hit model is valid\n")
cat("  - Higher lines = higher stem cell frequency\n\n")

# -----------------------------------------------------------------------------
# ANALYSIS COMPLETE
# -----------------------------------------------------------------------------

cat("===============================================================================\n")
cat("                    ANALYSIS COMPLETE - SUMMARY\n")
cat("===============================================================================\n\n")

cat("Results Summary:\n")
cat("  1. Stem cell frequencies calculated for", nlevels(elda_data$group), "groups\n")
cat("  2. Overall comparison:", 
    ifelse(overall_p < 0.05, "SIGNIFICANT differences found", "No significant differences"), "\n")
cat("  3. Model validation:", 
    ifelse(all(c(lr_test$`Pr(>Chi)`[2], p_gof, p_pearson) > 0.05), 
           "ALL TESTS PASSED", "Some warnings present"), "\n")
cat("  4. Pairwise comparisons: Completed for all group pairs\n")
cat("  5. Visualization: Plot generated\n\n")

cat("Analysis finished successfully!\n")
cat("===============================================================================\n\n")