# Extreme Limiting Dilution Analysis (ELDA) - Complete Statistical Framework

A comprehensive R script for performing Extreme Limiting Dilution Analysis (ELDA) to estimate stem cell frequencies in cancer research, with complete model validation and statistical testing.

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Understanding the Analysis](#understanding-the-analysis)
- [Output Interpretation](#output-interpretation)
- [Example Data](#example-data)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

## 🔬 Overview

Extreme Limiting Dilution Analysis (ELDA) is a statistical method used to estimate the frequency of stem cells or cancer stem cells in a population. This script provides:

- **Stem cell frequency estimation** with 95% confidence intervals
- **Statistical comparison** between treatment groups
- **Model validation tests** to ensure assumptions are met
- **Pairwise comparisons** for multiple treatment groups
- **Visualization** of dose-response curves

## ✨ Features

### Core Functionality
- ✅ Bias-reduced maximum likelihood estimation
- ✅ Complementary log-log (cloglog) transformation for single-hit kinetics
- ✅ Multiple group comparison with automatic pairwise testing
- ✅ Comprehensive goodness-of-fit validation

### Statistical Tests Included
1. **Overall Group Comparison**: Tests if any groups differ significantly
2. **Likelihood Ratio Test**: Validates single-hit model assumption (parallel slopes)
3. **Goodness of Fit Test**: Ensures model adequately fits observed data
4. **Overdispersion Test**: Detects technical variability beyond expected

### Output Features
- 📊 Formatted tables for all statistical tests
- 📈 Log-fraction dose-response plot
- 📝 Automatic interpretation of results
- 💾 Easy-to-export results for publication

## 📦 Requirements

### R Version
- R ≥ 4.0.0 (recommended)

### Required Packages
```r
install.packages("statmod")
```

### Optional Packages (for enhanced functionality)
```r
install.packages(c("ggplot2", "knitr", "kableExtra"))
```

## 🚀 Installation

### Option 1: Clone the Repository
```bash
git clone https://github.com/yourusername/elda-analysis.git
cd elda-analysis
```

### Option 2: Download the Script
Download `elda_analysis.R` directly and save it to your working directory.

### Install Dependencies
```r
# In R console
install.packages("statmod")
library(statmod)
```

## 📖 Usage

### Basic Usage

1. **Prepare your data** in the following format:

```r
elda_data <- data.frame(
  cells = c(50, 40, 30, 20, 10, 50, 40, 30, 20, 10),  # Number of cells plated
  positive = c(10, 10, 10, 8, 6, 4, 3, 1, 0, 0),      # Wells with colonies
  tested = c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10), # Total wells tested
  group = c(rep("Control", 5), rep("Treatment", 5))    # Treatment groups
)
```

2. **Run the analysis**:

```r
source("elda_analysis.R")
```

### Data Format Requirements

Your dataframe must contain these columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `cells` | numeric | Number of cells plated per well | 10, 20, 30, 40, 50 |
| `positive` | numeric | Number of wells showing colony growth | 2, 4, 7, 9, 10 |
| `tested` | numeric | Total number of wells tested | 10, 10, 10, 10, 10 |
| `group` | character/factor | Treatment or condition name | "Control", "Drug_A" |

### Example Workflow

```r
# Load library
library(statmod)

# Load your data
elda_data <- read.csv("my_elda_data.csv")

# Ensure group is a factor
elda_data$group <- factor(elda_data$group)

# Run the complete analysis
source("elda_analysis.R")

# Results will be printed to console
# Plot will be generated automatically
```

## 🧮 Understanding the Analysis

### The ELDA Method

ELDA assumes **single-hit kinetics**: one stem cell is sufficient to generate a colony. The probability of observing a positive well follows:

```
P(positive) = 1 - exp(-dose / frequency)
```

Where:
- **dose** = number of cells plated
- **frequency** = number of cells needed to contain one stem cell

### Statistical Models

The script fits three nested GLM models:

#### 1. **Null Model** (Simplest)
```r
y ~ log(cells)
```
- Assumes all groups are identical
- One slope, one intercept

#### 2. **Single-Hit Model** (Standard ELDA)
```r
y ~ log(cells) + group
```
- Different frequencies (intercepts) for each group
- **Same slope** (parallel lines) - assumes single-hit kinetics
- **This is the model ELDA uses**

#### 3. **Full Model** (Most Complex)
```r
y ~ log(cells) * group
```
- Different slopes AND intercepts for each group
- Allows groups to violate single-hit assumption

### Model Validation Tests

#### Test 1: Likelihood Ratio Test
**Question**: Do we need different slopes, or are parallel slopes sufficient?

```
H0: Parallel slopes are adequate (single-hit model is valid)
H1: Groups have different slopes (single-hit model violated)
```

- **p > 0.05**: ✅ Single-hit model is appropriate
- **p < 0.05**: ❌ Single-hit assumption violated

#### Test 2: Goodness of Fit Test
**Question**: Does the model fit the observed data well?

```
Uses deviance statistic to compare observed vs predicted values
```

- **p > 0.05**: ✅ Model fits well
- **p < 0.05**: ❌ Poor fit - check for outliers or data issues

#### Test 3: Overdispersion Test
**Question**: Is there more variability than expected?

```
Checks if observed variance exceeds binomial expectation
```

- **p > 0.05**: ✅ No overdispersion
- **p < 0.05**: ❌ Extra variability detected (possible technical issues)

**Common causes of overdispersion**:
- Well-to-well contamination
- Plate edge effects
- Batch effects
- Technical variability

## 📊 Output Interpretation

### 1. Stem Cell Frequency Estimates

```
                      Group Estimate Lower_CI Upper_CI
1      GSC8-11_mock_DMSOctrl     5.79     3.46     9.99
2 GSC8-11_2Gyx3_DMSOctrl    63.18    39.83   100.39
3  GSC8-11_2Gyx3_ICM2.5uM   165.13    83.32   327.76
```

**Interpretation**:
- **Estimate**: Number of cells needed to contain one stem cell
- Lower values = **Higher stem cell frequency**
- Example: "1 in 5.79" means ~17% of cells are stem cells

### 2. Overall Test Results

```
     Test Chisq DF   P_Value
1 Overall 92.53  2 8.09e-21
```

**Interpretation**:
- **Chisq**: Chi-square test statistic
- **p < 0.001**: Highly significant difference between groups

### 3. Goodness of Fit Results

```
              Test  Chisq DF   P_Value
1 Likelihood_Ratio   2.19  2 3.35e-01
2  Goodness_of_Fit  15.57 11 1.58e-01
3  Overdispersion   15.57 11 1.58e-01
```

**All tests PASSED** (p > 0.05):
- ✅ Single-hit model is valid
- ✅ Model fits data well
- ✅ No technical issues detected

### 4. Pairwise Comparisons

```
                    Group1                   Group2 Chisq DF   P_Value
1 GSC8-11_2Gyx3_DMSOctrl GSC8-11_2Gyx3_ICM2.5uM  5.56  1 1.84e-02
2 GSC8-11_2Gyx3_DMSOctrl   GSC8-11_mock_DMSOctrl 55.05  1 1.18e-13
3 GSC8-11_2Gyx3_ICM2.5uM   GSC8-11_mock_DMSOctrl 87.35  1 9.08e-21
```

**Interpretation**:
- Each row compares two groups
- **p < 0.05**: Significant difference between this pair

## 📝 Example Data

The script includes example data from a cancer stem cell study:

### Experimental Design
- **Groups**: 3 treatment conditions
  - Mock DMSO control
  - 2Gy×3 radiation + DMSO
  - 2Gy×3 radiation + ICM 2.5µM
- **Cell doses**: 10, 20, 30, 40, 50 cells per well
- **Replicates**: 10 wells per dose

### Expected Results
1. **Radiation reduces stem cell frequency** (vs mock)
2. **Drug + radiation reduces it further** (vs radiation alone)
3. **All model validation tests pass**

## 🔧 Troubleshooting

### Common Issues

#### Issue 1: "Error in elda(): object not found"
**Solution**: Make sure `statmod` package is loaded
```r
library(statmod)
```

#### Issue 2: "Groups not found" or "Factor level missing"
**Solution**: Convert group column to factor
```r
elda_data$group <- factor(elda_data$group)
```

#### Issue 3: Model validation tests fail (p < 0.05)

**Likelihood Ratio Test fails**:
- Different groups may have different dose-response patterns
- Consider analyzing groups separately
- Check for biological reasons (e.g., different mechanisms)

**Goodness of Fit fails**:
- Check for outliers in your data
- Verify data entry is correct
- Consider if single-hit model is appropriate

**Overdispersion detected**:
- Possible well contamination
- Check for plate effects (edge wells vs center)
- Consider technical replicates on different days
- May need to adjust standard errors

#### Issue 4: Very wide confidence intervals
**Possible causes**:
- Too few replicates at each dose
- Doses not spanning the optimal range
- Very low response rates

**Solutions**:
- Increase number of replicates (≥8 per dose recommended)
- Test more cell doses
- Ensure doses bracket the expected frequency

### Data Quality Checks

Before running analysis, verify:

```r
# Check for missing values
any(is.na(elda_data))

# Verify positive ≤ tested for all rows
all(elda_data$positive <= elda_data$tested)

# Check if you have multiple doses per group
table(elda_data$group, elda_data$cells)

# Verify adequate replicates
table(elda_data$group)
```

## 📚 Citation

### For the ELDA Method
```
Hu, Y. and Smyth, G.K. (2009). 
ELDA: Extreme limiting dilution analysis for comparing depleted and enriched populations in stem cell and other assays. 
Journal of Immunological Methods, 347, 70-78.
```

### For the statmod Package
```
Giner, G. and Smyth, G.K. (2016). 
statmod: probability calculations for the inverse Gaussian distribution.
R Journal, 8, 339-351.
```

### BibTeX
```bibtex
@article{hu2009elda,
  title={ELDA: extreme limiting dilution analysis for comparing depleted and enriched populations in stem cell and other assays},
  author={Hu, Yifang and Smyth, Gordon K},
  journal={Journal of immunological methods},
  volume={347},
  number={1-2},
  pages={70--78},
  year={2009},
  publisher={Elsevier}
}
```

## 🔗 Additional Resources

### Online ELDA Calculator
- [http://bioinf.wehi.edu.au/software/elda/](http://bioinf.wehi.edu.au/software/elda/)
- Web-based interface for simple analyses

### Documentation
- [statmod Package Documentation](https://cran.r-project.org/package=statmod)
- [ELDA User Guide](http://bioinf.wehi.edu.au/software/elda/index.html)

### Related Methods
- **Flow cytometry** for prospective isolation
- **Single-cell RNA sequencing** for stem cell identification
- **In vivo limiting dilution** for functional validation

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### How to Contribute
1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Reporting Issues
Please use the GitHub issue tracker to report bugs or suggest enhancements.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Authors

- **Your Name** - *Initial work* - [YourGitHub](https://github.com/yourusername)

## 🙏 Acknowledgments

- Gordon Smyth and Yifang Hu for developing the ELDA method
- The R Core Team and statmod package maintainers
- Cancer research community for validation and feedback

## 📞 Contact

For questions or support:
- **Email**: khoih2023@gmail.com

---

**Note**: This script is for research purposes only. Always validate results with appropriate biological experiments and consult with a statistician for complex experimental designs.

## 🔖 Version History

### v1.0.0 (2025-01-13)
- Initial release
- Complete ELDA workflow
- Model validation tests
- Pairwise comparisons
- Automated plotting

---

*Last updated: January 2025*
