# Extreme Limiting Dilution Analysis (ELDA) - R Implementation

An R script that **exactly replicates** the statistical analysis and output formatting from the official ELDA web tool ([bioinf.wehi.edu.au/software/elda](http://bioinf.wehi.edu.au/software/elda/)), enabling reproducible stem cell frequency estimation with professional PDF reports.

## 📋 Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Why Use This Script?](#why-use-this-script)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Output Format](#output-format)
- [Understanding the Results](#understanding-the-results)
- [Validation](#validation)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

## 🔬 Overview

This R script performs Extreme Limiting Dilution Analysis (ELDA) to estimate stem cell frequencies in cancer research. **All calculations and statistical tests produce identical results to the official ELDA web tool**, with the added benefit of:

- ✅ **Automated PDF report generation** with publication-ready tables
- ✅ **Batch processing** for multiple datasets
- ✅ **Reproducible workflows** via R scripts
- ✅ **Integration** with existing R analysis pipelines

### What is ELDA?

ELDA is a statistical method that estimates the frequency of stem cells (or cancer stem cells) in a population by testing multiple cell doses and counting how many wells form colonies. It assumes **single-hit kinetics**: one stem cell is sufficient to generate a colony.

## ✨ Key Features

### Perfect Replication of Web Tool
- ✅ Identical stem cell frequency estimates
- ✅ Same confidence intervals (95% CI)
- ✅ Matching chi-square test statistics
- ✅ Equivalent p-values for all comparisons
- ✅ Identical goodness-of-fit results

### Professional Output
- 📊 **Formatted tables** matching web tool layout
- 📈 **Log-fraction dose-response plot** 
- 📝 **Publication-ready PDF reports**
- 💾 **Easy export** for manuscripts

### Comprehensive Statistical Analysis
1. **Stem cell frequency estimation** with 95% confidence intervals
2. **Overall group comparison** (chi-square test)
3. **Pairwise comparisons** between all treatment groups
4. **Goodness-of-fit tests** to validate model assumptions:
   - Likelihood ratio test of single-hit model
   - Score test of heterogeneity

## 🎯 Why Use This Script?

### Advantages Over Web Tool

| Feature | Web Tool | This Script |
|---------|----------|-------------|
| Results accuracy | ✅ Reference standard | ✅ **Identical results** |
| Batch processing | ❌ One dataset at a time | ✅ Multiple datasets |
| Reproducibility | ⚠️ Manual entry required | ✅ Fully scripted |
| PDF reports | ✅ Download available | ✅ **Auto-generated** |
| Integration | ❌ Standalone | ✅ Part of R workflow |
| Offline use | ❌ Requires internet | ✅ Works offline |
| Version control | ❌ Not applicable | ✅ Git-trackable |

### Perfect for:
- 🔬 **Researchers** needing reproducible analysis
- 📊 **Labs** processing multiple ELDA experiments
- 📝 **Publications** requiring documented methods
- 🔄 **Longitudinal studies** with repeated measurements
- 🤝 **Collaborations** sharing standardized protocols

## 📦 Requirements

### R Version
- R ≥ 4.0.0 (recommended)

### Required Packages
```r
