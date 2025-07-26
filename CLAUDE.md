# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains a comprehensive R framework for analyzing thermal tolerance in marine organisms through LD50/LT50 experiments. The framework processes survival data from organisms exposed to different temperatures for various durations, analyzing survival rates to determine lethal temperatures (LD50) and lethal times (LT50).

## Key Commands

### Running Analysis Scripts
```r
# Source and run full analysis (use Ctrl+Shift+S in RStudio)
source("Unified_Data_Analysis_2h_comparison.R")

# For different analysis types:
source("Unified_Data_Analysis_FR.R")              # Full-range with 3D visualization
source("Unified_Data_Analysis_FR_anova.R")        # ANOVA comparisons
source("Unified_Data_Analysis_FR_single-csv.R")   # Single file analysis
```

### Installing Required Packages
The scripts auto-install missing packages, but manually install with:
```r
install.packages(c("ggplot2", "Cairo", "lmtest"))
```

## Architecture and Code Organization

### Core Module Structure
The framework follows a modular design with three main analysis modules:

1. **plot_utilities.R** - Cross-platform plotting utilities
   - `check_os_func()`: OS detection and graphics driver setup
   - `chart_subtitle_func()`: Standardized chart labeling
   - `save_plot_func()`: Consistent plot saving (300 DPI PNG)

2. **non_linear_regression.R** - Sigmoid survival curve modeling
   - `non_linear_regression_file_func()`: Single CSV file processing
   - `non_linear_regression_dir_func()`: Directory-level analysis with visualization
   - Uses 3-parameter sigmoid model: `survival = p1 / (1 + (temperature / p2)^p3)`

3. **linear_regression.R** - Thermal Death Time (TDT) analysis
   - `linear_regression_func()`: TDT curve generation and statistics
   - `anova_analysis_func()`: Multi-condition comparison
   - `t_test_func()`: Two-condition statistical testing
   - `verify_normality()` and `verify_homogeneity()`: Statistical assumption validation

### Analysis Workflow
The typical analysis follows this sequence:
1. **Data Loading**: CSV files with Temperature, Alive, Dead columns
2. **Non-linear Fitting**: Sigmoid curves to extract LD50 values per time point
3. **Linear Regression**: TDT analysis on LD50 vs log(time) relationship
4. **Statistical Testing**: ANOVA or t-tests for condition comparisons
5. **Visualization**: LT50 curves, TDT plots, and comparative charts

### Data Structure Requirements
- CSV files must follow naming convention: `Dataxx.csv` (where xx = 01, 02, etc.)
- Directory structure represents experimental conditions
- Required CSV columns: Temperature, Alive, Dead
- Files organized by exposure duration (extracted from filename numbers)

### Statistical Analysis Features
- **Normality Testing**: Shapiro-Wilk test with Q-Q plots and histograms
- **Variance Homogeneity**: Breusch-Pagan test for heteroscedasticity
- **Robust Comparisons**: Welch's t-test for unequal variances
- **Non-parametric Options**: Mann-Whitney U test for non-normal data
- **Model Validation**: R-squared calculations for goodness-of-fit

### Output Generation
All visualizations are automatically saved to `charts/` directory with:
- High resolution (300 DPI) PNG format
- Standardized naming based on condition and chart type
- Consistent styling and color schemes
- Three visualization types:
  - LT50 Survival Curves (`_lt50.png`)
  - Thermal Death Time curves (`_tdt.png`) 
  - Combined comparative plots

### Mathematical Models
1. **Survival Model**: `Survival = p1 / (1 + ((Temperature - (z * log10(Time / t0))) / p2)^p3)`
2. **TDT Model**: `LD50 = intercept + slope * log10(Time)`

Where p1 = max survival, p2 = LD50 temperature, p3 = curve steepness, z = temperature coefficient.

## Configuration and Customization

### Modifying Analysis Parameters
Edit the main analysis script parameters:
```r
main_title <- "Your experiment description"
dir_paths <- list("path/to/your/data")
nls_param_list <- list(list(100, 30, 4))  # [max_survival, LD50_guess, slope]
```

### Graphics Configuration
The framework automatically handles Windows-specific Cairo graphics requirements. For other platforms, the `check_os_func()` ensures appropriate graphics drivers are loaded.

### Statistical Thresholds
Default significance level is α = 0.05. Modify in `linear_regression.R`:
```r
t_test_threshold <- 0.05
```