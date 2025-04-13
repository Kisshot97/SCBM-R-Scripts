# SCBM R Scripts

This repository contains a collection of R scripts for the SCBM (Structured Causal Boosting with Modularity) project. The scripts support data simulation, causal inference modeling, subgroup effect analysis, model evaluation, and visualization in both randomized and observational data settings.

## Project Structure

All scripts are modularized and annotated in English. File paths are dynamically constructed using `getwd()` to ensure compatibility across environments.

## Included Scripts

| Script Name                  | Description |
|-----------------------------|-------------|
| `SCBM_datagen.R`            | Functions for simulating treatment assignment, covariates, and heterogeneous treatment effects (HTEs). |
| `SCBM_simulation.R`         | Main simulation pipeline for training/testing HTE estimators across multiple scenarios. |
| `SCBM_func_gl_ps.R`         | Group lasso + MARS model for HTE estimation, with IPTW and transformed outcome variants. |
| `SCBM_func_gl_greed.R`      | A variant of the above model with fixed resampling strategies. |
| `SCBM_make_plot_data.R`     | Aggregates MSE and bias values from model outputs across simulation runs. |
| `SCBM_plot_box.R`           | Generates boxplots for MSE or bias comparison across multiple methods. |
| `SCBM_plot_greedy.R`        | Visualizes CV error heatmaps for tuning parameters (e.g., sampling rate, boost rounds). |
| `SCBM_trend_plot.R`         | Shows subgroup-wise trends of true vs. predicted treatment effects (HTEs). |
| `SCBM_summary_table.R`      | Computes average bias and MSE per method, exports LaTeX-ready tables. |
| `SCBM_real_pd.R`            | Partial dependence plots (PDPs) for a real dataset using group lasso + MARS. |
| `SCBM_real_vimp.R`          | Variable importance analysis for real-world data. |

## Requirements

These scripts rely on the following R packages:

- `grf`
- `earth`
- `glmnet`
- `grpreg`
- `BART`
- `randomForest`
- `rpart`
- `ggplot2`
- `reshape2`
- `doParallel`
- `foreach`
- `cowplot`, `gridExtra`

Install all dependencies with:

```r
install.packages(c("grf", "earth", "glmnet", "grpreg", "BART", "randomForest",
                   "rpart", "ggplot2", "reshape2", "doParallel", "foreach",
                   "cowplot", "gridExtra"))
