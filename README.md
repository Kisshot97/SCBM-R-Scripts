# SCBM R Scripts

This repository contains a collection of R scripts for the SCBM (Structured Causal Boosting with Modularity) project. The scripts support data simulation, causal inference modeling, subgroup effect analysis, model evaluation, and visualization in both randomized and observational data settings.

## Project Structure

All scripts are modularized and annotated in English. File paths are dynamically constructed using `getwd()` to ensure compatibility across environments.

## Included Scripts

| Script Name                  | Description |
|-----------------------------|-------------|
| `fun_GL_PS.R`                | Group lasso + bootstrap bagging MARS model for HTE estimation, with IPTW and transformed outcome variants.  |
| `fun_GL_PS_stratify.R`         | Group lasso + bootstrap bagging MARS model for HTE estimation, with stratification variant. |
| `he.heatmapdata.R`         | Group lasso + bootstrap bagging MARS model for HTE estimation, with IPTW and transformed outcome variants. |

| `he.heatmapplot.R`       | for the heatmapplot |
| `he.simulaitonplot.R`           | Generates boxplots for MSE or bias comparison across multiple methods. (plot)|
| `he.simulaitonplot.R`        | Generates boxplots for MSE or bias comparison across multiple methods. (data)|
| `he.trendplot.R`         | Shows subgroup-wise trends of true vs. predicted treatment effects (HTEs). |
| `SCBM_summary_table.R`      | Computes average bias and MSE per method, exports LaTeX-ready tables. |
| `he.pd plot.R`            | Partial dependence plots (PDPs) for a real dataset using group lasso + MARS. |
| `SCBM_real_vimp.R`          | Variable importance analysis for real-world data. |
| `SCBM_realdata_greedy.R`    | SEARCH the best parameters in real dataset. |
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
