# Linear-Mix-Model

## Overview
This script performs linear mixed model (LMM) analysis on a dataset using various statistical methods. It loads the required libraries, processes data, fits LMMs, and generates comprehensive reports in HTML format. The models investigate the interaction between **Plasma_AB42_40_V0_dic** and **Time** while adjusting for other covariates.

## Dependencies
The script requires the following R libraries:
```r
library(readxl)
library(nlme)
library(lme4)
library(lmerTest)
library(Matrix)
library(interactions)
library(plyr)
library(jtools)
library(ggplot2)
library(ggeffects)
library(dplyr)
library(emmeans)
library(performance)
library(sjPlot)
library(knitr)
library(rmarkdown)
library(kableExtra)
```
Ensure these packages are installed before running the script.

## Data
The script reads two tab-separated value files:
- `FACEHBI_BD_long_Araclon.ver2.0.txt`
- `FACEHBI_BD_long_Araclon.ver2.0.long.txt`

These files contain clinical and cognitive data required for model fitting.

## Variables Used
The script analyzes multiple dependent variables, including:
- **Neuroimaging metrics**: `fbb_pet_centilod`, `aHipocampal_vol`, `eTIV`, `aVentricular_vol`, `Cortical_vol`
- **Cognitive scores**: `mmse_fac`, `fname_total_fac`, `tmta_time_fac`, `tmtb_time_fac`, `composite_global_cognition`, etc.

## Function: `mlms()`
The core function, `mlms()`, performs the following:
1. Iterates through a list of variables.
2. Removes rows with missing values in the relevant columns.
3. Fits an LMM using the `lme()` function:
   ```r
   modelo <- lme(fixed = variable_actual ~ basal_age + apoe4_al + Plasma_AB42_40_V0_dic*Time,
                 random = ~ Time | facehbi_id_fac,
                 correlation = corCompSymm(form = ~ Time | facehbi_id_fac),
                 data = datos_no_na,
                 control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100))
   ```
4. Extracts results, including fixed effects, simple slopes, and model performance metrics.
5. Generates plots using `ggpredict()`.
6. Saves results as an RDS file for reporting.
7. Renders two RMarkdown reports:
   - `Reporte_completo_modelos_lineales_mixtos.html`: Full report with model summaries and plots.
   - `Resumen_reporte_modelos_lineales_mixtos.html`: Summary report with key statistics (AIC, interaction coefficients).

## Output
- **Full Report (`.html`)**: Detailed results for each model, including coefficients, statistical tests, and visualizations.
- **Summary Report (`.html`)**: Key model parameters such as AIC and interaction effects.
- **RDS File (`resultados_modelos.rds`)**: Saved model results for further analysis.

## Usage
Run the script in an R environment with the required datasets available:
```r
mlms(lista_variables, dfl, main_output_file = "Reporte_completo_modelos_lineales_mixtos.html",
     supplementary_output_file = "Resumen_reporte_modelos_lineales_mixtos.html")
```
This function will process the data, fit the models, and generate the reports automatically.


