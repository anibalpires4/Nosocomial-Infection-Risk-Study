# Nosocomial Infection Risk Study

This repository analyses hospital-acquired (nosocomial) infection risk using data from 113 U.S. hospitals. The study employs Bayesian statistical methods and Markov Chain Monte Carlo (MCMC) techniques to evaluate infection control programs and identify significant factors influencing infection risk.

---

## Overview

Nosocomial infections are a critical healthcare concern. This project aims to:
1. Evaluate the effectiveness of infection surveillance and control programs in reducing infection rates.
2. Identify factors contributing to increased infection risk.
3. Develop a generalized linear model to accurately predict infection risk.

### Key Features
- **Data Source**: Collected during 1975–1976 from a random sample of U.S. hospitals.
- **Target Variable**: Infection Risk (average probability of acquiring an infection in the hospital).
- **Methodology**:
  - Initial exploratory analysis to identify data issues (e.g., outliers, correlations).
  - Bayesian modelling using MCMC to infer posterior distributions of model parameters.

---

## Data and Preprocessing

### Dataset Description
The dataset consists of 113 observations and 11 variables:
- **Identification Number**: Unique ID for each hospital.
- **Length of Stay (L.Stay)**: Average stay duration (days).
- **Age**: Average patient age (years).
- **Infection Risk (Inf.Risk)**: Probability of acquiring infection (%).
- **Routine Culturing Ratio (R.C.Ratio)**: Cultures performed per non-infected patient.
- **Routine Chest X-ray Ratio (R.C.Xray)**: X-rays performed per non-infected pneumonia patient.
- **Number of Beds**: Average number of beds.
- **Medical School Affiliation**: 1 = Yes, 2 = No.
- **Region**: Geographic region (1 = NE, 2 = NC, 3 = S, 4 = W).
- **Average Daily Census (ADC)**: Average number of daily patients.
- **Number of Nurses**: Full-time equivalent nurses.

### Preprocessing Steps
1. **Data Normalization**:
   - Ensures uniform scaling for Bayesian analysis and reduces MCMC autocorrelation.
2. **Variable Transformations**:
   - Categorical variables (Region, Medical School Affiliation) converted to factors.
3. **Outlier Detection**:
   - Identified using histograms, boxplots, and Q-Q plots.
4. **Correlation Analysis**:
   - High correlations (>80%) between **Beds**, **ADC**, and **Nurses** led to removal of redundant variables in model testing.

---

## Bayesian Modeling

### Generalized Linear Models (GLM)
Six models were tested by varying the inclusion of highly correlated variables:
1. **Model 1**: Full model.
2. **Model 2**: Removed Average Daily Census.
3. **Model 3**: Removed Number of Beds.
4. **Model 4**: Removed both Beds and ADC (Best Model).
5. **Model 5**: Removed Beds and Nurses.
6. **Model 6**: Removed ADC and Nurses.

### Model Selection
- **MCMC Configuration**:
  - Iterations: 10,000 (Burn-in: 1,000).
  - Chains: 2 parallel chains.
- **Prior Distributions**:
  - Coefficients (\( \beta \)): Normal(0, 0.0001).
  - Precision (\( \tau \)): Gamma(0.001, 0.001).
- **Selection Metric**: Deviance Information Criterion (DIC).

### Results
Model 4 had the lowest DIC (247.08), indicating the best fit.

---

## Results and Conclusions

### Key Findings
- Hospitals with medical school affiliations had lower infection risks.
- Infection risks varied significantly by region (highest in Region 2, lowest in Region 4).
- Variables **Beds**, **ADC**, and **Nurses** were highly correlated, reducing model performance.

### Residual Analysis
- The chosen model's residual showed normal distribution and low autocorrelation, confirming its suitability.

### Prediction
- Predictions from Model 4 closely matched observed infection risks (Mean: 4.355%, Variance: 1.037%).

---

## Technologies
- **Programming Language**: R
- **Libraries**:
  - `rjags` (Bayesian modelling with MCMC)
  - `coda` (MCMC diagnostics)
  - `ggplot2` (Visualization)

---

## Usage

1. **Clone the repository**:
   ```bash
   git clone https://github.com/anibalpires4/Nosocomial-Infection-Risk-Study.git

