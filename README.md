# Youth Employment Gaps in Austria

This project examines differences in early-career labour market outcomes between young migrants and native-born individuals in Austria.

The analysis uses Austrian Microcensus Labour Force Survey data for 2022–2024 and focuses on individuals aged 15–34.

## Research Questions

- Are young migrants less likely to be employed than native-born youth?
- How do employment outcomes differ across migrant groups?
- Do the results differ in Vienna?
- How much of the employment gap can be explained by observable characteristics?

## Methods

- Data cleaning and variable construction in R
- Descriptive analysis
- Binary logistic regression
- Multinomial logistic regression
- Interaction effects
- Train-test model evaluation
- ROC curves and AUC
- Oaxaca–Blinder decomposition

## Main Findings

The results indicate that young people with a migration background have a lower probability of employment than native-born youth, after controlling for demographic, educational, household, regional, and time characteristics.

The analysis also finds substantial differences across origin groups and evidence that migrants are more likely to report wanting to work or wanting additional working hours. Oaxaca–Blinder decomposition suggests that only part of the employment gap can be attributed to observable characteristics.

## Data

The project uses Austrian Microcensus Labour Force Survey Scientific Use Files for 2022–2024. The original data are not included in this repository due to redistribution restrictions.

The analytical sample consists of individuals aged 15–34.

## Repository Structure

- `youth_employment.R` — data preparation, regression models, robustness checks, and decomposition analysis
- `tables/` — model evaluation and descriptive figures
- `data/README.md` — information about the data source

## Tools

R, tidyverse, dplyr, VGAM, nnet, pROC, margins, oaxaca, stargazer
