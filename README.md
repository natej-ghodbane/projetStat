# Association Between Dietary Patterns and Blood Lipid Profiles

## 📌 Project Overview

This project replicates and extends the methodology of the study:

**“Association between dietary patterns and blood lipid profiles among Chinese women”**

The objective is to:
- Identify dietary patterns using **Principal Component Analysis (PCA)**
- Examine associations between dietary patterns and blood lipid biomarkers:
  - HDL-C
  - LDL-C
  - Triglycerides (TG)
  - Total Cholesterol (TC)
- Demonstrate correct use of **parametric and non-parametric statistical tests** based on assumption checking

---

## 📊 Data Description

### Dietary Data (`c12diet.csv`)
- 3-day nutrient intake
- Variables used in PCA:
  - `d3kcal` (Energy)
  - `d3carbo` (Carbohydrates)
  - `d3fat` (Fat)
  - `d3protn` (Protein)

### Biomarker Data (`biomarker.csv`)
- Blood lipid and biochemical measurements:
  - HDL-C, LDL-C, TG, TC
  - Other clinical biomarkers

Only **wave 2009** observations were retained.

---

## 🧪 Methodology

### 1️⃣ Dietary Pattern Identification
- PCA with **varimax rotation**
- Selection based on eigenvalues > 1
- Three dietary patterns identified:
  - **Pattern 1**: Carbohydrate / energy-rich
  - **Pattern 2**: Fat-rich
  - **Pattern 3**: Protein-rich

### 2️⃣ Statistical Analyses
- Descriptive statistics and histograms
- Normality assessment:
  - QQ-plots
  - Shapiro–Wilk test on subsamples
- Simple linear regression
- Multiple linear regression with diagnostic checks:
  - Residual normality
  - Homoscedasticity (Breusch–Pagan)
  - Independence (Durdurbin–Watson)
- ANOVA across dietary pattern quartiles
- Non-parametric alternatives:
  - Kruskal–Wallis test
  - Spearman correlation

---

## 📈 Main Results (Summary)

- **Carbohydrate-rich pattern** → lower LDL-C and TC
- **Fat-rich pattern** → higher TG and TC
- **Protein-rich pattern** → higher TG and TC
- HDL-C showed no strong association with dietary patterns
- TG distributions were highly skewed → log-transformation improved model fit

---

## 📂 Repository Structure

```
.
├── biomarker.csv
├── c12diet.csv
├── test.Rmd            # Main analysis (R Markdown)
├── test.pdf
├── script.R            # Optional R script
├── README.md
├── .gitignore
└── projet2.Rproj
```

---

## 🛠️ Requirements

- R (≥ 4.2)
- Packages:
  - `psych`
  - `ggplot2`
  - `car`
  - `lmtest`
  - `reshape2`

Install packages if needed:

```r
install.packages(c("psych", "ggplot2", "car", "lmtest", "reshape2"))
```

---

## 📄 Report

The full statistical report is available as:
- `test.html`
- `test.pdf` 

---

## 👨‍🎓 Author

Natej Ghodbane
Academic project — Statistical Analysis / Nutritional Epidemiology
