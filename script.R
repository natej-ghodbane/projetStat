############################################################
# PROJECT: Dietary Patterns & Blood Lipids
############################################################


############################################################
# 1. LOAD PACKAGES
############################################################

library(psych)
library(car)
library(lmtest)
library(ggplot2)


############################################################
# 2. LOAD & MERGE DATA
############################################################

bio  <- read.csv("biomarker.csv")
diet <- read.csv("c12diet.csv")

dim(bio)
dim(diet)

diet <- subset(diet, wave == 2009)
data <- merge(diet, bio, by = "IDind")

dim(data)
str(data)

############################################################
# 2bis. DISTRIBUTIONS INITIALES DES VARIABLES BRUTES
############################################################

cat("\n\n")
cat("############################################################\n")
cat("INITIAL DISTRIBUTIONS OF RAW VARIABLES (PRE-PCA)\n")
cat("############################################################\n\n")

# Afficher les statistiques descriptives
cat("Descriptive Statistics of Dietary Variables:\n")
cat("-------------------------------------------\n")
diet_vars_desc <- data[, c("d3kcal", "d3carbo", "d3fat", "d3protn")]
print(summary(diet_vars_desc))

cat("\nDescriptive Statistics of Biomarkers:\n")
cat("-------------------------------------\n")
bio_vars_desc <- data[, c("HDL_C", "LDL_C", "TG", "TC")]
print(summary(bio_vars_desc))

## HISTOGRAMMES DES VARIABLES BRUTES
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))

# Variables diététiques
hist(data$d3kcal, main = "Calories (d3kcal)", col = "lightblue", 
     xlab = "Calories", breaks = 30)
abline(v = mean(data$d3kcal, na.rm = TRUE), col = "red", lwd = 2)
abline(v = median(data$d3kcal, na.rm = TRUE), col = "blue", lwd = 2, lty = 2)

hist(data$d3carbo, main = "Carbohydrates (d3carbo)", col = "lightgreen",
     xlab = "Carbohydrates", breaks = 30)
abline(v = mean(data$d3carbo, na.rm = TRUE), col = "red", lwd = 2)
abline(v = median(data$d3carbo, na.rm = TRUE), col = "blue", lwd = 2, lty = 2)

hist(data$d3fat, main = "Fat (d3fat)", col = "lightcoral",
     xlab = "Fat", breaks = 30)
abline(v = mean(data$d3fat, na.rm = TRUE), col = "red", lwd = 2)
abline(v = median(data$d3fat, na.rm = TRUE), col = "blue", lwd = 2, lty = 2)

hist(data$d3protn, main = "Protein (d3protn)", col = "lightgoldenrod",
     xlab = "Protein", breaks = 30)
abline(v = mean(data$d3protn, na.rm = TRUE), col = "red", lwd = 2)
abline(v = median(data$d3protn, na.rm = TRUE), col = "blue", lwd = 2, lty = 2)

# Biomarqueurs
hist(data$HDL_C, main = "HDL Cholesterol", col = "skyblue",
     xlab = "HDL-C (mg/dL)", breaks = 30)
abline(v = mean(data$HDL_C, na.rm = TRUE), col = "red", lwd = 2)
abline(v = median(data$HDL_C, na.rm = TRUE), col = "blue", lwd = 2, lty = 2)

hist(data$LDL_C, main = "LDL Cholesterol", col = "salmon",
     xlab = "LDL-C (mg/dL)", breaks = 30)
abline(v = mean(data$LDL_C, na.rm = TRUE), col = "red", lwd = 2)
abline(v = median(data$LDL_C, na.rm = TRUE), col = "blue", lwd = 2, lty = 2)

hist(data$TG, main = "Triglycerides", col = "lightgreen",
     xlab = "TG (mg/dL)", breaks = 30)
abline(v = mean(data$TG, na.rm = TRUE), col = "red", lwd = 2)
abline(v = median(data$TG, na.rm = TRUE), col = "blue", lwd = 2, lty = 2)

hist(data$TC, main = "Total Cholesterol", col = "plum",
     xlab = "TC (mg/dL)", breaks = 30)
abline(v = mean(data$TC, na.rm = TRUE), col = "red", lwd = 2)
abline(v = median(data$TC, na.rm = TRUE), col = "blue", lwd = 2, lty = 2)

## BOXPLOTS DES VARIABLES BRUTES
par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))

# Boxplots des variables diététiques
boxplot(data[, c("d3kcal", "d3carbo", "d3fat", "d3protn")],
        main = "Boxplots of Dietary Variables (Raw)",
        col = c("lightblue", "lightgreen", "lightcoral", "lightgoldenrod"),
        las = 2,  # Étiquettes verticales
        ylab = "Amount",
        cex.axis = 0.8)

# Boxplots des biomarqueurs
boxplot(data[, c("HDL_C", "LDL_C", "TG", "TC")],
        main = "Boxplots of Blood Lipids (Raw)",
        col = c("skyblue", "salmon", "lightgreen", "plum"),
        las = 2,  # Étiquettes verticales
        ylab = "Concentration (mg/dL)",
        cex.axis = 0.8)

## QQ-PLOTS DES VARIABLES BRUTES
par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))

# Variables diététiques
qqnorm(data$d3kcal, main = "QQ-plot d3kcal"); qqline(data$d3kcal)
qqnorm(data$d3carbo, main = "QQ-plot d3carbo"); qqline(data$d3carbo)
qqnorm(data$d3fat, main = "QQ-plot d3fat"); qqline(data$d3fat)
qqnorm(data$d3protn, main = "QQ-plot d3protn"); qqline(data$d3protn)

# Biomarqueurs
qqnorm(data$HDL_C, main = "QQ-plot HDL-C"); qqline(data$HDL_C)
qqnorm(data$LDL_C, main = "QQ-plot LDL-C"); qqline(data$LDL_C)
qqnorm(data$TG, main = "QQ-plot TG"); qqline(data$TG)
qqnorm(data$TC, main = "QQ-plot TC"); qqline(data$TC)

## MATRICE DE CORRÉLATION DES VARIABLES BRUTES
cat("\n\nCorrelation Matrix of Raw Dietary Variables:\n")
cat("--------------------------------------------\n")
cor_matrix_raw <- cor(diet_vars_desc, use = "complete.obs")
print(round(cor_matrix_raw, 3))

# Visualisation de la matrice de corrélation
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
corrplot::corrplot(cor_matrix_raw, 
                   method = "color",
                   type = "upper",
                   tl.col = "black",
                   tl.srt = 45,
                   addCoef.col = "black",
                   number.cex = 0.8,
                   title = "Correlation Matrix of Raw Dietary Variables",
                   mar = c(0, 0, 2, 0))

############################################################
# 3. PCA — DIETARY PATTERNS
############################################################

diet_vars <- data[, c("d3kcal", "d3carbo", "d3fat", "d3protn")]
diet_scaled <- scale(diet_vars)

pca <- principal(diet_scaled, nfactors = 3, rotate = "varimax", scores = TRUE)
print(pca, digits = 3)

scores <- as.data.frame(pca$scores)
colnames(scores) <- c("Pattern1", "Pattern2", "Pattern3")

data <- cbind(data, scores)


############################################################
# 4. CREATE QUARTILES (FOR ANOVA)
############################################################

data$Pattern1_Q <- cut(
  data$Pattern1,
  quantile(data$Pattern1, probs = seq(0,1,0.25), na.rm = TRUE),
  include.lowest = TRUE
)

data$Pattern2_Q <- cut(
  data$Pattern2,
  quantile(data$Pattern2, probs = seq(0,1,0.25), na.rm = TRUE),
  include.lowest = TRUE
)

data$Pattern3_Q <- cut(
  data$Pattern3,
  quantile(data$Pattern3, probs = seq(0,1,0.25), na.rm = TRUE),
  include.lowest = TRUE
)


############################################################
# 5. DISTRIBUTIONS & NORMALITY
############################################################

## Histograms
par(mfrow=c(2,3))

hist(data$HDL_C, main="HDL-C", col="lightblue")
hist(data$LDL_C, main="LDL-C", col="lightblue")
hist(data$TG,    main="TG",    col="lightblue")
hist(data$TC,    main="TC",    col="lightblue")

hist(data$Pattern1, main="Pattern1", col="pink")
hist(data$Pattern2, main="Pattern2", col="pink")
hist(data$Pattern3, main="Pattern3", col="pink")


############################################################
# 5bis. SHAPIRO–WILK ON RANDOM SUBSAMPLE (n ≤ 5000)
############################################################

set.seed(123)

hdl_sample <- sample(na.omit(data$HDL_C), 5000)
ldl_sample <- sample(na.omit(data$LDL_C), 5000)
tg_sample  <- sample(na.omit(data$TG),    5000)
tc_sample  <- sample(na.omit(data$TC),    5000)

shapiro.test(hdl_sample)
shapiro.test(ldl_sample)
shapiro.test(tg_sample)
shapiro.test(tc_sample)


############################################################
# 5ter. QQ-PLOTS FOR NORMALITY (FULL SAMPLE)
############################################################

par(mfrow = c(2,2))

qqnorm(data$HDL_C, main="QQ-plot HDL-C"); qqline(data$HDL_C)
qqnorm(data$LDL_C, main="QQ-plot LDL-C"); qqline(data$LDL_C)
qqnorm(data$TG,    main="QQ-plot TG");    qqline(data$TG)
qqnorm(data$TC,    main="QQ-plot TC");    qqline(data$TC)


############################################################
# 6. SIMPLE REGRESSION (PARAMETRIC)
############################################################

model_simple <- lm(TG ~ d3fat, data = data)
summary(model_simple)

par(mfrow=c(2,2))
plot(model_simple)

############################################################
# Residual normality (Shapiro on subsample)
############################################################

set.seed(123)
resid_sample <- sample(residuals(model_simple), 5000)
shapiro.test(resid_sample)

qqnorm(residuals(model_simple), main = "QQ-plot Residuals (TG ~ d3fat)")
qqline(residuals(model_simple))

############################################################
# 7. MULTIPLE LINEAR REGRESSIONS
############################################################

check_assumptions <- function(model, seed = 123, n_shapiro = 5000) {
  
  par(mfrow = c(2,2))
  plot(model)
  
  cat("\n--- Shapiro (Residuals, subsample) ---\n")
  set.seed(seed)
  res <- residuals(model)
  res <- res[!is.na(res)]
  
  if (length(res) > n_shapiro) {
    res <- sample(res, n_shapiro)
  }
  print(shapiro.test(res))
  
  cat("\n--- Breusch-Pagan (Homoscedasticity) ---\n")
  print(bptest(model))
  
  cat("\n--- Durbin-Watson (Independence) ---\n")
  print(dwtest(model))
  
  cat("\n--- Top Cook's Distance ---\n")
  cooks <- cooks.distance(model)
  print(sort(cooks, decreasing = TRUE)[1:5])
}

# HDL
model_hdl <- lm(HDL_C ~ Pattern1 + Pattern2 + Pattern3, data = data)
summary(model_hdl)
check_assumptions(model_hdl)

# LDL
model_ldl <- lm(LDL_C ~ Pattern1 + Pattern2 + Pattern3, data = data)
summary(model_ldl)
check_assumptions(model_ldl)

# TG
model_tg <- lm(TG ~ Pattern1 + Pattern2 + Pattern3, data = data)
summary(model_tg)
check_assumptions(model_tg)

# TC
model_tc <- lm(TC ~ Pattern1 + Pattern2 + Pattern3, data = data)
summary(model_tc)
check_assumptions(model_tc)


############################################################
# 8. LOG-TRANSFORMATION FOR SKEWED TG
############################################################

data$log_TG <- log(data$TG)

hist(data$log_TG, main="Log(TG)", col="orange")

model_log_tg <- lm(log_TG ~ Pattern1 + Pattern2 + Pattern3, data = data)
summary(model_log_tg)
plot(model_log_tg)


############################################################
# 9. ANOVA (PARAMETRIC)
############################################################
extract_anova <- function(aov_model, pattern_name) {
  s <- summary(aov_model)[[1]]
  data.frame(
    Pattern = pattern_name,
    Df = s[1, "Df"],
    F_value = s[1, "F value"],
    p_value = s[1, "Pr(>F)"]
  )
}
run_anova_patterns <- function(response, biomarker_name, data) {
  
  a1 <- aov(as.formula(paste(response, "~ Pattern1_Q")), data = data)
  a2 <- aov(as.formula(paste(response, "~ Pattern2_Q")), data = data)
  a3 <- aov(as.formula(paste(response, "~ Pattern3_Q")), data = data)
  
  r1 <- extract_anova(a1, "Pattern1_Q")
  r2 <- extract_anova(a2, "Pattern2_Q")
  r3 <- extract_anova(a3, "Pattern3_Q")
  
  out <- rbind(r1, r2, r3)
  out$Biomarker <- biomarker_name
  out <- out[, c("Biomarker", "Pattern", "Df", "F_value", "p_value")]
  out
}


anova_tc  <- run_anova_patterns("TC",     "TC",     data)
anova_tg  <- run_anova_patterns("TG",     "TG",     data)
anova_hdl <- run_anova_patterns("HDL_C",  "HDL-C",  data)
anova_ldl <- run_anova_patterns("LDL_C",  "LDL-C",  data)

anova_summary_all <- rbind(
  anova_tc,
  anova_tg,
  anova_hdl,
  anova_ldl
)

anova_summary_all


############################################################
# 10. NON-PARAMETRIC ALTERNATIVE (Kruskal-Wallis)
############################################################

kruskal.test(TC ~ Pattern1_Q, data = data)
kruskal.test(LDL_C ~ Pattern1_Q, data = data)
kruskal.test(HDL_C ~ Pattern1_Q, data = data)
kruskal.test(TG ~ Pattern1_Q,  data = data)

kruskal.test(TC ~ Pattern2_Q, data = data)
kruskal.test(LDL_C ~ Pattern2_Q, data = data)
kruskal.test(HDL_C ~ Pattern2_Q, data = data)
kruskal.test(TG ~ Pattern2_Q,  data = data)

kruskal.test(TC ~ Pattern3_Q, data = data)
kruskal.test(LDL_C ~ Pattern3_Q, data = data)
kruskal.test(HDL_C ~ Pattern3_Q, data = data)
kruskal.test(TG ~ Pattern3_Q,  data = data)

############################################################
# 11. CORRELATION: PEARSON vs SPEARMAN
############################################################

cor.test(data$TG, data$d3fat, method = "pearson")

extract_spearman <- function(x, y, var_x, var_y) {
  test <- suppressWarnings(cor.test(x, y, method = "spearman"))
  data.frame(
    Variable_X = var_x,
    Variable_Y = var_y,
    rho = as.numeric(test$estimate),
    p_value = test$p.value
  )
}
cor.test(data$TG, data$d3fat, method = "spearman")
cor.test(data$TC, data$d3fat, method = "spearman")
cor.test(data$LDL_C, data$d3fat, method = "spearman")
cor.test(data$HDL_C, data$d3fat, method = "spearman")

cor.test(data$TG, data$d3carbo, method = "spearman")
cor.test(data$TC, data$d3carbo, method = "spearman")
cor.test(data$LDL_C, data$d3carbo, method = "spearman")
cor.test(data$HDL_C, data$d3carbo, method = "spearman")

cor.test(data$TG, data$d3protn, method = "spearman")
cor.test(data$TC, data$d3protn, method = "spearman")
cor.test(data$LDL_C, data$d3protn, method = "spearman")
cor.test(data$HDL_C, data$d3protn, method = "spearman")


############################################################
# END OF FILE
############################################################
