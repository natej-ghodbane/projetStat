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

diet <- subset(diet, wave == 2009)
data <- merge(diet, bio, by = "IDind")

str(data)


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

anova_tc <- aov(TC ~ Pattern1_Q, data = data)
summary(anova_tc)

# Homogeneity of variances
bartlett.test(TC ~ Pattern1_Q, data = data)
leveneTest(TC ~ Pattern1_Q, data = data)

# Visualization
boxplot(TC ~ Pattern1_Q, data = data,
        col = "lightgreen",
        main = "TC across Pattern1 Quartiles",
        ylab = "Total Cholesterol")


############################################################
# 10. NON-PARAMETRIC ALTERNATIVE (Kruskal-Wallis)
############################################################

kruskal.test(TC ~ Pattern1_Q, data = data)
kruskal.test(LDL_C ~ Pattern1_Q, data = data)
kruskal.test(HDL_C ~ Pattern1_Q, data = data)
kruskal.test(TG ~ Pattern1_Q,  data = data)


############################################################
# 11. CORRELATION: PEARSON vs SPEARMAN
############################################################

cor.test(data$TG, data$d3fat, method = "pearson")
cor.test(data$TG, data$d3fat, method = "spearman")


############################################################
# END OF FILE
############################################################
