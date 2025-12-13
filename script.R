############################################################
# 1. LOAD DATA
############################################################

bio  <- read.csv("biomarker.csv")
diet <- read.csv("c12diet.csv")

diet <- subset(diet, wave == 2009)
data <- merge(diet, bio, by = "IDind")
str(data)

############################################################
# 2. PCA — Dietary Patterns
############################################################

diet_vars <- data[, c("d3kcal", "d3carbo", "d3fat", "d3protn")]
diet_scaled <- scale(diet_vars)

library(psych)
pca <- principal(diet_scaled, nfactors = 3, rotate = "varimax")
print(pca, digits = 3)

scores <- as.data.frame(pca$scores)
colnames(scores) <- c("Pattern1", "Pattern2", "Pattern3")
data <- cbind(data, scores)

############################################################
# 3. QUARTILES NEST7A9HOM LEL ANOVA
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
# 4. HISTOGRAMS 
############################################################

par(mfrow=c(2,3))

hist(data$HDL_C, main="HDL Distribution", col="lightblue")
hist(data$LDL_C, main="LDL Distribution", col="lightblue")
hist(data$TG, main="TG Distribution", col="lightblue")
hist(data$TC, main="TC Distribution", col="lightblue")

hist(data$Pattern1, main="Pattern1 (Carbs)", col="pink")
hist(data$Pattern2, main="Pattern2 (Fat)", col="pink")
hist(data$Pattern3, main="Pattern3 (Protein)", col="pink")

############################################################
# 5. SIMPLE REGRESSION 
############################################################

model_simple <- lm(TG ~ d3fat, data = data)
summary(model_simple)

par(mfrow=c(2,2))
plot(model_simple)  # Residuals, QQ, Scale-Location, Leverage

############################################################
# 6. MULTIPLE REGRESSION 
############################################################

check_assumptions <- function(model) {
  
  par(mfrow = c(2, 2))
  plot(model)
  
  cat("\n--- Breusch-Pagan Test (Homoscedasticity) ---\n")
  library(lmtest)
  print(bptest(model))
  
  cat("\n--- Durbin-Watson Test (Autocorrelation) ---\n")
  print(dwtest(model))
  
  cat("\n--- Cook's Distance (Influential points) ---\n")
  cooks <- cooks.distance(model)
  print(sort(cooks, decreasing = TRUE)[1:10])
}

### HDL Model
model_hdl <- lm(HDL_C ~ Pattern1 + Pattern2 + Pattern3, data = data)
summary(model_hdl)
check_assumptions(model_hdl)

### LDL Model
model_ldl <- lm(LDL_C ~ Pattern1 + Pattern2 + Pattern3, data = data)
summary(model_ldl)
check_assumptions(model_ldl)

### TG Model
model_tg <- lm(TG ~ Pattern1 + Pattern2 + Pattern3, data = data)
summary(model_tg)
check_assumptions(model_tg)

### TC Model
model_tc <- lm(TC ~ Pattern1 + Pattern2 + Pattern3, data = data)
summary(model_tc)
check_assumptions(model_tc)

############################################################
# 7. ANOVA 
############################################################

### Example: Does TC differ across Pattern1 quartiles?

cat("\n========== ANOVA: TC ~ Pattern1_Q ==========\n")
anova_tc_p1 <- aov(TC ~ Pattern1_Q, data = data)
summary(anova_tc_p1)

### Bartlett test for equality of variances 
cat("\n--- Bartlett Test (Equal variances) ---\n")
bartlett.test(TC ~ Pattern1_Q, data=data)

### Boxplot for ANOVA visualization
boxplot(TC ~ Pattern1_Q, data=data, col="lightgreen",
        main="TC Differences Across Pattern 1 Quartiles",
        ylab="Total Cholesterol", xlab="Pattern 1 Quartile")

### ANOVA for LDL, HDL, TG 
anova_ldl_p1 <- aov(LDL_C ~ Pattern1_Q, data=data)
anova_hdl_p1 <- aov(HDL_C ~ Pattern1_Q, data=data)
anova_tg_p1  <- aov(TG ~ Pattern1_Q,  data=data)

summary(anova_ldl_p1)
summary(anova_hdl_p1)
summary(anova_tg_p1)


