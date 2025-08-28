library(stats)
library(knitr)
library(ggplot2)
library(magick)
library(forecast)
library(fpp2)
library(GGally)
library(gridExtra)
library(knitr)
library(patchwork)
library(BSDA)
library(dplyr)
library(GLMsData)
library(MedDataSets)
library(NHSRdatasets)
library(medicaldata)
library(predictmeans)
library(tidyverse) 
library(MASS)
library(nnet)
library(survival)
library(VGAM) 
library(mlbench)
library(ResourceSelection)
library(pROC)
library(caret)
library(datasets) 
library(foreign)
library(brant)
library(svyVGAM)
library(pscl)
library(tibble)
library(cards)
library(cardx)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(survminer)
library(pwr)
library(blockrand)
library(randomizeR)


data(package = "MedDataSets")
library(MedDataSets)
?Melanoma_df
view(Melanoma_df)


summary(Melanoma_df)
# This code computes summary statistics for the entire data, indicating spread and more

cor(Melanoma_df)
# This code computes a Pearson correlation coefficient matrix of all variables

ggplot(Melanoma_df, aes(x = factor(sex, levels = c(0, 1), 
                                   labels = c("Female", "Male")))) +
  geom_bar(fill = "grey") +
  geom_text(stat = "count", aes(label = ..count..), vjust = 1.5, color = "black") +
  labs(title = "Frequency Bar Plot for Patient Sex",
       x = "Patient Sex",
       y = "Frequency")
ggplot(Melanoma_df, aes(x = factor(status, levels = c(1, 2, 3), 
                                   labels = c("Died from Melanoma", "Alive", "Died from Other Causes")))) +
  geom_bar(fill = "grey") +
  geom_text(stat = "count", aes(label = ..count..), vjust = 1.5, color = "black") +
  labs(title = "Frequency Bar Plot for Patient Status",
       x = "Patient Status",
       y = "Frequency")
ggplot(Melanoma_df, aes(x = factor(ulcer, levels = c(1, 0), 
                                   labels = c("Presence of Ulceration", "No Presence of Ulceration")))) +
  geom_bar(fill = "grey") +
  geom_text(stat = "count", aes(label = ..count..), vjust = 1.5, color = "black") +
  labs(title = "Frequency Bar Plot for Ulceration Presence",
       x = "Ulceration Presence",
       y = "Frequency")
# This code generates bar plots and visualisations for the three categorical variables

ggplot(Melanoma_df, aes(x = year, y = time)) +
  geom_point(alpha = 0.7) + 
  labs(title = "Scatter Plot", x = "Year of Diagnosis", y = "Survival Time of Patient") +
  theme_minimal()
ggplot(Melanoma_df, aes(x = age, y = time)) +
  geom_point(alpha = 0.7) + 
  labs(title = "Scatter Plot", x = "Age of Patient", y = "Survival Time of Patient") +
  theme_minimal()
ggplot(Melanoma_df, aes(x = thickness, y = time)) +
  geom_point(alpha = 0.7) + 
  labs(title = "Scatter Plot", x = "Thickness of Melanoma", y = "Survival Time of Patient") +
  theme_minimal()
# This code generates scatter plots and visualisations of correlation between some variables

ggplot(Melanoma_df, aes(y = time)) +
  geom_boxplot() +
  labs(title = "Boxplot of Survival Time", y = "Survival Time of Patients") +
  theme_minimal()
ggplot(Melanoma_df, aes(y = age)) +
  geom_boxplot() +
  labs(title = "Boxplot of Patient Age", y = "Patient Age") +
  theme_minimal()
ggplot(Melanoma_df, aes(y = year)) +
  geom_boxplot() +
  labs(title = "Boxplot of Year of Diagnosis", y = "Year of Diagnosis") +
  theme_minimal()
ggplot(Melanoma_df, aes(y = thickness)) +
  geom_boxplot() +
  labs(title = "Boxplot for Thickness of Melanoma", y = "Thickness of Melanoma") +
  theme_minimal()
# This code generates box plots and visualisations of spread of continuous variables

ggplot(Melanoma_df, aes(x = time)) +
  geom_histogram(binwidth = 80, color = "black", fill = "grey") +
  labs(title = "Histogram of Survival Time",
       x = "Survival Time of Patients",
       y = "Frequency") +
  theme_minimal()
ggplot(Melanoma_df, aes(x = age)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey") +
  labs(title = "Histogram of of Patient Age",
       x = "Patient Age",
       y = "Frequency") +
  theme_minimal()
ggplot(Melanoma_df, aes(x = thickness)) +
  geom_histogram(binwidth = 0.5, color = "black", fill = "grey") +
  labs(title = "Histogram of Thickness of Melanoma",
       x = "Thickness of Melanoma",
       y = "Frequency") +
  theme_minimal()
# This code generates histograms and visualisations of spread of continuous variables

q3time <- quantile(Melanoma_df$time, 0.75)
q1time <- quantile(Melanoma_df$time, 0.25)
iqrtime <- q3time - q1time
iqrtime
upperboundtime <- q3time + 1.5*iqrtime
upperboundtime
outliers_uppertime <- Melanoma_df$time[Melanoma_df$time > upperboundtime]
outliers_uppertime
lowerboundtime <- q1time - 1.5*iqrtime
lowerboundtime
outliers_lowertime <- Melanoma_df$time[Melanoma_df$time < lowerboundtime]
outliers_lowertime
Melanoma_df_nooutliers1 <- subset(Melanoma_df, Melanoma_df$time > (q1time - 1.5*iqrtime) & Melanoma_df$time < (q3time + 1.5*iqrtime))
dim(Melanoma_df_nooutliers1)
count(Melanoma_df_nooutliers1)
# Upon visual inspection, outliers are suspected for 'time'
# Therefore, observations above or below 1.5*IQR of upper or lower bound are removed from data
# Result of removal of outliers is stored in new dataframe, titled 'Melanoma_df_nooutliers1'
# 'count(Melanoma_df_nooutliers1)' confirms removal of 1, leading to reduced 204 observations

q3age <- quantile(Melanoma_df_nooutliers1$age, 0.75)
q1age <- quantile(Melanoma_df_nooutliers1$age, 0.25)
iqrage <- q3age - q1age
upperboundage <- q3age + 1.5*iqrage
upperboundage
outliers_upperage <- Melanoma_df$age[Melanoma_df_nooutliers1$age > upperboundage]
outliers_upperage
lowerboundage <- q1age - 1.5*iqrage
lowerboundage
outliers_lowerage <- Melanoma_df_nooutliers1$age[Melanoma_df_nooutliers1$age < lowerboundage]
outliers_lowerage
Melanoma_df_nooutliers2 <- subset(Melanoma_df_nooutliers1, Melanoma_df_nooutliers1$age > (q1age - 1.5*iqrage) & Melanoma_df_nooutliers1$age < (q3age + 1.5*iqrage))
dim(Melanoma_df_nooutliers2)
count(Melanoma_df_nooutliers2)
# Upon visual inspection, outliers are suspected for 'age'
# Therefore, observations above or below 1.5*IQR of upper or lower bound are removed from data
# Result of removal of outliers is stored in new dataframe, titled 'Melanoma_df_nooutliers2'
# 'count(Melanoma_df_nooutliers2)' confirms removal of 1, leading to reduced 203 observations

q3thickness <- quantile(Melanoma_df_nooutliers2$thickness, 0.75)
q1thickness <- quantile(Melanoma_df_nooutliers2$thickness, 0.25)
iqrthickness <- q3thickness - q1thickness
upperboundthickness <- q3thickness + 1.5*iqrthickness
upperboundthickness
outliers_upperthickness <- Melanoma_df_nooutliers2$thickness[Melanoma_df_nooutliers2$thickness > upperboundthickness]
outliers_upperthickness
lowerboundthickness <- q1thickness - 1.5*iqrthickness
lowerboundthickness
outliers_lowerthickness <- Melanoma_df_nooutliers2$thickness[Melanoma_df_nooutliers2$thickness < lowerboundthickness]
outliers_lowerthickness
Melanoma_df_cleaned <- subset(Melanoma_df_nooutliers2, Melanoma_df_nooutliers2$thickness > (q1thickness - 1.5*iqrthickness) & Melanoma_df_nooutliers2$thickness < (q3thickness + 1.5*iqrthickness))
dim(Melanoma_df_cleaned)
count(Melanoma_df_cleaned)
# Upon visual inspection, outliers are suspected for 'time'
# Therefore, observations above or below 1.5*IQR of upper or lower bound are removed from data
# Result of removal of outliers is stored in new dataframe, titled 'Melanoma_df_cleaned'
# 'count(Melanoma_df_cleaned)' confirms removal of 1, leading to reduced 190 observations

Melanoma_df_cleaned$status <- ifelse(Melanoma_df_cleaned$status == 1, 1, 0)
view(Melanoma_df_cleaned)
# For the analysis, one is interested in modelling death from a specific cause (melanoma)
# Therefore, 'status' column is adjusted; 1 is 'Death from Melanoma' and 0 is all other (censored)

ggplot(Melanoma_df_cleaned, aes(x = factor(status, levels = c(1, 0), 
                                   labels = c("Died from Melanoma", "Alive or Died from Other Causes")))) +
  geom_bar(fill = "grey") +
  geom_text(stat = "count", aes(label = ..count..), vjust = 1.5, color = "black") +
  labs(title = "Frequency Bar Plot for Patient Status",
       x = "Patient Status",
       y = "Frequency")
summary(Melanoma_df_cleaned)
# This code confirms the update of all variables, alongside bi-segmentation of 'status' 

cor(Melanoma_df_cleaned)
ggplot(Melanoma_df_cleaned, aes(x = year, y = time)) +
  geom_point(alpha = 0.7) + 
  labs(title = "Scatter Plot", x = "Year of Diagnosis", y = "Survival Time of Patient") +
  theme_minimal()
ggplot(Melanoma_df_cleaned, aes(x = age, y = time)) +
  geom_point(alpha = 0.7) + 
  labs(title = "Scatter Plot", x = "Age of Patient", y = "Survival Time of Patient") +
  theme_minimal()
ggplot(Melanoma_df_cleaned, aes(x = thickness, y = time)) +
  geom_point(alpha = 0.7) + 
  labs(title = "Scatter Plot", x = "Thickness of Melanoma", y = "Survival Time of Patient") +
  theme_minimal()
# This code re-generates scatter plots and visualisations of correlation between some variables



surv_object <- Surv(time = Melanoma_df_cleaned$time, event = Melanoma_df_cleaned$status == 1)
# First, one creates a survival object
# Here, 'time' is survival time, and 'status' is event indicator (1 = death, 0 = censored)

log_rank_testsex <- survdiff(surv_object ~ sex, data = Melanoma_df_cleaned)
log_rank_testsex
# Next, survival times can be compared across 'sex', where 1 is male and 0 is female
# Here, a log-rank test is performed for survival differences according to 'sex'
log_rank_testulcer <- survdiff(surv_object ~ ulcer, data = Melanoma_df_cleaned)
log_rank_testulcer
# Next, survival times can be compared across 'ulcer', where 1 is presence of ulceration
# Here, a log-rank test is performed for survival differences according to 'ulcer'

fit1 <- coxph(Surv(time = Melanoma_df_cleaned$time, 
                   event = Melanoma_df_cleaned$status) ~ sex + ulcer, 
              data = Melanoma_df_cleaned)
summary(fit1)

fit2 <- coxph(Surv(time = Melanoma_df_cleaned$time, 
                   event = Melanoma_df_cleaned$status) ~ sex + ulcer + age, 
              data = Melanoma_df_cleaned)
summary(fit2)

fit3 <- coxph(Surv(time = Melanoma_df_cleaned$time, 
                   event = Melanoma_df_cleaned$status) ~ sex + ulcer + age + thickness, 
              data = Melanoma_df_cleaned)
summary(fit3)
# Next, three Cox proportional hazards models are fitted to the cleaned dataframe

anova(fit1, fit2, fit3, test = 'LRT')
# This performs a Likelihood Ratio Test (LRT) to compare the three nested Cox PH models
# Pr(>|Chi|) = 0.0102; result supports fit3, containing all predictors excluding 'year'   

cox.zph(fit3)
ggcoxzph(cox.zph(fit3)) 
# This runs the Schoenfeld residuals test for the proportional hazards assumption
# The Schoenfeld residuals test tests whether effect of each covariate changes over time
# p-values smaller than 0.05 suggest violation of the proportional hazards assumption
# Global p-value of 0.0708 suggests no evidence of PH violation, but p = 0.0065 for 'thickness'

ggsurvplot(
  survfit(Surv(time, status) ~ sex, data = Melanoma_df_cleaned),
  fun = "cloglog", 
  palette = "Dark2",
  legend.title = "Patient Sex",
  legend.labs = c("Male", "Female"),
  xlab = "Time (Log-Scale)",
  ylab = "log(-log(Survival))",
  ggtheme = theme_minimal()
)
ggsurvplot(
  survfit(Surv(time, status) ~ ulcer, data = Melanoma_df_cleaned),
  fun = "cloglog",  
  palette = "Dark2",
  legend.title = "Patient Ulceration Presence",
  legend.labs = c("Ulceration Presence", "No Ulceration Presence"),
  xlab = "Time (Log-Scale)",
  ylab = "log(-log(Survival))",
  ggtheme = theme_minimal()
)
# Next, log-log survival plots (cumulative hazard) is computed for categorical variables

coxfit3 <- coxph(Surv(time = Melanoma_df_cleaned$time, 
                   event = Melanoma_df_cleaned$status) ~ sex + ulcer + age + thickness, 
                 method = "breslow", data = Melanoma_df_cleaned)
beta <- coef(coxfit3)   
base_hazard <- basehaz(coxfit3, centered = FALSE) 
base_hazard
head(base_hazard$time)
beta
sex_value1 <- 1
ulcer_value1 <- 1
age_value1 <- 50 
thickness_value1 <- 4.0
hazard_ratio1 <- unname(exp(
  beta["sex"] * sex_value1 +
    beta["ulcer"] * ulcer_value1 +
    beta["age"] * age_value1 +
    beta["thickness"] * thickness_value1
))
hazard_ratio1
summary(Melanoma_df_cleaned)
t = 2166
baseline_hazard_2166 <- base_hazard$hazard[which.min(abs(base_hazard$time - t))]
baseline_hazard_2166
predicted_hazard <- baseline_hazard_2166 * hazard_ratio1
survival_probability <- exp(-predicted_hazard)
survival_probability
print(predicted_hazard)
# Next, absolute risk prediction and relative risk can be computed
# Hazard ratio and absolute risk prediction is computed for 50-year old male, plus other parameters
# 'breslow' is the applied method for handling ties in specific event times
# 'coxfit3$means' yields mean values of coviarates in Cox model, needed for baseline hazards
# 'hazard_ratio1' is the exponentiated linear predictor for a specific individual
# The hazard ratio is relative to the baseline individuals, assuming all covariates are 0 
# The data's average survival time (2166) is used in the baseline hazard calculation
# The absolute hazard is calculated by multiplying baseline hazard (t = 2166) by hazard ratio

coxfit3 <- coxph(Surv(time, status) ~ sex + ulcer + age + thickness, 
                 data = Melanoma_df_cleaned)
pred_dat <- data.frame(sex = c(0, 1), ulcer = c(1, 0), age = c(35, 52), thickness = c(4.5, 6.1))
surv_obj <- survfit(coxfit3, newdata = pred_dat)
summary(surv_obj, times = 2166)
# To conceptualise survival analysis and prediction, a dataframe of two individuals is created
# Individual 1; sex = 0, ulcer = 1, age = 35, thickness = 4.5
# Individual 2; sex = 1, ulcer = 0, age = 52, thickness = 6.1
# The code computes survival curves for both individuals using the Cox model
# Here, survival probabilities at time (month) 2166 are extracted for each individual
# 'n.risk' is the number of individuals at risk at time (month) 2166 in the overall data
# 'n.event' is the number of events (deaths) at time (month) 2166
# 'survival1' is the estimated survival probability at time (month) 2166 for individual 1
# 'survival2' is the estimated survival probability at time (month) 2166 for individual 2

set.seed(123)
pred_dat2 <- data.frame(
  sex = sample(0:1, 100, replace = TRUE),         
  ulcer = sample(0:1, 100, replace = TRUE),     
  age = round(runif(100, 20, 80)),                
  thickness = round(runif(100, 0.5, 10.0), 1)     
)
surv_obj <- survfit(coxfit3, newdata = pred_dat2)
surv_summary_final <- summary(surv_obj, times = 2166)
surv_probs_final <- surv_summary_final$surv
head(surv_probs_final)
view(pred_dat2)
view(surv_probs_final) 
pred_dat2$survivalprobability <- as.numeric(t(surv_probs_final))
pred_dat2$surv_2166 <- NULL
view(pred_dat2)
# Here, a random selection of 100 individuals, with varying values of predictors, is generated
# Sex, ulcer, age and thickness could take particular values as defined in the code
# Then, 'survfit' applies the fitted Cox PH model to the randomised 100 individuals
# Survival probabilities are stored in 'surv_probs_final' and merged with dataframe 'pred_dat2'

subset_fit <- survfit(coxfit3, newdata = pred_dat2[1:10, ])
plot(subset_fit, col = 1:10, lty = 1:10,
     xlab = "Time (Months)", ylab = "Survival Probability",
     main = "Survival Curves for 10 Individuals of Simulation")
abline(v = 2166, col = "blue", lty = 2)

fit_complete <- survfit(coxfit3, newdata = pred_dat2)
plot(fit_complete, col = rainbow(100), lty = 1, xlab = "Time (Months)", ylab = "Survival Probability",
     main = "Survival Curves for 100 Simulated Individuals")
abline(v = 2166, col = "blue", lty = 2)
# This plots survival curves (decreasing in probability) for 10 or 100 (random) individuals
# The downward slope suggests individuals have decreasing survival from melanoma over time
# Some curves drop steeply, indicating higher risk, worse prognosis, and higher hazard
# There is wide spread among the curves, indicating covariates have strong effect on survival

