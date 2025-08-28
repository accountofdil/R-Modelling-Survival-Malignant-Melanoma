# MSMM - Data Loading & Integration

options(warn = -1)
suppressPackageStartupMessages({
  library(stats)
  library(knitr)
  library(ggplot2)
  library(magick)
  library(forecast)
  library(fpp2)
  library(GGally)
  library(gridExtra)
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
  library(randomizeR)})

data(package = "MedDataSets")
library(MedDataSets)
?Melanoma_df
view(Melanoma_df)



# MSMM - Data Preprocessing & Exploratory Data Analysis

summary(Melanoma_df)
cor(Melanoma_df)

library(ggplot2)

ggplot(Melanoma_df, aes(x = factor(sex, levels = c(0, 1), 
                                   labels = c("Female", "Male")))) +
  geom_bar(fill = "grey") +
  geom_text(stat = "count", aes(label = ..count..), 
            vjust = 1.5, color = "black", family = "Garamond") +
  labs(title = "Distribution of Patient Sex (Bar Chart)",
       x = "Patient Sex",
       y = "Frequency") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(x = factor(status, levels = c(1, 2, 3), 
                                   labels = c("Died from Melanoma", "Alive", "Died from Other Causes")))) +
  geom_bar(fill = "grey") +
  geom_text(stat = "count", aes(label = ..count..), 
            vjust = 1.5, color = "black", family = "Garamond") +
  labs(title = "Distribution of Patient Status (Bar Chart)",
       x = "Patient Status",
       y = "Frequency") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),  
    axis.title = element_text(size = 13, face = "bold"),  
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(x = factor(ulcer, levels = c(1, 0), 
                                   labels = c("Presence of Ulceration", "No Presence of Ulceration")))) +
  geom_bar(fill = "grey") +
  geom_text(stat = "count", aes(label = ..count..), 
            vjust = 1.5, color = "black", family = "Garamond") +
  labs(title = "Distribution of Ulceration Presence (Bar Chart)",
       x = "Ulceration Presence",
       y = "Frequency") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),   
    axis.title = element_text(size = 13, face = "bold"),  
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(x = year, y = time)) +
  geom_point(alpha = 0.7) + 
  labs(title = "Year of Diagnosis vs Survival Time (Scatter Plot)",
       x = "Year of Diagnosis",
       y = "Survival Time of Patient") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(x = age, y = time)) +
  geom_point(alpha = 0.7) + 
  labs(title = "Age vs Survival Time (Scatter Plot)",
       x = "Age of Patient",
       y = "Survival Time of Patient") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(x = thickness, y = time)) +
  geom_point(alpha = 0.7) + 
  labs(title = "Melanoma Thickness vs Survival Time (Scatter Plot)",
       x = "Thickness of Melanoma",
       y = "Survival Time of Patient") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(y = time)) +
  geom_boxplot() +
  labs(title = "Distribution of Survival Time (Box Plot)",
       y = "Survival Time of Patients") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(y = age)) +
  geom_boxplot() +
  labs(title = "Distribution of Age (Box Plot)",
       y = "Patient Age") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(y = year)) +
  geom_boxplot() +
  labs(title = "Distribution of Diagnosis Year (Box Plot)",
       y = "Year of Diagnosis") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(y = thickness)) +
  geom_boxplot() +
  labs(title = "Distribution of Melanoma Thickness (Box Plot)",
       y = "Thickness of Melanoma") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(x = time)) +
  geom_histogram(binwidth = 80, color = "black", fill = "grey") +
  labs(title = "Distribution of Survival Time (Histogram)",
       x = "Survival Time of Patients",
       y = "Frequency") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(x = age)) +
  geom_histogram(binwidth = 1, color = "black", fill = "grey") +
  labs(title = "Distribution of Age (Histogram)",
       x = "Patient Age",
       y = "Frequency") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

ggplot(Melanoma_df, aes(x = thickness)) +
  geom_histogram(binwidth = 0.5, color = "black", fill = "grey") +
  labs(title = "Distribution of Melanoma Thickness (Histogram)",
       x = "Thickness of Melanoma",
       y = "Frequency") +
  theme_minimal(base_family = "Garamond") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text  = element_text(size = 12, family = "Garamond"))

q3time <- quantile(Melanoma_df$time, 0.75)
q1time <- quantile(Melanoma_df$time, 0.25)
iqrtime <- q3time - q1time
iqrtime
upperboundtime <- q3time + 1.5 * iqrtime
upperboundtime
outliers_uppertime <- Melanoma_df$time[Melanoma_df$time > upperboundtime]
outliers_uppertime
lowerboundtime <- q1time - 1.5 * iqrtime
lowerboundtime
outliers_lowertime <- Melanoma_df$time[Melanoma_df$time < lowerboundtime]
outliers_lowertime
Melanoma_df_nooutliers1 <- subset(
  Melanoma_df,
  Melanoma_df$time > (q1time - 1.5 * iqrtime) &
    Melanoma_df$time < (q3time + 1.5 * iqrtime))
dim(Melanoma_df_nooutliers1)
count(Melanoma_df_nooutliers1)

q3age <- quantile(Melanoma_df_nooutliers1$age, 0.75)
q1age <- quantile(Melanoma_df_nooutliers1$age, 0.25)
iqrage <- q3age - q1age
upperboundage <- q3age + 1.5 * iqrage
upperboundage
outliers_upperage <- Melanoma_df$age[Melanoma_df_nooutliers1$age > upperboundage]
outliers_upperage
lowerboundage <- q1age - 1.5 * iqrage
lowerboundage
outliers_lowerage <- Melanoma_df_nooutliers1$age[Melanoma_df_nooutliers1$age < lowerboundage]
outliers_lowerage
Melanoma_df_nooutliers2 <- subset(Melanoma_df_nooutliers1, Melanoma_df_nooutliers1$age > (q1age - 1.5 * iqrage) & Melanoma_df_nooutliers1$age < (q3age + 1.5 * iqrage))
dim(Melanoma_df_nooutliers2)
count(Melanoma_df_nooutliers2)

q3thickness <- quantile(Melanoma_df_nooutliers2$thickness, 0.75)
q1thickness <- quantile(Melanoma_df_nooutliers2$thickness, 0.25)
iqrthickness <- q3thickness - q1thickness
upperboundthickness <- q3thickness + 1.5 * iqrthickness
upperboundthickness
outliers_upperthickness <- Melanoma_df_nooutliers2$thickness[Melanoma_df_nooutliers2$thickness > upperboundthickness]
outliers_upperthickness
lowerboundthickness <- q1thickness - 1.5 * iqrthickness
lowerboundthickness
outliers_lowerthickness <- Melanoma_df_nooutliers2$thickness[Melanoma_df_nooutliers2$thickness < lowerboundthickness]
outliers_lowerthickness
Melanoma_df_cleaned <- subset(Melanoma_df_nooutliers2, Melanoma_df_nooutliers2$thickness > (q1thickness - 1.5 * iqrthickness) & Melanoma_df_nooutliers2$thickness < (q3thickness + 1.5 * iqrthickness))
dim(Melanoma_df_cleaned)
count(Melanoma_df_cleaned)



# MSMM - Methodology

Melanoma_df_cleaned$status <- ifelse(Melanoma_df_cleaned$status == 1, 1, 0)
view(Melanoma_df_cleaned)
surv_object <- Surv(time = Melanoma_df_cleaned$time, event = Melanoma_df_cleaned$status == 1)


## Methodology - Log-Rank Test

log_rank_testsex <- survdiff(surv_object ~ sex, data = Melanoma_df_cleaned)
log_rank_testsex
log_rank_testulcer <- survdiff(surv_object ~ ulcer, data = Melanoma_df_cleaned)
log_rank_testulcer


## Methodology - Cox Proportional Hazards Regression Model

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


## Methodology - Proportional Hazards Assumption & Schoenfeld Residual Test

library(survival)
zph <- cox.zph(fit3)
print(zph)


## Methodology - Survival Analysis

library(showtext)
font_add_google("EB Garamond", "ebgaramond")
showtext_auto(TRUE)
fit <- survfit(Surv(time, status) ~ sex, data = Melanoma_df_cleaned)
p <- ggsurvplot(
  fit,
  fun = "cloglog",
  palette = "Dark2",
  legend.title = "Patient Sex",
  legend.labs  = c("Male", "Female"),
  xlab = "Time (Log-Scale)",
  ylab = "log(-log(Survival))",
  title = "Log-Log Survival Curves per Patient Sex",
  legend = "bottom",
  ggtheme = theme_minimal(base_family = "ebgaramond"))
p$plot <- p$plot +
  theme(
    plot.title   = element_text(size = 17, face = "bold", family = "ebgaramond"),
    axis.title   = element_text(size = 13, face = "bold", family = "ebgaramond"),
    axis.text    = element_text(size = 11, family = "ebgaramond"),
    legend.title = element_text(size = 12, face = "bold", family = "ebgaramond"),
    legend.text  = element_text(size = 11, family = "ebgaramond"))
print(p)

library(survival)
library(survminer)
library(ggplot2)
library(showtext)
font_add_google("EB Garamond", "ebgaramond")
showtext_auto(TRUE)
fit_ulcer <- survfit(Surv(time, status) ~ ulcer, data = Melanoma_df_cleaned)
p_ulcer <- ggsurvplot(
  fit_ulcer,
  fun = "cloglog",
  palette = "Dark2",
  legend.title = "Patient Ulceration Presence",
  legend.labs  = c("Ulceration Presence", "No Ulceration Presence"),
  xlab = "Time (Log-Scale)",
  ylab = "log(-log(Survival))",
  title = "Log-Log Survival Curves per Ulceration Presence",
  legend = "bottom",
  ggtheme = theme_minimal(base_family = "ebgaramond"))
p_ulcer$plot <- p_ulcer$plot +
  theme(
    plot.title   = element_text(size = 17, face = "bold", family = "ebgaramond"),
    axis.title   = element_text(size = 13, face = "bold", family = "ebgaramond"),
    axis.text    = element_text(size = 11, family = "ebgaramond"),
    legend.title = element_text(size = 12, face = "bold", family = "ebgaramond"),
    legend.text  = element_text(size = 11, family = "ebgaramond"))
print(p_ulcer)

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
    beta["thickness"] * thickness_value1))
hazard_ratio1
summary(Melanoma_df_cleaned)
t = 2166
baseline_hazard_2166 <- base_hazard$hazard[which.min(abs(base_hazard$time - t))]
baseline_hazard_2166
predicted_hazard <- baseline_hazard_2166 * hazard_ratio1
survival_probability <- exp(-predicted_hazard)
survival_probability
print(predicted_hazard)

coxfit3 <- coxph(Surv(time, status) ~ sex + ulcer + age + thickness, 
                 data = Melanoma_df_cleaned)
pred_dat <- data.frame(sex = c(0, 1), ulcer = c(1, 0), age = c(35, 52), thickness = c(4.5, 6.1))
surv_obj <- survfit(coxfit3, newdata = pred_dat)
summary(surv_obj, times = 2166)

set.seed(123)
pred_dat2 <- data.frame(
  sex = sample(0:1, 100, replace = TRUE),         
  ulcer = sample(0:1, 100, replace = TRUE),     
  age = round(runif(100, 20, 80)),                
  thickness = round(runif(100, 0.5, 10.0), 1))
surv_obj <- survfit(coxfit3, newdata = pred_dat2)
surv_summary_final <- summary(surv_obj, times = 2166)
surv_probs_final <- surv_summary_final$surv
head(surv_probs_final)
view(pred_dat2)
view(surv_probs_final) 
pred_dat2$survivalprobability <- as.numeric(t(surv_probs_final))
pred_dat2$surv_2166 <- NULL
view(pred_dat2)

library(showtext)
font_add_google("EB Garamond", "ebgaramond")
showtext_auto(TRUE)
par(
  family = "ebgaramond",  
  font.main = 2,         
  font.lab = 2,          
  cex.main = 1.5,          
  cex.lab = 1.1,           
  cex.axis = 0.9)
fit_complete <- survfit(coxfit3, newdata = pred_dat2)
plot(
  fit_complete,
  col = rainbow(100),
  lty = 1,
  xlab = "Time (Months)",
  ylab = "Survival Probability",
  main = "Survival Curves for 100 Simulated Individuals")
abline(v = 2166, col = "blue", lty = 2)
par(family = "", font.main = 1, font.lab = 1, cex.main = 1, cex.lab = 1, cex.axis = 1)