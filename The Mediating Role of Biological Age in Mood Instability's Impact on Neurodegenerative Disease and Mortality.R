# ANALYSIS 1: Association of Mood Instability with Adverse Outcomes
# =====================================================================
library(haven)
library(survival)

# Load data
data1 <- read_dta("D:/C-copy/data/data2/final11.21.dta")
data1$ct <- as.factor(data1$category)
data1$ct <- relevel(data1$ct, ref = "3")  

# 1.1 Dementia Analysis --------------------------------------------------
# Model 1
fit1 <- coxph(Surv(survival_time1, dementia_occurred == 1) ~ ct + age + sex + ethni + edu + TOW, data = data1)
summary(fit1)
cox.zph(fit1)
# Model 2
fit2 <- coxph(Surv(survival_time1, dementia_occurred == 1) ~ ct + age + sex + ethni + edu + TOW + 
                sleep + smoke + alcoholstatus + MET_group + BMI + family_history, data = data1)
summary(fit2)
cox.zph(fit2)
# Model 3
fit3 <- coxph(Surv(survival_time1, dementia_occurred == 1) ~ ct + age + sex + ethni + edu + TOW + 
                sleep + smoke + alcoholstatus + BMI + MET_group + hypertension + diabetes + depressed, 
              data = data1)  
summary(fit3)
cox.zph(fit3)
# 1.2 Parkinson's Disease Analysis ---------------------------------------
# Model 1: Basic adjustment
fit1 <- coxph(Surv(time1, parkinson_occurred == 1) ~ ct + age + sex + ethni + edu + TOW, data = data1)
summary(fit1)
cox.zph(fit1)
# Model 2: Lifestyle adjustments
fit2 <- coxph(Surv(time1, parkinson_occurred == 1) ~ ct + age + sex + ethni + edu + TOW + 
                sleep + smoke + alcoholstatus + MET_group + BMI + family_history, data = data1)
summary(fit2)
cox.zph(fit2)
# Model 3: Full clinical adjustments
fit3 <- coxph(Surv(time1, parkinson_occurred == 1) ~ ct + age + sex + ethni + TOW + 
                sleep + smoke + alcoholstatus + BMI + family_parkin + diabetes + MET_group + hypertension + depressed,
              data = data1)  # Fixed variable name
summary(fit3)
cox.zph(fit3)
# 1.3 All-Cause Mortality Analysis ---------------------------------------
# Use consistent dataset (data1 instead of data4)
# Model 1: Basic adjustment
fit1 <- coxph(Surv(survival_time4, status == 1) ~ ct + age + sex + ethni + edu + TOW, data = data1)
summary(fit1)
cox.zph(fit1)
# Model 2: Lifestyle adjustments
fit2 <- coxph(Surv(survival_time4, status == 1) ~ ct + age + sex + ethni + edu + TOW + 
                sleep + smoke + alcoholstatus + MET_group + BMI, data = data1)
summary(fit2)
cox.zph(fit2)
# Model 3: Full clinical adjustments
fit3 <- coxph(Surv(survival_time4, status == 1) ~ ct + age + sex + ethni + TOW + 
                sleep + smoke + alcoholstatus + BMI + diabetes + MET_group + hypertension,  # Fixed duplicate
              data = data1)
summary(fit3)
cox.zph(fit3)
# =====================================================================
# ANALYSIS 2: Association of Age Acceleration with Adverse Outcomes
# =====================================================================
# 2.1 KDM Acceleration --------------------------------------------------
# Dementia
fit1 <- coxph(Surv(survival_time1, dementia_occurred == 1) ~ kdm_res + age + sex + ethni + 
                sleep + smoke + alcoholstatus + edu + TOW + BMI + MET_group + hypertension + diabetes + depressed, 
              data = data1)
summary(fit1)
cox.zph(fit1)
# Parkinson's
fit2 <- coxph(Surv(time1, parkinson_occurred == 1) ~ kdm_res + age + sex + ethni + 
                sleep + smoke + alcoholstatus + edu + TOW + BMI + MET_group + hypertension + diabetes + depressed, 
              data = data1)
summary(fit2)
cox.zph(fit2)
# Mortality (fixed syntax error)
fit3 <- coxph(Surv(survival_time4, status == 1) ~ kdm_res + age + sex + ethni + TOW + 
                smoke + alcoholstatus + BMI + diabetes + MET_group + hypertension + depressed + sleep,  # Moved sleep
              data = data1)
summary(fit3)
cox.zph(fit3)
# 2.2 PhenoAge Acceleration ---------------------------------------------
# Dementia
fit1 <- coxph(Surv(survival_time1, dementia_occurred == 1) ~ phenoage_res + age + sex + ethni + 
                sleep + smoke + alcoholstatus + edu + TOW + BMI + MET_group + hypertension + diabetes + depressed, 
              data = data1)
summary(fit1)
cox.zph(fit1)

# Parkinson's
fit2 <- coxph(Surv(time1, parkinson_occurred == 1) ~ phenoage_res + age + sex + ethni + 
                sleep + smoke + alcoholstatus + edu + TOW + BMI + MET_group + hypertension + diabetes + depressed, 
              data = data1)
summary(fit2)
cox.zph(fit2)

# Mortality (fixed syntax error)
fit3 <- coxph(Surv(survival_time4, status == 1) ~ phenoage_res + age + sex + ethni + TOW + 
                smoke + alcoholstatus + BMI + diabetes + MET_group + hypertension + depressed + sleep,
              data = data1)
summary(fit3)
cox.zph(fit3)
# =====================================================================
# ANALYSIS 3: Mediation Analysis 
# =====================================================================
library(CMAverse)
set.seed(111)
# Prepare analysis dataset with consistent naming
analysis_data <- data1[, c(
  "survival_time1", "dementia_occurred", "ct", "phenoage_res", "kdm_res", "age",
  "sex", "ethni", "TOW", "sleep", "smoke", "alcoholstatus",   # Renamed gender->sex
  "BMI", "family_parkin", "family_dementia", "diabetes", "MET_group", 
  "hypertension", "depressed", "time1", "parkinson_occurred",  # Standardized names
  "survival_time4", "status"
)]

# Create exposure indicators
analysis_data$D1 <- ifelse(analysis_data$ct == "1", 1, 0)  # Severe instability
analysis_data$D2 <- ifelse(analysis_data$ct == "2", 1, 0)  # Mild instability
analysis_data$D3 <- ifelse(analysis_data$ct == "3", 1, 0)  # Reference

# Define common covariates
base_covariates <- c("age", "sex", "ethni", "TOW", "sleep", "smoke", 
                     "alcoholstatus", "BMI", "diabetes", "MET_group", 
                     "hypertension", "depressed")

# 3.1-3.3: SEVERE MOOD INSTABILITY (D1) ---------------------------------
data_severe <- subset(analysis_data, ct %in% c("1", "3"))

# 3.1: Dementia Outcome (Severe)
med_dementia_phenoage_sev <- cmest(
  data = data_severe, 
  model = "rb", 
  outcome = "survival_time1",
  event = "dementia_occurred",
  exposure = "D1",
  mediator = "phenoage_res", 
  basec = c(base_covariates, "family_dementia"),
  mreg = list("linear"), 
  yreg = "coxph", 
  EMint = FALSE, 
  mval = list(0),
  inference = "bootstrap",
  nboot = 1000, 
  boot.ci.type = "bca"
)

# KDM version (similar structure, change mediator)

# 3.2: Parkinson's Outcome (Severe)
med_parkinson_phenoage_sev <- cmest(
  data = data_severe, 
  model = "rb", 
  outcome = "time1",
  event = "parkinson_occurred",
  exposure = "D1",
  mediator = "phenoage_res", 
  basec = c(base_covariates, "family_parkin"),
  mreg = list("linear"), 
  yreg = "coxph", 
  EMint = FALSE, 
  mval = list(0), 
  inference = "bootstrap",
  nboot = 1000, 
  boot.ci.type = "bca"
)

# 3.3: Mortality Outcome (Severe)
med_mortality_phenoage_sev <- cmest(
  data = data_severe, 
  model = "rb", 
  outcome = "survival_time4",
  event = "status",
  exposure = "D1",
  mediator = "phenoage_res", 
  basec = base_covariates,  
  mreg = list("linear"), 
  yreg = "coxph", 
  EMint = FALSE, 
  mval = list(0), 
  inference = "bootstrap",
  nboot = 1000, 
  boot.ci.type = "bca"
)

# 3.4-3.6: MILD MOOD INSTABILITY (D2) -----------------------------------
data_mild <- subset(analysis_data, ct %in% c("2", "3"))

# 3.4: Dementia Outcome (Mild)
med_dementia_phenoage_mild <- cmest(
  data = data_mild,  
  model = "rb", 
  outcome = "survival_time1",
  event = "dementia_occurred",
  exposure = "D2",    
  mediator = "phenoage_res", 
  basec = c(base_covariates, "family_dementia"),
  mreg = list("linear"), 
  yreg = "coxph", 
  EMint = FALSE, 
  mval = list(0),
  inference = "bootstrap",
  nboot = 1000, 
  boot.ci.type = "bca"
)

# 3.5: Parkinson's Outcome (Mild)
med_parkinson_phenoage_mild <- cmest(
  data = data_mild,  
  model = "rb", 
  outcome = "time1",
  event = "parkinson_occurred",
  exposure = "D2",   
  mediator = "phenoage_res", 
  basec = c(base_covariates, "family_parkin"),
  mreg = list("linear"), 
  yreg = "coxph", 
  EMint = FALSE, 
  mval = list(0), 
  inference = "bootstrap",
  nboot = 1000, 
  boot.ci.type = "bca"
)

# 3.6: Mortality Outcome (Mild)
med_mortality_phenoage_mild <- cmest(
  data = data_mild,   
  model = "rb", 
  outcome = "survival_time4",
  event = "status",
  exposure = "D2",    # Mild exposure
  mediator = "phenoage_res", 
  basec = base_covariates,
  mreg = list("linear"), 
  yreg = "coxph", 
  EMint = FALSE, 
  mval = list(0), 
  inference = "bootstrap",
  nboot = 1000, 
  boot.ci.type = "bca"
)