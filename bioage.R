#安装R包
install.packages("devtools")
devtools::install_github("dayoonkwon/BioAge")
install.packages("ggeffects") 
library(haven)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2); 
library(patchwork);
library(grid); 
library(gridExtra)
library(psych)
library(splines) 
library(ggeffects)
library(haven)
names(data1)
ukb <- read_dta("D:/C-copy/data2/bage.dta")
### Change units of biomarkers
ukb$lnalp <- log(ukb$alp)
ukb$bun <- ukb$bun_mmol*2.8
ukb$lnbun <- log(ukb$bun)
ukb$lncreat_umol <- log(ukb$creat_umol)
ukb$wbc <- ukb$wbc_L
ukb$crp <- ukb$crp_mgL/10
ukb$lncrp <- log(ukb$crp)
ukb$hba1c <- ukb$hba1c_mmol*0.0915 + 2.15
ukb$lnhba1c <- log(ukb$hba1c)
ukb$totchol <- ukb$totchol_mmol*38.67
ukb$bun_mmol <- ukb$bun*2.8
ukb$albumin<-ukb$albumin_gL / 10

# 读取CSV文件
library(haven)
names(data1)
data1 <- read_dta("D:/C-copy/data2/bioage2.dta")

NHANES3 %>% 
  filter(age >= 30 & age <= 75 & pregnant == 0) %>% 
  select(sampleID, year, wave, gender, age, fev, sbp, totchol, hba1c, albumin,
         creat, lncrp, alp, crp) %>%
  na.omit() %>%describe()
#kdm using NHANES (separate training for men and women)
train1 = kdm_calc(NHANES3,biomarkers=c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","albumin_gL","alp","crp"))

#KDM using
#female
kdm_fem = kdm_calc(data = ukb %>%  
                   filter(gender == 0),  
                   biomarkers = c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","albumin_gL","alp","crp"),                   
                   fit = train1$fit$female,                   
                   s_ba2 = train1$fit$female$s_b2)
                 
# male
kdm_male = kdm_calc(data = ukb %>%                      
                      filter(gender == 1),                   
                    biomarkers = c("fev_1000","sbp","bun","hba1c","totchol","creat_umol","albumin_gL","alp","crp"),                   
                    fit = train1$fit$male,                    
                    s_ba2 = train1$fit$male$s_b2)

#merge
kdm_data = rbind(kdm_fem$data, kdm_male$data)

#phenoage using NHANES
# Descriptive statistics of the dataset used for training of PhenoAge
NHANES3 %>% 
  filter(age >= 20 & age <= 84) %>%
  select(sampleID, year, wave, gender, age, albumin_gL,lymph,mcv,glucose_mmol,
         rdw,creat_umol,lncrp,alp,wbc)%>%
  na.omit() %>%
  describe()

train2 = phenoage_calc(NHANES3,biomarkers=c("albumin_gL","lymph","mcv","glucose_mmol", "rdw","creat_umol","crp","alp","wbc"))
#phenoage计算
  phenoage_test = phenoage_calc(data = ukb,                             
  biomarkers = c("albumin_gL","lymph","mcv","glucose_mmol",
                 "rdw","creat_umol","crp","alp","wbc"),                            
  fit = train2$fit)
 phenoage_data = phenoage_test$data
 
newdata = left_join(ukb, kdm_data[, c("n_eid", "kdm", "kdm_advance")], by = "n_eid") %>%  
left_join(., phenoage_data[, c("n_eid", "phenoage","phenoage_advance")], by = "n_eid")
summary(newdata %>% select(age,kdm,phenoage,kdm_advance, phenoage_advance))

### Create BA residuals using regression model of a natural spline of CA with 3 degrees of freedom
get_BA_resids <- function(BA){
  data = newdata %>% drop_na(BA)
  model <- parse(text = sprintf("lm(%s ~ ns(age, df = 3), data = newdata)", BA)) %>% eval()
  model_predict <- ggpredict(model, terms = c("age"))
  data[,"BA_res"] <- NA
  data[!is.na(data[BA]),"BA_res"] <- resid(model)
  return(residuals(model))
}
for(BA in c("kdm","phenoage")){
  BA_res <- paste0(BA,"_res")
  newdata[,BA_res] = NA
  newdata[!is.na(newdata[BA]),BA_res] <- get_BA_resids(BA)
}
rm(list=c("BA","BA_res"))

summary(newdata %>% 
          select("kdm_res","phenoage_res", ))

write.csv(newdata, file = "newdata.csv", row.names = FALSE)
