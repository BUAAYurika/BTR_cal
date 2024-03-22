library(lavaan)

# code for Figure 3d
vardata <- read.csv('BTRData_PRECISE.csv')
model <-'
VRF ~ Age 
AS ~ Age 
BTR ~ VRF + AS + Age
volume ~ VRF + AS + Age
MoCA ~ BTR + AS + VRF + Age + volume
'
fit <- sem(model, data=vardata, se="bootstrap", bootstrap = 1000)
fitdata <- summary(fit,fit.measures=TRUE)
pardata <- parameterEstimates(fit,standardized = TRUE)
stdbeta <- standardizedsolution(fit)
write.csv(fitdata$fit,file = "finalRe_fit.csv")
write.csv(pardata,file = "finalRe_coef.csv")
write.csv(stdbeta,file = "finalRe_stdcoef.csv")

####################################################
# code for Figure S8
vardata <- read.csv('BTRData_PRECISE.csv')
model <-'
VRF ~ Age 
AS ~ Age 
BTR ~ VRF + AS + Age
volume ~ VRF + AS + Age
MoCA ~ BTR + AS + VRF + Age + volume + Sex
'
fit <- sem(model, data=vardata, se="bootstrap", bootstrap = 1000)
fitdata <- summary(fit,fit.measures=TRUE)
pardata <- parameterEstimates(fit,standardized = TRUE)
stdbeta <- standardizedsolution(fit)
write.csv(fitdata$fit,file = "finalRe_fit.csv")
write.csv(pardata,file = "finalRe_coef.csv")
write.csv(stdbeta,file = "finalRe_stdcoef.csv")

####################################################
# code for Figure S9a
vardata <- read.csv('BTRData_PRECISE.csv')
model <-'
VRF ~ Age 
BTR ~ VRF + Age
volume ~ VRF + Age
MoCA ~ BTR + VRF + Age + volume 
'
fit <- sem(model, data=vardata, se="bootstrap", bootstrap = 1000)
fitdata <- summary(fit,fit.measures=TRUE)
pardata <- parameterEstimates(fit,standardized = TRUE)
stdbeta <- standardizedsolution(fit)
write.csv(fitdata$fit,file = "finalRe_fit.csv")
write.csv(pardata,file = "finalRe_coef.csv")
write.csv(stdbeta,file = "finalRe_stdcoef.csv")

####################################################
# code for Figure S9b
vardata <- read.csv('BTRData_MAS.csv')
model <-'
VRF ~ Age 
BTR ~ VRF + Age
volume ~ VRF + Age
MMSE ~ BTR + VRF + Age + volume 
'
fit <- sem(model, data=vardata, se="bootstrap", bootstrap = 1000)
fitdata <- summary(fit,fit.measures=TRUE)
pardata <- parameterEstimates(fit,standardized = TRUE)
stdbeta <- standardizedsolution(fit)
write.csv(fitdata$fit,file = "finalRe_fit.csv")
write.csv(pardata,file = "finalRe_coef.csv")
write.csv(stdbeta,file = "finalRe_stdcoef.csv")
