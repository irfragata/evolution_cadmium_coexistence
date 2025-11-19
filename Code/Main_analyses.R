#' ---
#' title: "R Notebook"
#' output:
#'   html_document:
#'     df_print: paged
#' ---
rm(list=ls())

# Note the script should be run from the main directory of the repository.

# Checking if we are inside the correct folder, otherwise changing directory to the previous folder
file_exists<-"Session_Info"

if(file.exists(file_exists)){
  setwd("../")
  if(!file.exists("Code")){
    print("Please ensure that the script is run from the main repository folder")
  }
}else{
  if(!file.exists("Code")){
    print("Please ensure that the script is run from the main repository folder")
  }
}

# Functions, general information and packages------------------------------------

if(!require("plyr")){
  install.packages("plyr")
}
if(!require("ggplot2")){
  install.packages("ggplot2")
}
if(!require("dplyr")){
  install.packages("dplyr")
}
if(!require("car")){
  install.packages("car")
}
if(!require("fitdistrplus")){
  install.packages("fitdistrplus")
}
if(!require("tidyr")){
  install.packages("tidyr")
}
if(!require("cxr")){
  install.packages("cxr")
}
if(!require("MASS")){
  install.packages("MASS")
}
if(!require("mvtnorm")){
  install.packages("mvtnorm")
}
if(!require("lme4")){
  install.packages("lme4")
}
if(!require("lmerTest")){
  install.packages("lmerTest")
}
if(!require("emmeans")){
  install.packages("emmeans")
}
if(!require("glmmTMB")){
  install.packages("glmmTMB")
}
if(!require("DescTools")){
  install.packages("DescTools")
}
if(!require("performance")){
  install.packages("performance")
}
if(!require("DHARMa")){
  install.packages("DHARMa")
}
if(!require("effects")){
  install.packages("effects")
}
if(!require("marginaleffects")){
  install.packages("marginaleffects")
}
if(!require("LSAfun")){
  install.packages("LSAfun")
}
if(!require("arm")){
  install.packages("arm")
}
if(!require("cowplot")){
  install.packages("cowplot")
}
if(!require("grid")){
  install.packages("grid")
}
if(!require("gridExtra")){
  install.packages("gridExtra")
}

library(plyr)
library(dplyr)
library(ggplot2)
library(car)
library(fitdistrplus)
library(tidyr)
library(cxr)
library(MASS)
library(mvtnorm)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(DescTools)
library(performance)
library(DHARMa)
library(effects)
library(marginaleffects)
library(LSAfun)
library(arm)
library(cowplot)
library(grid)
library(gridExtra)

# Creating vectors with regime names to use for plots
regimeTu<-c("Tu \ncontrol", "Tu evolved \n in cadmium")
names(regimeTu)<-c("SR1", "SR2")

regimeTe<-c("Te \n control", "Te evolved \n in cadmium")
names(regimeTe)<-c("SR4", "SR5")

#' 
#' # Evaluation
#' Some pieces of the code take a lot of time to run. We provide the intermediate results to speed up the process. To run the whole code change eval from FALSE to TRUE.
#' 
# Evaluation------------------------------------
evaluation<-TRUE
if(!file.exists("./Analyses")){
  dir.create("./Analyses/")
}

#' 
#' # Running cxr
#' 
#' This code takes a bit of time to run. So the intermediate files are provided in the Analyses folder.
#' 
#' To run the code to estimate the parameters with cxr use the scripts Running_cxr_PooledData.R (to obtain the estimates for pooled data) and Running_cxr_Replicates.R (to obtain the estimates for each replicate).
#' 
# Running cxr------------------------------------
if(evaluation){
  print("Running cxr")
  source("./Code/Running_cxr_pooledData.R")
  source("./Code/Running_cxr_Replicates.R")
}

print("Importing parameter estimation")
#' 
#' 
#' # Importing parameter estimates
#' 
#' 
#' ### Estimated parameters per replicate
# Importing parameter estimates------------------------------------
# Importing the mean parameters
param_all_w0<-read.csv("./Analyses/cxr_normal_Replicates/parameters_cxr_normal.csv")
# Importing the upper values of the parameters
param_all_w0_upper<-read.csv("./Analyses/cxr_normal_Replicates/parameters_cxr_normal_upper.csv")
# Importing the lower values of the parameters
param_all_w0_lower<-read.csv( "./Analyses/cxr_normal_Replicates/parameters_cxr_normal_lower.csv")

# Removing the first column
param_all_w0<-param_all_w0[,-1]
param_all_w0_upper<-param_all_w0_upper[,-1]
param_all_w0_lower<-param_all_w0_lower[,-1]

#' 
#' ### Estimated parameters for pooled data

# Importing the mean parameters
param_all_REP<-read.csv("./Analyses/cxr_normal_Pooled/parameters_cxr_normal_Pooled.csv")
# Importing the upper values of the parameters
  param_all_REP_upper<-read.csv("./Analyses/cxr_normal_Pooled/parameters_cxr_normal_Pooled_upper.csv")
# Importing the lower values of the parameters
  param_all_REP_lower<-read.csv( "./Analyses/cxr_normal_Pooled/parameters_cxr_normal_Pooled_lower.csv")

param_all_REP<-param_all_REP[,-1]
param_all_REP_upper<-param_all_REP_upper[,-1]
param_all_REP_lower<-param_all_REP_lower[,-1]

print("Statistical tests")
#' 
#' # Testing differences in estimated parameters
#' 
#' First thing is to test the distribution of each trait
# Statistical analyses------------------------------------
#Using descdist to check potential distributions to test in the model

descdist(param_all_w0$Tu_lambda, discrete=FALSE, boot=1000)
descdist(param_all_w0$Te_lambda, discrete=FALSE, boot=1000)

descdist(param_all_w0$Tu_intra, discrete=FALSE, boot=1000)
descdist(param_all_w0$Te_intra, discrete=FALSE, boot=1000)

descdist(param_all_w0$Tu_inter, discrete=FALSE, boot=1000)
descdist(param_all_w0$Te_inter, discrete=FALSE, boot=1000)



#' This next section corresponds to a series of statistical tests to evaluate changes in growth rate, intraspecific competition and interspecific competition. For each trait we perform a test for Tu and Te species.
#' In the end of each section we provide a summary of the results. These results are cited directly in text or in supplementary tables.
## Does cadmium change parameters? (no evolution)
#' 
### Growth rate-------
#### Tu-------
##
# Testing the effect of cadmium on growth rate using different distributions  
gr_tu_cd_1<-glmmTMB(Tu_lambda~Environment, data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ))
gr_tu_cd_2<-glmmTMB(Tu_lambda~Environment, data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=Gamma(link="log"))
gr_tu_cd_3<-glmmTMB(Tu_lambda~Environment, data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=gaussian(link="log"))

AIC(gr_tu_cd_1,gr_tu_cd_2,gr_tu_cd_3)
# best model gr_tu_cd_2

summary(gr_tu_cd_2)

# Testing residual of the best model
simulationOutput <- simulateResiduals(fittedModel = gr_tu_cd_2, plot = F)
plot(simulationOutput)
#No problems

#' 
#' Cadmium significantly decreases growth rate for Tu.
#' 
#### Te----------
## ------------------------------------
# Testing the effect of cadmium on growth rate using different distributions  
gr_te_cd_1<-glmmTMB(Te_lambda~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ))
gr_te_cd_2<-glmmTMB(Te_lambda~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=Gamma(link="log"))
gr_te_cd_3<-glmmTMB(Te_lambda~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=gaussian(link="log"))

AIC(gr_te_cd_1,gr_te_cd_2,gr_te_cd_3)
#Best model is gr_te_cd_2

summary(gr_te_cd_2)

# Testing model residuals
simulationOutput <- simulateResiduals(fittedModel = gr_te_cd_2, plot = F)
plot(simulationOutput)# No problems

#' Cadmium significantly decreases growth rate for Te. Best model is with gamma distribution.
#' 
### Intraspecific competition----------
#' 
#### Tu----------

# Testing the effect of cadmium on intraspecific competition using different distributions  
intra_tu_cd_1<-glmmTMB(Tu_intra~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ))
intra_tu_cd_2<-glmmTMB(Tu_intra+1~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=Gamma(link="log"))
intra_tu_cd_3<-glmmTMB(Tu_intra+1~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=gaussian(link="log"))

AIC(intra_tu_cd_1,intra_tu_cd_2,intra_tu_cd_3)
#No difference between the models. In this case we should use the simplest model, which will be the gaussian distribution

summary(intra_tu_cd_1) 

# Testing model residuals
simulationOutput <- simulateResiduals(fittedModel = intra_tu_cd_1, plot = F)
plot(simulationOutput)# No problems

#' 
#' Presence of cadmium decreases intraspecific competition.
#' 
#### Te----------
#' 

# Testing the effect of cadmium on intraspecific competition using different distributions
intra_te_cd_1<-glmmTMB(Te_intra~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ))
intra_te_cd_2<-glmmTMB(Te_intra+1~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=Gamma(link="log"))
intra_te_cd_3<-glmmTMB(Te_intra+1~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=gaussian(link="log"))

AIC(intra_te_cd_1,intra_te_cd_2,intra_te_cd_3)
#No difference between the models. In this case we should use the simplest model, which will be the gaussian distribution

summary(intra_te_cd_1) #Again very similar estimates

# Testing the model residuals
simulationOutput <- simulateResiduals(fittedModel = intra_te_cd_1, plot = F)
plot(simulationOutput) #No problems

#' 
#' Presence of cadmium decreases intraspecific competition.
#' 
### Interspecific competition----------
#' 
#### Tu----------
#' 
# Testing the effect of cadmium on interspecific competition using different distributions. Here we only use the data from the regimes that evolved in the no-cadmium environment
inter_tu_cd_1<-glmmTMB(Tu_inter~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")))
inter_tu_cd_2<-glmmTMB(Tu_inter+1~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")), family=Gamma(link="log"))
inter_tu_cd_3<-glmmTMB(Tu_inter+1~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")), family=gaussian(link="log"))

AIC(inter_tu_cd_1,inter_tu_cd_2,inter_tu_cd_3)
# Same as above, similar AIC so we use the simplest model

summary(inter_tu_cd_1)

# Testing model residuals
simulationOutput <- simulateResiduals(fittedModel = inter_tu_cd_1, plot = F)
plot(simulationOutput) #no problems

#'  
#' Presence of cadmium decreases intraspecific competition for Tu
#' 
#### Te----------
#' 
# Testing the effect of cadmium on interspecific competition using different distributions. Here we only use the data from the regimes that evolved in the no-cadmium environment
inter_te_cd_1<-glmmTMB(Te_inter~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")))
inter_te_cd_2<-glmmTMB(Te_inter+1~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")), family=Gamma(link="log"))
inter_te_cd_3<-glmmTMB(Te_inter+1~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")), family=gaussian(link="log"))

#testing which is the best model
AIC(inter_te_cd_1,inter_te_cd_2,inter_te_cd_3)
# Same as above, similar AIC so we use the simplest model

summary(inter_te_cd_1)

# Testing residuals
simulationOutput <- simulateResiduals(fittedModel = inter_te_cd_1, plot = F)
plot(simulationOutput) #No problems

#' 
#' Presence of cadmium decreases interspecific competition for Te
#' 
#' 
## Summary----------
#' Results are cited in Table S5. The estimate, z-value and p-value correspond to the row indicated by EnvironmentN, which correspond to the differences between the no cadmium and cadmium environment.
#' 

summary(gr_tu_cd_2)

summary(gr_te_cd_2)

summary(intra_tu_cd_1) 

summary(intra_te_cd_1)

summary(inter_tu_cd_1)

summary(inter_te_cd_1)

Anova(gr_tu_cd_2)

Anova(gr_te_cd_2)

Anova(intra_tu_cd_1) 

Anova(intra_te_cd_1)

Anova(inter_tu_cd_1)

Anova(inter_te_cd_1)

#' 
## Does evolution change the performance in cadmium?----------
#' 
### Growth rate----------
#' 
#### Tu----------

# Testing the effect of evolution on cadmium on growth rate on cadmium watered plants using different distributions  
gr_tu_ev_1<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="Cd" & Te_Regime=="SR4"))
gr_tu_ev_2<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"), family=Gamma(link="log"))
gr_tu_ev_3<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"), family=gaussian(link="log"))

AIC(gr_tu_ev_1,gr_tu_ev_2,gr_tu_ev_3)
# There are no differences between AIC. In this case we will use the gamma, because that was the best distribution for the other tests. The results are very similar between normal and gamma distributions

summary(gr_tu_ev_2)

# Testing the residuals
simulationOutput <- simulateResiduals(fittedModel = gr_tu_ev_2, plot = F)
plot(simulationOutput) #no problems

#' No differences between models.
#' Evolution in cadmium increases the growth rate in cadmium for Tu
#' 
#### Te----------
#' 

# Testing the effect of evolution on cadmium on growth rate on cadmium watered plants using different distributions  
gr_te_ev_1<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"))
gr_te_ev_2<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"), family=Gamma(link="log"))
gr_te_ev_3<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"), family=gaussian(link="log"))

AIC(gr_te_ev_1,gr_te_ev_2,gr_te_ev_3)
# There are no differences between AIC. In this case we will use the gamma, because that was the best distribution for the other tests. The results are very similar between normal and gamma distributions

summary(gr_te_ev_2)

# Testing residuals
simulationOutput <- simulateResiduals(fittedModel = gr_te_ev_2, plot = F)
plot(simulationOutput)# No problems

#' No difference between models
#' Evolution in cadmium increases growth rate in cadmium, but only marginally
#' 
#' 
### Intraspecif competition----------
#' 
#### Tu----------
#' 

# Testing the effect of evolution on cadmium on intraspecific competition on cadmium watered plants using different distributions
intra_tu_ev_1<-glmmTMB(Tu_intra~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"))
intra_tu_ev_2<-glmmTMB(Tu_intra+1~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"), family=Gamma(link="log"))
intra_tu_ev_3<-glmmTMB(Tu_intra+1~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"), family=gaussian(link="log"))

AIC(intra_tu_ev_1,intra_tu_ev_2,intra_tu_ev_3)
# No differences between AICs, so the best model is the simplest (normal)

summary(intra_tu_ev_1)

# Testing the residuals
simulationOutput <- simulateResiduals(fittedModel = intra_tu_ev_1, plot = F)
plot(simulationOutput) #no problems

#' Similar models.
#' Evolution in cadmium does not affect intraspecific competition in Tu
#' 
#### Te----------

# Testing the effect of evolution on cadmium on intraspecific competition on cadmium watered plants using different distributions
intra_te_ev_1<-glmmTMB(Te_intra~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"))
intra_te_ev_2<-glmmTMB(Te_intra+1~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"), family=Gamma(link="log"))
intra_te_ev_3<-glmmTMB(Te_intra+1~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"), family=gaussian(link="log"))

AIC(intra_te_ev_1,intra_te_ev_2,intra_te_ev_3)
# No differences between AICs, so the best model is the simplest (normal)  

summary(intra_te_ev_1)

# Testing residuals
simulationOutput <- simulateResiduals(fittedModel = intra_te_ev_1, plot = F)
plot(simulationOutput) #no problems

#' 
#' Evolution in cadmium increases intraspecific competition for Te
#' 
#' 
### Interspecific competition----------
#' 
#### Tu----------
#' 

# Testing the effect of evolution on cadmium on interspecific competition on cadmium watered plants using different distributions
inter_tu_ev_1<-glmmTMB(Tu_inter~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="Cd"))
inter_tu_ev_2<-glmmTMB(Tu_inter+1~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="Cd"), family=Gamma(link="log"))
inter_tu_ev_3<-glmmTMB(Tu_inter+1~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="Cd"), family=gaussian(link="log"))

AIC(inter_tu_ev_1,inter_tu_ev_2,inter_tu_ev_3)
# No differences between AICs, so the best model is the simplest (normal)  

summary(inter_tu_ev_1)

# Testing residuals
simulationOutput <- simulateResiduals(fittedModel = inter_tu_ev_1, plot = F)
plot(simulationOutput) #no problem

#' There is no difference between AIC for the models
#' 
#### Te----------
#' 

# Testing the effect of evolution on cadmium on interspecific competition on cadmium watered plants using different distributions
inter_te_ev_1<-glmmTMB(Te_inter~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="Cd"))
inter_te_ev_2<-glmmTMB(Te_inter+1~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="Cd"), family=Gamma(link="log"))
inter_te_ev_3<-glmmTMB(Te_inter+1~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="Cd"), family=gaussian(link="log"))

AIC(inter_te_ev_1,inter_te_ev_2,inter_te_ev_3)
# No differences between models so we will use the simpler one

summary(inter_te_ev_1)

# Testing residuals
simulationOutput <- simulateResiduals(fittedModel = inter_te_ev_1, plot = F, )

plotResiduals(simulationOutput, subset(param_all_w0, Environment=="Cd")$Te_Regime)
plotResiduals(simulationOutput, subset(param_all_w0, Environment=="Cd")$Tu_Regime)
plot(simulationOutput) #no problem

# emmeans to check which comparisons are significant
emmeans(inter_te_ev_1, pairwise~Te_Regime+Tu_Regime, type="response")

#' No difference between AIC
#' 
#' Evolution in cadmium affects interspecific competition of Te and this effect depends on the evolution of Tu. Namely against Tu that did not evolve in cadmium the populations of Te that evolved in cadmium suffer more from interspecific competition. However, this does not happen for Tu that evolved in cadmium, where competition slightly decreases when Te evolved in cadmium.
#' 
#' 
## Summary----------
#' Results of the ANOVA for growth rate and intraspecific competition are cited in text. The ANOVA results for interspecific competition are cited in Table S3A. The  results from emmeans for Te are cited in Table S3B.

summary(gr_tu_ev_2)

summary(gr_te_ev_2)

## intra

summary(intra_tu_ev_1)

summary(intra_te_ev_1)

## inter
summary(inter_tu_ev_1)

summary(inter_te_ev_1)

# Anova

Anova(gr_tu_ev_2)

Anova(gr_te_ev_2)

## intra

Anova(intra_tu_ev_1)

Anova(intra_te_ev_1)

## Inter
Anova(inter_tu_ev_1, type=3)

Anova(inter_te_ev_1, type=3)

# emmeans to check which comparisons are significant
emmeans(inter_te_ev_1, pairwise~Te_Regime+Tu_Regime, type="response")

#' 
## Does evolution change the performance in the ancestral environment?----------
#' 
### Growth rate----------
#' 
#### Tu----------
#' 
# Testing the effect of evolution on cadmium on growth rate on plants without cadmium using different distributions 
gr_tu_an_1<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="N" & Te_Regime=="SR4"))
gr_tu_an_2<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"), family=Gamma(link="log"))
gr_tu_an_3<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"), family=gaussian(link="log"))

AIC(gr_tu_an_1,gr_tu_an_2,gr_tu_an_3)

summary(gr_tu_an_2)

# Test residuals
simulationOutput <- simulateResiduals(fittedModel = gr_tu_an_2, plot = F, )
plot(simulationOutput)# no problem

#' Gamma shows slightly lower AIC
#' No effect, although evolution in cadmium slightly reduces growth rate.
#' 
#' #### Te----------
#' 
# Testing the effect of evolution on cadmium on growth rate on plants without cadmium using different distributions 
gr_te_an_1<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"))
gr_te_an_2<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"), family=Gamma(link="log"))
gr_te_an_3<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"), family=gaussian(link="log"))

AIC(gr_te_an_1,gr_te_an_2,gr_te_an_3)
# Best model is with the gamma distribution

summary(gr_te_an_2)

simulationOutput <- simulateResiduals(fittedModel = gr_te_an_2, plot = F, )
plot(simulationOutput)

#' 
#' No effect of evolution in cadmium,  although evolution in cadmium slightly reduces growth rate.
#' 
### Intra----------
#' 
#### Tu----------
#' 

# Testing the effect of evolution on cadmium on intraspecific competition on plants without cadmium using different distributions 
intra_tu_an_1<-glmmTMB(Tu_intra~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"))
intra_tu_an_2<-glmmTMB(Tu_intra+1~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"), family=Gamma(link="log"))
intra_tu_an_3<-glmmTMB(Tu_intra+1~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"), family=gaussian(link="log"))

AIC(intra_tu_an_1,intra_tu_an_2,intra_tu_an_3)
# No differences in AIC, best model is the simplest

summary(intra_tu_an_1)

# Testing the residuals
simulationOutput <- simulateResiduals(fittedModel = intra_tu_an_1, plot = F, )
plot(simulationOutput) #no problem

#' 
#' No effect of evolution in cadmium, although it slighlty decreases the intraspecific competition
#' 
#### Te----------
#' 

# Testing the effect of evolution on cadmium on intraspecific competition on plants without cadmium using different distributions 
intra_te_an_1<-glmmTMB(Te_intra~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"))
intra_te_an_2<-glmmTMB(Te_intra+1~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"), family=Gamma(link="log"))
intra_te_an_3<-glmmTMB(Te_intra+1~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"), family=gaussian(link="log"))

AIC(intra_te_an_1,intra_te_an_2,intra_te_an_3)
# No differences in AIC, best model is the simplest

summary(intra_te_an_1)

# Test residuals
simulationOutput <- simulateResiduals(fittedModel = intra_te_an_1, plot = F, )
plot(simulationOutput)# no problem

#' 
#' Evolution in cadmium does not affect Te intraspecific competition in the ancestral environment
#' 
### Interspecific competition----------
#' 
#### Tu----------

# Testing the effect of evolution on cadmium on interspecific competition on plants without cadmium using different distributions 
inter_tu_an_1<-glmmTMB(Tu_inter~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="N"))
inter_tu_an_2<-glmmTMB(Tu_inter+1~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="N"), family=Gamma(link="log"))
inter_tu_an_3<-glmmTMB(Tu_inter+1~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="N"), family=gaussian(link="log"))

AIC(inter_tu_an_1,inter_tu_an_2,inter_tu_an_3)
# No differences in AIC

summary(inter_tu_an_1)

# Test residuals
simulationOutput <- simulateResiduals(fittedModel = inter_tu_an_1, plot = F, )
plot(simulationOutput) # no problem

#' 
#' Evolution in cadmium does not affect sensitivity to interspecific competition for Tu in the ancestral environment
#' 
#### Te----------
#' 

# Testing the effect of evolution on cadmium on interspecific competition on plants without cadmium using different distributions 
inter_te_an_1<-glmmTMB(Te_inter~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="N"))
inter_te_an_2<-glmmTMB(Te_inter+1~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="N"), family=Gamma(link="log"))
inter_te_an_3<-glmmTMB(Te_inter+1~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="N"), family=gaussian(link="log"))

AIC(inter_te_an_1,inter_te_an_2,inter_te_an_3)
# No differences in AIC

summary(inter_te_an_1)

# Test residuals
simulationOutput <- simulateResiduals(fittedModel = inter_te_an_1, plot = F, )
plot(simulationOutput) # no problems


#' No difference in model AIC
#' No changes in interspecific competition for Te in the ancestral environment, after evolution in cadmium
#' 
### Summary----------
#' Results of the ANOVA for growth rate and intraspecific competition are cited in text. The ANOVA results for interspecific competition are cited in Table S4.

summary(gr_tu_an_2)

summary(gr_te_an_2)

summary(intra_tu_an_1)

summary(intra_te_an_1)

summary(inter_tu_an_1)

summary(inter_te_an_1)


Anova(gr_tu_an_2)
Anova(gr_te_an_2)

Anova(intra_tu_an_1)
Anova(intra_te_an_1)

Anova(inter_tu_an_1, type=3)
Anova(inter_te_an_1, type=3)

## Bootstrap analyses ----------
#' 
#' ## Bootstrap differences between selection regimes
#' 
#' Beware this script takes a long time to run.
#' For each question we will randomize the replicates between selection regimes, using 10000 bootstrap samples and test how many times do we recover significant differences between the different regimes just by chance.
#' 
set.seed(1809) # Setting the seed to ensure reproducibility
nboot<-10000 # Number of bootstrap samples

if(!evaluation){
  load(file="./Analyses/Bootstrap_stats_tests.RData")
}

#' 
#' 
#' ### Does cadmium change parameters?
#' #### Tu
#' 
#+ eval=FALSE
if(evaluation){
#Bootstrap to reestimate the p-value obtained for growth rate and intraspecific competition.
boot_tu_gr_intra_env<-as.data.frame(t(sapply(c(1:nboot),function(x){
  # Printing every 10 samples
  if(x%%10 ==0){
    print(x)
  }
  
  # Subset data to use only the Tu and Te regimes
  auxi<-subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" )
  # randomize the rows
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  # replace the lambda, intra and inter with the randomized data
  auxi$Tu_lambda<-auxi[rand_numb,"Tu_lambda"] # randomizing the trais
  auxi$Tu_intra<-auxi[rand_numb,"Tu_intra"]
  auxi$Tu_inter<-auxi[rand_numb,"Tu_inter"]
  
  #Testing the impact of environment with the new randomized data for the three traits
  gr <-glmmTMB(Tu_lambda~Environment, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Tu_intra~Environment, data=auxi)
  inter<-glmmTMB(Tu_inter~Environment, data=auxi)
  
  # Estimate mean per environment for each trait
  sum_auxi<-auxi %>% group_by(Environment)%>% summarize(meanGr=mean(Tu_lambda, na.rm=TRUE), meanIntra=mean(Tu_intra, na.rm=TRUE), meaninter=mean(Tu_inter, na.rm=TRUE)) %>% as.data.frame()
  
  # Estimate difference between no cadmium and cadmium
  # N - Cd
  diff<-sum_auxi[2,c(2:4)]-sum_auxi[1,(2:4)]
  
  # Store the p-value
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  
  # Vector with p-value and differences
  c(gr_p, intra_p, inter_p, diff[1,1], diff[1,2], diff[1,3])
  
} )))}

#' 
#'
#' 
#+ eval=FALSE
if(evaluation){
#Bootstrap to reestimate the p-value obtained for growth rate and intraspecific competition.
boot_te_gr_intra_env<-as.data.frame(t(sapply(c(1:nboot),function(x){
  
  # Printing every 10 samples
  if(x%%10 ==0){
    print(x)
  }
    # Subset data to use only the Tu and Te regimes that evolved in no cadmium
  auxi<-subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" )
    # randomize the rows
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  
  # replace the lambda, intra and inter with the randomized data
  auxi$Te_lambda<-auxi[rand_numb,"Te_lambda"] # randomizing the trais
  auxi$Te_intra<-auxi[rand_numb,"Te_intra"]
  auxi$Te_inter<-auxi[rand_numb,"Te_inter"]
  
    #Testing the impact of environment with the new randomized data for the three traits
  gr <-glmmTMB(Te_lambda~Environment, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Te_intra~Environment, data=auxi)
  inter <-glmmTMB(Te_inter~Environment, data=auxi)
  
    # Estimate mean per environment for each trait
  sum_auxi<-auxi %>% group_by(Environment)%>% summarize(meanGr=mean(Te_lambda, na.rm=TRUE), meanIntra=mean(Te_intra, na.rm=TRUE), meaninter=mean(Te_inter, na.rm=TRUE)) %>% as.data.frame()
  
  # Estimate difference between no cadmium and cadmium
  # N - Cd
  diff<-sum_auxi[2,c(2:4)]-sum_auxi[1,(2:4)]
  
  # Store the p-value
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  
  # Vector with p-value and differences
  c(gr_p, intra_p, inter_p,diff[1,1], diff[1,2], diff[1,3])
  
} )))

# Adding column names to the bootstrap
colnames(boot_tu_gr_intra_env)<-c("lambda_p","intra_p","inter_p", "lambda_diff", "intra_diff", "inter_diff")
colnames(boot_te_gr_intra_env)<-c("lambda_p","intra_p","inter_p", "lambda_diff", "intra_diff", "inter_diff")
str(boot_te_gr_intra_env)
}

#' 
#' #### Printing the p-value
#' 

#Number of times that p-value was lower than 0.05
length(which(boot_tu_gr_intra_env$lambda_p<=0.05))/(nboot)

length(which(boot_tu_gr_intra_env$intra_p<=0.05))/(nboot)

length(which(boot_tu_gr_intra_env$inter_p<=0.05))/(nboot)

length(which(boot_te_gr_intra_env$lambda_p<=0.05))/(nboot)

length(which(boot_te_gr_intra_env$intra_p<=0.05))/(nboot)

length(which(boot_te_gr_intra_env$inter_p<=0.05))/(nboot)

#' 

print("Boot p-values for tests between environments")
length(which(boot_tu_gr_intra_env$lambda_p<=as.data.frame(Anova(gr_tu_cd_2))[,3]))/(nboot+1)

length(which(boot_tu_gr_intra_env$intra_p<=as.data.frame(Anova(intra_tu_cd_1))[,3]))/(nboot+1)

length(which(boot_tu_gr_intra_env$inter_p<=as.data.frame(Anova(inter_tu_cd_1))[,3]))/(nboot+1)

length(which(boot_te_gr_intra_env$lambda_p<=as.data.frame(Anova(gr_te_cd_2))[,3]))/(nboot+1)

length(which(boot_te_gr_intra_env$intra_p<=as.data.frame(Anova(intra_te_cd_1))[,3]))/(nboot+1)

length(which(boot_te_gr_intra_env$inter_p<=as.data.frame(Anova(inter_te_cd_1))[,3]))/(nboot+1)

#' 
#' 
#' ### Does evolution change the performance in cadmium?
#' 
#' 
#+ eval=FALSE
if(evaluation){
#Bootstrap to reestimate the p-value obtained for growth rate and intraspecific competition to test if evoution affected the performance in cadmium
boot_tu_evolcd<-as.data.frame(t(sapply(c(1:nboot),function(x){
    # Printing every 10 samples
  if(x%%10 ==0){
    print(x)
  }
    
  # Subset data to get both Tu regimes (that is why we select the SR4 which is a Te regime, this way we guarantee that we have the two Tu regimes).
    
  auxi<-subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4")
  # randomize the rows (to randomize the data)
    rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
    # replace the lambda and intra with  the randomized data
  auxi$Tu_lambda<-auxi[rand_numb,"Tu_lambda"] # randomizing the trais
  auxi$Tu_intra<-auxi[rand_numb,"Tu_intra"]
  
  # For interspecific competition we need both regimes for Te and Tu, ence we need to subset only for the environment
  # Here we do the same thing, reshuffle the row numbers and then replace the data
  auxi2<-subset(param_all_w0, Environment=="Cd")
  rand_numb2<-sample(c(1:dim(auxi2)[1]), dim(auxi2)[1], replace = TRUE)
  auxi2$Tu_inter<-auxi2[rand_numb2,"Tu_inter"]
  
   #Testing the impact of evolution in cadmium in the performance  in cadium with the new randomized data for the three traits
  gr <-glmmTMB(Tu_lambda~Tu_Regime, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Tu_intra~Tu_Regime, data=auxi)
  inter<-glmmTMB(Tu_inter~Tu_Regime*Te_Regime, data=auxi2)
  
   # Estimate mean per selection regime for each trait
  sum_auxi<-auxi %>% group_by(Tu_Regime)%>% summarize(meanGr=mean(Tu_lambda, na.rm=TRUE), meanIntra=mean(Tu_intra, na.rm=TRUE)) %>% as.data.frame()
  
  sum_auxi2<-auxi2 %>% group_by(Tu_Regime, Te_Regime)%>% summarize( meanInter=mean(Tu_inter, na.rm=TRUE)) %>% as.data.frame()
  
    # Estimate difference between evolved and non-evolved
  diff<-sum_auxi[2,c(2:3)]-sum_auxi[1,(2:3)]
  # For inter it is more complicated because we have several pairwises
  diff2<-sum_auxi2[1,3]-sum_auxi2[2,3] # Evolution of the competitor with control focal
  diff3<-sum_auxi2[1,3]-sum_auxi2[3,3] # Evolution of the focal with control competitor
  diff4<-sum_auxi2[2,3]-sum_auxi2[4,3] # Evolution of the focal with evolved competitor
  diff5<-sum_auxi2[3,3]-sum_auxi2[4,3] # Evolution of the competitor with evolved focal
  
    # Store the p-value
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  inter_p2<-as.data.frame(Anova(inter))[2,3]
  inter_p3<-as.data.frame(Anova(inter))[3,3]
  
  # Vector with p-value and differences
  c(gr_p, intra_p, inter_p,inter_p2, inter_p3, diff[1,1], diff[1,2], diff2, diff3, diff4,diff5)
  
} )))
}

#' 
#'
#' 
#+ eval=FALSE
if(evaluation){
#Bootstrap to reestimate the p-value obtained for growth rate and intraspecific competition to test if evoution affected the performance in cadmium. For T. evansi
boot_te_evolcd<-as.data.frame(t(sapply(c(1:nboot),function(x){
  # Printing every 10 samples
  if(x%%10 ==0){
    print(x)
  }
  # Subset data to get both Te regimes (that is why we select the SR1 which is a Tu regime, this way we guarantee that we have the two Te regimes).
  auxi<-subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1")
  # randomize the rows
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  # replace the lambda and intra with  the randomized data
  auxi$Te_lambda<-auxi[rand_numb,"Te_lambda"] # randomizing the trais
  auxi$Te_intra<-auxi[rand_numb,"Te_intra"]
  
  # For interspecific competition we need both regimes for Te and Tu, hence we need to subset only for the environment
  # Here we do the same thing, reshuffle the row numbers and then replace the data
  auxi2<-subset(param_all_w0, Environment=="Cd")
  # reshuffle rows
  rand_numb2<-sample(c(1:dim(auxi2)[1]), dim(auxi2)[1], replace = TRUE)
  auxi2$Te_inter<-auxi2[rand_numb2,"Te_inter"]
  
  #Testing the impact of evolution in cadmium in the performance  in cadium with the new randomized data for the three traits
  gr <-glmmTMB(Te_lambda~Te_Regime, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Te_intra~Te_Regime, data=auxi)
  inter<-glmmTMB(Te_inter~Tu_Regime*Te_Regime, data=auxi2)
 
    # Estimate mean per selection regime for each trait
  sum_auxi<-auxi %>% group_by(Te_Regime)%>% summarize(meanGr=mean(Te_lambda, na.rm=TRUE), meanIntra=mean(Te_intra, na.rm=TRUE)) %>% as.data.frame()
  
  sum_auxi2<-auxi2 %>% group_by(Tu_Regime, Te_Regime)%>% summarize( meanInter=mean(Te_inter, na.rm=TRUE)) %>% as.data.frame()
  
   # Estimate the difference between evolved in cadmium and evolved in no cadmium
  diff<-sum_auxi[2,c(2:3)]-sum_auxi[1,(2:3)]
  diff2<-sum_auxi2[1,3]-sum_auxi2[2,3] # Evolution of the competitor with control focal
  diff3<-sum_auxi2[1,3]-sum_auxi2[3,3] # Evolution of the focal with control competitor
  diff4<-sum_auxi2[2,3]-sum_auxi2[4,3] # Evolution of the focal with evolved competitor
  diff5<-sum_auxi2[3,3]-sum_auxi2[4,3] # Evolution of the competitor with evolved focal
  
    # Store the p-value
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  inter_p2<-as.data.frame(Anova(inter))[2,3]
  inter_p3<-as.data.frame(Anova(inter))[3,3]
    # Vector with p-value and differences
  c(gr_p, intra_p, inter_p,inter_p2, inter_p3, diff[1,1], diff[1,2], diff2, diff3, diff4,diff5)
  
} )))

# Add column names
colnames(boot_tu_evolcd)<-c("lambda_p","intra_p","inter_p_TuReg","inter_p_TeReg","inter_p_int", "lambda_diff", "intra_diff", "inter_diffEvolComp_focalControl","inter_diffEvolFocal_CompControl","inter_diffEvolFocal_EvolComp","inter_diffEvolComp_focalEvol" )
colnames(boot_te_evolcd)<-c("lambda_p","intra_p","inter_p_TuReg","inter_p_TeReg","inter_p_int", "lambda_diff", "intra_diff", "inter_diffEvolComp_focalControl","inter_diffEvolFocal_CompControl","inter_diffEvolFocal_EvolComp","inter_diffEvolComp_focalEvol" )
}

#' 
#' #### Printing the p-value
#' 

#Number of times that p-value was lower than 0.05
length(which(boot_tu_evolcd$lambda_p<=0.05))/(nboot)
length(which(boot_tu_evolcd$intra_p<=0.05))/(nboot)
length(which(boot_tu_evolcd$inter_p_TuReg<=0.05))/(nboot)
length(which(boot_tu_evolcd$inter_p_TeReg<=0.05))/(nboot)
length(which(boot_tu_evolcd$inter_p_int<=0.05))/(nboot)
length(which(boot_te_evolcd$lambda_p<=0.05))/(nboot)
length(which(boot_te_evolcd$intra_p<=0.05))/(nboot)
length(which(boot_te_evolcd$inter_p_TuReg<=0.05))/(nboot)
length(which(boot_te_evolcd$inter_p_TeReg<=0.05))/(nboot)
length(which(boot_te_evolcd$inter_p_int<=0.05))/(nboot)

#' 


print("Boot p-values for tests for evolution in cadmium")
length(which(boot_tu_evolcd$lambda_p<=as.data.frame(Anova(gr_tu_ev_2))[,3]))/(nboot+1)

length(which(boot_tu_evolcd$intra_p<=as.data.frame(Anova(intra_tu_ev_1))[,3]))/(nboot+1)

length(which(boot_tu_evolcd$inter_p_TuReg<=as.data.frame(Anova(inter_tu_ev_1))[1,3]))/(nboot+1)

length(which(boot_tu_evolcd$inter_p_TeReg<=as.data.frame(Anova(inter_tu_ev_1))[2,3]))/(nboot+1)

length(which(boot_tu_evolcd$inter_p_int<=as.data.frame(Anova(inter_tu_ev_1))[3,3]))/(nboot+1)

length(which(boot_te_evolcd$lambda_p<=as.data.frame(Anova(gr_te_ev_2))[,3]))/(nboot+1)

length(which(boot_te_evolcd$intra_p<=as.data.frame(Anova(intra_te_ev_1))[,3]))/(nboot+1)

length(which(boot_te_evolcd$inter_p_TuReg<=as.data.frame(Anova(inter_te_ev_1))[1,3]))/(nboot+1)

length(which(boot_te_evolcd$inter_p_TeReg<=as.data.frame(Anova(inter_te_ev_1))[2,3]))/(nboot+1)

length(which(boot_te_evolcd$inter_p_int<=as.data.frame(Anova(inter_te_ev_1))[3,3]))/(nboot+1)

#' 
#' 
#' ### Does evolution change the ancestral?
#' 
#'
#' 
#+ eval=FALSE
if(evaluation){
#Bootstrap to reestimate the p-value obtained for growth rate and intraspecific competition to test if evoution affected the performance in the ancestral environment
boot_tu_evolN<-as.data.frame(t(sapply(c(1:nboot),function(x){
  
    # Printing every 10 samples
  if(x%%10 ==0){
    print(x)
  }
  
  # Subset data to get both Tu regimes (that is why we select the SR4 which is a Te regime, this way we guarantee that we have the two Tu regimes).
  auxi<-subset(param_all_w0, Environment=="N"& Te_Regime=="SR4")
  
  # randomize the rows
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  # replace the lambda, intra and inter with  the randomized data
  auxi$Tu_lambda<-auxi[rand_numb,"Tu_lambda"] # randomizing the trais
  auxi$Tu_intra<-auxi[rand_numb,"Tu_intra"]
  
  # Subset data to get both Tu regimes (that is why we select the SR4 which is a Te regime, this way we guarantee that we have the two Tu regimes).
  auxi2<-subset(param_all_w0, Environment=="N")
  # randomized rows
  rand_numb2<-sample(c(1:dim(auxi2)[1]), dim(auxi2)[1], replace = TRUE)
  auxi2$Tu_inter<-auxi2[rand_numb2,"Tu_inter"]
  
  #Testing the impact of evolution in cadmium in the performance  in the ancestral environment with the new randomized data for the three traits
  gr <-glmmTMB(Tu_lambda~Tu_Regime, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Tu_intra~Tu_Regime, data=auxi)
  inter<-glmmTMB(Tu_inter~Tu_Regime*Te_Regime, data=auxi2)
  
    # Estimate mean per selection regime for each trait
  sum_auxi<-auxi %>% group_by(Tu_Regime)%>% summarize(meanGr=mean(Tu_lambda, na.rm=TRUE), meanIntra=mean(Tu_intra, na.rm=TRUE)) %>% as.data.frame()
  
  sum_auxi2<-auxi2 %>% group_by(Tu_Regime, Te_Regime)%>% summarize( meanInter=mean(Tu_inter, na.rm=TRUE)) %>% as.data.frame()
  
  # Estimate the difference between evolved in cadmium and evolved in no cadmium
  diff<-sum_auxi[2,c(2:3)]-sum_auxi[1,(2:3)]
  diff2<-sum_auxi2[1,3]-sum_auxi2[2,3] # Evolution of the competitor with control focal
  diff3<-sum_auxi2[1,3]-sum_auxi2[3,3] # Evolution of the focal with control competitor
  diff4<-sum_auxi2[2,3]-sum_auxi2[4,3] # Evolution of the focal with evolved competitor
  diff5<-sum_auxi2[3,3]-sum_auxi2[4,3] # Evolution of the competitor with evolved focal
  
    # Store the p-value
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  inter_p2<-as.data.frame(Anova(inter))[2,3]
  inter_p3<-as.data.frame(Anova(inter))[3,3]
  
  # Vector with p-value and differences
  c(gr_p, intra_p, inter_p,inter_p2, inter_p3, diff[1,1], diff[1,2], diff2, diff3, diff4,diff5)
  
} )))
}

#' 
#'
#' 
#+ eval=FALSE
if(evaluation){
#Bootstrap to reestimate the p-value obtained for growth rate and intraspecific competition to test if evoution affected the performance in the ancestral environment. for Te
boot_te_evolN<-as.data.frame(t(sapply(c(1:nboot),function(x){
  # Printing every 10 samples
  if(x%%10 ==0){
    print(x)
  }
  # Subset data to get both Te regimes (that is why we select the SR1 which is a Tu regime, this way we guarantee that we have the two Te regimes).
  auxi<-subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1")
  # randomize the rows
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
   # replace the lambda, intra and inter with  the randomized data
  auxi$Te_lambda<-auxi[rand_numb,"Te_lambda"] # randomizing the trais
  auxi$Te_intra<-auxi[rand_numb,"Te_intra"]
  
  # print(x)
  
  auxi2<-subset(param_all_w0, Environment=="N")
  rand_numb2<-sample(c(1:dim(auxi2)[1]), dim(auxi2)[1], replace = TRUE)
  auxi2$Te_inter<-auxi2[rand_numb2,"Te_inter"]
  
    #Testing the impact of evolution in cadmium on the performance in the ancestral environment with the new randomized data for the three traits
  gr <-glmmTMB(Te_lambda~Te_Regime, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Te_intra~Te_Regime, data=auxi)
  inter<-glmmTMB(Te_inter~Tu_Regime*Te_Regime, data=auxi2)
  
    # Estimate mean per selection regime for each trait
  sum_auxi<-auxi %>% group_by(Te_Regime)%>% summarize(meanGr=mean(Te_lambda, na.rm=TRUE), meanIntra=mean(Te_intra, na.rm=TRUE)) %>% as.data.frame()
  
  sum_auxi2<-auxi2 %>% group_by(Tu_Regime, Te_Regime)%>% summarize( meanInter=mean(Te_inter, na.rm=TRUE)) %>% as.data.frame()
  
    # Estimate the difference between evolved in cadmium and evolved in no cadmium
  diff<-sum_auxi[2,c(2:3)]-sum_auxi[1,(2:3)]
  diff2<-sum_auxi2[1,3]-sum_auxi2[2,3] # Evolution of the competitor with control focal
  diff3<-sum_auxi2[1,3]-sum_auxi2[3,3] # Evolution of the focal with control competitor
  diff4<-sum_auxi2[2,3]-sum_auxi2[4,3] # Evolution of the focal with evolved competitor
  diff5<-sum_auxi2[3,3]-sum_auxi2[4,3] # Evolution of the competitor with evolved focal
  
  # Store the p-value
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  inter_p2<-as.data.frame(Anova(inter))[2,3]
  inter_p3<-as.data.frame(Anova(inter))[3,3]
  
  # Vector with p-value and differences
  c(gr_p, intra_p, inter_p,inter_p2, inter_p3, diff[1,1], diff[1,2], diff2, diff3, diff4,diff5)
  
} )))

# adding column names
colnames(boot_tu_evolN)<-c("lambda_p","intra_p","inter_p_TuReg","inter_p_TeReg","inter_p_int", "lambda_diff", "intra_diff", "inter_diffEvolComp_focalControl","inter_diffEvolFocal_CompControl","inter_diffEvolFocal_EvolComp","inter_diffEvolComp_focalEvol" )

colnames(boot_te_evolN)<-c("lambda_p","intra_p","inter_p_TuReg","inter_p_TeReg","inter_p_int", "lambda_diff", "intra_diff", "inter_diffEvolComp_focalControl","inter_diffEvolFocal_CompControl","inter_diffEvolFocal_EvolComp","inter_diffEvolComp_focalEvol" )
}


#' 
#' ### Printing the p-value

#Number of times that p-value was lower than 0.05
length(which(boot_tu_evolN$lambda_p<=0.05))/(nboot)
length(which(boot_tu_evolN$intra_p<=0.05))/(nboot)
length(which(boot_tu_evolN$inter_p_TuReg<=0.05))/(nboot)
length(which(boot_tu_evolN$inter_p_TeReg<=0.05))/(nboot)
length(which(boot_tu_evolN$inter_p_int<=0.05))/(nboot)
length(which(boot_te_evolN$lambda_p<=0.05))/(nboot)
length(which(boot_te_evolN$intra_p<=0.05))/(nboot)
length(which(boot_te_evolN$inter_p_TuReg<=0.05))/(nboot)
length(which(boot_te_evolN$inter_p_TeReg<=0.05))/(nboot)
length(which(boot_te_evolN$inter_p_int<=0.05))/(nboot)


#' 

print("Boot p-values for tests for evolution in cadmium")
length(which(boot_tu_evolN$lambda_p<=as.data.frame(Anova(gr_tu_an_2))[,3]))/(nboot+1)

length(which(boot_tu_evolN$intra_p<=as.data.frame(Anova(intra_tu_an_1))[,3]))/(nboot+1)

length(which(boot_tu_evolN$inter_p_TuReg<=as.data.frame(Anova(inter_tu_an_1))[1,3]))/(nboot+1)

length(which(boot_tu_evolN$inter_p_TeReg<=as.data.frame(Anova(inter_tu_an_1))[2,3]))/(nboot+1)

length(which(boot_tu_evolN$inter_p_int<=as.data.frame(Anova(inter_tu_an_1))[3,3]))/(nboot+1)

length(which(boot_te_evolN$lambda_p<=as.data.frame(Anova(gr_te_an_2))[,3]))/(nboot+1)

length(which(boot_te_evolN$intra_p<=as.data.frame(Anova(intra_te_an_1))[,3]))/(nboot+1)

length(which(boot_te_evolN$inter_p_TuReg<=as.data.frame(Anova(inter_te_an_1))[1,3]))/(nboot+1)

length(which(boot_te_evolN$inter_p_TeReg<=as.data.frame(Anova(inter_te_an_1))[2,3]))/(nboot+1)

length(which(boot_te_evolN$inter_p_int<=as.data.frame(Anova(inter_te_an_1))[3,3]))/(nboot+1)

#' 
#' 
#' #### Save the results
# save the results so we can reuse it later
if(evaluation){
save.image(file="./Analyses/Bootstrap_stats_tests.RData")
  }
print("Structural stability")
#' 
# Predicting coexistence using structural stability----------
#' 
#' ## Defining functions

#input parameters:
#alpha = competition strength matrix 
#r = vector of intrinsic growth rates

# This functions were obtained from Saavedra et al. 2017, 2020 
#Function to estimate structural niche difference (output on a log scale)
Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha, tol = 1e-40)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}

#Function to estimate the vector defining the centroid of the feasibility domain
r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}

#Function to estimate structural fitness difference (in degree)
theta <- function(alpha,r){
  r_c <- r_centroid(alpha)
  out <- acos(sum(r_c*r, na.rm = TRUE)/(sqrt(sum(r^2, na.rm = TRUE))*sqrt(sum(r_c^2, na.rm = TRUE))))*90/pi
  return(out)
}

#Function to test if a system (alpha and r) is feasible (output 1 = feasible, 0 = not feasible)
test_feasibility <- function(alpha,r){
  out <- prod(solve(alpha,r)>0)
  return(out)
}

#' 
## Estimating structural stability for pooled replicates----------
#' 
#' Here we use the parameter estimates to calculate niche and fitness differences to predict coexistence.

# Creates a loop to predict coexistence for each combination of selection regimes and environments
struct_mat_REP<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_REP[,1])), function(x){
  #print(x)
  # Create the matrix with the intra and interspecific competition
  aux_alpha<-matrix(c(param_all_REP$Tu_intra[x],param_all_REP$Tu_inter[x],param_all_REP$Te_inter[x], param_all_REP$Te_intra[x]), ncol=2, byrow=TRUE)
  # Vector with growth rates
  aux_lambda<-c(param_all_REP$Tu_lambda[x],param_all_REP$Te_lambda[x])
  
  # Estimate niche differences
  om<-Omega(aux_alpha)
  # Estimate fitness fitness differences
  tta<-theta(aux_alpha, aux_lambda)
  # Predict coexistence (i.e. feasibility)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  # vector with niche and fitness differences and to predict coexistence
  c(om, tta, feas)
}))))

# Add column names
colnames(struct_mat_REP)<-c("ND", "FD", "Feasibility")

# Creates a loop to predict coexistence for each combination of selection regimes and environments but for the UPPER values
struct_mat_REP_U<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_REP_upper[,1])), function(x){
  
  # competitive matrix
  aux_alpha<-matrix(c(param_all_REP_upper$Te_intra[x],param_all_REP_upper$Te_inter[x], param_all_REP_upper$Tu_inter[x],param_all_REP_upper$Tu_intra[x]), ncol=2, byrow=TRUE)
  # growth rate vector
  aux_lambda<-c(param_all_REP_upper$Te_lambda[x],param_all_REP_upper$Tu_lambda[x] )
  
  # estimate niche, fitness differences and feasibility
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

# Creates a loop to predict coexistence for each combination of selection regimes and environments but for the LOWER values
struct_mat_REP_L<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_REP_lower[,1])), function(x){
  
  # Competitive matrix
  aux_alpha<-matrix(c(param_all_REP_lower$Te_intra[x],param_all_REP_lower$Te_inter[x],param_all_REP_lower$Tu_inter[x],param_all_REP_lower$Tu_intra[x]), ncol=2, byrow=TRUE)
  # vector of growth rates
  aux_lambda<-c(param_all_REP_lower$Te_lambda[x],param_all_REP_lower$Tu_lambda[x] )
  
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  # estimate niche, fitness differences and feasibility
  c(om, tta, feas)
}))))

# Add column names
colnames(struct_mat_REP_U)<-c("ND_U", "FD_U", "Feasibility_U")
colnames(struct_mat_REP_L)<-c("ND_L", "FD_L", "Feasibility_L")

# Join parameter and structural stability estimates (mean, upper and lower) 
struct_mat_REP<-cbind(param_all_REP, struct_mat_REP,struct_mat_REP_L,struct_mat_REP_U)


#' 
## Estimating structural stability per replicate----------
#' Here we apply the same loop as above, but doing the estimates per replicate.
#' 

# Creates a loop to predict coexistence for each combination of selection regimes and environments
struct_mat_w0<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_w0[,1])), function(x){
  
  # Competition matrix
  aux_alpha<-matrix(c(param_all_w0$Te_intra[x],param_all_w0$Te_inter[x],param_all_w0$Tu_inter[x], param_all_w0$Tu_intra[x]), ncol=2, byrow=TRUE)
  # vector of growth rates
  aux_lambda<-c(param_all_w0$Te_lambda[x],param_all_w0$Tu_lambda[x] )
  
  #Estimate niche and fitness differences and feasibility domain
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))
# add column names
colnames(struct_mat_w0)<-c("ND", "FD", "Feasibility")

# Creates a loop to predict coexistence for each combination of selection regimes and environments but for LOWER estimates
struct_mat_w0_L<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_w0_lower[,1])), function(x){
  
  # Competition matrix
  aux_alpha<-matrix(c(param_all_w0_lower$Te_intra[x],param_all_w0_lower$Te_inter[x],param_all_w0_lower$Tu_inter[x], param_all_w0_lower$Tu_intra[x]), ncol=2, byrow=TRUE)
  # Vector of growth rates
  aux_lambda<-c(param_all_w0_lower$Te_lambda[x],param_all_w0_lower$Tu_lambda[x] )
  
  #Estimate niche and fitness differences and feasibility domain
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

# Add column names
colnames(struct_mat_w0_L)<-c("ND_L", "FD_L", "Feasibility_L")

# Creates a loop to predict coexistence for each combination of selection regimes and environments but for UPPER estimates
struct_mat_w0_U<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_w0_upper[,1])), function(x){
  # Competition matrix
  aux_alpha<-matrix(c(param_all_w0_upper$Te_intra[x],param_all_w0_upper$Te_inter[x],param_all_w0_upper$Tu_inter[x], param_all_w0_upper$Tu_intra[x]), ncol=2, byrow=TRUE)
  # vector of growth rates
  aux_lambda<-c(param_all_w0_upper$Te_lambda[x],param_all_w0_upper$Tu_lambda[x] )
  
  #Estimate niche and fitness differences and feasibility domain
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

# Add column names
colnames(struct_mat_w0_U)<-c("ND_U", "FD_U", "Feasibility_U")

# Join parameter and structural stability estimates (mean, upper and lower) 
struct_mat_w0<-cbind(param_all_w0, struct_mat_w0,struct_mat_w0_L,struct_mat_w0_U)


#' 
#' 
## Estimate distance to the edge----------

### Pooled data----------
#' 
#' We calculate the distance to the edge and the species who will win, following Allen-Perkins et al. 2023; Medeiros et al. 2021.

# Calculate the slope that define the edges of the feasibility domain for each combination of selection regimes 
struct_mat_REP$a21_a11<-struct_mat_REP$Te_inter/struct_mat_REP$Tu_intra
struct_mat_REP$a22_a12<-struct_mat_REP$Te_intra/struct_mat_REP$Tu_inter

## Same as above but for the lower and upper estimates
struct_mat_REP$a21_a11_lower<-param_all_REP_lower$Te_inter/param_all_REP_lower$Tu_intra
struct_mat_REP$a22_a12_lower<-param_all_REP_lower$Te_intra/param_all_REP_lower$Tu_inter

struct_mat_REP$a21_a11_upper<-param_all_REP_upper$Te_inter/param_all_REP_upper$Tu_intra
struct_mat_REP$a22_a12_upper<-param_all_REP_upper$Te_intra/param_all_REP_upper$Tu_inter

# Adding the upper and lower estimates of the intrinsic growth rate of both species
struct_mat_REP$Tu_lambda_lower<-param_all_REP_lower$Tu_lambda
struct_mat_REP$Te_lambda_lower<-param_all_REP_lower$Te_lambda
struct_mat_REP$Tu_lambda_upper<-param_all_REP_upper$Tu_lambda
struct_mat_REP$Te_lambda_upper<-param_all_REP_upper$Te_lambda

##Because of facilitation we need to check that the order is maintained (facilitation sometimes switches the signs)
struct_mat_REP$min_a21_a11<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  # checking the minimum
  min(c(abs(struct_mat_REP$a21_a11_lower[x]),abs(struct_mat_REP$a21_a11_upper[x])))})

struct_mat_REP$min_a22_a12<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  # checking the minimum
  min(c(abs(struct_mat_REP$a22_a12_lower[x]),abs(struct_mat_REP$a22_a12_upper[x])))})

struct_mat_REP$max_a21_a11<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  # checking the maximum
  max(c(abs(struct_mat_REP$a21_a11_lower[x]),abs(struct_mat_REP$a21_a11_upper[x])))})

struct_mat_REP$max_a22_a12<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  #checking the maximum
  max(c(abs(struct_mat_REP$a22_a12_lower[x]),abs(struct_mat_REP$a22_a12_upper[x])))})

#' 
#' 
### Replicate data----------
### Doing the same thing per replicate
# Calculate the slope that define the edges of the feasibility domain for each combination of selection regimes 
struct_mat_w0$a21_a11<-struct_mat_w0$Te_inter/struct_mat_w0$Tu_intra
struct_mat_w0$a22_a12<-struct_mat_w0$Te_intra/struct_mat_w0$Tu_inter
## Same as above but for the lower and upper estimates
struct_mat_w0$a21_a11_lower<-param_all_w0_lower$Te_inter/param_all_w0_lower$Tu_intra
struct_mat_w0$a22_a12_lower<-param_all_w0_lower$Te_intra/param_all_w0_lower$Tu_inter

struct_mat_w0$a21_a11_upper<-param_all_w0_upper$Te_inter/param_all_w0_upper$Tu_intra
struct_mat_w0$a22_a12_upper<-param_all_w0_upper$Te_intra/param_all_w0_upper$Tu_inter

#' 
#' ## Estimating the distance for the pooled data
#' 
# CI is estimated using as reference the vector of the growth rate and then estimating the difference to the lower or upper edges. Negative distances indicate unfeasible systems

# Distance to the edge of Turticae
struct_mat_REP$distanceTu<-sapply(c(1:length(struct_mat_REP$ND)), function(x){
  # testing if its the system is unfeasible. If so then we multiply the results -1
  if(struct_mat_REP$Feasibility[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x],struct_mat_REP$Tu_lambda[x]*struct_mat_REP$a21_a11[x]))))*-1
  }else{ # Otherwise its just the distance
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x],struct_mat_REP$Tu_lambda[x]*struct_mat_REP$a21_a11[x]))))
  }
  a
})

# estimating distance to the Te edge
struct_mat_REP$distanceTe<-sapply(c(1:length(struct_mat_REP$ND)), function(x){ 
  if(struct_mat_REP$Feasibility[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Te_lambda[x]/struct_mat_REP$a22_a12[x],struct_mat_REP$Te_lambda[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Te_lambda[x]/struct_mat_REP$a22_a12[x],struct_mat_REP$Te_lambda[x]))))
  }
  a
})

#If the growth rate vector is within the feasibility domain, it will be closer to the smaller feasibility domain, if the growth rate vector is outside of the feasibility domain, the smallest distance will be to the wider feasibility 

# Estimating distance to the lower bound for the edge of Tu
struct_mat_REP$distanceTu_lower<-sapply(c(1:length(struct_mat_REP$ND)), function(x){
  if(struct_mat_REP$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x],struct_mat_REP$Tu_lambda[x]*struct_mat_REP$a21_a11_lower[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x],struct_mat_REP$Tu_lambda[x]*struct_mat_REP$a21_a11_lower[x]))))
  }
  a
})
# Estimating distance to the upper bound for the edge of Tu
struct_mat_REP$distanceTu_upper<-sapply(c(1:length(struct_mat_REP$ND)), function(x){
  if(struct_mat_REP$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x],struct_mat_REP$Tu_lambda[x]*struct_mat_REP$a21_a11_upper[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x],struct_mat_REP$Tu_lambda[x]*struct_mat_REP$a21_a11_upper[x]))))
  }
  a
})
# Estimating distance to the lower bound for the edge of Te
struct_mat_REP$distanceTe_lower<-sapply(c(1:length(struct_mat_REP$ND)), function(x){
  if(struct_mat_REP$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Te_lambda[x]/struct_mat_REP$a22_a12_lower[x],struct_mat_REP$Te_lambda[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Te_lambda[x]/struct_mat_REP$a22_a12_lower[x],struct_mat_REP$Te_lambda[x]))))
  }
  a
})

# Estimating distance to the upper bound for the edge of Te
struct_mat_REP$distanceTe_upper<-sapply(c(1:length(struct_mat_REP$ND)), function(x){
  if(struct_mat_REP$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Te_lambda[x]/struct_mat_REP$a22_a12_upper[x],struct_mat_REP$Te_lambda[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP$Tu_lambda[x], struct_mat_REP$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP$Te_lambda[x]/struct_mat_REP$a22_a12_upper[x],struct_mat_REP$Te_lambda[x]))))
  }
  a
})


#' 
#' ## Estimating distance to the edge per replicate
#' 
# Creating a loop to estimate the distance to the Tu edge per replicate (same code as above, just applied to a different matrix)
struct_mat_w0$distanceTu<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11[x]))))
  }
  a
})
# Creating a loop to estimate the distance to the Te edge per replicate (same code as above, just applied to a different matrix)
struct_mat_w0$distanceTe<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12[x],struct_mat_w0$Te_lambda[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12[x],struct_mat_w0$Te_lambda[x]))))
  }
  a
})
# Code to estimate distance to lower bounds for Tu
struct_mat_w0$distanceTu_lower<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_lower[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_lower[x]))))
  }
  a
})

# Code to estimate distance to upper bounds for Tu
struct_mat_w0$distanceTu_upper<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_upper[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_upper[x]))))
  }
  a
} )
# Code to estimate distance to lower bounds for Te
struct_mat_w0$distanceTe_lower<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12_lower[x],struct_mat_w0$Te_lambda[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12_lower[x],struct_mat_w0$Te_lambda[x]))))
  }
  
  a
} )

# Code to estimate distance to upper bounds for Te
struct_mat_w0$distanceTe_upper<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12_upper[x],struct_mat_w0$Te_lambda[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12_upper[x],struct_mat_w0$Te_lambda[x]))))
  }
})


#' 
#' 
## Saving----------
#' 
## 
# If you are running the code to evaluate and store data 
if(evaluation){
  # write the results for the pooled data
  write.csv(struct_mat_REP, "./Analyses/structural_Pooled.csv")
  # write the results for the replicate data
  write.csv(struct_mat_w0, "./Analyses/structural_Replicates.csv")
}

# # Because of facilitation (specially when there is double facilitation estimated) the distance of the upper and lower estimates is switched. 
struct_mat_REP_final<-struct_mat_REP
# loop to order the results
aux_Tu_orders<-as.data.frame(t(sapply(c(1:nrow(struct_mat_REP_final)), function(x){
  # vector of distances
  vect<-c(struct_mat_REP_final$distanceTu[x],struct_mat_REP_final$distanceTu_lower[x], struct_mat_REP_final$distanceTu_upper[x])
  ord<-order(vect,decreasing = FALSE)
  
  # ordered vector
  vect[ord]
})))

# adding column names
colnames(aux_Tu_orders)<-c("Tu_distance_lower","Tu_distance", "Tu_distance_upper" )

# same thing as above but now for distance of Te
aux_Te_orders<-as.data.frame(t(sapply(c(1:nrow(struct_mat_REP_final)), function(x){
  vect<-c(struct_mat_REP_final$distanceTe[x],struct_mat_REP_final$distanceTe_lower[x], struct_mat_REP_final$distanceTe_upper[x])
  ord<-order(vect,decreasing = FALSE)

  vect[ord]
})))

# adding column names
colnames(aux_Te_orders)<-c("Te_distance_lower","Te_distance", "Te_distance_upper" )

# Ordering but now for the replicates. For Tu
aux_Tu_orders_w0<-as.data.frame(t(sapply(c(1:nrow(struct_mat_w0)), function(x){
  vect<-c(struct_mat_w0$distanceTu[x],struct_mat_w0$distanceTu_lower[x], struct_mat_w0$distanceTu_upper[x])
  ord<-order(vect,decreasing = FALSE)

  vect[ord]
})))

colnames(aux_Tu_orders_w0)<-c("Tu_distance_lower","Tu_distance", "Tu_distance_upper" )

# Ordering but now for the replicates. For Te
aux_Te_orders_w0<-as.data.frame(t(sapply(c(1:nrow(struct_mat_w0)), function(x){
  vect<-c(struct_mat_w0$distanceTe[x],struct_mat_w0$distanceTe_lower[x], struct_mat_w0$distanceTe_upper[x])
  ord<-order(vect,decreasing = FALSE)

  vect[ord]
})))

colnames(aux_Te_orders_w0)<-c("Te_distance_lower","Te_distance", "Te_distance_upper" )

# Creating a data structure with the information of the distances (only)

struct_mat_REP_final<-as.data.frame(cbind(struct_mat_REP[,c(1,2,3)], aux_Tu_orders, aux_Te_orders))

struct_mat_w0_final<-as.data.frame(cbind(struct_mat_w0[,c(1,2,3,4)], aux_Tu_orders_w0, aux_Te_orders_w0))

#' 
## Estimating the minimum distance to the edge----------
#' To know what is the likelihood for a system to leave (or enter the feasibility domain) we will estimate the minimum distance to the edges of feasibility domain. This could be the distance from the edge that corresponds to the competitive interactions for Te or Tu.
#' 
### minimum distance
struct_mat_REP_final2<-struct_mat_REP_final

# Loop to test which is the minimum distance. To maintain consistency the upper and lower estimates will always correspond to the same edge (otherwise it would not be consistent)
aux_min_distance<-as.data.frame(t(sapply(c(1:nrow(struct_mat_REP_final2)), function(x){
  # Checking which is the minimum distance
  a<-min(c(abs(struct_mat_REP_final2$Tu_distance[x]),abs(struct_mat_REP_final2$Te_distance[x])))
  
  # store the name of the column that has the lower distance
  a_name<-colnames(struct_mat_REP_final2)[which(struct_mat_REP_final2[x,]==a | struct_mat_REP_final2[x,]==-a)]
  
  # If its Tu then we store the upper, lower and mean estimated distance of Tu
  if(a_name=="Tu_distance"){
    vecF<-c(struct_mat_REP_final2$Tu_distance_lower[x], struct_mat_REP_final2$Tu_distance[x], struct_mat_REP_final2$Tu_distance_upper[x])
  }else if(a_name=="Te_distance"){ # same as above, but for Te
    vecF<- c(struct_mat_REP_final2$Te_distance_lower[x], struct_mat_REP_final2$Te_distance[x], struct_mat_REP_final2$Te_distance_upper[x])
  }
  
  vecF
})))

#adding columns
colnames(aux_min_distance)<-c("minDistance_L", "minDistance", "minDistance_U")

# Loop to test which is the minimum distance for REPLICATES. To maintain consistency the upper and lower estimates will always correspond to the same edge (otherwise it would not be consistent). The code is similar, just applied to the replicate data frame
aux_min_distance_w0<-as.data.frame(t(sapply(c(1:nrow(struct_mat_w0_final)), function(x){
  a<-min(c(abs(struct_mat_w0_final$Tu_distance[x]),abs(struct_mat_w0_final$Te_distance[x])))
  
  a_name<-colnames(struct_mat_w0_final)[which(struct_mat_w0_final[x,]==a | struct_mat_w0_final[x,]==-a)][1]
  
  if(a_name=="Tu_distance" |a_name=="Tu_distance_lower" | a_name=="Tu_distance_upper" ){
    vecF<-c(struct_mat_w0_final$Tu_distance_lower[x],struct_mat_w0_final$Tu_distance[x] , struct_mat_w0_final$Tu_distance_upper[x])
  }else if(a_name=="Te_distance" |a_name=="Te_distance_lower" | a_name=="Te_distance_upper"){
    vecF<- c(struct_mat_w0_final$Te_distance_lower[x], struct_mat_w0_final$Te_distance[x], struct_mat_w0_final$Te_distance_upper[x])
  }
  
  vecF
})))
# Adding column names
colnames(aux_min_distance_w0)<-c("minDistance_L", "minDistance", "minDistance_U")

# Joining the data to a single data frame
struct_mat_REP_final<-as.data.frame(cbind(struct_mat_REP_final, aux_min_distance))

struct_mat_w0_final<-as.data.frame(cbind(struct_mat_w0_final, aux_min_distance_w0))

if(evaluation){
# Write for the figure
write.csv(struct_mat_REP_final,"./Analyses/min_distance_pooled.csv")
write.csv(struct_mat_w0_final,"./Analyses/min_distance_per_replicate.csv")
}
print("Population growth experiments")
# Population growth experiments ----------
#' 
#' 
#' # Testing fit of predictions
#' We will test our ability to predict the population dynamics data using the estimated parameters.
#' 
#' ## Importing data
#' 
coex_g42<-read.csv("./Data/Coexistence_Cd_G42_submit.csv", header=TRUE) # Data from the population dynamics experiment

# Transforming columns into factors
coex_g42$Rep2<-as.factor(coex_g42$Rep)
coex_g42$X1st.pair<-as.factor(coex_g42$X1st_pair)
coex_g42$X2nd.pair<-as.factor(coex_g42$X2nd_pair)
coex_g42$SRTu<-as.factor(coex_g42$SRTu)
coex_g42$SRTe<-as.factor(coex_g42$SRTe)
coex_g42$Box2<-as.factor(coex_g42$Box)

### summary data per leaf (because the leaflets are not attributable)
coex_g42_res<-gather(coex_g42, leaf, females, Leaf_2_Up_Tu:Leaf_5_Down_Te, factor_key=TRUE)
str(coex_g42_res)

# transform the leaf into character so we can subset it
coex_g42_res$char<-as.character(coex_g42_res$leaf)

# Loop to create a data frame to separate the leaf into the leaf direction and species counted.
aux4<-as.data.frame(t(as.data.frame(sapply(c(1:length(coex_g42_res$Rep)), function(x){
  a<-strsplit(coex_g42_res$char[x], split="_")[[1]]
  
  c(a[2:4])
}))))
# adding column names
colnames(aux4)<-c("Leaf2", "Direction", "Species")

# Joining the two data frames so we have all information to then summarize
coex_g42_res<-cbind(coex_g42_res, aux4)

# Summarizing data by selection regimes, replicate, environment, box, leaf, direction and environment. This sums the number of individuals found in the leaf
sum_coex_g42<-coex_g42_res %>%
  group_by(Rep2, SRTu, SRTe, Box2, Leaf2,X1st.pair,X2nd.pair, Direction, Species, Env) %>%
  summarize(av_females=sum(females, na.rm=TRUE))

# Separate into columns the information from different species 
sum_coex_g42_res<-as.data.frame(spread(sum_coex_g42, key=Species, value=av_females))

# Calculate the sum of females of each species per leaf
coex_g42_rep<-sum_coex_g42_res %>%
  group_by(Rep2, Leaf2, SRTu, SRTe, Env, Box2) %>%
  summarize( sum_Te=sum(Te, na.rm=TRUE), sum_Tu=mean(Tu, na.rm=TRUE)) %>% as.data.frame()

# Ensure that the Cd is well written
coex_g42_rep$Env[which(coex_g42_rep$Env=="Cd")]<-"Cd"

#' 
#' ### Can we predict the outcome of species interactions?
#' 
### Predicting first generation pooled data----------
#' 
#' Here we will use the parameter estimates and predict the number of offspring produced after one generation.
#' The code applied will always be the same to Te and Tu and lower and upper bounds, the only thing that changes are the estimates used (that match the species and interval considered).
#' 
# creating the data frame to predict
pred_coex_RK_REP<-expand_grid(Te=c("SR4","SR5"), Tu=c("SR1", "SR2"), Environment= c("N", "Cd"))

# Doing the first generation of predictions for Tu
pred_coex_RK_REP$predTu1<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  # Subset data
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  # Predict offspring production
  bl<-aux_alphas$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6 - aux_alphas$Tu_inter[1]*6)
  
  bl
})

# Predicting first generation for Te
pred_coex_RK_REP$predTe1<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})

#' 
### Predicting second generation pooled data----------
#' 
#' Here we will use the parameter estimates and predict the number of offspring produced after two generations. For that we use the number of individuals predicted in the previous generation and use that as competitors.
#' The code applied will always be the same to Te and Tu and lower and upper bounds, the only thing that changes are the estimates used (that match the species and interval considered).

# Predicting Tu
pred_coex_RK_REP$predTu2<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  # subset data to get parameters
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  # Use the riker model with the number of offspring predicted for the previous generation
  bl<-aux_alphas$Tu_lambda[1]*pred_coex_RK_REP$predTu1[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_REP$predTu1[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_REP$predTe1[x])
  
  bl
})
# predicting second generation for Te
pred_coex_RK_REP$predTe2<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  # subset data to get parameters
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
   # Use the riker model with the number of offspring predicted for the previous generation
  bl<-aux_alphas$Te_lambda[1]*pred_coex_RK_REP$predTe1[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_REP$predTe1[x] - aux_alphas$Te_inter[1]*pred_coex_RK_REP$predTu1[x])
  
  bl
})

#' 
### Predicting for lower and upper boundaries----------
#' 
#lower - stronger alpha and lower lambda
# Tu populations first generation
pred_coex_RK_REP$predTu1_L<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6- aux_alphas$Tu_inter[1]*6)
  
  bl
})
# Te populations first generation
pred_coex_RK_REP$predTe1_L<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})

# Tu populations second generation
pred_coex_RK_REP$predTu2_L<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Tu_lambda[1]*pred_coex_RK_REP$predTu1_L[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_REP$predTu1_L[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_REP$predTe1_L[x])
  
  bl
})

# Te populations second generation
pred_coex_RK_REP$predTe2_L<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Te_lambda[1]*pred_coex_RK_REP$predTe1_L[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_REP$predTe1_L[x] - aux_alphas$Te_inter[1]*pred_coex_RK_REP$predTu1_L[x])
  
  bl
})

#' 
#' Upper estimates
##
# upper bounds
# Tu populations first generation
pred_coex_RK_REP$predTu1_U<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6- aux_alphas$Tu_inter[1]*6)
  
  bl
})
# Te populations first generation
pred_coex_RK_REP$predTe1_U<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})

# Tu populations second generation
pred_coex_RK_REP$predTu2_U<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Tu_lambda[1]* pred_coex_RK_REP$predTu1_U[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_REP$predTu1_U[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_REP$predTe1_U[x])
  
  bl
})
# Te populations second generation
pred_coex_RK_REP$predTe2_U<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Te_lambda[1]*pred_coex_RK_REP$predTe1_U[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_REP$predTe1_U[x] - aux_alphas$Te_inter[1]*pred_coex_RK_REP$predTu1_U[x])
  
  bl
})

# rename to match the other data frame
names(pred_coex_RK_REP)[1:3]<-c("SRTe", "SRTu", "Env")
# transform to a data frame
pred_coex_RK_REP<-as.data.frame(pred_coex_RK_REP)

#' 
### Predicting first generation per replicate----------
#' We will use the same approach as we did for the pooled data
#' 

# creating the data frame to store data
pred_coex_RK_w0<-as.data.frame(expand_grid(Te=c("SR4","SR5"), Tu=c("SR1", "SR2"), Environment= c("N", "Cd"), Replicate=c(1,2,3,4,5)))

# Removing replicate two for the SR2
pred_coex_RK_w0<- pred_coex_RK_w0[- which(pred_coex_RK_w0$Replicate==2 & pred_coex_RK_w0$Tu=="SR2" & pred_coex_RK_w0$Environment=="Cd"),]

# predicting Tu for first generation
pred_coex_RK_w0$predTu1<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_alphas$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6 - aux_alphas$Tu_inter[1]*6)
  
  bl
})
# predicting Te for first generation
pred_coex_RK_w0$predTe1<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_alphas$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})

# predicting Tu for second generation
pred_coex_RK_w0$predTu2<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_alphas$Tu_lambda[1]*pred_coex_RK_w0$predTu1[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_w0$predTu1[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_w0$predTe1[x])
  
  bl
})
# predicting Te for second generation
pred_coex_RK_w0$predTe2<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_alphas$Te_lambda[1]*pred_coex_RK_w0$predTe1[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_w0$predTe1[x] - aux_alphas$Te_inter[1]*pred_coex_RK_w0$predTu1[x])
  
  bl
})

#' 
### Predicting for lower and upper boundaries----------
#' 

# lower - stronger alpha and lower lambda
# predicting Tu for first generation
pred_coex_RK_w0$predTu1_L<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6- aux_alphas$Tu_inter[1]*6)
  
  bl
})
# predicting Te for first generation
pred_coex_RK_w0$predTe1_L<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})

# predicting Tu for second generation
pred_coex_RK_w0$predTu2_L<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Tu_lambda[1]*pred_coex_RK_w0$predTu1_L[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_w0$predTu1_L[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_w0$predTe1_L[x])
  
  bl
})
# predicting Te for second generation
pred_coex_RK_w0$predTe2_L<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Te_lambda[1]*pred_coex_RK_w0$predTe1_L[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_w0$predTe1_L[x] - aux_alphas$Te_inter[1]*pred_coex_RK_w0$predTu1_L[x])
  
  bl
})

#' 
#' Upper estimates

# upper
# predicting Tu for first generation
pred_coex_RK_w0$predTu1_U<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6- aux_alphas$Tu_inter[1]*6)
  
  bl
})
# predicting Te for first generation
pred_coex_RK_w0$predTe1_U<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})

# predicting Tu for second generation
pred_coex_RK_w0$predTu2_U<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Tu_lambda[1]* pred_coex_RK_w0$predTu1_U[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_w0$predTu1_U[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_w0$predTe1_U[x])
  
  bl
})
# predicting Te for second generation
pred_coex_RK_w0$predTe2_U<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Te_lambda[1]*pred_coex_RK_w0$predTe1_U[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_w0$predTe1_U[x] - aux_alphas$Te_inter[1]*pred_coex_RK_w0$predTu1_U[x])
  
  bl
})

# replacing names so they match the data frame for the empirical data
names(pred_coex_RK_w0)[1:3]<-c("SRTe", "SRTu", "Env")

pred_coex_RK_w0<-as.data.frame(pred_coex_RK_w0)

#' 
## Estimating number of observed females----------
#' Using the data from the population dynamics

# Summarizing per box
coex_g42_rep2<-coex_g42_rep %>%
  group_by(SRTu, SRTe, Rep2, Box2, Env) %>%
  summarize(sumTe=sum(sum_Te, na.rm=TRUE), sumTu=sum(sum_Tu, na.rm=TRUE))

# Estimating the proportion of Te in the box
coex_g42_rep2$Te_ratio<-sapply(c(1:dim(coex_g42_rep2)[1]), function(x) coex_g42_rep2$sumTe[x]/sum(coex_g42_rep2$sumTu[x],coex_g42_rep2$sumTe[x]))

# Estimating mean ratio per box, or summing the number of individuals in all the boxes to then calculate the proportion (because mean of ratios sometimes don't behave well)
coex_g42_rep3<-coex_g42_rep2 %>%
  group_by(SRTu, SRTe, Rep2, Env) %>%
  summarize(sum_Te=sum(sumTe, na.rm=TRUE), sum_Tu=sum(sumTu, na.rm=TRUE), sdTe=sd(sumTe, na.rm=TRUE), sdTu=sd(sumTu, na.rm=TRUE), meanTeRatio=mean(Te_ratio, na.rm=TRUE))

# estimating the proportion of Te females
coex_g42_rep3$Te_ratio<-sapply(c(1:dim(coex_g42_rep3)[1]), function(x) coex_g42_rep3$sum_Te[x]/sum(coex_g42_rep3$sum_Tu[x],coex_g42_rep3$sum_Te[x]))

# Changing the factor levels so they match between empirical and predicted data frames
coex_g42_rep3$Env<-plyr::mapvalues(coex_g42_rep3$Env, c("Cd","water"), c("Cd", "N"))
coex_g42_rep3$SRTu2<-plyr::mapvalues(coex_g42_rep3$SRTu, c("Tu1","Tu2"), c("SR1", "SR2"))
coex_g42_rep3$SRTe2<-plyr::mapvalues(coex_g42_rep3$SRTe, c("Te4","Te5"), c("SR4", "SR5"))

# This is just a new variable created to make the code cleaner
sum_observed_coex<-coex_g42_rep3


#' 
# Summarizing data at the level of the Regime
sum_observed_coex2<-sum_observed_coex %>%
  group_by(SRTu2, SRTe2, Env)%>%
  summarise(sumTe=sum(sum_Te, na.rm=TRUE),sumTu=mean(sum_Tu, na.rm=TRUE), sdTe2=sd(sum_Te, na.rm=TRUE)/sqrt(5), sdTu2=sd(sum_Tu, na.rm=TRUE)/sqrt(5), TeRatio=mean(meanTeRatio, na.rm=TRUE), sdTeRatio2=sd(meanTeRatio, na.rm=TRUE)/sqrt(5)) %>% as.data.frame()

# estimating the standard deviation surrounding the ratio
sum_observed_coex2$TeRatio_L<-sum_observed_coex2$TeRatio-sum_observed_coex2$sdTeRatio2
sum_observed_coex2$TeRatio_U<-sum_observed_coex2$TeRatio+sum_observed_coex2$sdTeRatio2

#' 
#' Setting up data frames to join data
## Estimating predicted proportion of Te for the pooled data
pred_coex_RK_REP$TeRatio<-sapply(c(1:dim(pred_coex_RK_REP)[1]), function(x){
  pred_coex_RK_REP$predTe2[x]/(pred_coex_RK_REP$predTe2[x]+pred_coex_RK_REP$predTu2[x])
})
## Estimating predicted proportion of Te for each replicate
pred_coex_RK_w0$TeRatio<-sapply(c(1:dim(pred_coex_RK_w0)[1]), function(x){
  pred_coex_RK_w0$predTe2[x]/(pred_coex_RK_w0$predTe2[x]+pred_coex_RK_w0$predTu2[x])
})

# Same as above but for the lower bounds - pooled data
pred_coex_RK_REP$TeRatio_L<-sapply(c(1:dim(pred_coex_RK_REP)[1]), function(x){
  pred_coex_RK_REP$predTe2_L[x]/(pred_coex_RK_REP$predTe2_L[x]+pred_coex_RK_REP$predTu2_L[x])
})
# Same as above but for the lower bounds - replicate data
pred_coex_RK_w0$TeRatio_L<-sapply(c(1:dim(pred_coex_RK_w0)[1]), function(x){
  pred_coex_RK_w0$predTe2_L[x]/(pred_coex_RK_w0$predTe2_L[x]+pred_coex_RK_w0$predTu2_L[x])
})

# Same as above but for the upper bounds - pooled data
pred_coex_RK_REP$TeRatio_U<-sapply(c(1:dim(pred_coex_RK_REP)[1]), function(x){
  pred_coex_RK_REP$predTe2_U[x]/(pred_coex_RK_REP$predTe2_U[x]+pred_coex_RK_REP$predTu2_U[x])
})
# Same as above but for the upper bounds - replicate data
pred_coex_RK_w0$TeRatio_U<-sapply(c(1:dim(pred_coex_RK_w0)[1]), function(x){
  pred_coex_RK_w0$predTe2_U[x]/(pred_coex_RK_w0$predTe2_U[x]+pred_coex_RK_w0$predTu2_U[x])
})

#' 
#' 
# Joining empirical and predicted data----------
#Predicted data
sum_pred_coex_RK_REP<-pred_coex_RK_REP %>%
  group_by(SRTu, SRTe, Env)%>%
  summarise(predTe=mean(predTe2, na.rm=TRUE),predTu=mean(predTu2, na.rm=TRUE), sumTeRatio=(sum(predTe2, na.rm=TRUE)/(sum(predTe2, na.rm=TRUE)+sum(predTu2, na.rm=TRUE)))) %>% as.data.frame()

## Empirical data
sum_observed_coex_rep2<-sum_observed_coex %>%
  group_by(SRTe, SRTu, Env) %>%
  summarize(obs_TeRatio=mean(meanTeRatio), SE_obs=sd(meanTeRatio)/sqrt(5), meanTe=mean(sum_Te, na.rm=TRUE), meanTu=mean(sum_Tu, na.rm=TRUE)) %>% as.data.frame()

# Creating columns to match the levels of the factors in the two data frames
sum_observed_coex_rep2$SRTe2<-(plyr::mapvalues(as.character(sum_observed_coex_rep2$SRTe), c("Te4","Te5"), c("SR4", "SR5")))
sum_observed_coex_rep2$SRTu2<-(plyr::mapvalues(as.character(sum_observed_coex_rep2$SRTu), c("Tu1","Tu2"), c("SR1", "SR2")))
colnames(sum_observed_coex_rep2)[c(1:2, 8,9)]<-c("SRTe2", "SRTu2","SRTe", "SRTu" )

# Renaming column names so they match
sum_observed_coex_rep2<-sum_observed_coex_rep2[,c(8,9,3:7)]
colnames(sum_observed_coex)[c(1,2,3,11:12)]<-c("SRTu2","SRTe2","Replicate","SRTu","SRTe")

#Join the two data frames at the pooled data level (i.e. per regime)
sum_observed_coex_rep<-inner_join(sum_observed_coex_rep2, pred_coex_RK_REP, by=c("SRTe", "SRTu", "Env"))

# thinning down to remove unnecessary variables
sum_observed_coex_rep<-as.data.frame(sum_observed_coex_rep[,c("SRTe", "SRTu", "Env", "obs_TeRatio","SE_obs", "TeRatio", "TeRatio_L", "TeRatio_U", "meanTe","meanTu", "predTu2","predTe2")])

colnames(sum_observed_coex_rep)[6:8]<-c("pred_T1", "T1_L","T1_U")

# Transform into a factor
pred_coex_RK_w0$Replicate<-as.factor(pred_coex_RK_w0$Replicate)

# Join everything at the replicate level
sum_observed_coex_ALL<-inner_join(sum_observed_coex, pred_coex_RK_w0, by=c("SRTe", "SRTu", "Env", "Replicate"))

# thinning down to remove unnecessary variables
sum2_observed_coex_ALL<-as.data.frame(sum_observed_coex_ALL[,c(11,12,3,4,5,6,9,15,16)])

# Estimating proportion of T evansi
sum2_observed_coex_ALL$predTeRatio<-sapply(c(1:dim(sum2_observed_coex_ALL)[1]), function(x){sum2_observed_coex_ALL$predTe2[x]/sum(sum2_observed_coex_ALL$predTe2[x],sum2_observed_coex_ALL$predTu2[x])})



## Comparing to data---------
#' Here we will compare the predicted and observed proportion of Te females for each combination of regimes, at the replicate and pooled data
#' 
#' #### Testing pooled replicates
#' 
#This forces the line to pass by the 0,0, as the reviewer suggested

m7<-glmmTMB(obs_TeRatio~0+pred_T1, data=sum_observed_coex_rep, family="Gamma")
summary(m7)
emtrends(m7, var="pred_T1", type="response")
test(emtrends(m7, var="pred_T1", type="response"))

# writing file
if(evaluation){
write.csv(sum_observed_coex_rep, file="./Analyses/popDyn_pooledData.csv")
}
print(paste("Slope pooled:", round(log(coefficients(summary(m7))$cond[1]),3), sep=" "))

#' 
#' #### Testing per replicate
#' 

sum_observed_coex_ALL2<-sum2_observed_coex_ALL %>%
  group_by(SRTe, SRTu, Env) %>%
  summarize(meanRatio=mean(meanTeRatio, na.rm=TRUE), mean_pred=mean(predTeRatio, na.rm=TRUE), sdRatio=sd(meanTeRatio, na.rm=TRUE)/sqrt(5), sdPred=sd(predTeRatio, na.rm=TRUE)/sqrt(5), meanTe=mean(sum_Te, na.rm=TRUE), meanTu=mean(sum_Tu, na.rm=TRUE))

if(evaluation){
#writing file
write.csv(sum_observed_coex_ALL2, file="./Analyses/popDyn_replicates.csv")
}

m7_rep2<-glmmTMB(meanRatio~0+mean_pred, data=sum_observed_coex_ALL2, family="Gamma")
summary(m7_rep2)
emtrends(m7_rep2, var="mean_pred", type="response")
test(emtrends(m7_rep2, var="mean_pred", type="response"))

print(paste("Slope replicates:", round(log(coefficients(summary(m7_rep2))$cond[1]),3), sep=" "))
#' 
#' 
#' 
#' 
# Figures --------
# Script to run and get all figures, except figure S2
source("./Code/Figures.R")

