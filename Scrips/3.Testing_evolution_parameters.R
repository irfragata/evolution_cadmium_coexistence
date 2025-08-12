rm(list=ls())

library(plyr)
library(tidyverse)
library(car)
library(fitdistrplus)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(MASS)
library(DescTools)
library(performance)
library(DHARMa)
library(effects)
library(marginaleffects)
library(LSAfun)
library(arm)


###### Importing data files for replicate estimation

## Importing
param_all_w0<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal.csv")
param_all_w0_upper<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal_upper.csv")
param_all_w0_lower<-read.csv( "./Analyses/cxr_normal/parameters_cxr_normal_lower.csv")

param_all_w0<-param_all_w0[,-1]
param_all_w0_upper<-param_all_w0_upper[,-1]
param_all_w0_lower<-param_all_w0_lower[,-1]


## Testing differences in parameters (estimated)

#### Distribution
str(param_all_w0)

descdist(param_all_w0$Tu_lambda, discrete=FALSE, boot=1000)
descdist(param_all_w0$Te_lambda, discrete=FALSE, boot=1000)

descdist(param_all_w0$Tu_intra, discrete=FALSE, boot=1000)
descdist(param_all_w0$Te_intra, discrete=FALSE, boot=1000)

descdist(param_all_w0$Tu_inter, discrete=FALSE, boot=1000)
descdist(param_all_w0$Te_inter, discrete=FALSE, boot=1000)

hist(param_all_w0$Te_lambda)

hist(param_all_w0$Tu_intra)
hist(param_all_w0$Te_intra)

hist(param_all_w0$Tu_inter)
hist(param_all_w0$Te_inter)

#### Does cadmium change parameters? (no evolution)

gr_tu_cd_1<-glmmTMB(Tu_lambda~Environment, data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ))
gr_tu_cd_2<-glmmTMB(Tu_lambda~Environment, data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=Gamma(link="log"))
gr_tu_cd_3<-glmmTMB(Tu_lambda~Environment, data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=gaussian(link="log"))
anova(gr_tu_cd_1,gr_tu_cd_2,gr_tu_cd_3)

summary(gr_tu_cd_2)

simulationOutput <- simulateResiduals(fittedModel = gr_tu_cd_2, plot = F)
plot(simulationOutput)
#No problems

gr_te_cd_1<-glmmTMB(Te_lambda~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ))
gr_te_cd_2<-glmmTMB(Te_lambda~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=Gamma(link="log"))
gr_te_cd_3<-glmmTMB(Te_lambda~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=gaussian(link="log"))
anova(gr_te_cd_1,gr_te_cd_2,gr_te_cd_3)

summary(gr_te_cd_2)
simulationOutput <- simulateResiduals(fittedModel = gr_te_cd_2, plot = F)
plot(simulationOutput)# No problems

########

intra_tu_cd_1<-glmmTMB(Tu_intra~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ))
intra_tu_cd_2<-glmmTMB(Tu_intra+1~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=Gamma(link="log"))
intra_tu_cd_3<-glmmTMB(Tu_intra+1~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=gaussian(link="log"))
anova(intra_tu_cd_1,intra_tu_cd_2,intra_tu_cd_3)

summary(intra_tu_cd_1) 

simulationOutput <- simulateResiduals(fittedModel = intra_tu_cd_1, plot = F)
plot(simulationOutput)# No problems



#Gamma and gaussian give very similar estimates and values
intra_te_cd_1<-glmmTMB(Te_intra~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ))
intra_te_cd_2<-glmmTMB(Te_intra+1~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=Gamma(link="log"))
intra_te_cd_3<-glmmTMB(Te_intra+1~Environment,  data=subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" ), family=gaussian(link="log"))
anova(intra_te_cd_1,intra_te_cd_2,intra_te_cd_3)

summary(intra_te_cd_2)
summary(intra_te_cd_1) #Again very similar estimates

simulationOutput <- simulateResiduals(fittedModel = intra_te_cd_1, plot = F)
plot(simulationOutput) #No problems

######
inter_tu_cd_1<-glmmTMB(Tu_inter~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")))
inter_tu_cd_2<-glmmTMB(Tu_inter+1~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")), family=Gamma(link="log"))
inter_tu_cd_3<-glmmTMB(Tu_inter+1~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")), family=gaussian(link="log"))
anova(inter_tu_cd_1,inter_tu_cd_2,inter_tu_cd_3)

summary(inter_tu_cd_2)
summary(inter_tu_cd_1)

simulationOutput <- simulateResiduals(fittedModel = inter_tu_cd_1, plot = F)
plot(simulationOutput) #no problems

inter_te_cd_1<-glmmTMB(Te_inter~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")))
inter_te_cd_2<-glmmTMB(Te_inter+1~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")), family=Gamma(link="log"))
inter_te_cd_3<-glmmTMB(Te_inter+1~Environment, data=subset(param_all_w0, (Tu_Regime=="SR1" & Te_Regime=="SR4")), family=gaussian(link="log"))
anova(inter_te_cd_1,inter_te_cd_2,inter_te_cd_3)

summary(inter_te_cd_2)
summary(inter_te_cd_1)

simulationOutput <- simulateResiduals(fittedModel = inter_te_cd_1, plot = F)
plot(simulationOutput) #No problems

###### Summary
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


#### Does evolution change the performance?

gr_tu_ev_1<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="Cd" & Te_Regime=="SR4"))
gr_tu_ev_2<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"), family=Gamma(link="log"))
gr_tu_ev_3<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"), family=gaussian(link="log"))
anova(gr_tu_ev_1,gr_tu_ev_2,gr_tu_ev_3)

summary(gr_tu_ev_2)
summary(gr_tu_ev_1)

simulationOutput <- simulateResiduals(fittedModel = gr_tu_ev_2, plot = F)
plot(simulationOutput) #no problems

gr_te_ev_1<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"))
gr_te_ev_2<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"), family=Gamma(link="log"))
gr_te_ev_3<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"), family=gaussian(link="log"))
anova(gr_te_ev_1,gr_te_ev_2,gr_te_ev_3)

summary(gr_te_ev_2)
summary(gr_te_ev_1)

simulationOutput <- simulateResiduals(fittedModel = gr_te_ev_2, plot = F)
plot(simulationOutput)# No problems

## intra

intra_tu_ev_1<-glmmTMB(Tu_intra~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"))
intra_tu_ev_2<-glmmTMB(Tu_intra+1~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"), family=Gamma(link="log"))
intra_tu_ev_3<-glmmTMB(Tu_intra+1~Tu_Regime, data=subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4"), family=gaussian(link="log"))
anova(intra_tu_ev_1,intra_tu_ev_2,intra_tu_ev_3)

summary(intra_tu_ev_2)
summary(intra_tu_ev_1)

simulationOutput <- simulateResiduals(fittedModel = intra_tu_ev_1, plot = F)
plot(simulationOutput) #no problems

intra_te_ev_1<-glmmTMB(Te_intra~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"))
intra_te_ev_2<-glmmTMB(Te_intra+1~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"), family=Gamma(link="log"))
intra_te_ev_3<-glmmTMB(Te_intra+1~Te_Regime, data=subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1"), family=gaussian(link="log"))
anova(intra_te_ev_1,intra_te_ev_2,intra_te_ev_3)

summary(intra_te_ev_3)
summary(intra_te_ev_1)

simulationOutput <- simulateResiduals(fittedModel = intra_te_ev_1, plot = F)
plot(simulationOutput) #no problems

## inter

inter_tu_ev_1<-glmmTMB(Tu_inter~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="Cd"))
inter_tu_ev_2<-glmmTMB(Tu_inter+1~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="Cd"), family=Gamma(link="log"))
inter_tu_ev_3<-glmmTMB(Tu_inter+1~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="Cd"), family=gaussian(link="log"))
anova(inter_tu_ev_1,inter_tu_ev_2,inter_tu_ev_3)

summary(inter_tu_ev_3)
summary(inter_tu_ev_1)

simulationOutput <- simulateResiduals(fittedModel = inter_tu_ev_1, plot = F)
plot(simulationOutput) #no problem

inter_te_ev_1<-glmmTMB(Te_inter~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="Cd"))
inter_te_ev_2<-glmmTMB(Te_inter+1~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="Cd"), family=Gamma(link="log"))
inter_te_ev_3<-glmmTMB(Te_inter+1~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="Cd"), family=gaussian(link="log"))
anova(inter_te_ev_1,inter_te_ev_2,inter_te_ev_3)

summary(inter_te_ev_3)
summary(inter_te_ev_1)

simulationOutput <- simulateResiduals(fittedModel = inter_te_ev_1, plot = F, )
#plot(simulationOutput)
plotResiduals(simulationOutput, subset(param_all_w0, Environment=="Cd")$Te_Regime)
plotResiduals(simulationOutput, subset(param_all_w0, Environment=="Cd")$Tu_Regime)

# No problem


###### Summary

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


#### 2.3.4 - Does evolution change the performance in the ancestral environment?

gr_tu_an_1<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="N" & Te_Regime=="SR4"))
gr_tu_an_2<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"), family=Gamma(link="log"))
gr_tu_an_3<-glmmTMB(Tu_lambda~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"), family=gaussian(link="log"))
anova(gr_tu_an_1,gr_tu_an_2,gr_tu_an_3)

summary(gr_tu_an_2)
summary(gr_tu_an_1)

simulationOutput <- simulateResiduals(fittedModel = gr_tu_an_2, plot = F, )
plot(simulationOutput)# no problem

gr_te_an_1<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"))
gr_te_an_2<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"), family=Gamma(link="log"))
gr_te_an_3<-glmmTMB(Te_lambda~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"), family=gaussian(link="log"))
anova(gr_te_an_1,gr_te_an_2,gr_te_an_3)

summary(gr_te_an_2)

simulationOutput <- simulateResiduals(fittedModel = gr_te_an_2, plot = F, )
plot(simulationOutput)
## intra

intra_tu_an_1<-glmmTMB(Tu_intra~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"))
intra_tu_an_2<-glmmTMB(Tu_intra+1~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"), family=Gamma(link="log"))
intra_tu_an_3<-glmmTMB(Tu_intra+1~Tu_Regime, data=subset(param_all_w0, Environment=="N"& Te_Regime=="SR4"), family=gaussian(link="log"))
anova(intra_tu_an_1,intra_tu_an_2,intra_tu_an_3)

summary(intra_tu_an_3)
summary(intra_tu_an_1)

simulationOutput <- simulateResiduals(fittedModel = intra_tu_an_1, plot = F, )
plot(simulationOutput) #no problem

intra_te_an_1<-glmmTMB(Te_intra~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"))
intra_te_an_2<-glmmTMB(Te_intra+1~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"), family=Gamma(link="log"))
intra_te_an_3<-glmmTMB(Te_intra+1~Te_Regime, data=subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1"), family=gaussian(link="log"))
anova(intra_te_an_1,intra_te_an_2,intra_te_an_3)

summary(intra_te_an_3)
summary(intra_te_an_1)

simulationOutput <- simulateResiduals(fittedModel = intra_te_an_1, plot = F, )
plot(simulationOutput)# no problem


## inter

inter_tu_an_1<-glmmTMB(Tu_inter~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="N"))
inter_tu_an_2<-glmmTMB(Tu_inter+1~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="N"), family=Gamma(link="log"))
inter_tu_an_3<-glmmTMB(Tu_inter+1~Tu_Regime*Te_Regime, data=subset(param_all_w0, Environment=="N"), family=gaussian(link="log"))
anova(inter_tu_an_1,inter_tu_an_2,inter_tu_an_3)

summary(inter_tu_an_3)
summary(inter_tu_an_1)

simulationOutput <- simulateResiduals(fittedModel = inter_tu_an_1, plot = F, )
plot(simulationOutput) # no problem


inter_te_an_1<-glmmTMB(Te_inter~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="N"))
inter_te_an_2<-glmmTMB(Te_inter+1~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="N"), family=Gamma(link="log"))
inter_te_an_3<-glmmTMB(Te_inter+1~Te_Regime*Tu_Regime, data=subset(param_all_w0, Environment=="N"), family=gaussian(link="log"))
anova(inter_te_an_1,inter_te_an_2,inter_te_an_3)

summary(inter_te_an_2)
summary(inter_te_an_1)


simulationOutput <- simulateResiduals(fittedModel = inter_te_an_1, plot = F, )
plot(simulationOutput) # no problems

inter_te_an_1_2<-glmmTMB(Te_inter~Te_Regime+Tu_Regime, data=subset(param_all_w0, Environment=="N"))

summary(inter_te_an_1_2)

###### Summary

summary(gr_tu_an_2)

summary(gr_te_an_2)

summary(intra_tu_an_1)

summary(intra_te_an_1)


emmeans(inter_tu_an_1, pairwise~Te_Regime:Tu_Regime, adjust="none")

emmeans(inter_te_an_1, pairwise~Te_Regime:Tu_Regime, adjust="none")


Anova(gr_tu_an_2)
Anova(gr_te_an_2)

Anova(intra_tu_an_1)
Anova(intra_te_an_1)

Anova(inter_tu_an_1, type=3)
Anova(inter_te_an_1, type=3)

emmeans(inter_te_an_1, pairwise~Tu_Regime| Te_Regime, adjust="none")


# Bootstrap differences between selection regimes
#For each question we will randomize the replicates between selection regimes.

#### Does cadmium change parameters?
nboot<-1000

#Bootstrap to reestimate the p-value obtained for growth rate and intraspecific competition.
boot_tu_gr_intra_env<-as.data.frame(t(sapply(c(1:nboot),function(x){
  
  if(x%%10 ==0){
    print(x)
  }
  
  auxi<-subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" )
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  auxi$Tu_lambda<-auxi[rand_numb,"Tu_lambda"] # randomizing the trais
  auxi$Tu_intra<-auxi[rand_numb,"Tu_intra"]
  auxi$Tu_inter<-auxi[rand_numb,"Tu_inter"]
  
  gr <-glmmTMB(Tu_lambda~Environment, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Tu_intra~Environment, data=auxi)
  inter<-glmmTMB(Tu_inter~Environment, data=auxi)
  
  sum_auxi<-auxi %>% group_by(Environment)%>% summarize(meanGr=mean(Tu_lambda, na.rm=TRUE), meanIntra=mean(Tu_intra, na.rm=TRUE), meaninter=mean(Tu_inter, na.rm=TRUE)) %>% as.data.frame()
  
  # N - Cd
  diff<-sum_auxi[2,c(2:4)]-sum_auxi[1,(2:4)]
  
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  
  c(gr_p, intra_p, inter_p, diff[1,1], diff[1,2], diff[1,3])
  
} )))


boot_te_gr_intra_env<-as.data.frame(t(sapply(c(1:nboot),function(x){
  
  if(x%%10 ==0){
    print(x)
  }
  
  auxi<-subset(param_all_w0, Tu_Regime=="SR1" & Te_Regime=="SR4" )
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  auxi$Te_lambda<-auxi[rand_numb,"Te_lambda"] # randomizing the trais
  auxi$Te_intra<-auxi[rand_numb,"Te_intra"]
  auxi$Te_inter<-auxi[rand_numb,"Te_inter"]
  
  gr <-glmmTMB(Te_lambda~Environment, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Te_intra~Environment, data=auxi)
  inter <-glmmTMB(Te_inter~Environment, data=auxi)
  
  sum_auxi<-auxi %>% group_by(Environment)%>% summarize(meanGr=mean(Te_lambda, na.rm=TRUE), meanIntra=mean(Te_intra, na.rm=TRUE), meaninter=mean(Te_inter, na.rm=TRUE)) %>% as.data.frame()
  
  # N - Cd
  diff<-sum_auxi[2,c(2:4)]-sum_auxi[1,(2:4)]
  
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  
  c(gr_p, intra_p, inter_p,diff[1,1], diff[1,2], diff[1,3])
  
} )))
colnames(boot_tu_gr_intra_env)<-c("lambda_p","intra_p","inter_p", "lambda_diff", "intra_diff", "inter_diff")
colnames(boot_te_gr_intra_env)<-c("lambda_p","intra_p","inter_p", "lambda_diff", "intra_diff", "inter_diff")
str(boot_te_gr_intra_env)

ggplot(boot_tu_gr_intra_env, aes(x=lambda_p))+
  geom_histogram()+
  geom_vline(data=as.data.frame(Anova(gr_tu_cd_2)), aes_string(xintercept=as.data.frame(Anova(gr_tu_cd_2))[,3]))


print("Boot p-values for tests between environments")
length(which(boot_tu_gr_intra_env$lambda_p<=as.data.frame(Anova(gr_tu_cd_2))[,3]))/(nboot+1)

length(which(boot_tu_gr_intra_env$intra_p<=as.data.frame(Anova(intra_tu_cd_1))[,3]))/(nboot+1)

length(which(boot_tu_gr_intra_env$inter_p<=as.data.frame(Anova(inter_tu_cd_1))[,3]))/(nboot+1)

length(which(boot_te_gr_intra_env$lambda_p<=as.data.frame(Anova(gr_te_cd_2))[,3]))/(nboot+1)

length(which(boot_te_gr_intra_env$intra_p<=as.data.frame(Anova(intra_te_cd_1))[,3]))/(nboot+1)

length(which(boot_te_gr_intra_env$inter_p<=as.data.frame(Anova(inter_te_cd_1))[,3]))/(nboot+1)

as.data.frame(Anova(gr_tu_cd_2))[,3]
as.data.frame(Anova(intra_tu_cd_1))[,3]
as.data.frame(Anova(inter_tu_cd_1))[,3]
as.data.frame(Anova(gr_te_cd_2))[,3]
as.data.frame(Anova(intra_te_cd_1))[,3]
as.data.frame(Anova(inter_te_cd_1))[,3]


length(which(boot_tu_gr_intra_env$lambda_p<=0.05))/(nboot)

length(which(boot_tu_gr_intra_env$intra_p<=0.05))/(nboot)

length(which(boot_tu_gr_intra_env$inter_p<=0.05))/(nboot)

length(which(boot_te_gr_intra_env$lambda_p<=0.05))/(nboot)

length(which(boot_te_gr_intra_env$intra_p<=0.05))/(nboot)

length(which(boot_te_gr_intra_env$inter_p<=0.05))/(nboot)


####Does evolution change the performance in cadmium?

#Bootstrap to reestimate the p-value obtained for growth rate and intraspecific competition.
boot_tu_evolcd<-as.data.frame(t(sapply(c(1:nboot),function(x){
  
  if(x%%10 ==0){
    print(x)
  }
  
  auxi<-subset(param_all_w0, Environment=="Cd"& Te_Regime=="SR4")
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  auxi$Tu_lambda<-auxi[rand_numb,"Tu_lambda"] # randomizing the trais
  auxi$Tu_intra<-auxi[rand_numb,"Tu_intra"]
  
  auxi2<-subset(param_all_w0, Environment=="Cd")
  rand_numb2<-sample(c(1:dim(auxi2)[1]), dim(auxi2)[1], replace = TRUE)
  auxi2$Tu_inter<-auxi2[rand_numb2,"Tu_inter"]
  
  gr <-glmmTMB(Tu_lambda~Tu_Regime, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Tu_intra~Tu_Regime, data=auxi)
  inter<-glmmTMB(Tu_inter~Tu_Regime*Te_Regime, data=auxi2)
  
  #inter<-glmmTMB(Tu_inter~Environment, data=auxi)
  
  sum_auxi<-auxi %>% group_by(Tu_Regime)%>% summarize(meanGr=mean(Tu_lambda, na.rm=TRUE), meanIntra=mean(Tu_intra, na.rm=TRUE)) %>% as.data.frame()
  
  sum_auxi2<-auxi2 %>% group_by(Tu_Regime, Te_Regime)%>% summarize( meanInter=mean(Tu_inter, na.rm=TRUE)) %>% as.data.frame()
  
  # N - Cd
  diff<-sum_auxi[2,c(2:3)]-sum_auxi[1,(2:3)]
  diff2<-sum_auxi2[1,3]-sum_auxi2[2,3] # Evolution of the competitor with control focal
  diff3<-sum_auxi2[1,3]-sum_auxi2[3,3] # Evolution of the focal with control competitor
  diff4<-sum_auxi2[2,3]-sum_auxi2[4,3] # Evolution of the focal with evolved competitor
  diff5<-sum_auxi2[3,3]-sum_auxi2[4,3] # Evolution of the competitor with evolved focal
  
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  inter_p2<-as.data.frame(Anova(inter))[2,3]
  inter_p3<-as.data.frame(Anova(inter))[3,3]
  
  c(gr_p, intra_p, inter_p,inter_p2, inter_p3, diff[1,1], diff[1,2], diff2, diff3, diff4,diff5)
  
} )))


boot_te_evolcd<-as.data.frame(t(sapply(c(1:nboot),function(x){
  
  if(x%%10 ==0){
    print(x)
  }
  
  auxi<-subset(param_all_w0, Environment=="Cd"& Tu_Regime=="SR1")
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  auxi$Te_lambda<-auxi[rand_numb,"Te_lambda"] # randomizing the trais
  auxi$Te_intra<-auxi[rand_numb,"Te_intra"]
  
  # print(x)
  
  auxi2<-subset(param_all_w0, Environment=="Cd")
  rand_numb2<-sample(c(1:dim(auxi2)[1]), dim(auxi2)[1], replace = TRUE)
  auxi2$Te_inter<-auxi2[rand_numb2,"Te_inter"]
  
  gr <-glmmTMB(Te_lambda~Te_Regime, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Te_intra~Te_Regime, data=auxi)
  inter<-glmmTMB(Te_inter~Tu_Regime*Te_Regime, data=auxi2)
  
  #inter<-glmmTMB(Tu_inter~Environment, data=auxi)
  
  sum_auxi<-auxi %>% group_by(Te_Regime)%>% summarize(meanGr=mean(Te_lambda, na.rm=TRUE), meanIntra=mean(Te_intra, na.rm=TRUE)) %>% as.data.frame()
  
  sum_auxi2<-auxi2 %>% group_by(Tu_Regime, Te_Regime)%>% summarize( meanInter=mean(Te_inter, na.rm=TRUE)) %>% as.data.frame()
  
  # N - Cd
  diff<-sum_auxi[2,c(2:3)]-sum_auxi[1,(2:3)]
  diff2<-sum_auxi2[1,3]-sum_auxi2[2,3] # Evolution of the competitor with control focal
  diff3<-sum_auxi2[1,3]-sum_auxi2[3,3] # Evolution of the focal with control competitor
  diff4<-sum_auxi2[2,3]-sum_auxi2[4,3] # Evolution of the focal with evolved competitor
  diff5<-sum_auxi2[3,3]-sum_auxi2[4,3] # Evolution of the competitor with evolved focal
  
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  inter_p2<-as.data.frame(Anova(inter))[2,3]
  inter_p3<-as.data.frame(Anova(inter))[3,3]
  
  c(gr_p, intra_p, inter_p,inter_p2, inter_p3, diff[1,1], diff[1,2], diff2, diff3, diff4,diff5)
  
} )))


colnames(boot_tu_evolcd)<-c("lambda_p","intra_p","inter_p_TuReg","inter_p_TeReg","inter_p_int", "lambda_diff", "intra_diff", "inter_diffEvolComp_focalControl","inter_diffEvolFocal_CompControl","inter_diffEvolFocal_EvolComp","inter_diffEvolComp_focalEvol" )
colnames(boot_te_evolcd)<-c("lambda_p","intra_p","inter_p_TuReg","inter_p_TeReg","inter_p_int", "lambda_diff", "intra_diff", "inter_diffEvolComp_focalControl","inter_diffEvolFocal_CompControl","inter_diffEvolFocal_EvolComp","inter_diffEvolComp_focalEvol" )

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

hist(boot_tu_evolcd$lambda_p)
abline(v=as.data.frame(Anova(gr_tu_ev_2))[,3], col="red")

hist(boot_tu_evolcd$intra_p)
abline(v=as.data.frame(Anova(intra_tu_ev_1))[,3], col="red")

hist(boot_tu_evolcd$inter_p_TuReg)
abline(v=as.data.frame(Anova(inter_tu_ev_1))[1,3], col="red")

hist(boot_tu_evolcd$inter_p_TeReg)
abline(v=as.data.frame(Anova(inter_tu_ev_1))[2,3], col="red")

hist(boot_tu_evolcd$inter_p_int)
abline(v=as.data.frame(Anova(inter_tu_ev_1))[3,3], col="red")


hist(boot_te_evolcd$lambda_p)
abline(v=as.data.frame(Anova(gr_te_ev_2))[,3], col="red")

hist(boot_te_evolcd$intra_p)
abline(v=as.data.frame(Anova(intra_te_ev_1))[,3], col="red")

hist(boot_te_evolcd$inter_p_TuReg)
abline(v=as.data.frame(Anova(inter_te_ev_1))[1,3], col="red")

hist(boot_te_evolcd$inter_p_TeReg)
abline(v=as.data.frame(Anova(inter_te_ev_1))[2,3], col="red")

hist(boot_te_evolcd$inter_p_int)
abline(v=as.data.frame(Anova(inter_te_ev_1))[3,3], col="red")


as.data.frame(Anova(gr_tu_ev_2))[,3]
as.data.frame(Anova(intra_tu_ev_1))[,3]
as.data.frame(Anova(inter_tu_ev_1))[1,3]
as.data.frame(Anova(inter_tu_ev_1))[2,3]
as.data.frame(Anova(inter_tu_ev_1))[3,3]
as.data.frame(Anova(gr_te_ev_2))[,3]
as.data.frame(Anova(intra_te_ev_1))[,3]
as.data.frame(Anova(inter_te_ev_1))[1,3]
as.data.frame(Anova(inter_te_ev_1))[2,3]
as.data.frame(Anova(inter_te_ev_1))[3,3]


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

#### Does evolution change the ancestral?
#Bootstrap to reestimate the p-value obtained for growth rate and intraspecific competition.
boot_tu_evolN<-as.data.frame(t(sapply(c(1:nboot),function(x){
  
  if(x%%10 ==0){
    print(x)
  }
  
  auxi<-subset(param_all_w0, Environment=="N"& Te_Regime=="SR4")
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  auxi$Tu_lambda<-auxi[rand_numb,"Tu_lambda"] # randomizing the trais
  auxi$Tu_intra<-auxi[rand_numb,"Tu_intra"]
  
  auxi2<-subset(param_all_w0, Environment=="N")
  rand_numb2<-sample(c(1:dim(auxi2)[1]), dim(auxi2)[1], replace = TRUE)
  auxi2$Tu_inter<-auxi2[rand_numb2,"Tu_inter"]
  
  gr <-glmmTMB(Tu_lambda~Tu_Regime, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Tu_intra~Tu_Regime, data=auxi)
  inter<-glmmTMB(Tu_inter~Tu_Regime*Te_Regime, data=auxi2)
  
  #inter<-glmmTMB(Tu_inter~Environment, data=auxi)
  
  sum_auxi<-auxi %>% group_by(Tu_Regime)%>% summarize(meanGr=mean(Tu_lambda, na.rm=TRUE), meanIntra=mean(Tu_intra, na.rm=TRUE)) %>% as.data.frame()
  
  sum_auxi2<-auxi2 %>% group_by(Tu_Regime, Te_Regime)%>% summarize( meanInter=mean(Tu_inter, na.rm=TRUE)) %>% as.data.frame()
  
  # N - N
  diff<-sum_auxi[2,c(2:3)]-sum_auxi[1,(2:3)]
  diff2<-sum_auxi2[1,3]-sum_auxi2[2,3] # Evolution of the competitor with control focal
  diff3<-sum_auxi2[1,3]-sum_auxi2[3,3] # Evolution of the focal with control competitor
  diff4<-sum_auxi2[2,3]-sum_auxi2[4,3] # Evolution of the focal with evolved competitor
  diff5<-sum_auxi2[3,3]-sum_auxi2[4,3] # Evolution of the competitor with evolved focal
  
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  inter_p2<-as.data.frame(Anova(inter))[2,3]
  inter_p3<-as.data.frame(Anova(inter))[3,3]
  
  c(gr_p, intra_p, inter_p,inter_p2, inter_p3, diff[1,1], diff[1,2], diff2, diff3, diff4,diff5)
  
} )))


boot_te_evolN<-as.data.frame(t(sapply(c(1:nboot),function(x){
  
  if(x%%10 ==0){
    print(x)
  }
  
  auxi<-subset(param_all_w0, Environment=="N"& Tu_Regime=="SR1")
  rand_numb<-sample(c(1:dim(auxi)[1]), dim(auxi)[1], replace = TRUE)
  auxi$Te_lambda<-auxi[rand_numb,"Te_lambda"] # randomizing the trais
  auxi$Te_intra<-auxi[rand_numb,"Te_intra"]
  
  # print(x)
  
  auxi2<-subset(param_all_w0, Environment=="N")
  rand_numb2<-sample(c(1:dim(auxi2)[1]), dim(auxi2)[1], replace = TRUE)
  auxi2$Te_inter<-auxi2[rand_numb2,"Te_inter"]
  
  gr <-glmmTMB(Te_lambda~Te_Regime, data=auxi, family=Gamma(link="log"))
  intra <-glmmTMB(Te_intra~Te_Regime, data=auxi)
  inter<-glmmTMB(Te_inter~Tu_Regime*Te_Regime, data=auxi2)
  
  #inter<-glmmTMB(Tu_inter~Environment, data=auxi)
  
  sum_auxi<-auxi %>% group_by(Te_Regime)%>% summarize(meanGr=mean(Te_lambda, na.rm=TRUE), meanIntra=mean(Te_intra, na.rm=TRUE)) %>% as.data.frame()
  
  sum_auxi2<-auxi2 %>% group_by(Tu_Regime, Te_Regime)%>% summarize( meanInter=mean(Te_inter, na.rm=TRUE)) %>% as.data.frame()
  
  # N - N
  diff<-sum_auxi[2,c(2:3)]-sum_auxi[1,(2:3)]
  diff2<-sum_auxi2[1,3]-sum_auxi2[2,3] # Evolution of the competitor with control focal
  diff3<-sum_auxi2[1,3]-sum_auxi2[3,3] # Evolution of the focal with control competitor
  diff4<-sum_auxi2[2,3]-sum_auxi2[4,3] # Evolution of the focal with evolved competitor
  diff5<-sum_auxi2[3,3]-sum_auxi2[4,3] # Evolution of the competitor with evolved focal
  
  gr_p<-as.data.frame(Anova(gr))[1,3]
  intra_p<-as.data.frame(Anova(intra))[1,3]
  inter_p<-as.data.frame(Anova(inter))[1,3]
  inter_p2<-as.data.frame(Anova(inter))[2,3]
  inter_p3<-as.data.frame(Anova(inter))[3,3]
  
  c(gr_p, intra_p, inter_p,inter_p2, inter_p3, diff[1,1], diff[1,2], diff2, diff3, diff4,diff5)
  
} )))


colnames(boot_tu_evolN)<-c("lambda_p","intra_p","inter_p_TuReg","inter_p_TeReg","inter_p_int", "lambda_diff", "intra_diff", "inter_diffEvolComp_focalControl","inter_diffEvolFocal_CompControl","inter_diffEvolFocal_EvolComp","inter_diffEvolComp_focalEvol" )
colnames(boot_te_evolN)<-c("lambda_p","intra_p","inter_p_TuReg","inter_p_TeReg","inter_p_int", "lambda_diff", "intra_diff", "inter_diffEvolComp_focalControl","inter_diffEvolFocal_CompControl","inter_diffEvolFocal_EvolComp","inter_diffEvolComp_focalEvol" )

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

hist(boot_tu_evolN$lambda_p)
abline(v=as.data.frame(Anova(gr_tu_an_2))[,3], col="red")

hist(boot_tu_evolN$intra_p)
abline(v=as.data.frame(Anova(intra_tu_an_1))[,3], col="red")

hist(boot_tu_evolN$inter_p_TuReg)
abline(v=as.data.frame(Anova(inter_tu_an_1))[1,3], col="red")

hist(boot_tu_evolN$inter_p_TeReg)
abline(v=as.data.frame(Anova(inter_tu_an_1))[2,3], col="red")

hist(boot_tu_evolN$inter_p_int)
abline(v=as.data.frame(Anova(inter_tu_an_1))[3,3], col="red")


hist(boot_te_evolN$lambda_p)
abline(v=as.data.frame(Anova(gr_te_an_2))[,3], col="red")

hist(boot_te_evolN$intra_p)
abline(v=as.data.frame(Anova(intra_te_an_1))[,3], col="red")

hist(boot_te_evolN$inter_p_TuReg)
abline(v=as.data.frame(Anova(inter_te_an_1))[1,3], col="red")

hist(boot_te_evolN$inter_p_TeReg)
abline(v=as.data.frame(Anova(inter_te_an_1))[2,3], col="red")

hist(boot_te_evolN$inter_p_int)
abline(v=as.data.frame(Anova(inter_te_an_1))[3,3], col="red")

as.data.frame(Anova(gr_tu_an_2))[,3]
as.data.frame(Anova(intra_tu_an_1))[,3]
as.data.frame(Anova(inter_tu_an_1))[1,3]
as.data.frame(Anova(inter_tu_an_1))[2,3]
as.data.frame(Anova(inter_tu_an_1))[3,3]
as.data.frame(Anova(gr_te_an_2))[,3]
as.data.frame(Anova(intra_te_an_1))[,3]
as.data.frame(Anova(inter_te_an_1))[1,3]
as.data.frame(Anova(inter_te_an_1))[2,3]
as.data.frame(Anova(inter_te_an_1))[3,3]


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
