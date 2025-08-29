rm(list=ls())
# Packages and functions
library(ggplot2)
library(plyr)
library(dplyr)
library(car)
library(fitdistrplus)
library(tidyr)
library(tidyverse)
library(ggtext)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(ggbreak)
library(nlme)
library(cxr)
library(MASS)
library(mvtnorm)
library(DescTools)
library(phia)
library(performance)
library(DHARMa)
library(effects)
library(cowplot)

theme_plot<-theme(axis.text = element_text(size=14), axis.title = element_text(size=14, face="bold"), legend.text = element_text(size=12), strip.text = element_text(size=14), plot.title = element_text(size=14, face="bold"), panel.grid=element_line(colour="white"), panel.background = element_rect(fill="white") , axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"), strip.background = element_rect(fill="white"))

save_plot<-function(dir, width=15, height=10, ...){
  ggsave(dir, width = width, height = height, units = c("cm"))
}

Env<-c("Water", "Cadmium")
names(Env)<-c("N", "Cd")

## Warning! This code takes a very long time to run!
# 1 - Importing data and checking it

coex<-read.csv("./Data/Coexistence_Cd_G42_submit.csv", header=TRUE) # Data from the coexistence experiment
ca_raw<-read.csv(file = "./Data/CompetitiveAbility_Cd_G40_submit.csv", header=TRUE) # cdata from the competitive ability

str(ca_raw) 
# Summary of the data to be sure that everything is ok!
summary(as.factor(ca_raw$Foca_rawlSR))

ca_raw$Block2<-as.factor(ca_raw$Block)
ca_raw$Rep2<-as.factor(ca_raw$Rep)
ca_raw$Disk2<-as.factor(ca_raw$Disk)
ca_raw$Leaf2<-as.factor(ca_raw$Leaf)
ca_raw$Env2<-as.factor(ca_raw$Env)
ca_raw$FocalSR2<-as.factor(ca_raw$FocalSR)
ca_raw$CompSR2<-as.factor(ca_raw$CompSR)
ca_raw$Type2<-as.factor(ca_raw$Type)
ca_raw$Focal_Female2<-as.factor(ca_raw$Focalfemale)


regimeTu<-c("Tu \ncontrol", "Tu evolved \n in cadmium")
names(regimeTu)<-c("SR1", "SR2")

regimeTe<-c("Te \n control", "Te evolved \n in cadmium")
names(regimeTe)<-c("SR4", "SR5")

#Creating columns that are needed
ca_raw$Nr_Focal_Females_Tu_Alive_G0<-sapply(c(1:length(ca_raw$Block)), function(x){
  if(ca_raw$Focalfemale[x]=="Tu"){
    if(ca_raw$Type[x]=="INTRA"){
      a<-ca_raw$Dens[x]-ca_raw$FocalDead[x]-ca_raw$FocalDrowned[x]-ca_raw$FocalMissing[x]
    }else
      a<-1-ca_raw$FocalDead[x]-ca_raw$FocalDrowned[x]-ca_raw$FocalMissing[x]
    
  }else
    a<-NA
})

ca_raw$Nr_Focal_Females_Te_Alive_G0<-sapply(c(1:length(ca_raw$Block)), function(x){
  if(ca_raw$Focalfemale[x]=="Te"){
    if(ca_raw$Type[x]=="INTRA"){
      a<-ca_raw$Dens[x]-ca_raw$FocalDead[x]-ca_raw$FocalDrowned[x]-ca_raw$FocalMissing[x]
    }else
      a<-1-ca_raw$FocalDead[x]-ca_raw$FocalDrowned[x]-ca_raw$FocalMissing[x]
    
  }else
    a<-NA
})


ca_raw$Num_Comp_Tu_Alive_G0<-sapply(c(1:length(ca_raw$Block)), function(x){
  if(ca_raw$Focalfemale[x]=="Te"){
    if(ca_raw$Type[x]=="INTER"){
      a<-ca_raw$Dens[x]-ca_raw$NumbDeadComp[x]-1
    }else
      a<-NA
    
  }else
    a<-NA
})


ca_raw$Num_Comp_Te_Alive_G0<-sapply(c(1:length(ca_raw$Block)), function(x){
  if(ca_raw$Focalfemale[x]=="Tu"){
    if(ca_raw$Type[x]=="INTER"){
      a<-ca_raw$Dens[x]-ca_raw$NumbDeadComp[x]-1
    }else
      a<-NA
    
  }else
    a<-NA
})

ca_raw$Nr_Focal_Females_G0<-sapply(c(1:length(ca_raw$Block)), function(x){
  if(ca_raw$Type[x]=="INTRA"){
    a<-ca_raw$Dens[x]
  }else
    a<-1
  
})

ca_raw$Nr_Comp_Females_G0<-sapply(c(1:length(ca_raw$Block)), function(x){
  a<-ca_raw$Dens[x]-1
  
  a
  
})

# Removing rows where there were less than 0 females
ca_raw<-ca_raw[-c(which(ca_raw$Num_Comp_Te_Alive_G0<0),which(ca_raw$Num_Comp_Tu_Alive_G0<0), which(ca_raw$Nr_Focal_Females_Te_Alive_G0<0),which(ca_raw$Nr_Focal_Females_Tu_Alive_G0<0) ),]


# # Removing virgin females
# ca<-ca_raw[-c(which(ca_raw$TeFemales==0 &ca_raw$TeMales>0 & ca_raw$Focalfemale=="Te" ),which(ca_raw$TuFemales==0 &ca_raw$TuMales>0 & ca_raw$Focalfemale=="Tu" )),]

# 2 - Estimate growth rate per generation

ca_raw$GrowthRateOA<-sapply(c(1:length(ca_raw[,1])), function(x){
  #print(x)
  if(ca_raw$Focal_Female[x]=="Tu"){
    a<-ca_raw$TuFemales[x]/ca_raw$Nr_Focal_Females_G0[x]
  }else if(ca_raw$Focal_Female[x]=="Te"){
    a<-ca_raw$TeFemales[x]/ca_raw$Nr_Focal_Females_G0[x]
  }else
    a<-NA
  
  a
})

ca<-ca_raw

# 3 - Estimate competitive ability & predict data

#Here we have two differeny methods, using the cxr package or with the optim. We will also vary if we estimate lambda from the data or from the model and if using cxr with the nested approach is better or not. So the different hypothesis are

#A - CXR normal: using cxr with the normal approach
#B - CXR lambda fixed: using cxr but lambda comes from the data
#C - CXR nested: lambda comes the data, and we use the same nested approach as the optim --> for that we can put intra as another species (column)
#D - optim normal: the same approach as used in Fragata 2022
#E - optim lambda fixed: using optim, but lambda is fixed 

#In all the models we will use density -1 for the intra, which basically corresponds to the number of competitors.


### A - CXR normal

#cxr accepts a data frame with a first column called fitness with positive values and numeric columns with number of individuals. Each row is one individual. For multiple species the easier is to create a list, each with a data frame that has in the first column number of individuals produced and then the number of neighbours
#this case we transformed all 0s into 1 (so that the log is 0) For that we need to add +1 to all data so that the variance is not changed

##### normal

dir.create("./Analyses/MethodComparison/cxr_normal", showWarnings = FALSE)

# modifying data frame to fit the type of setup that is need for CXR
forCXR_N<-subset(ca, Env=="N")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "TuFemales")]

forCXR_N$Focal<-mapvalues(forCXR_N$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
forCXR_N$CompSR2<-mapvalues(forCXR_N$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

forCXR_N$Comp<-sapply(c(1:length(forCXR_N[,1])), function(x){
  if(is.na(forCXR_N$CompSR2[x])){
    a<- forCXR_N$Focal[x]
  }else{
    a<-forCXR_N$CompSR2[x]
  }
  
  a
})

aux<-data.frame(SR1=rep(0, length(forCXR_N[,1])), SR2=rep(0, length(forCXR_N[,1])), SR4=rep(0, length(forCXR_N[,1])), SR5=rep(0, length(forCXR_N[,1])))

for(i in 1:length(forCXR_N[,1])){
  #coluna onde por focais
  colunaF<-which(colnames(aux)==forCXR_N$Focal[i])
  #coluna onde por competidors
  colunaC<-which(colnames(aux)==forCXR_N$Comp[i])
  
  #if its the same regime
  if(forCXR_N$Focal[i]==forCXR_N$Comp[i] & forCXR_N$Dens[i]==1){
    aux[i,colunaF]<-forCXR_N$Dens[i]-1
    
  }else if(forCXR_N$Focal[i]==forCXR_N$Comp[i]){
    aux[i,colunaF]<-forCXR_N$Dens[i]-1
  }else{ #if it is heterospecific then its -1 for the competitors (because of the focal) and its one for the focal
    aux[i,colunaC]<-forCXR_N$Dens[i]-1
    aux[i, colunaF]<-1
  }
  
}

forCXR_N<-cbind(forCXR_N, aux)

forCXR_N$fitness<-sapply(c(1:length(forCXR_N[,1])), function(x){
  colF<-which(colnames(forCXR_N)==forCXR_N$Focal[x])
  
  if(forCXR_N$Focal[x]=="SR1"){
    a<-forCXR_N$TuFemales[x]/forCXR_N$SR1[x]
  } else if(forCXR_N$Focal[x]=="SR2"){
    a<-forCXR_N$TuFemales[x]/forCXR_N$SR2[x]
  } else if(forCXR_N$Focal[x]=="SR4"){
    a<-forCXR_N$TeFemales[x]/forCXR_N$SR4[x]
  } else if(forCXR_N$Focal[x]=="SR5"){
    a<-forCXR_N$TeFemales[x]/forCXR_N$SR5[x]
  }
  
  a
})

#removing rows for which there is no data for fitness
forCXR_N<-forCXR_N[-which(is.na(forCXR_N$fitness)),]

# adding +1 to all data
#forCXR_N$fitness<-forCXR_N$fitness+1

forCXR_N[which(forCXR_N$fitness=="-Inf" | forCXR_N$fitness=="Inf"),"fitness"]<-0


# all data gets +1 because of the 0 problem
forCXR_N$fitness<-forCXR_N$fitness+1

# vector that tells which are the selection regimes, the columns have to have the same name
my.reg <- c("SR1", "SR2","SR4","SR5")

# Do list per replicate and environment
R1<-list(SR1= subset(forCXR_N, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R2<-list(SR1= subset(forCXR_N, Rep==2 & Focal=="SR1")[,c("fitness", "SR1", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR4", "SR5")])

R3<-list(SR1= subset(forCXR_N, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R4<-list(SR1= subset(forCXR_N, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R5<-list(SR1= subset(forCXR_N, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])


obs.R1_w0<-cxr_pm_multifit(data = R1,
                           focal_column = my.reg,
                           model_family = "RK",
                           covariates = NULL,
                           optimization_method = "Nelder-Mead",
                           alpha_form = "pairwise",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none",
                           initial_values = list(lambda = 1,
                                                 alpha_intra = 0.1,
                                                 alpha_inter = 0.1),
                           fixed_terms = NULL,
                           # no standard errors
                           bootstrap_samples = 200)

str(obs.R1_w0)
obs.R1_w0$lambda_standard_error
obs.R1_w0$alpha_matrix_standard_error



obs.R3_w0<-cxr_pm_multifit(data = R3,
                           focal_column = my.reg,
                           model_family = "RK",
                           covariates = NULL,
                           optimization_method = "Nelder-Mead",
                           alpha_form = "pairwise",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none",
                           initial_values = list(lambda = 1,
                                                 alpha_intra = 0.1,
                                                 alpha_inter = 0.1),
                           fixed_terms = NULL,
                           # no standard errors
                           bootstrap_samples = 200)

obs.R4_w0<-cxr_pm_multifit(data = R4,
                           focal_column = my.reg,
                           model_family = "RK",
                           covariates = NULL,
                           optimization_method = "Nelder-Mead",
                           alpha_form = "pairwise",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none",
                           initial_values = list(lambda = 1,
                                                 alpha_intra = 0.1,
                                                 alpha_inter = 0.1),
                           fixed_terms = NULL,
                           # no standard errors
                           bootstrap_samples = 200)

obs.R5_w0<-cxr_pm_multifit(data = R5,
                           focal_column = my.reg,
                           model_family = "RK",
                           covariates = NULL,
                           optimization_method = "Nelder-Mead",
                           alpha_form = "pairwise",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none",
                           initial_values = list(lambda = 1,
                                                 alpha_intra = 0.1,
                                                 alpha_inter = 0.1),
                           fixed_terms = NULL,
                           # no standard errors
                           bootstrap_samples = 200)

summary(obs.R1_w0)
summary(obs.R3_w0)
summary(obs.R4_w0)
summary(obs.R5_w0)


# For replicate 2 we need to do it differently


obs.R2_w0_sr1<-cxr_pm_fit(data = R2[[1]],
                          focal_column = my.reg[1],
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(lambda = 1,
                                                alpha_intra = 0.1,
                                                alpha_inter = 0.1),
                          fixed_terms = NULL,
                          # no standard errors
                          bootstrap_samples = 200)

obs.R2_w0_sr4<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]==0),c("fitness", "SR4")],
                          focal_column =NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "global",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(lambda = 1,
                                                alpha_inter = 0.1),
                          fixed_terms = NULL,
                          # no standard errors
                          bootstrap_samples = 200)

obs.R2_w0_sr4_inter<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]!=0),c("fitness", "SR4")],
                                focal_column =NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = 0.1),
                                fixed_terms = list(lambda=obs.R2_w0_sr4$lambda),
                                # no standard errors
                                bootstrap_samples = 200)

obs.R2_w0_sr5<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]==0),c("fitness", "SR5")],
                          focal_column =NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "global",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(lambda = 1,
                                                alpha_inter = 0.1),
                          fixed_terms = NULL,
                          # no standard errors
                          bootstrap_samples = 200)

obs.R2_w0_sr5_inter<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]!=0),c("fitness", "SR5")],
                                focal_column =NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = 0.1),
                                fixed_terms = list(lambda=obs.R2_w0_sr5$lambda),
                                # no standard errors
                                bootstrap_samples = 200)


#rows in the alpha element of the returning list correspond to species i and columns to species j for each αij coefficient.

###### data table summary
cxr_param_w0<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_w0$Tu_lambda<-0
cxr_param_w0$Te_lambda<-0
cxr_param_w0$Tu_intra<-0
cxr_param_w0$Te_intra<-0
cxr_param_w0$Tu_inter<-0
cxr_param_w0$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0<-cxr_param_w0[-which(cxr_param_w0$Replicate==2 & cxr_param_w0$Tu_Regime=="SR2"),]


cxr_param_w0[which(cxr_param_w0$Replicate==1),"Tu_lambda"]<-obs.R1_w0$lambda[1:2]
cxr_param_w0[which(cxr_param_w0$Replicate==1),"Te_lambda"]<-obs.R1_w0$lambda[c(3,3,4,4)]

cxr_param_w0[which(cxr_param_w0$Replicate==2),"Tu_lambda"]<-obs.R2_w0_sr1$lambda
cxr_param_w0[which(cxr_param_w0$Replicate==2),"Te_lambda"]<-c(obs.R2_w0_sr4$lambda,obs.R2_w0_sr5$lambda)

cxr_param_w0[which(cxr_param_w0$Replicate==3),"Tu_lambda"]<-obs.R3_w0$lambda[1:2]
cxr_param_w0[which(cxr_param_w0$Replicate==3),"Te_lambda"]<-obs.R3_w0$lambda[c(3,3,4,4)]

cxr_param_w0[which(cxr_param_w0$Replicate==4),"Tu_lambda"]<-obs.R4_w0$lambda[1:2]
cxr_param_w0[which(cxr_param_w0$Replicate==4),"Te_lambda"]<-obs.R4_w0$lambda[c(3,3,4,4)]

cxr_param_w0[which(cxr_param_w0$Replicate==5),"Tu_lambda"]<-obs.R5_w0$lambda[1:2]
cxr_param_w0[which(cxr_param_w0$Replicate==5),"Te_lambda"]<-obs.R5_w0$lambda[c(3,3,4,4)]


cxr_param_w0[which(cxr_param_w0$Replicate==1),"Tu_intra"]<-rep(c(obs.R1_w0$alpha_matrix[1,1], obs.R1_w0$alpha_matrix[2,2]), 2)
cxr_param_w0[which(cxr_param_w0$Replicate==1),"Te_intra"]<-rep(c(obs.R1_w0$alpha_matrix[3,3], obs.R1_w0$alpha_matrix[4,4]), each=2)

cxr_param_w0[which(cxr_param_w0$Replicate==2),"Tu_intra"]<-obs.R2_w0_sr1$alpha_intra
cxr_param_w0[which(cxr_param_w0$Replicate==2),"Te_intra"]<-c(obs.R2_w0_sr4$alpha_inter, obs.R2_w0_sr5$alpha_inter)

cxr_param_w0[which(cxr_param_w0$Replicate==3),"Tu_intra"]<-rep(c(obs.R3_w0$alpha_matrix[1,1], obs.R3_w0$alpha_matrix[2,2]), 2)
cxr_param_w0[which(cxr_param_w0$Replicate==3),"Te_intra"]<-rep(c(obs.R3_w0$alpha_matrix[3,3], obs.R3_w0$alpha_matrix[4,4]), each=2)

cxr_param_w0[which(cxr_param_w0$Replicate==4),"Tu_intra"]<-rep(c(obs.R4_w0$alpha_matrix[1,1], obs.R4_w0$alpha_matrix[2,2]), 2)
cxr_param_w0[which(cxr_param_w0$Replicate==4),"Te_intra"]<-rep(c(obs.R4_w0$alpha_matrix[3,3], obs.R4_w0$alpha_matrix[4,4]), each=2)

cxr_param_w0[which(cxr_param_w0$Replicate==5),"Tu_intra"]<-rep(c(obs.R5_w0$alpha_matrix[1,1], obs.R5_w0$alpha_matrix[2,2]), 2)
cxr_param_w0[which(cxr_param_w0$Replicate==5),"Te_intra"]<-rep(c(obs.R5_w0$alpha_matrix[3,3], obs.R5_w0$alpha_matrix[4,4]), each=2)


cxr_param_w0[which(cxr_param_w0$Replicate==1),"Tu_inter"]<-c(obs.R1_w0$alpha_matrix[1,3], obs.R1_w0$alpha_matrix[2,3],obs.R1_w0$alpha_matrix[1,4], obs.R1_w0$alpha_matrix[2,4])
cxr_param_w0[which(cxr_param_w0$Replicate==1),"Te_inter"]<-c(obs.R1_w0$alpha_matrix[3,1], obs.R1_w0$alpha_matrix[3,2],obs.R1_w0$alpha_matrix[4,1], obs.R1_w0$alpha_matrix[4,2])

cxr_param_w0[which(cxr_param_w0$Replicate==2),"Tu_inter"]<-obs.R2_w0_sr1$alpha_inter
cxr_param_w0[which(cxr_param_w0$Replicate==2),"Te_inter"]<-c(obs.R2_w0_sr4_inter$alpha_inter, obs.R2_w0_sr5_inter$alpha_inter)

cxr_param_w0[which(cxr_param_w0$Replicate==3),"Tu_inter"]<-c(obs.R3_w0$alpha_matrix[1,3], obs.R3_w0$alpha_matrix[2,3],obs.R3_w0$alpha_matrix[1,4], obs.R3_w0$alpha_matrix[2,4])
cxr_param_w0[which(cxr_param_w0$Replicate==3),"Te_inter"]<-c(obs.R3_w0$alpha_matrix[3,1], obs.R3_w0$alpha_matrix[3,2],obs.R3_w0$alpha_matrix[4,1], obs.R3_w0$alpha_matrix[4,2])

cxr_param_w0[which(cxr_param_w0$Replicate==4),"Tu_inter"]<-c(obs.R4_w0$alpha_matrix[1,3], obs.R4_w0$alpha_matrix[2,3],obs.R4_w0$alpha_matrix[1,4], obs.R4_w0$alpha_matrix[2,4])
cxr_param_w0[which(cxr_param_w0$Replicate==4),"Te_inter"]<-c(obs.R4_w0$alpha_matrix[3,1], obs.R4_w0$alpha_matrix[3,2],obs.R4_w0$alpha_matrix[4,1], obs.R4_w0$alpha_matrix[4,2])

cxr_param_w0[which(cxr_param_w0$Replicate==5),"Tu_inter"]<-c(obs.R5_w0$alpha_matrix[1,3], obs.R5_w0$alpha_matrix[2,3],obs.R5_w0$alpha_matrix[1,4], obs.R5_w0$alpha_matrix[2,4])
cxr_param_w0[which(cxr_param_w0$Replicate==5),"Te_inter"]<-c(obs.R5_w0$alpha_matrix[3,1], obs.R5_w0$alpha_matrix[3,2],obs.R5_w0$alpha_matrix[4,1], obs.R5_w0$alpha_matrix[4,2])

### Lower

cxr_param_w0_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_w0_lower$Tu_lambda<-0
cxr_param_w0_lower$Te_lambda<-0
cxr_param_w0_lower$Tu_intra<-0
cxr_param_w0_lower$Te_intra<-0
cxr_param_w0_lower$Tu_inter<-0
cxr_param_w0_lower$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0_lower<-cxr_param_w0_lower[-which(cxr_param_w0_lower$Replicate==2 & cxr_param_w0_lower$Tu_Regime=="SR2"),]


cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==1),"Tu_lambda"]<-rep(c(obs.R1_w0$lambda[1]-obs.R1_w0$lambda_standard_error[1], obs.R1_w0$lambda[2]-obs.R1_w0$lambda_standard_error[2]), 2)
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==1),"Te_lambda"]<-rep(c(obs.R1_w0$lambda[3]-obs.R1_w0$lambda_standard_error[3], obs.R1_w0$lambda[4]-obs.R1_w0$lambda_standard_error[4]), each=2)

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==2),"Tu_lambda"]<-obs.R2_w0_sr1$lambda-obs.R2_w0_sr1$lambda_standard_error
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==2),"Te_lambda"]<-c(obs.R2_w0_sr4$lambda-obs.R2_w0_sr4$lambda_standard_error,obs.R2_w0_sr5$lambda-obs.R2_w0_sr5$lambda_standard_error)

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==3),"Tu_lambda"]<-rep(c(obs.R3_w0$lambda[1]-obs.R3_w0$lambda_standard_error[1], obs.R3_w0$lambda[2]-obs.R3_w0$lambda_standard_error[2]), 2)
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==3),"Te_lambda"]<-rep(c(obs.R3_w0$lambda[3]-obs.R3_w0$lambda_standard_error[3], obs.R3_w0$lambda[4]-obs.R3_w0$lambda_standard_error[4]), each=2)

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==4),"Tu_lambda"]<-rep(c(obs.R4_w0$lambda[1]-obs.R4_w0$lambda_standard_error[1], obs.R4_w0$lambda[2]-obs.R4_w0$lambda_standard_error[2]), 2)
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==4),"Te_lambda"]<-rep(c(obs.R4_w0$lambda[3]-obs.R4_w0$lambda_standard_error[3], obs.R4_w0$lambda[4]-obs.R4_w0$lambda_standard_error[4]), each=2)

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==5),"Tu_lambda"]<-rep(c(obs.R5_w0$lambda[1]-obs.R5_w0$lambda_standard_error[1], obs.R5_w0$lambda[2]-obs.R5_w0$lambda_standard_error[2]), 2)
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==5),"Te_lambda"]<-rep(c(obs.R5_w0$lambda[3]-obs.R5_w0$lambda_standard_error[3], obs.R5_w0$lambda[4]-obs.R5_w0$lambda_standard_error[4]), each=2)


cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==1),"Tu_intra"]<-rep(c(obs.R1_w0$alpha_matrix[1,1]-obs.R1_w0$alpha_matrix_standard_error[1,1], obs.R1_w0$alpha_matrix[2,2]-obs.R1_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==1),"Te_intra"]<-rep(c(obs.R1_w0$alpha_matrix[3,3]-obs.R1_w0$alpha_matrix_standard_error[3,3], obs.R1_w0$alpha_matrix[4,4]-obs.R1_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==2),"Tu_intra"]<-obs.R2_w0_sr1$alpha_intra-obs.R2_w0_sr1$alpha_intra_standard_error
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==2),"Te_intra"]<-c(obs.R2_w0_sr4$alpha_inter-obs.R2_w0_sr4$alpha_inter_standard_error, obs.R2_w0_sr5$alpha_inter-obs.R2_w0_sr5$alpha_inter_standard_error)

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==3),"Tu_intra"]<-rep(c(obs.R3_w0$alpha_matrix[1,1]-obs.R3_w0$alpha_matrix_standard_error[1,1], obs.R3_w0$alpha_matrix[2,2]-obs.R3_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==3),"Te_intra"]<-rep(c(obs.R3_w0$alpha_matrix[3,3]-obs.R3_w0$alpha_matrix_standard_error[3,3], obs.R3_w0$alpha_matrix[4,4]-obs.R3_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==4),"Tu_intra"]<-rep(c(obs.R4_w0$alpha_matrix[1,1]-obs.R4_w0$alpha_matrix_standard_error[1,1], obs.R4_w0$alpha_matrix[2,2]-obs.R4_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==4),"Te_intra"]<-rep(c(obs.R4_w0$alpha_matrix[3,3]-obs.R4_w0$alpha_matrix_standard_error[3,3], obs.R4_w0$alpha_matrix[4,4]-obs.R4_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==5),"Tu_intra"]<-rep(c(obs.R5_w0$alpha_matrix[1,1]-obs.R5_w0$alpha_matrix_standard_error[1,1], obs.R5_w0$alpha_matrix[2,2]-obs.R5_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==5),"Te_intra"]<-rep(c(obs.R5_w0$alpha_matrix[3,3]-obs.R5_w0$alpha_matrix_standard_error[3,3], obs.R5_w0$alpha_matrix[4,4]-obs.R5_w0$alpha_matrix_standard_error[4,4]), each=2)


cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==1),"Tu_inter"]<-c(obs.R1_w0$alpha_matrix[1,3]-obs.R1_w0$alpha_matrix_standard_error[1,3], obs.R1_w0$alpha_matrix[2,3]-obs.R1_w0$alpha_matrix_standard_error[2,3],obs.R1_w0$alpha_matrix[1,4]-obs.R1_w0$alpha_matrix_standard_error[1,4], obs.R1_w0$alpha_matrix[2,4]-obs.R1_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==1),"Te_inter"]<-c(obs.R1_w0$alpha_matrix[3,1]-obs.R1_w0$alpha_matrix_standard_error[3,1], obs.R1_w0$alpha_matrix[3,2]-obs.R1_w0$alpha_matrix_standard_error[3,2],obs.R1_w0$alpha_matrix[4,1]-obs.R1_w0$alpha_matrix_standard_error[4,1], obs.R1_w0$alpha_matrix[4,2]-obs.R1_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==2),"Tu_inter"]<-obs.R2_w0_sr1$alpha_inter-obs.R2_w0_sr1$alpha_inter_standard_error
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==2),"Te_inter"]<-c(obs.R2_w0_sr4_inter$alpha_inter-obs.R2_w0_sr4_inter$alpha_inter_standard_error, obs.R2_w0_sr5_inter$alpha_inter-obs.R2_w0_sr5_inter$alpha_inter_standard_error)

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==3),"Tu_inter"]<-c(obs.R3_w0$alpha_matrix[1,3]-obs.R3_w0$alpha_matrix_standard_error[1,3], obs.R3_w0$alpha_matrix[2,3]-obs.R3_w0$alpha_matrix_standard_error[2,3],obs.R3_w0$alpha_matrix[1,4]-obs.R3_w0$alpha_matrix_standard_error[1,4], obs.R3_w0$alpha_matrix[2,4]-obs.R3_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==3),"Te_inter"]<-c(obs.R3_w0$alpha_matrix[3,1]-obs.R3_w0$alpha_matrix_standard_error[3,1], obs.R3_w0$alpha_matrix[3,2]-obs.R3_w0$alpha_matrix_standard_error[3,2],obs.R3_w0$alpha_matrix[4,1]-obs.R3_w0$alpha_matrix_standard_error[4,1], obs.R3_w0$alpha_matrix[4,2]-obs.R3_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==4),"Tu_inter"]<-c(obs.R4_w0$alpha_matrix[1,3]-obs.R4_w0$alpha_matrix_standard_error[1,3], obs.R4_w0$alpha_matrix[2,3]-obs.R4_w0$alpha_matrix_standard_error[2,3],obs.R4_w0$alpha_matrix[1,4]-obs.R4_w0$alpha_matrix_standard_error[1,4], obs.R4_w0$alpha_matrix[2,4]-obs.R4_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==4),"Te_inter"]<-c(obs.R4_w0$alpha_matrix[3,1]-obs.R4_w0$alpha_matrix_standard_error[3,1], obs.R4_w0$alpha_matrix[3,2]-obs.R4_w0$alpha_matrix_standard_error[3,2],obs.R4_w0$alpha_matrix[4,1]-obs.R4_w0$alpha_matrix_standard_error[4,1], obs.R4_w0$alpha_matrix[4,2]-obs.R4_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==5),"Tu_inter"]<-c(obs.R5_w0$alpha_matrix[1,3]-obs.R5_w0$alpha_matrix_standard_error[1,3], obs.R5_w0$alpha_matrix[2,3]-obs.R5_w0$alpha_matrix_standard_error[2,3],obs.R5_w0$alpha_matrix[1,4]-obs.R5_w0$alpha_matrix_standard_error[1,4], obs.R5_w0$alpha_matrix[2,4]-obs.R5_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0_lower[which(cxr_param_w0_lower$Replicate==5),"Te_inter"]<-c(obs.R5_w0$alpha_matrix[3,1]-obs.R5_w0$alpha_matrix_standard_error[3,1], obs.R5_w0$alpha_matrix[3,2]-obs.R5_w0$alpha_matrix_standard_error[3,2],obs.R5_w0$alpha_matrix[4,1]-obs.R5_w0$alpha_matrix_standard_error[4,1], obs.R5_w0$alpha_matrix[4,2]-obs.R5_w0$alpha_matrix_standard_error[4,2])

### upper

cxr_param_w0_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_w0_upper$Tu_lambda<-0
cxr_param_w0_upper$Te_lambda<-0
cxr_param_w0_upper$Tu_intra<-0
cxr_param_w0_upper$Te_intra<-0
cxr_param_w0_upper$Tu_inter<-0
cxr_param_w0_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0_upper<-cxr_param_w0_upper[-which(cxr_param_w0_upper$Replicate==2 & cxr_param_w0_upper$Tu_Regime=="SR2"),]


cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==1),"Tu_lambda"]<-rep(c(obs.R1_w0$lambda[1]+obs.R1_w0$lambda_standard_error[1], obs.R1_w0$lambda[2]+obs.R1_w0$lambda_standard_error[2]), 2)
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==1),"Te_lambda"]<-rep(c(obs.R1_w0$lambda[3]+obs.R1_w0$lambda_standard_error[3], obs.R1_w0$lambda[4]+obs.R1_w0$lambda_standard_error[4]), each=2)

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==2),"Tu_lambda"]<-obs.R2_w0_sr1$lambda+ obs.R2_w0_sr1$lambda_standard_error
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==2),"Te_lambda"]<-c(obs.R2_w0_sr4$lambda+obs.R2_w0_sr4$lambda_standard_error, obs.R2_w0_sr5$lambda+obs.R2_w0_sr5$lambda_standard_error)

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==3),"Tu_lambda"]<-rep(c(obs.R3_w0$lambda[1]+obs.R3_w0$lambda_standard_error[1], obs.R3_w0$lambda[2]+obs.R3_w0$lambda_standard_error[2]), 2)
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==3),"Te_lambda"]<-rep(c(obs.R3_w0$lambda[3]+obs.R3_w0$lambda_standard_error[3], obs.R3_w0$lambda[4]+obs.R3_w0$lambda_standard_error[4]), each=2)

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==4),"Tu_lambda"]<-rep(c(obs.R4_w0$lambda[1]+obs.R4_w0$lambda_standard_error[1], obs.R4_w0$lambda[2]+obs.R4_w0$lambda_standard_error[2]), 2)
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==4),"Te_lambda"]<-rep(c(obs.R4_w0$lambda[3]+obs.R4_w0$lambda_standard_error[3], obs.R4_w0$lambda[4]+obs.R4_w0$lambda_standard_error[4]), each=2)

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==5),"Tu_lambda"]<-rep(c(obs.R5_w0$lambda[1]+obs.R5_w0$lambda_standard_error[1], obs.R5_w0$lambda[2]+obs.R5_w0$lambda_standard_error[2]), 2)
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==5),"Te_lambda"]<-rep(c(obs.R5_w0$lambda[3]+obs.R5_w0$lambda_standard_error[3], obs.R5_w0$lambda[4]+obs.R5_w0$lambda_standard_error[4]), each=2)


cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==1),"Tu_intra"]<-rep(c(obs.R1_w0$alpha_matrix[1,1]+obs.R1_w0$alpha_matrix_standard_error[1,1], obs.R1_w0$alpha_matrix[2,2]+obs.R1_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==1),"Te_intra"]<-rep(c(obs.R1_w0$alpha_matrix[3,3]+obs.R1_w0$alpha_matrix_standard_error[3,3], obs.R1_w0$alpha_matrix[4,4]+obs.R1_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==2),"Tu_intra"]<-obs.R2_w0_sr1$alpha_intra+obs.R2_w0_sr1$alpha_intra_standard_error
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==2),"Te_intra"]<-c(obs.R2_w0_sr4$alpha_inter+obs.R2_w0_sr4$alpha_inter_standard_error,  obs.R2_w0_sr5$alpha_inter+obs.R2_w0_sr5$alpha_inter_standard_error)

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==3),"Tu_intra"]<-rep(c(obs.R3_w0$alpha_matrix[1,1]+obs.R3_w0$alpha_matrix_standard_error[1,1], obs.R3_w0$alpha_matrix[2,2]+obs.R3_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==3),"Te_intra"]<-rep(c(obs.R3_w0$alpha_matrix[3,3]+obs.R3_w0$alpha_matrix_standard_error[3,3], obs.R3_w0$alpha_matrix[4,4]+obs.R3_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==4),"Tu_intra"]<-rep(c(obs.R4_w0$alpha_matrix[1,1]+obs.R4_w0$alpha_matrix_standard_error[1,1], obs.R4_w0$alpha_matrix[2,2]+obs.R4_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==4),"Te_intra"]<-rep(c(obs.R4_w0$alpha_matrix[3,3]+obs.R4_w0$alpha_matrix_standard_error[3,3], obs.R4_w0$alpha_matrix[4,4]+obs.R4_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==5),"Tu_intra"]<-rep(c(obs.R5_w0$alpha_matrix[1,1]+obs.R5_w0$alpha_matrix_standard_error[1,1], obs.R5_w0$alpha_matrix[2,2]+obs.R5_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==5),"Te_intra"]<-rep(c(obs.R5_w0$alpha_matrix[3,3]+obs.R5_w0$alpha_matrix_standard_error[3,3], obs.R5_w0$alpha_matrix[4,4]+obs.R5_w0$alpha_matrix_standard_error[4,4]), each=2)


cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==1),"Tu_inter"]<-c(obs.R1_w0$alpha_matrix[1,3]+obs.R1_w0$alpha_matrix_standard_error[1,3], obs.R1_w0$alpha_matrix[2,3]+obs.R1_w0$alpha_matrix_standard_error[2,3],obs.R1_w0$alpha_matrix[1,4]+obs.R1_w0$alpha_matrix_standard_error[1,4], obs.R1_w0$alpha_matrix[2,4]+obs.R1_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==1),"Te_inter"]<-c(obs.R1_w0$alpha_matrix[3,1]+obs.R1_w0$alpha_matrix_standard_error[3,1], obs.R1_w0$alpha_matrix[3,2]+obs.R1_w0$alpha_matrix_standard_error[3,2],obs.R1_w0$alpha_matrix[4,1]+obs.R1_w0$alpha_matrix_standard_error[4,1], obs.R1_w0$alpha_matrix[4,2]+obs.R1_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==2),"Tu_inter"]<-c(obs.R2_w0_sr1$alpha_inter+obs.R2_w0_sr1$alpha_inter_standard_error)
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==2),"Te_inter"]<-c(obs.R2_w0_sr4_inter$alpha_inter+obs.R2_w0_sr4_inter$alpha_inter_standard_error, obs.R2_w0_sr5_inter$alpha_inter+obs.R2_w0_sr5_inter$alpha_inter_standard_error)

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==3),"Tu_inter"]<-c(obs.R3_w0$alpha_matrix[1,3]+obs.R3_w0$alpha_matrix_standard_error[1,3], obs.R3_w0$alpha_matrix[2,3]+obs.R3_w0$alpha_matrix_standard_error[2,3],obs.R3_w0$alpha_matrix[1,4]+obs.R3_w0$alpha_matrix_standard_error[1,4], obs.R3_w0$alpha_matrix[2,4]+obs.R3_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==3),"Te_inter"]<-c(obs.R3_w0$alpha_matrix[3,1]+obs.R3_w0$alpha_matrix_standard_error[3,1], obs.R3_w0$alpha_matrix[3,2]+obs.R3_w0$alpha_matrix_standard_error[3,2],obs.R3_w0$alpha_matrix[4,1]+obs.R3_w0$alpha_matrix_standard_error[4,1], obs.R3_w0$alpha_matrix[4,2]+obs.R3_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==4),"Tu_inter"]<-c(obs.R4_w0$alpha_matrix[1,3]+obs.R4_w0$alpha_matrix_standard_error[1,3], obs.R4_w0$alpha_matrix[2,3]+obs.R4_w0$alpha_matrix_standard_error[2,3],obs.R4_w0$alpha_matrix[1,4]+obs.R4_w0$alpha_matrix_standard_error[1,4], obs.R4_w0$alpha_matrix[2,4]+obs.R4_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==4),"Te_inter"]<-c(obs.R4_w0$alpha_matrix[3,1]+obs.R4_w0$alpha_matrix_standard_error[3,1], obs.R4_w0$alpha_matrix[3,2]+obs.R4_w0$alpha_matrix_standard_error[3,2],obs.R4_w0$alpha_matrix[4,1]+obs.R4_w0$alpha_matrix_standard_error[4,1], obs.R4_w0$alpha_matrix[4,2]+obs.R4_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==5),"Tu_inter"]<-c(obs.R5_w0$alpha_matrix[1,3]+obs.R5_w0$alpha_matrix_standard_error[1,3], obs.R5_w0$alpha_matrix[2,3]+obs.R5_w0$alpha_matrix_standard_error[2,3],obs.R5_w0$alpha_matrix[1,4]+obs.R5_w0$alpha_matrix_standard_error[1,4], obs.R5_w0$alpha_matrix[2,4]+obs.R5_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0_upper[which(cxr_param_w0_upper$Replicate==5),"Te_inter"]<-c(obs.R5_w0$alpha_matrix[3,1]+obs.R5_w0$alpha_matrix_standard_error[3,1], obs.R5_w0$alpha_matrix[3,2]+obs.R5_w0$alpha_matrix_standard_error[3,2],obs.R5_w0$alpha_matrix[4,1]+obs.R5_w0$alpha_matrix_standard_error[4,1], obs.R5_w0$alpha_matrix[4,2]+obs.R5_w0$alpha_matrix_standard_error[4,2])


##### Cadmium
# modifying data frame to fit the type of setup that is need for CXR
forCXR_Cd<-subset(ca, Env=="Cd")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "TuFemales")]

forCXR_Cd$Focal<-mapvalues(forCXR_Cd$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
forCXR_Cd$CompSR2<-mapvalues(forCXR_Cd$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

forCXR_Cd$Comp<-sapply(c(1:length(forCXR_Cd[,1])), function(x){
  if(is.na(forCXR_Cd$CompSR2[x])){
    a<- forCXR_Cd$Focal[x]
  }else{
    a<-forCXR_Cd$CompSR2[x]
  }
  
  a
})

aux<-data.frame(SR1=rep(0, length(forCXR_Cd[,1])), SR2=rep(0, length(forCXR_Cd[,1])), SR4=rep(0, length(forCXR_Cd[,1])), SR5=rep(0, length(forCXR_Cd[,1])))

for(i in 1:length(forCXR_Cd[,1])){
  #coluna onde por focais
  colunaF<-which(colnames(aux)==forCXR_Cd$Focal[i])
  #coluna onde por competidors
  colunaC<-which(colnames(aux)==forCXR_Cd$Comp[i])
  
  #if its the same regime
  if(forCXR_Cd$Focal[i]==forCXR_Cd$Comp[i] & forCXR_Cd$Dens[i]==1){
    aux[i,colunaF]<-forCXR_Cd$Dens[i]-1
    
  }else if(forCXR_Cd$Focal[i]==forCXR_Cd$Comp[i]){
    aux[i,colunaF]<-forCXR_Cd$Dens[i]-1
  }else{ #if it is heterospecific then its -1 for the competitors (because of the focal) and its one for the focal
    aux[i,colunaC]<-forCXR_Cd$Dens[i]-1
    aux[i, colunaF]<-1
  }
  
}

forCXR_Cd<-cbind(forCXR_Cd, aux)

forCXR_Cd$fitness<-sapply(c(1:length(forCXR_Cd[,1])), function(x){
  colF<-which(colnames(forCXR_Cd)==forCXR_Cd$Focal[x])
  
  if(forCXR_Cd$Focal[x]=="SR1"){
    a<-forCXR_Cd$TuFemales[x]/forCXR_Cd$SR1[x]
  } else if(forCXR_Cd$Focal[x]=="SR2"){
    a<-forCXR_Cd$TuFemales[x]/forCXR_Cd$SR2[x]
  } else if(forCXR_Cd$Focal[x]=="SR4"){
    a<-forCXR_Cd$TeFemales[x]/forCXR_Cd$SR4[x]
  } else if(forCXR_Cd$Focal[x]=="SR5"){
    a<-forCXR_Cd$TeFemales[x]/forCXR_Cd$SR5[x]
  }
  
  a
})

subset(ca, Env=="Cd" & Rep=="2" & FocalSR==5 &Type=="INTER")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "Block")]

#removing rows for which there is no data for fitness
#forCXR_Cd<-forCXR_Cd[-which(is.na(forCXR_Cd$fitness)),]
#forCXR_Cd$fitness<-forCXR_Cd$fitness+1

forCXR_Cd[which(forCXR_Cd$fitness=="-Inf" | forCXR_Cd$fitness=="Inf"),"fitness"]<-0

#0 to 1 to mainrain data
forCXR_Cd<-forCXR_Cd[-which(is.na(forCXR_Cd$fitness)),]
forCXR_Cd$fitness<-forCXR_Cd$fitness+1



# vector that tells which are the selection regimes, the columns have to have the same name
my.reg <- c("SR1", "SR2","SR4","SR5")

# Do list per replicate and environment
R1_Cd<-list(SR1= subset(forCXR_Cd, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R2_Cd<-list(SR1= subset(forCXR_Cd, Rep==2 & Focal=="SR1")[,c("fitness", "SR1", "SR2","SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R3_Cd<-list(SR1= subset(forCXR_Cd, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R4_Cd<-list(SR1= subset(forCXR_Cd, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R5_Cd<-list(SR1= subset(forCXR_Cd, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])


obs.R1_Cd_w0<-cxr_pm_multifit(data = R1_Cd,
                              focal_column = my.reg,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = 1,
                                                    alpha_intra = 0.1,
                                                    alpha_inter = 0.1),
                              fixed_terms = NULL,
                              # no standard errors
                              bootstrap_samples = 200)

# replicate 2 below


obs.R3_Cd_w0<-cxr_pm_multifit(data = R3_Cd,
                              focal_column = my.reg,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = 1,
                                                    alpha_intra = 0.1,
                                                    alpha_inter = 0.1),
                              fixed_terms = NULL,
                              # no standard errors
                              bootstrap_samples =10)

obs.R4_Cd_w0<-cxr_pm_multifit(data = R4_Cd,
                              focal_column = my.reg,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = 1,
                                                    alpha_intra = 0.1,
                                                    alpha_inter = 0.1),
                              fixed_terms = NULL,
                              # no standard errors
                              bootstrap_samples = 200)

obs.R5_Cd_w0<-cxr_pm_multifit(data = R5_Cd,
                              focal_column = my.reg,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = 1,
                                                    alpha_intra = 0.1,
                                                    alpha_inter = 0.1),
                              fixed_terms = NULL,
                              # no standard errors
                              bootstrap_samples = 200)

summary(obs.R1_Cd_w0)
#summary(obs.R2_Cd_w0)
summary(obs.R3_Cd_w0)
summary(obs.R4_Cd_w0)
summary(obs.R5_Cd_w0)

# This one works well
obs.R2_Cd_w0_sr1<-cxr_pm_fit(data = R2_Cd[[1]],
                             focal_column = my.reg[1],
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = 1,
                                                   alpha_intra = 0.1,
                                                   alpha_inter = 0.1),
                             fixed_terms = NULL,
                             # no standard errors
                             bootstrap_samples = 1000)

obs.R2_Cd_w0_sr4<-cxr_pm_fit(data = R2_Cd[[2]][which(R2_Cd[[2]][,"SR1"]==0),c("fitness", "SR4")],
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "global",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = 1,
                                                   alpha_inter = 0.1),
                             fixed_terms = NULL,
                             # no standard errors
                             bootstrap_samples = 1000)

obs.R2_Cd_w0_sr5<-cxr_pm_fit(data = R2_Cd[[3]][which(R2_Cd[[3]][,"SR1"]==0),c("fitness", "SR5")],
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "global",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = 1,
                                                   alpha_inter = 0.1),
                             fixed_terms = NULL,
                             # no standard errors
                             bootstrap_samples = 1000)

obs.R2_Cd_w0_sr4_inter<-cxr_pm_fit(data = R2_Cd[[2]][which(R2_Cd[[2]][,"SR1"]!=0),c("fitness", "SR1")],
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list( alpha_inter = 0.1),
                                   fixed_terms = list(lambda=obs.R2_Cd_w0_sr4$lambda),
                                   # no standard errors
                                   bootstrap_samples = 1000)

obs.R2_Cd_w0_sr5_inter<-cxr_pm_fit(data = R2_Cd[[3]][which(R2_Cd[[3]][,"SR1"]!=0),c("fitness", "SR1")],
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list( alpha_inter = 0.1),
                                   fixed_terms = list(lambda=obs.R2_Cd_w0_sr5$lambda),
                                   # no standard errors
                                   bootstrap_samples = 1000)

###### data table summary

cxr_param_w0C<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_w0C$Tu_lambda<-0
cxr_param_w0C$Te_lambda<-0
cxr_param_w0C$Tu_intra<-0
cxr_param_w0C$Te_intra<-0
cxr_param_w0C$Tu_inter<-0
cxr_param_w0C$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0C<-cxr_param_w0C[-which(cxr_param_w0C$Replicate==2 & cxr_param_w0C$Tu_Regime=="SR2"),]


cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Tu_lambda"]<-obs.R1_Cd_w0$lambda[1:2]
cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Te_lambda"]<-obs.R1_Cd_w0$lambda[c(3,3,4,4)]

cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Tu_lambda"]<-obs.R2_Cd_w0_sr1$lambda
cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Te_lambda"]<-c(obs.R2_Cd_w0_sr4$lambda, obs.R2_Cd_w0_sr5$lambda)

cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Tu_lambda"]<-obs.R3_Cd_w0$lambda[1:2]
cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Te_lambda"]<-obs.R3_Cd_w0$lambda[c(3,3,4,4)]

cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Tu_lambda"]<-obs.R4_Cd_w0$lambda[1:2]
cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Te_lambda"]<-obs.R4_Cd_w0$lambda[c(3,3,4,4)]

cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Tu_lambda"]<-obs.R5_Cd_w0$lambda[1:2]
cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Te_lambda"]<-obs.R5_Cd_w0$lambda[c(3,3,4,4)]


cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Tu_intra"]<-rep(c(obs.R1_Cd_w0$alpha_matrix[1,1], obs.R1_Cd_w0$alpha_matrix[2,2]), 2)
cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Te_intra"]<-rep(c(obs.R1_Cd_w0$alpha_matrix[3,3], obs.R1_Cd_w0$alpha_matrix[4,4]), each=2)

cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Tu_intra"]<-obs.R2_Cd_w0_sr1$alpha_intra
cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Te_intra"]<-c(obs.R2_Cd_w0_sr4$alpha_inter, obs.R2_Cd_w0_sr5$alpha_inter)

cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Tu_intra"]<-rep(c(obs.R3_Cd_w0$alpha_matrix[1,1], obs.R3_Cd_w0$alpha_matrix[2,2]), 2)
cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Te_intra"]<-rep(c(obs.R3_Cd_w0$alpha_matrix[3,3], obs.R3_Cd_w0$alpha_matrix[4,4]), each=2)

cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Tu_intra"]<-rep(c(obs.R4_Cd_w0$alpha_matrix[1,1], obs.R4_Cd_w0$alpha_matrix[2,2]), 2)
cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Te_intra"]<-rep(c(obs.R4_Cd_w0$alpha_matrix[3,3], obs.R4_Cd_w0$alpha_matrix[4,4]), each=2)

cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Tu_intra"]<-rep(c(obs.R5_Cd_w0$alpha_matrix[1,1], obs.R5_Cd_w0$alpha_matrix[2,2]), 2)
cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Te_intra"]<-rep(c(obs.R5_Cd_w0$alpha_matrix[3,3], obs.R5_Cd_w0$alpha_matrix[4,4]), each=2)


cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Tu_inter"]<-c(obs.R1_Cd_w0$alpha_matrix[1,3], obs.R1_Cd_w0$alpha_matrix[2,3],obs.R1_Cd_w0$alpha_matrix[1,4], obs.R1_Cd_w0$alpha_matrix[2,4])
cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Te_inter"]<-c(obs.R1_Cd_w0$alpha_matrix[3,1], obs.R1_Cd_w0$alpha_matrix[3,2],obs.R1_Cd_w0$alpha_matrix[4,1], obs.R1_Cd_w0$alpha_matrix[4,2])

cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Tu_inter"]<-obs.R2_Cd_w0_sr1$alpha_inter[2:3]
cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Te_inter"]<-c(obs.R2_Cd_w0_sr4_inter$alpha_inter,  obs.R2_Cd_w0_sr5_inter$alpha_inter)

cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Tu_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[1,3], obs.R3_Cd_w0$alpha_matrix[2,3],obs.R3_Cd_w0$alpha_matrix[1,4], obs.R3_Cd_w0$alpha_matrix[2,4])
cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Te_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[3,1], obs.R3_Cd_w0$alpha_matrix[3,2],obs.R3_Cd_w0$alpha_matrix[4,1], obs.R3_Cd_w0$alpha_matrix[4,2])

cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Tu_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[1,3], obs.R4_Cd_w0$alpha_matrix[2,3],obs.R4_Cd_w0$alpha_matrix[1,4], obs.R4_Cd_w0$alpha_matrix[2,4])
cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Te_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[3,1], obs.R4_Cd_w0$alpha_matrix[3,2],obs.R4_Cd_w0$alpha_matrix[4,1], obs.R4_Cd_w0$alpha_matrix[4,2])

cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Tu_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[1,3], obs.R5_Cd_w0$alpha_matrix[2,3],obs.R5_Cd_w0$alpha_matrix[1,4], obs.R5_Cd_w0$alpha_matrix[2,4])
cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Te_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[3,1], obs.R5_Cd_w0$alpha_matrix[3,2],obs.R5_Cd_w0$alpha_matrix[4,1], obs.R5_Cd_w0$alpha_matrix[4,2])


### Lower

cxr_param_w0C_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_w0C_lower$Tu_lambda<-0
cxr_param_w0C_lower$Te_lambda<-0
cxr_param_w0C_lower$Tu_intra<-0
cxr_param_w0C_lower$Te_intra<-0
cxr_param_w0C_lower$Tu_inter<-0
cxr_param_w0C_lower$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0C_lower<-cxr_param_w0C_lower[-which(cxr_param_w0C_lower$Replicate==2 & cxr_param_w0C_lower$Tu_Regime=="SR2"),]


cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==1),"Tu_lambda"]<-rep(c(obs.R1_Cd_w0$lambda[1]-obs.R1_Cd_w0$lambda_standard_error[1], obs.R1_Cd_w0$lambda[2]-obs.R1_Cd_w0$lambda_standard_error[2]), 2)
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==1),"Te_lambda"]<-rep(c(obs.R1_Cd_w0$lambda[3]-obs.R1_Cd_w0$lambda_standard_error[3], obs.R1_Cd_w0$lambda[4]-obs.R1_Cd_w0$lambda_standard_error[4]), each=2)

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==2),"Tu_lambda"]<-c(obs.R2_Cd_w0_sr1$lambda-obs.R2_Cd_w0_sr1$lambda_standard_error)
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==2),"Te_lambda"]<-c(obs.R2_Cd_w0_sr4$lambda-obs.R2_Cd_w0_sr4$lambda_standard_error, obs.R2_Cd_w0_sr5$lambda-obs.R2_Cd_w0_sr5$lambda_standard_error)

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==3),"Tu_lambda"]<-rep(c(obs.R3_Cd_w0$lambda[1]-obs.R3_Cd_w0$lambda_standard_error[1], obs.R3_Cd_w0$lambda[2]-obs.R3_Cd_w0$lambda_standard_error[2]), 2)
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==3),"Te_lambda"]<-rep(c(obs.R3_Cd_w0$lambda[3]-obs.R3_Cd_w0$lambda_standard_error[3], obs.R3_Cd_w0$lambda[4]-obs.R3_Cd_w0$lambda_standard_error[4]), each=2)

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==4),"Tu_lambda"]<-rep(c(obs.R4_Cd_w0$lambda[1]-obs.R4_Cd_w0$lambda_standard_error[1], obs.R4_Cd_w0$lambda[2]-obs.R4_Cd_w0$lambda_standard_error[2]), 2)
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==4),"Te_lambda"]<-rep(c(obs.R4_Cd_w0$lambda[3]-obs.R4_Cd_w0$lambda_standard_error[3], obs.R4_Cd_w0$lambda[4]-obs.R4_Cd_w0$lambda_standard_error[4]), each=2)

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==5),"Tu_lambda"]<-rep(c(obs.R5_Cd_w0$lambda[1]-obs.R5_Cd_w0$lambda_standard_error[1], obs.R5_Cd_w0$lambda[2]-obs.R5_Cd_w0$lambda_standard_error[2]), 2)
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==5),"Te_lambda"]<-rep(c(obs.R5_Cd_w0$lambda[3]-obs.R5_Cd_w0$lambda_standard_error[3], obs.R5_Cd_w0$lambda[4]-obs.R5_Cd_w0$lambda_standard_error[4]), each=2)


cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==1),"Tu_intra"]<-rep(c(obs.R1_Cd_w0$alpha_matrix[1,1]-obs.R1_Cd_w0$alpha_matrix_standard_error[1,1], obs.R1_Cd_w0$alpha_matrix[2,2]-obs.R1_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==1),"Te_intra"]<-rep(c(obs.R1_Cd_w0$alpha_matrix[3,3]-obs.R1_Cd_w0$alpha_matrix_standard_error[3,3], obs.R1_Cd_w0$alpha_matrix[4,4]-obs.R1_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==2),"Tu_intra"]<-obs.R2_Cd_w0_sr1$alpha_intra-obs.R2_Cd_w0_sr1$alpha_intra_standard_error
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==2),"Te_intra"]<-c(obs.R2_Cd_w0_sr4$alpha_inter-obs.R2_Cd_w0_sr4$alpha_inter_standard_error, obs.R2_Cd_w0_sr5$alpha_inter-obs.R2_Cd_w0_sr5$alpha_inter_standard_error)

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==3),"Tu_intra"]<-rep(c(obs.R3_Cd_w0$alpha_matrix[1,1]-obs.R3_Cd_w0$alpha_matrix_standard_error[1,1], obs.R3_Cd_w0$alpha_matrix[2,2]-obs.R3_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==3),"Te_intra"]<-rep(c(obs.R3_Cd_w0$alpha_matrix[3,3]-obs.R3_Cd_w0$alpha_matrix_standard_error[3,3], obs.R3_Cd_w0$alpha_matrix[4,4]-obs.R3_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==4),"Tu_intra"]<-rep(c(obs.R4_Cd_w0$alpha_matrix[1,1]-obs.R4_Cd_w0$alpha_matrix_standard_error[1,1], obs.R4_Cd_w0$alpha_matrix[2,2]-obs.R4_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==4),"Te_intra"]<-rep(c(obs.R4_Cd_w0$alpha_matrix[3,3]-obs.R4_Cd_w0$alpha_matrix_standard_error[3,3], obs.R4_Cd_w0$alpha_matrix[4,4]-obs.R4_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==5),"Tu_intra"]<-rep(c(obs.R5_Cd_w0$alpha_matrix[1,1]-obs.R5_Cd_w0$alpha_matrix_standard_error[1,1], obs.R5_Cd_w0$alpha_matrix[2,2]-obs.R5_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==5),"Te_intra"]<-rep(c(obs.R5_Cd_w0$alpha_matrix[3,3]-obs.R5_Cd_w0$alpha_matrix_standard_error[3,3], obs.R5_Cd_w0$alpha_matrix[4,4]-obs.R5_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)


cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==1),"Tu_inter"]<-c(obs.R1_Cd_w0$alpha_matrix[1,3]-obs.R1_Cd_w0$alpha_matrix_standard_error[1,3], obs.R1_Cd_w0$alpha_matrix[2,3]-obs.R1_Cd_w0$alpha_matrix_standard_error[2,3],obs.R1_Cd_w0$alpha_matrix[1,4]-obs.R1_Cd_w0$alpha_matrix_standard_error[1,4], obs.R1_Cd_w0$alpha_matrix[2,4]-obs.R1_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==1),"Te_inter"]<-c(obs.R1_Cd_w0$alpha_matrix[3,1]-obs.R1_Cd_w0$alpha_matrix_standard_error[3,1], obs.R1_Cd_w0$alpha_matrix[3,2]-obs.R1_Cd_w0$alpha_matrix_standard_error[3,2],obs.R1_Cd_w0$alpha_matrix[4,1]-obs.R1_Cd_w0$alpha_matrix_standard_error[4,1], obs.R1_Cd_w0$alpha_matrix[4,2]-obs.R1_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==2),"Tu_inter"]<-obs.R2_Cd_w0_sr1$alpha_inter[2:3]-obs.R2_Cd_w0_sr1$alpha_inter_standard_error[2:3]
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==2),"Te_inter"]<-c(obs.R2_Cd_w0_sr4_inter$alpha_inter-obs.R2_Cd_w0_sr4_inter$alpha_inter_standard_error, obs.R2_Cd_w0_sr5_inter$alpha_inter-obs.R2_Cd_w0_sr5_inter$alpha_inter_standard_error)

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==3),"Tu_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[1,3]-obs.R3_Cd_w0$alpha_matrix_standard_error[1,3], obs.R3_Cd_w0$alpha_matrix[2,3]-obs.R3_Cd_w0$alpha_matrix_standard_error[2,3],obs.R3_Cd_w0$alpha_matrix[1,4]-obs.R3_Cd_w0$alpha_matrix_standard_error[1,4], obs.R3_Cd_w0$alpha_matrix[2,4]-obs.R3_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==3),"Te_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[3,1]-obs.R3_Cd_w0$alpha_matrix_standard_error[3,1], obs.R3_Cd_w0$alpha_matrix[3,2]-obs.R3_Cd_w0$alpha_matrix_standard_error[3,2],obs.R3_Cd_w0$alpha_matrix[4,1]-obs.R3_Cd_w0$alpha_matrix_standard_error[4,1], obs.R3_Cd_w0$alpha_matrix[4,2]-obs.R3_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==4),"Tu_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[1,3]-obs.R4_Cd_w0$alpha_matrix_standard_error[1,3], obs.R4_Cd_w0$alpha_matrix[2,3]-obs.R4_Cd_w0$alpha_matrix_standard_error[2,3],obs.R4_Cd_w0$alpha_matrix[1,4]-obs.R4_Cd_w0$alpha_matrix_standard_error[1,4], obs.R4_Cd_w0$alpha_matrix[2,4]-obs.R4_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==4),"Te_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[3,1]-obs.R4_Cd_w0$alpha_matrix_standard_error[3,1], obs.R4_Cd_w0$alpha_matrix[3,2]-obs.R4_Cd_w0$alpha_matrix_standard_error[3,2],obs.R4_Cd_w0$alpha_matrix[4,1]-obs.R4_Cd_w0$alpha_matrix_standard_error[4,1], obs.R4_Cd_w0$alpha_matrix[4,2]-obs.R4_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==5),"Tu_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[1,3]-obs.R5_Cd_w0$alpha_matrix_standard_error[1,3], obs.R5_Cd_w0$alpha_matrix[2,3]-obs.R5_Cd_w0$alpha_matrix_standard_error[2,3],obs.R5_Cd_w0$alpha_matrix[1,4]-obs.R5_Cd_w0$alpha_matrix_standard_error[1,4], obs.R5_Cd_w0$alpha_matrix[2,4]-obs.R5_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==5),"Te_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[3,1]-obs.R5_Cd_w0$alpha_matrix_standard_error[3,1], obs.R5_Cd_w0$alpha_matrix[3,2]-obs.R5_Cd_w0$alpha_matrix_standard_error[3,2],obs.R5_Cd_w0$alpha_matrix[4,1]-obs.R5_Cd_w0$alpha_matrix_standard_error[4,1], obs.R5_Cd_w0$alpha_matrix[4,2]-obs.R5_Cd_w0$alpha_matrix_standard_error[4,2])

### upper

cxr_param_w0C_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_w0C_upper$Tu_lambda<-0
cxr_param_w0C_upper$Te_lambda<-0
cxr_param_w0C_upper$Tu_intra<-0
cxr_param_w0C_upper$Te_intra<-0
cxr_param_w0C_upper$Tu_inter<-0
cxr_param_w0C_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0C_upper<-cxr_param_w0C_upper[-which(cxr_param_w0C_upper$Replicate==2 & cxr_param_w0C_upper$Tu_Regime=="SR2"),]


cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==1),"Tu_lambda"]<-rep(c(obs.R1_Cd_w0$lambda[1]+obs.R1_Cd_w0$lambda_standard_error[1], obs.R1_Cd_w0$lambda[2]+obs.R1_Cd_w0$lambda_standard_error[2]), 2)
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==1),"Te_lambda"]<-rep(c(obs.R1_Cd_w0$lambda[3]+obs.R1_Cd_w0$lambda_standard_error[3], obs.R1_Cd_w0$lambda[4]+obs.R1_Cd_w0$lambda_standard_error[4]), each=2)

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==2),"Tu_lambda"]<-c(obs.R2_Cd_w0_sr1$lambda+obs.R2_Cd_w0_sr1$lambda_standard_error)
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==2),"Te_lambda"]<-c(obs.R2_Cd_w0_sr4$lambda+obs.R2_Cd_w0_sr4$lambda_standard_error, obs.R2_Cd_w0_sr5$lambda+obs.R2_Cd_w0_sr5$lambda_standard_error)

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==3),"Tu_lambda"]<-rep(c(obs.R3_Cd_w0$lambda[1]+obs.R3_Cd_w0$lambda_standard_error[1], obs.R3_Cd_w0$lambda[2]+obs.R3_Cd_w0$lambda_standard_error[2]), 2)
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==3),"Te_lambda"]<-rep(c(obs.R3_Cd_w0$lambda[3]+obs.R3_Cd_w0$lambda_standard_error[3], obs.R3_Cd_w0$lambda[4]+obs.R3_Cd_w0$lambda_standard_error[4]), each=2)

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==4),"Tu_lambda"]<-rep(c(obs.R4_Cd_w0$lambda[1]+obs.R4_Cd_w0$lambda_standard_error[1], obs.R4_Cd_w0$lambda[2]+obs.R4_Cd_w0$lambda_standard_error[2]), 2)
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==4),"Te_lambda"]<-rep(c(obs.R4_Cd_w0$lambda[3]+obs.R4_Cd_w0$lambda_standard_error[3], obs.R4_Cd_w0$lambda[4]+obs.R4_Cd_w0$lambda_standard_error[4]), each=2)

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==5),"Tu_lambda"]<-rep(c(obs.R5_Cd_w0$lambda[1]+obs.R5_Cd_w0$lambda_standard_error[1], obs.R5_Cd_w0$lambda[2]+obs.R5_Cd_w0$lambda_standard_error[2]), 2)
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==5),"Te_lambda"]<-rep(c(obs.R5_Cd_w0$lambda[3]+obs.R5_Cd_w0$lambda_standard_error[3], obs.R5_Cd_w0$lambda[4]+obs.R5_Cd_w0$lambda_standard_error[4]), each=2)


cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==1),"Tu_intra"]<-rep(c(obs.R1_Cd_w0$alpha_matrix[1,1]+obs.R1_Cd_w0$alpha_matrix_standard_error[1,1], obs.R1_Cd_w0$alpha_matrix[2,2]+obs.R1_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==1),"Te_intra"]<-rep(c(obs.R1_Cd_w0$alpha_matrix[3,3]+obs.R1_Cd_w0$alpha_matrix_standard_error[3,3], obs.R1_Cd_w0$alpha_matrix[4,4]+obs.R1_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==2),"Tu_intra"]<-obs.R2_Cd_w0_sr1$alpha_intra + obs.R2_Cd_w0_sr1$alpha_intra_standard_error
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==2),"Te_intra"]<-c(obs.R2_Cd_w0_sr4$alpha_inter+obs.R2_Cd_w0_sr4$alpha_inter_standard_error, obs.R2_Cd_w0_sr5$alpha_inter+obs.R2_Cd_w0_sr5$alpha_inter_standard_error)


cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==3),"Tu_intra"]<-rep(c(obs.R3_Cd_w0$alpha_matrix[1,1]+obs.R3_Cd_w0$alpha_matrix_standard_error[1,1], obs.R3_Cd_w0$alpha_matrix[2,2]+obs.R3_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==3),"Te_intra"]<-rep(c(obs.R3_Cd_w0$alpha_matrix[3,3]+obs.R3_Cd_w0$alpha_matrix_standard_error[3,3], obs.R3_Cd_w0$alpha_matrix[4,4]+obs.R3_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==4),"Tu_intra"]<-rep(c(obs.R4_Cd_w0$alpha_matrix[1,1]+obs.R4_Cd_w0$alpha_matrix_standard_error[1,1], obs.R4_Cd_w0$alpha_matrix[2,2]+obs.R4_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==4),"Te_intra"]<-rep(c(obs.R4_Cd_w0$alpha_matrix[3,3]+obs.R4_Cd_w0$alpha_matrix_standard_error[3,3], obs.R4_Cd_w0$alpha_matrix[4,4]+obs.R4_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==5),"Tu_intra"]<-rep(c(obs.R5_Cd_w0$alpha_matrix[1,1]+obs.R5_Cd_w0$alpha_matrix_standard_error[1,1], obs.R5_Cd_w0$alpha_matrix[2,2]+obs.R5_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==5),"Te_intra"]<-rep(c(obs.R5_Cd_w0$alpha_matrix[3,3]+obs.R5_Cd_w0$alpha_matrix_standard_error[3,3], obs.R5_Cd_w0$alpha_matrix[4,4]+obs.R5_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)


cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==1),"Tu_inter"]<-c(obs.R1_Cd_w0$alpha_matrix[1,3]+obs.R1_Cd_w0$alpha_matrix_standard_error[1,3], obs.R1_Cd_w0$alpha_matrix[2,3]+obs.R1_Cd_w0$alpha_matrix_standard_error[2,3],obs.R1_Cd_w0$alpha_matrix[1,4]+obs.R1_Cd_w0$alpha_matrix_standard_error[1,4], obs.R1_Cd_w0$alpha_matrix[2,4]+obs.R1_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==1),"Te_inter"]<-c(obs.R1_Cd_w0$alpha_matrix[3,1]+obs.R1_Cd_w0$alpha_matrix_standard_error[3,1], obs.R1_Cd_w0$alpha_matrix[3,2]+obs.R1_Cd_w0$alpha_matrix_standard_error[3,2],obs.R1_Cd_w0$alpha_matrix[4,1]+obs.R1_Cd_w0$alpha_matrix_standard_error[4,1], obs.R1_Cd_w0$alpha_matrix[4,2]+obs.R1_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==2),"Tu_inter"]<-obs.R2_Cd_w0_sr1$alpha_inter[2:3]+obs.R2_Cd_w0_sr1$alpha_inter_standard_error[2:3]
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==2),"Te_inter"]<-c(obs.R2_Cd_w0_sr4_inter$alpha_inter+obs.R2_Cd_w0_sr4_inter$alpha_inter_standard_error, obs.R2_Cd_w0_sr5_inter$alpha_inter+obs.R2_Cd_w0_sr5_inter$alpha_inter_standard_error)

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==3),"Tu_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[1,3]+obs.R3_Cd_w0$alpha_matrix_standard_error[1,3], obs.R3_Cd_w0$alpha_matrix[2,3]+obs.R3_Cd_w0$alpha_matrix_standard_error[2,3],obs.R3_Cd_w0$alpha_matrix[1,4]+obs.R3_Cd_w0$alpha_matrix_standard_error[1,4], obs.R3_Cd_w0$alpha_matrix[2,4]+obs.R3_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==3),"Te_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[3,1]+obs.R3_Cd_w0$alpha_matrix_standard_error[3,1], obs.R3_Cd_w0$alpha_matrix[3,2]+obs.R3_Cd_w0$alpha_matrix_standard_error[3,2],obs.R3_Cd_w0$alpha_matrix[4,1]+obs.R3_Cd_w0$alpha_matrix_standard_error[4,1], obs.R3_Cd_w0$alpha_matrix[4,2]+obs.R3_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==4),"Tu_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[1,3]+obs.R4_Cd_w0$alpha_matrix_standard_error[1,3], obs.R4_Cd_w0$alpha_matrix[2,3]+obs.R4_Cd_w0$alpha_matrix_standard_error[2,3],obs.R4_Cd_w0$alpha_matrix[1,4]+obs.R4_Cd_w0$alpha_matrix_standard_error[1,4], obs.R4_Cd_w0$alpha_matrix[2,4]+obs.R4_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==4),"Te_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[3,1]+obs.R4_Cd_w0$alpha_matrix_standard_error[3,1], obs.R4_Cd_w0$alpha_matrix[3,2]+obs.R4_Cd_w0$alpha_matrix_standard_error[3,2],obs.R4_Cd_w0$alpha_matrix[4,1]+obs.R4_Cd_w0$alpha_matrix_standard_error[4,1], obs.R4_Cd_w0$alpha_matrix[4,2]+obs.R4_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==5),"Tu_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[1,3]+obs.R5_Cd_w0$alpha_matrix_standard_error[1,3], obs.R5_Cd_w0$alpha_matrix[2,3]+obs.R5_Cd_w0$alpha_matrix_standard_error[2,3],obs.R5_Cd_w0$alpha_matrix[1,4]+obs.R5_Cd_w0$alpha_matrix_standard_error[1,4], obs.R5_Cd_w0$alpha_matrix[2,4]+obs.R5_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==5),"Te_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[3,1]+obs.R5_Cd_w0$alpha_matrix_standard_error[3,1], obs.R5_Cd_w0$alpha_matrix[3,2]+obs.R5_Cd_w0$alpha_matrix_standard_error[3,2],obs.R5_Cd_w0$alpha_matrix[4,1]+obs.R5_Cd_w0$alpha_matrix_standard_error[4,1], obs.R5_Cd_w0$alpha_matrix[4,2]+obs.R5_Cd_w0$alpha_matrix_standard_error[4,2])


##### joining data frame
param_all_w0<-as.data.frame(rbind(cxr_param_w0, cxr_param_w0C))

param_all_w0_lower<-as.data.frame(rbind(cxr_param_w0_lower, cxr_param_w0C_lower))
param_all_w0_upper<-as.data.frame(rbind(cxr_param_w0_upper, cxr_param_w0C_upper))


param_all_w0_lower
param_all_w0_upper

#write.csv(param_all_w0, "./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal.csv")
#write.csv(param_all_w0_upper, "./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal_upper.csv")
#write.csv(param_all_w0_lower, "./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal_lower.csv")


#### Importing parameters

param_all_w0<-read.csv("./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal.csv")
param_all_w0_upper<-read.csv("./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal_upper.csv")
param_all_w0_lower<-read.csv("./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal_lower.csv")

param_all_w0<-param_all_w0[,-1]
param_all_w0_upper<-param_all_w0_upper[,-1]
param_all_w0_lower<-param_all_w0_lower[,-1]

str(param_all_w0)
str(param_all_w0_upper)
str(param_all_w0_lower)

##### Plotting data

#### Predicting densities
density_aux<-seq(0, 10, by=(10/100))

pred_df_cxr<-as.data.frame(expand_grid(Density=density_aux, Tu_Regime=c("SR1","SR2"), Te_Regime=c("SR4","SR5"), Replicate=c(1:5), Environment=c("N", "Cd")))

pred_df_cxr$Tu_mean_intra<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Tu_mean_inter<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr$Density[x])
  
  pred
})


pred_df_cxr$Tu_intra_L<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Tu_inter_L<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Tu_intra_U<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Tu_inter_U<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Te_mean_intra<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Te_mean_inter<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_w0, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Te_intra_L<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Te_inter_L<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_w0_lower, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Te_intra_U<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr$Density[x])
  
  pred
})

pred_df_cxr$Te_inter_U<-sapply(c(1:length(pred_df_cxr[,1])), function(x){
  alpha_i<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Tu_Regime==pred_df_cxr$Tu_Regime[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_w0_upper, Environment==pred_df_cxr$Environment[x] & Te_Regime==pred_df_cxr$Te_Regime[x] & Replicate==pred_df_cxr$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr$Density[x])
  
  pred
})

# Removing Tu evolved replicate 2 because there is no data
pred_df_cxr<-pred_df_cxr[-which(pred_df_cxr$Tu_Regime=="SR2" & pred_df_cxr$Replicate==2),]



# Transforming everything bellow 0 into 0 for the lower interval

pred_df_cxr$Te_inter_L[which(pred_df_cxr$Te_inter_L<0)]<-0
pred_df_cxr$Te_intra_L[which(pred_df_cxr$Te_intra_L<0)]<-0
pred_df_cxr$Tu_inter_L[which(pred_df_cxr$Tu_inter_L<0)]<-0
pred_df_cxr$Tu_intra_L[which(pred_df_cxr$Tu_intra_L<0)]<-0

##### Predicted vs observed

rk_func<- function(lambda, alpha_ii, alpha_ij, dens_i, dens_j, ...){
  gr<-lambda*exp(-alpha_ii*dens_i - alpha_ij*dens_j)
  
  return(gr)
}

red_ca<-ca[,c("Env", "Rep", "FocalSR", "CompSR", "Dens", "Type", "TeFemales", "TuFemales", "GrowthRateOA")]

red_ca

red_ca$Dens_Focal<-sapply(c(1:length(red_ca[,1])), function(x){
  if(red_ca$Type[x]=="INTRA"){
    a<-red_ca$Dens[x]-1
  }else if(red_ca$Type[x]=="INTER"){
    a<-1
  }
  
  a
})

red_ca$Dens_Comp<-sapply(c(1:length(red_ca[,1])), function(x){
  if(red_ca$Type[x]=="INTRA"){
    a<-0
  }else if(red_ca$Type[x]=="INTER"){
    a<-red_ca$Dens[x]-1
  }
  
  a
})

red_ca$Focal<-mapvalues(red_ca$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4", "SR5"))
red_ca$Comp<-mapvalues(red_ca$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4", "SR5"))

red_ca$pred<-sapply(c(1:length(red_ca[,1])), function(x){
  
  if(red_ca$Focal[x]=="SR1" | red_ca$Focal[x]=="SR2"){
    aux_data<-subset(param_all_w0, Environment==red_ca$Env[x] & Replicate== red_ca$Rep[x] & as.character(Tu_Regime)==red_ca$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Tu_lambda[1], alpha_ii =aux_data$Tu_intra[1], alpha_ij = aux_data$Tu_inter[1], dens_i = red_ca$Dens_Focal[x], dens_j =  red_ca$Dens_Comp[x])
    
  }else if(red_ca$Focal[x]=="SR4" | red_ca$Focal[x]=="SR5"){
    aux_data<-subset(param_all_w0, Environment==red_ca$Env[x] & Replicate== red_ca$Rep[x] & as.character(Te_Regime)==red_ca$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Te_lambda[1], alpha_ii =aux_data$Te_intra[1], alpha_ij = aux_data$Te_inter[1], dens_i = red_ca$Dens_Focal[x], dens_j =  red_ca$Dens_Comp[x])
  }
  
  aux_pred
})

red_ca$pred_L<-sapply(c(1:length(red_ca[,1])), function(x){
  
  if(red_ca$Focal[x]=="SR1" | red_ca$Focal[x]=="SR2"){
    aux_data<-subset(param_all_w0_lower, Environment==red_ca$Env[x] & Replicate== red_ca$Rep[x] & as.character(Tu_Regime)==red_ca$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Tu_lambda[1], alpha_ii =aux_data$Tu_intra[1], alpha_ij = aux_data$Tu_inter[1], dens_i = red_ca$Dens_Focal[x], dens_j =  red_ca$Dens_Comp[x])
    
  }else if(red_ca$Focal[x]=="SR4" | red_ca$Focal[x]=="SR5"){
    aux_data<-subset(param_all_w0_lower, Environment==red_ca$Env[x] & Replicate== red_ca$Rep[x] & as.character(Te_Regime)==red_ca$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Te_lambda[1], alpha_ii =aux_data$Te_intra[1], alpha_ij = aux_data$Te_inter[1], dens_i = red_ca$Dens_Focal[x], dens_j =  red_ca$Dens_Comp[x])
  }
  
  aux_pred
})

red_ca$pred_U<-sapply(c(1:length(red_ca[,1])), function(x){
  
  if(red_ca$Focal[x]=="SR1" | red_ca$Focal[x]=="SR2"){
    aux_data<-subset(param_all_w0_upper, Environment==red_ca$Env[x] & Replicate== red_ca$Rep[x] & as.character(Tu_Regime)==red_ca$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Tu_lambda[1], alpha_ii =aux_data$Tu_intra[1], alpha_ij = aux_data$Tu_inter[1], dens_i = red_ca$Dens_Focal[x], dens_j =  red_ca$Dens_Comp[x])
    
  }else if(red_ca$Focal[x]=="SR4" | red_ca$Focal[x]=="SR5"){
    aux_data<-subset(param_all_w0_upper, Environment==red_ca$Env[x] & Replicate== red_ca$Rep[x] & as.character(Te_Regime)==red_ca$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Te_lambda[1], alpha_ii =aux_data$Te_intra[1], alpha_ij = aux_data$Te_inter[1], dens_i = red_ca$Dens_Focal[x], dens_j =  red_ca$Dens_Comp[x])
  }
  
  aux_pred
})

red_ca$Replicate<-red_ca$Rep
str(red_ca)

### B - CXR lambda fixed

#cxr accepts a data frame with a first column called fitness with positive values and numeric columns with number of individuals. Each row is one individual. For multiple species the easier is to create a list, each with a data frame that has in the first column number of individuals produced and then the number of neighbours
#this case we transformed all 0s into 1 (so that the log is 0) For that we need to add +1 to all data so that the variance is not changed

##### normal
dir.create("./Analyses/MethodComparison/cxr_lambda_fixed_log", showWarnings = FALSE)

# modifying data frame to fit the type of setup that is need for CXR
CXR_B_N<-subset(ca, Env=="N")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "TuFemales")]

CXR_B_N$Focal<-mapvalues(CXR_B_N$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
CXR_B_N$CompSR2<-mapvalues(CXR_B_N$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

CXR_B_N$Comp<-sapply(c(1:length(CXR_B_N[,1])), function(x){
  if(is.na(CXR_B_N$CompSR2[x])){
    a<- CXR_B_N$Focal[x]
  }else{
    a<-CXR_B_N$CompSR2[x]
  }
  
  a
})

aux<-data.frame(SR1=rep(0, length(CXR_B_N[,1])), SR2=rep(0, length(CXR_B_N[,1])), SR4=rep(0, length(CXR_B_N[,1])), SR5=rep(0, length(CXR_B_N[,1])))

for(i in 1:length(CXR_B_N[,1])){
  #coluna onde por focais
  colunaF<-which(colnames(aux)==CXR_B_N$Focal[i])
  #coluna onde por competidors
  colunaC<-which(colnames(aux)==CXR_B_N$Comp[i])
  
  #if its the same regime
  if(CXR_B_N$Focal[i]==CXR_B_N$Comp[i] & CXR_B_N$Dens[i]==1){
    aux[i,colunaF]<-CXR_B_N$Dens[i]-1
    
  }else if(CXR_B_N$Focal[i]==CXR_B_N$Comp[i]){
    aux[i,colunaF]<-CXR_B_N$Dens[i]-1
  }else{ #if it is heterospecific then its -1 for the competitors (because of the focal) and its one for the focal
    aux[i,colunaC]<-CXR_B_N$Dens[i]-1
    aux[i, colunaF]<-1
  }
  
}

CXR_B_N<-cbind(CXR_B_N, aux)

CXR_B_N$fitness<-sapply(c(1:length(CXR_B_N[,1])), function(x){
  colF<-which(colnames(CXR_B_N)==CXR_B_N$Focal[x])
  
  if(CXR_B_N$Focal[x]=="SR1"){
    a<-CXR_B_N$TuFemales[x]/CXR_B_N$SR1[x]
  } else if(CXR_B_N$Focal[x]=="SR2"){
    a<-CXR_B_N$TuFemales[x]/CXR_B_N$SR2[x]
  } else if(CXR_B_N$Focal[x]=="SR4"){
    a<-CXR_B_N$TeFemales[x]/CXR_B_N$SR4[x]
  } else if(CXR_B_N$Focal[x]=="SR5"){
    a<-CXR_B_N$TeFemales[x]/CXR_B_N$SR5[x]
  }
  
  a
})

#removing rows for which there is no data for fitness
CXR_B_N<-CXR_B_N[-which(is.na(CXR_B_N$fitness)),]

# adding +1 to all data
#CXR_B_N$fitness<-CXR_B_N$fitness+1

CXR_B_N[which(CXR_B_N$fitness=="-Inf" | CXR_B_N$fitness=="Inf"),"fitness"]<-0


# all data gets +1 because of the 0 problem
CXR_B_N$fitness<-CXR_B_N$fitness+1

# vector that tells which are the selection regimes, the columns have to have the same name
my.reg <- c("SR1", "SR2","SR4","SR5")

# Do list per replicate and environment
R1<-list(SR1= subset(CXR_B_N, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(CXR_B_N, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(CXR_B_N, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(CXR_B_N, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R2<-list(SR1= subset(CXR_B_N, Rep==2 & Focal=="SR1")[,c("fitness", "SR1", "SR4", "SR5")], SR4= subset(CXR_B_N, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR4", "SR5")], SR5= subset(CXR_B_N, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR4", "SR5")])

R3<-list(SR1= subset(CXR_B_N, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(CXR_B_N, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(CXR_B_N, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(CXR_B_N, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R4<-list(SR1= subset(CXR_B_N, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(CXR_B_N, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(CXR_B_N, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(CXR_B_N, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R5<-list(SR1= subset(CXR_B_N, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(CXR_B_N, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(CXR_B_N, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(CXR_B_N, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

mean_dens1<-data.frame(SR=c(rep(1,10), rep(2,8), rep(4,10),rep(5,10)), Env=c(rep("N", 5),rep("Cd", 5), rep("N", 4),rep("Cd", 4),rep("N", 5),rep("Cd", 5),rep("N", 5),rep("Cd", 5)), Rep=c(rep(c(1,2,3,4,5),2),rep(c(1,3,4,5),2),rep(c(1,2,3,4,5),2),rep(c(1,2,3,4,5),2)))


#since in the model we use the log of data +1, here we also have to use the +1 to estimate the lambda
mean_dens1$lambda<-sapply(c(1:length(mean_dens1[,1])), function(x){
  mean(subset(ca, FocalSR==mean_dens1$SR[x] & Dens==1 & Env==mean_dens1$Env[x] & Rep==mean_dens1$Rep[x] )$GrowthRateOA+1, na.rm=TRUE)
})


mean_dens1$sd_lambda<-sapply(c(1:length(mean_dens1[,1])), function(x){
  sd(subset(ca, FocalSR==mean_dens1$SR[x] & Dens==1 & Env==mean_dens1$Env[x] & Rep==mean_dens1$Rep[x])$GrowthRateOA+1, na.rm=TRUE)
})

mean_dens1$sd_lambda[which(is.na(mean_dens1$sd_lambda))]<-0.01
mean_dens1$sd_lambda[which(mean_dens1$sd_lambda==0)]<-0.01

#### lambda

fixed_terms_1N <- list(list(lambda = subset(mean_dens1, Rep==1 & Env=="N" & SR==1)$lambda ), # focal sp 1
                       list(lambda = subset(mean_dens1, Rep==1 & Env=="N" & SR==2)$lambda), # focal sp 2
                       list(lambda = subset(mean_dens1, Rep==1 & Env=="N" & SR==4)$lambda),
                       list(lambda= subset(mean_dens1, Rep==1 & Env=="N" & SR==5)$lambda))

fixed_terms_2N <- list(list(lambda = subset(mean_dens1, Rep==2 & Env=="N" & SR==1)$lambda ), # focal sp 1
                       list(lambda = subset(mean_dens1, Rep==2 & Env=="N" & SR==4)$lambda),
                       list(lambda= subset(mean_dens1, Rep==2 & Env=="N" & SR==5)$lambda))

fixed_terms_3N <- list(list(lambda = subset(mean_dens1, Rep==3 & Env=="N" & SR==1)$lambda ), # focal sp 1
                       list(lambda = subset(mean_dens1, Rep==3 & Env=="N" & SR==2)$lambda), # focal sp 2
                       list(lambda = subset(mean_dens1, Rep==3 & Env=="N" & SR==4)$lambda),
                       list(lambda= subset(mean_dens1, Rep==3 & Env=="N" & SR==5)$lambda))

fixed_terms_4N <- list(list(lambda = subset(mean_dens1, Rep==4 & Env=="N" & SR==1)$lambda ), # focal sp 1
                       list(lambda = subset(mean_dens1, Rep==4 & Env=="N" & SR==2)$lambda), # focal sp 2
                       list(lambda = subset(mean_dens1, Rep==4 & Env=="N" & SR==4)$lambda),
                       list(lambda= subset(mean_dens1, Rep==4 & Env=="N" & SR==5)$lambda))

fixed_terms_5N <- list(list(lambda = subset(mean_dens1, Rep==5 & Env=="N" & SR==1)$lambda ), # focal sp 1
                       list(lambda = subset(mean_dens1, Rep==5 & Env=="N" & SR==2)$lambda), # focal sp 2
                       list(lambda = subset(mean_dens1, Rep==5 & Env=="N" & SR==4)$lambda),
                       list(lambda= subset(mean_dens1, Rep==5 & Env=="N" & SR==5)$lambda))


cxr_B.R1_w0<-cxr_pm_multifit(data = R1,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_intra = 0.1,
                                                   alpha_inter = 0.1),
                             fixed_terms = fixed_terms_1N,
                             # no standard errors
                             bootstrap_samples = 200)


cxr_B.R2_w0<-cxr_pm_multifit(data = R2,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_intra = 0.1,
                                                   alpha_inter = 0.1),
                             fixed_terms = fixed_terms_2N,
                             # no standard errors
                             bootstrap_samples = 200)

cxr_B.R3_w0<-cxr_pm_multifit(data = R3,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_intra = 0.1,
                                                   alpha_inter = 0.1),
                             fixed_terms = fixed_terms_3N,
                             # no standard errors
                             bootstrap_samples = 200)

cxr_B.R4_w0<-cxr_pm_multifit(data = R4,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_intra = 0.1,
                                                   alpha_inter = 0.1),
                             fixed_terms = fixed_terms_3N,
                             # no standard errors
                             bootstrap_samples = 200)

cxr_B.R5_w0<-cxr_pm_multifit(data = R5,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_intra = 0.1,
                                                   alpha_inter = 0.1),
                             fixed_terms = fixed_terms_5N,
                             # no standard errors
                             bootstrap_samples = 200)

summary(cxr_B.R1_w0)
summary(cxr_B.R2_w0)
summary(cxr_B.R3_w0)
summary(cxr_B.R4_w0)
summary(cxr_B.R5_w0)

#ab<-abundance_projection(cxr_B.R1_w0, timesteps = 1, initial_abundances = c(3,3,3,3))



# rows in the alpha element of the returning list correspond to species i and columns to species j for each αij coefficient.

###### data table summary



cxr_param_B<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_B$Tu_lambda<-0
cxr_param_B$Te_lambda<-0
cxr_param_B$Tu_intra<-0
cxr_param_B$Te_intra<-0
cxr_param_B$Tu_inter<-0
cxr_param_B$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_B<-cxr_param_B[-which(cxr_param_B$Replicate==2 & cxr_param_B$Tu_Regime=="SR2"),]


cxr_param_B[which(cxr_param_B$Replicate==1),"Tu_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[1]]$lambda,cxr_B.R1_w0$fixed_terms[[2]]$lambda)
cxr_param_B[which(cxr_param_B$Replicate==1),"Te_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[3]]$lambda,cxr_B.R1_w0$fixed_terms[[3]]$lambda, cxr_B.R1_w0$fixed_terms[[4]]$lambda,cxr_B.R1_w0$fixed_terms[[4]]$lambda)

cxr_param_B[which(cxr_param_B$Replicate==2),"Tu_lambda"]<-c(cxr_B.R2_w0$fixed_terms[[1]]$lambda,cxr_B.R2_w0$fixed_terms[[1]]$lambda)
cxr_param_B[which(cxr_param_B$Replicate==2),"Te_lambda"]<-c(cxr_B.R2_w0$fixed_terms[[2]]$lambda,cxr_B.R2_w0$fixed_terms[[3]]$lambda)

cxr_param_B[which(cxr_param_B$Replicate==3),"Tu_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[1]]$lambda,cxr_B.R3_w0$fixed_terms[[2]]$lambda)
cxr_param_B[which(cxr_param_B$Replicate==3),"Te_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[3]]$lambda,cxr_B.R3_w0$fixed_terms[[3]]$lambda, cxr_B.R3_w0$fixed_terms[[4]]$lambda,cxr_B.R3_w0$fixed_terms[[4]]$lambda)

cxr_param_B[which(cxr_param_B$Replicate==4),"Tu_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[1]]$lambda,cxr_B.R4_w0$fixed_terms[[2]]$lambda)
cxr_param_B[which(cxr_param_B$Replicate==4),"Te_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[3]]$lambda,cxr_B.R4_w0$fixed_terms[[3]]$lambda, cxr_B.R4_w0$fixed_terms[[4]]$lambda,cxr_B.R4_w0$fixed_terms[[4]]$lambda)

cxr_param_B[which(cxr_param_B$Replicate==5),"Tu_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[1]]$lambda,cxr_B.R5_w0$fixed_terms[[2]]$lambda)
cxr_param_B[which(cxr_param_B$Replicate==5),"Te_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[3]]$lambda,cxr_B.R5_w0$fixed_terms[[3]]$lambda, cxr_B.R5_w0$fixed_terms[[4]]$lambda,cxr_B.R5_w0$fixed_terms[[4]]$lambda)


cxr_param_B[which(cxr_param_B$Replicate==1),"Tu_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[1,1], cxr_B.R1_w0$alpha_matrix[2,2]), 2)
cxr_param_B[which(cxr_param_B$Replicate==1),"Te_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[3,3], cxr_B.R1_w0$alpha_matrix[4,4]), each=2)

cxr_param_B[which(cxr_param_B$Replicate==2),"Tu_intra"]<-rep(c(cxr_B.R2_w0$alpha_matrix[1,1]), 2)
cxr_param_B[which(cxr_param_B$Replicate==2),"Te_intra"]<-rep(c(cxr_B.R2_w0$alpha_matrix[2,2], cxr_B.R2_w0$alpha_matrix[3,3]))

cxr_param_B[which(cxr_param_B$Replicate==3),"Tu_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[1,1], cxr_B.R3_w0$alpha_matrix[2,2]), 2)
cxr_param_B[which(cxr_param_B$Replicate==3),"Te_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[3,3], cxr_B.R3_w0$alpha_matrix[4,4]), each=2)

cxr_param_B[which(cxr_param_B$Replicate==4),"Tu_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[1,1], cxr_B.R4_w0$alpha_matrix[2,2]), 2)
cxr_param_B[which(cxr_param_B$Replicate==4),"Te_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[3,3], cxr_B.R4_w0$alpha_matrix[4,4]), each=2)

cxr_param_B[which(cxr_param_B$Replicate==5),"Tu_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[1,1], cxr_B.R5_w0$alpha_matrix[2,2]), 2)
cxr_param_B[which(cxr_param_B$Replicate==5),"Te_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[3,3], cxr_B.R5_w0$alpha_matrix[4,4]), each=2)


cxr_param_B[which(cxr_param_B$Replicate==1),"Tu_inter"]<-c(cxr_B.R1_w0$alpha_matrix[1,3], cxr_B.R1_w0$alpha_matrix[2,3],cxr_B.R1_w0$alpha_matrix[1,4], cxr_B.R1_w0$alpha_matrix[2,4])
cxr_param_B[which(cxr_param_B$Replicate==1),"Te_inter"]<-c(cxr_B.R1_w0$alpha_matrix[3,1], cxr_B.R1_w0$alpha_matrix[3,2],cxr_B.R1_w0$alpha_matrix[4,1], cxr_B.R1_w0$alpha_matrix[4,2])

cxr_param_B[which(cxr_param_B$Replicate==2),"Tu_inter"]<-c(cxr_B.R2_w0$alpha_matrix[1,2], cxr_B.R2_w0$alpha_matrix[1,3])
cxr_param_B[which(cxr_param_B$Replicate==2),"Te_inter"]<-c(cxr_B.R2_w0$alpha_matrix[2,1],cxr_B.R2_w0$alpha_matrix[3,1])

cxr_param_B[which(cxr_param_B$Replicate==3),"Tu_inter"]<-c(cxr_B.R3_w0$alpha_matrix[1,3], cxr_B.R3_w0$alpha_matrix[2,3],cxr_B.R3_w0$alpha_matrix[1,4], cxr_B.R3_w0$alpha_matrix[2,4])
cxr_param_B[which(cxr_param_B$Replicate==3),"Te_inter"]<-c(cxr_B.R3_w0$alpha_matrix[3,1], cxr_B.R3_w0$alpha_matrix[3,2],cxr_B.R3_w0$alpha_matrix[4,1], cxr_B.R3_w0$alpha_matrix[4,2])

cxr_param_B[which(cxr_param_B$Replicate==4),"Tu_inter"]<-c(cxr_B.R4_w0$alpha_matrix[1,3], cxr_B.R4_w0$alpha_matrix[2,3],cxr_B.R4_w0$alpha_matrix[1,4], cxr_B.R4_w0$alpha_matrix[2,4])
cxr_param_B[which(cxr_param_B$Replicate==4),"Te_inter"]<-c(cxr_B.R4_w0$alpha_matrix[3,1], cxr_B.R4_w0$alpha_matrix[3,2],cxr_B.R4_w0$alpha_matrix[4,1], cxr_B.R4_w0$alpha_matrix[4,2])

cxr_param_B[which(cxr_param_B$Replicate==5),"Tu_inter"]<-c(cxr_B.R5_w0$alpha_matrix[1,3], cxr_B.R5_w0$alpha_matrix[2,3],cxr_B.R5_w0$alpha_matrix[1,4], cxr_B.R5_w0$alpha_matrix[2,4])
cxr_param_B[which(cxr_param_B$Replicate==5),"Te_inter"]<-c(cxr_B.R5_w0$alpha_matrix[3,1], cxr_B.R5_w0$alpha_matrix[3,2],cxr_B.R5_w0$alpha_matrix[4,1], cxr_B.R5_w0$alpha_matrix[4,2])

### Lower

cxr_param_B_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_B_lower$Tu_lambda<-0
cxr_param_B_lower$Te_lambda<-0
cxr_param_B_lower$Tu_intra<-0
cxr_param_B_lower$Te_intra<-0
cxr_param_B_lower$Tu_inter<-0
cxr_param_B_lower$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_B_lower<-cxr_param_B_lower[-which(cxr_param_B_lower$Replicate==2 & cxr_param_B_lower$Tu_Regime=="SR2"),]

#Since the error comes directly from the data we need to create some lists with that information
sd_1N <- list(list(lambda = subset(mean_dens1, Rep==1 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==1 & Env=="N" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==1 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==1 & Env=="N" & SR==5)$sd_lambda))

sd_2N <- list(list(lambda = subset(mean_dens1, Rep==2 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==2 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==2 & Env=="N" & SR==5)$sd_lambda))

sd_3N <- list(list(lambda = subset(mean_dens1, Rep==3 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==3 & Env=="N" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==3 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==3 & Env=="N" & SR==5)$sd_lambda))

sd_4N <- list(list(lambda = subset(mean_dens1, Rep==4 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==4 & Env=="N" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==4 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==4 & Env=="N" & SR==5)$sd_lambda))

sd_5N <- list(list(lambda = subset(mean_dens1, Rep==5 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==5 & Env=="N" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==5 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==5 & Env=="N" & SR==5)$sd_lambda))

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Tu_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[1]]$lambda-sd_1N[[1]]$lambda,cxr_B.R1_w0$fixed_terms[[2]]$lambda-sd_1N[[2]]$lambda)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Te_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[3]]$lambda-sd_1N[[3]]$lambda,cxr_B.R1_w0$fixed_terms[[3]]$lambda-sd_1N[[3]]$lambda, cxr_B.R1_w0$fixed_terms[[4]]$lambda-sd_1N[[4]]$lambda,cxr_B.R1_w0$fixed_terms[[4]]$lambda-sd_1N[[4]]$lambda)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Tu_lambda"]<-c(cxr_B.R2_w0$fixed_terms[[1]]$lambda-sd_2N[[1]]$lambda,cxr_B.R2_w0$fixed_terms[[1]]$lambda-sd_2N[[1]]$lambda)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Te_lambda"]<-c(cxr_B.R2_w0$fixed_terms[[2]]$lambda-sd_2N[[2]]$lambda,cxr_B.R2_w0$fixed_terms[[3]]$lambda-sd_2N[[3]]$lambda)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Tu_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[1]]$lambda-sd_3N[[1]]$lambda,cxr_B.R3_w0$fixed_terms[[2]]$lambda-sd_3N[[2]]$lambda)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Te_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[3]]$lambda-sd_3N[[3]]$lambda,cxr_B.R3_w0$fixed_terms[[3]]$lambda-sd_3N[[3]]$lambda, cxr_B.R3_w0$fixed_terms[[4]]$lambda-sd_3N[[4]]$lambda,cxr_B.R3_w0$fixed_terms[[4]]$lambda-sd_3N[[4]]$lambda)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Tu_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[1]]$lambda-sd_4N[[1]]$lambda,cxr_B.R4_w0$fixed_terms[[2]]$lambda-sd_4N[[2]]$lambda)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Te_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[3]]$lambda-sd_4N[[3]]$lambda,cxr_B.R4_w0$fixed_terms[[3]]$lambda-sd_4N[[3]]$lambda, cxr_B.R4_w0$fixed_terms[[4]]$lambda-sd_4N[[4]]$lambda,cxr_B.R4_w0$fixed_terms[[4]]$lambda-sd_4N[[4]]$lambda)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Tu_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[1]]$lambda-sd_5N[[1]]$lambda,cxr_B.R5_w0$fixed_terms[[2]]$lambda-sd_5N[[2]]$lambda)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Te_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[3]]$lambda-sd_5N[[3]]$lambda,cxr_B.R5_w0$fixed_terms[[3]]$lambda-sd_5N[[3]]$lambda, cxr_B.R5_w0$fixed_terms[[4]]$lambda-sd_5N[[4]]$lambda,cxr_B.R5_w0$fixed_terms[[4]]$lambda-sd_5N[[4]]$lambda)


cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Tu_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[1,1]-cxr_B.R1_w0$alpha_matrix_standard_error[1,1], cxr_B.R1_w0$alpha_matrix[2,2]-cxr_B.R1_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Te_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[3,3]-cxr_B.R1_w0$alpha_matrix_standard_error[3,3], cxr_B.R1_w0$alpha_matrix[4,4]-cxr_B.R1_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Tu_intra"]<-rep(c(cxr_B.R2_w0$alpha_matrix[1,1]-cxr_B.R2_w0$alpha_matrix_standard_error[1,1]), 2)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Te_intra"]<-c(cxr_B.R2_w0$alpha_matrix[2,2]-cxr_B.R2_w0$alpha_matrix_standard_error[2,2], cxr_B.R2_w0$alpha_matrix[3,3]-cxr_B.R2_w0$alpha_matrix_standard_error[3,3])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Tu_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[1,1]-cxr_B.R3_w0$alpha_matrix_standard_error[1,1], cxr_B.R3_w0$alpha_matrix[2,2]-cxr_B.R3_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Te_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[3,3]-cxr_B.R3_w0$alpha_matrix_standard_error[3,3], cxr_B.R3_w0$alpha_matrix[4,4]-cxr_B.R3_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Tu_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[1,1]-cxr_B.R4_w0$alpha_matrix_standard_error[1,1], cxr_B.R4_w0$alpha_matrix[2,2]-cxr_B.R4_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Te_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[3,3]-cxr_B.R4_w0$alpha_matrix_standard_error[3,3], cxr_B.R4_w0$alpha_matrix[4,4]-cxr_B.R4_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Tu_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[1,1]-cxr_B.R5_w0$alpha_matrix_standard_error[1,1], cxr_B.R5_w0$alpha_matrix[2,2]-cxr_B.R5_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Te_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[3,3]-cxr_B.R5_w0$alpha_matrix_standard_error[3,3], cxr_B.R5_w0$alpha_matrix[4,4]-cxr_B.R5_w0$alpha_matrix_standard_error[4,4]), each=2)


cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Tu_inter"]<-c(cxr_B.R1_w0$alpha_matrix[1,3]-cxr_B.R1_w0$alpha_matrix_standard_error[1,3], cxr_B.R1_w0$alpha_matrix[2,3]-cxr_B.R1_w0$alpha_matrix_standard_error[2,3],cxr_B.R1_w0$alpha_matrix[1,4]-cxr_B.R1_w0$alpha_matrix_standard_error[1,4], cxr_B.R1_w0$alpha_matrix[2,4]-cxr_B.R1_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Te_inter"]<-c(cxr_B.R1_w0$alpha_matrix[3,1]-cxr_B.R1_w0$alpha_matrix_standard_error[3,1], cxr_B.R1_w0$alpha_matrix[3,2]-cxr_B.R1_w0$alpha_matrix_standard_error[3,2],cxr_B.R1_w0$alpha_matrix[4,1]-cxr_B.R1_w0$alpha_matrix_standard_error[4,1], cxr_B.R1_w0$alpha_matrix[4,2]-cxr_B.R1_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Tu_inter"]<-c(cxr_B.R2_w0$alpha_matrix[1,2]-cxr_B.R2_w0$alpha_matrix_standard_error[1,2], cxr_B.R2_w0$alpha_matrix[1,3]-cxr_B.R2_w0$alpha_matrix_standard_error[1,3])
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Te_inter"]<-c(cxr_B.R2_w0$alpha_matrix[2,1]-cxr_B.R2_w0$alpha_matrix_standard_error[2,1], cxr_B.R2_w0$alpha_matrix[3,1]-cxr_B.R2_w0$alpha_matrix_standard_error[3,1])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Tu_inter"]<-c(cxr_B.R3_w0$alpha_matrix[1,3]-cxr_B.R3_w0$alpha_matrix_standard_error[1,3], cxr_B.R3_w0$alpha_matrix[2,3]-cxr_B.R3_w0$alpha_matrix_standard_error[2,3],cxr_B.R3_w0$alpha_matrix[1,4]-cxr_B.R3_w0$alpha_matrix_standard_error[1,4], cxr_B.R3_w0$alpha_matrix[2,4]-cxr_B.R3_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Te_inter"]<-c(cxr_B.R3_w0$alpha_matrix[3,1]-cxr_B.R3_w0$alpha_matrix_standard_error[3,1], cxr_B.R3_w0$alpha_matrix[3,2]-cxr_B.R3_w0$alpha_matrix_standard_error[3,2],cxr_B.R3_w0$alpha_matrix[4,1]-cxr_B.R3_w0$alpha_matrix_standard_error[4,1], cxr_B.R3_w0$alpha_matrix[4,2]-cxr_B.R3_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Tu_inter"]<-c(cxr_B.R4_w0$alpha_matrix[1,3]-cxr_B.R4_w0$alpha_matrix_standard_error[1,3], cxr_B.R4_w0$alpha_matrix[2,3]-cxr_B.R4_w0$alpha_matrix_standard_error[2,3],cxr_B.R4_w0$alpha_matrix[1,4]-cxr_B.R4_w0$alpha_matrix_standard_error[1,4], cxr_B.R4_w0$alpha_matrix[2,4]-cxr_B.R4_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Te_inter"]<-c(cxr_B.R4_w0$alpha_matrix[3,1]-cxr_B.R4_w0$alpha_matrix_standard_error[3,1], cxr_B.R4_w0$alpha_matrix[3,2]-cxr_B.R4_w0$alpha_matrix_standard_error[3,2],cxr_B.R4_w0$alpha_matrix[4,1]-cxr_B.R4_w0$alpha_matrix_standard_error[4,1], cxr_B.R4_w0$alpha_matrix[4,2]-cxr_B.R4_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Tu_inter"]<-c(cxr_B.R5_w0$alpha_matrix[1,3]-cxr_B.R5_w0$alpha_matrix_standard_error[1,3], cxr_B.R5_w0$alpha_matrix[2,3]-cxr_B.R5_w0$alpha_matrix_standard_error[2,3],cxr_B.R5_w0$alpha_matrix[1,4]-cxr_B.R5_w0$alpha_matrix_standard_error[1,4], cxr_B.R5_w0$alpha_matrix[2,4]-cxr_B.R5_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Te_inter"]<-c(cxr_B.R5_w0$alpha_matrix[3,1]-cxr_B.R5_w0$alpha_matrix_standard_error[3,1], cxr_B.R5_w0$alpha_matrix[3,2]-cxr_B.R5_w0$alpha_matrix_standard_error[3,2],cxr_B.R5_w0$alpha_matrix[4,1]-cxr_B.R5_w0$alpha_matrix_standard_error[4,1], cxr_B.R5_w0$alpha_matrix[4,2]-cxr_B.R5_w0$alpha_matrix_standard_error[4,2])

### upper

cxr_param_B_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_B_upper$Tu_lambda<-0
cxr_param_B_upper$Te_lambda<-0
cxr_param_B_upper$Tu_intra<-0
cxr_param_B_upper$Te_intra<-0
cxr_param_B_upper$Tu_inter<-0
cxr_param_B_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_B_upper<-cxr_param_B_upper[-which(cxr_param_B_upper$Replicate==2 & cxr_param_B_upper$Tu_Regime=="SR2"),]


cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Tu_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[1]]$lambda+sd_1N[[1]]$lambda,cxr_B.R1_w0$fixed_terms[[2]]$lambda+sd_1N[[2]]$lambda)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Te_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[3]]$lambda+sd_1N[[3]]$lambda,cxr_B.R1_w0$fixed_terms[[3]]$lambda+sd_1N[[3]]$lambda, cxr_B.R1_w0$fixed_terms[[4]]$lambda+sd_1N[[4]]$lambda,cxr_B.R1_w0$fixed_terms[[4]]$lambda+sd_1N[[4]]$lambda)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Tu_lambda"]<-c(cxr_B.R2_w0$fixed_terms[[1]]$lambda+sd_2N[[1]]$lambda,cxr_B.R2_w0$fixed_terms[[1]]$lambda+sd_2N[[1]]$lambda)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Te_lambda"]<-c(cxr_B.R2_w0$fixed_terms[[2]]$lambda+sd_2N[[2]]$lambda,cxr_B.R2_w0$fixed_terms[[3]]$lambda+sd_2N[[3]]$lambda)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Tu_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[1]]$lambda+sd_3N[[1]]$lambda,cxr_B.R3_w0$fixed_terms[[2]]$lambda+sd_3N[[2]]$lambda)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Te_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[3]]$lambda+sd_3N[[3]]$lambda,cxr_B.R3_w0$fixed_terms[[3]]$lambda+sd_3N[[3]]$lambda, cxr_B.R3_w0$fixed_terms[[4]]$lambda+sd_3N[[4]]$lambda,cxr_B.R3_w0$fixed_terms[[4]]$lambda+sd_3N[[4]]$lambda)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Tu_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[1]]$lambda+sd_4N[[1]]$lambda,cxr_B.R4_w0$fixed_terms[[2]]$lambda+sd_4N[[2]]$lambda)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Te_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[3]]$lambda+sd_4N[[3]]$lambda,cxr_B.R4_w0$fixed_terms[[3]]$lambda+sd_4N[[3]]$lambda, cxr_B.R4_w0$fixed_terms[[4]]$lambda+sd_4N[[4]]$lambda,cxr_B.R4_w0$fixed_terms[[4]]$lambda+sd_4N[[4]]$lambda)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Tu_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[1]]$lambda+sd_5N[[1]]$lambda,cxr_B.R5_w0$fixed_terms[[2]]$lambda+sd_5N[[2]]$lambda)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Te_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[3]]$lambda+sd_5N[[3]]$lambda,cxr_B.R5_w0$fixed_terms[[3]]$lambda+sd_5N[[3]]$lambda, cxr_B.R5_w0$fixed_terms[[4]]$lambda+sd_5N[[4]]$lambda,cxr_B.R5_w0$fixed_terms[[4]]$lambda+sd_5N[[4]]$lambda)


cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Tu_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[1,1]+cxr_B.R1_w0$alpha_matrix_standard_error[1,1], cxr_B.R1_w0$alpha_matrix[2,2]+cxr_B.R1_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Te_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[3,3]+cxr_B.R1_w0$alpha_matrix_standard_error[3,3], cxr_B.R1_w0$alpha_matrix[4,4]+cxr_B.R1_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Tu_intra"]<-rep(c(cxr_B.R2_w0$alpha_matrix[1,1]+cxr_B.R2_w0$alpha_matrix_standard_error[1,1]), 2)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Te_intra"]<-c(cxr_B.R2_w0$alpha_matrix[2,2]+cxr_B.R2_w0$alpha_matrix_standard_error[2,2], cxr_B.R2_w0$alpha_matrix[3,3]+cxr_B.R2_w0$alpha_matrix_standard_error[3,3])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Tu_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[1,1]+cxr_B.R3_w0$alpha_matrix_standard_error[1,1], cxr_B.R3_w0$alpha_matrix[2,2]+cxr_B.R3_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Te_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[3,3]+cxr_B.R3_w0$alpha_matrix_standard_error[3,3], cxr_B.R3_w0$alpha_matrix[4,4]+cxr_B.R3_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Tu_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[1,1]+cxr_B.R4_w0$alpha_matrix_standard_error[1,1], cxr_B.R4_w0$alpha_matrix[2,2]+cxr_B.R4_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Te_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[3,3]+cxr_B.R4_w0$alpha_matrix_standard_error[3,3], cxr_B.R4_w0$alpha_matrix[4,4]+cxr_B.R4_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Tu_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[1,1]+cxr_B.R5_w0$alpha_matrix_standard_error[1,1], cxr_B.R5_w0$alpha_matrix[2,2]+cxr_B.R5_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Te_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[3,3]+cxr_B.R5_w0$alpha_matrix_standard_error[3,3], cxr_B.R5_w0$alpha_matrix[4,4]+cxr_B.R5_w0$alpha_matrix_standard_error[4,4]), each=2)


cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Tu_inter"]<-c(cxr_B.R1_w0$alpha_matrix[1,3]+cxr_B.R1_w0$alpha_matrix_standard_error[1,3], cxr_B.R1_w0$alpha_matrix[2,3]+cxr_B.R1_w0$alpha_matrix_standard_error[2,3],cxr_B.R1_w0$alpha_matrix[1,4]+cxr_B.R1_w0$alpha_matrix_standard_error[1,4], cxr_B.R1_w0$alpha_matrix[2,4]+cxr_B.R1_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Te_inter"]<-c(cxr_B.R1_w0$alpha_matrix[3,1]+cxr_B.R1_w0$alpha_matrix_standard_error[3,1], cxr_B.R1_w0$alpha_matrix[3,2]+cxr_B.R1_w0$alpha_matrix_standard_error[3,2],cxr_B.R1_w0$alpha_matrix[4,1]+cxr_B.R1_w0$alpha_matrix_standard_error[4,1], cxr_B.R1_w0$alpha_matrix[4,2]+cxr_B.R1_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Tu_inter"]<-c(cxr_B.R2_w0$alpha_matrix[1,2]+cxr_B.R2_w0$alpha_matrix_standard_error[1,2], cxr_B.R2_w0$alpha_matrix[1,3]+cxr_B.R2_w0$alpha_matrix_standard_error[1,3])
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Te_inter"]<-c(cxr_B.R2_w0$alpha_matrix[2,1]+cxr_B.R2_w0$alpha_matrix_standard_error[2,1], cxr_B.R2_w0$alpha_matrix[3,1]+cxr_B.R2_w0$alpha_matrix_standard_error[3,1])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Tu_inter"]<-c(cxr_B.R3_w0$alpha_matrix[1,3]+cxr_B.R3_w0$alpha_matrix_standard_error[1,3], cxr_B.R3_w0$alpha_matrix[2,3]+cxr_B.R3_w0$alpha_matrix_standard_error[2,3],cxr_B.R3_w0$alpha_matrix[1,4]+cxr_B.R3_w0$alpha_matrix_standard_error[1,4], cxr_B.R3_w0$alpha_matrix[2,4]+cxr_B.R3_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Te_inter"]<-c(cxr_B.R3_w0$alpha_matrix[3,1]+cxr_B.R3_w0$alpha_matrix_standard_error[3,1], cxr_B.R3_w0$alpha_matrix[3,2]+cxr_B.R3_w0$alpha_matrix_standard_error[3,2],cxr_B.R3_w0$alpha_matrix[4,1]+cxr_B.R3_w0$alpha_matrix_standard_error[4,1], cxr_B.R3_w0$alpha_matrix[4,2]+cxr_B.R3_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Tu_inter"]<-c(cxr_B.R4_w0$alpha_matrix[1,3]+cxr_B.R4_w0$alpha_matrix_standard_error[1,3], cxr_B.R4_w0$alpha_matrix[2,3]+cxr_B.R4_w0$alpha_matrix_standard_error[2,3],cxr_B.R4_w0$alpha_matrix[1,4]+cxr_B.R4_w0$alpha_matrix_standard_error[1,4], cxr_B.R4_w0$alpha_matrix[2,4]+cxr_B.R4_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Te_inter"]<-c(cxr_B.R4_w0$alpha_matrix[3,1]+cxr_B.R4_w0$alpha_matrix_standard_error[3,1], cxr_B.R4_w0$alpha_matrix[3,2]+cxr_B.R4_w0$alpha_matrix_standard_error[3,2],cxr_B.R4_w0$alpha_matrix[4,1]+cxr_B.R4_w0$alpha_matrix_standard_error[4,1], cxr_B.R4_w0$alpha_matrix[4,2]+cxr_B.R4_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Tu_inter"]<-c(cxr_B.R5_w0$alpha_matrix[1,3]+cxr_B.R5_w0$alpha_matrix_standard_error[1,3], cxr_B.R5_w0$alpha_matrix[2,3]+cxr_B.R5_w0$alpha_matrix_standard_error[2,3],cxr_B.R5_w0$alpha_matrix[1,4]+cxr_B.R5_w0$alpha_matrix_standard_error[1,4], cxr_B.R5_w0$alpha_matrix[2,4]+cxr_B.R5_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Te_inter"]<-c(cxr_B.R5_w0$alpha_matrix[3,1]+cxr_B.R5_w0$alpha_matrix_standard_error[3,1], cxr_B.R5_w0$alpha_matrix[3,2]+cxr_B.R5_w0$alpha_matrix_standard_error[3,2],cxr_B.R5_w0$alpha_matrix[4,1]+cxr_B.R5_w0$alpha_matrix_standard_error[4,1], cxr_B.R5_w0$alpha_matrix[4,2]+cxr_B.R5_w0$alpha_matrix_standard_error[4,2])



##### Cadmium


# modifying data frame to fit the type of setup that is need for CXR
CXR_B_Cd<-subset(ca, Env=="Cd")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "TuFemales")]

CXR_B_Cd$Focal<-mapvalues(CXR_B_Cd$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
CXR_B_Cd$CompSR2<-mapvalues(CXR_B_Cd$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

CXR_B_Cd$Comp<-sapply(c(1:length(CXR_B_Cd[,1])), function(x){
  if(is.na(CXR_B_Cd$CompSR2[x])){
    a<- CXR_B_Cd$Focal[x]
  }else{
    a<-CXR_B_Cd$CompSR2[x]
  }
  
  a
})

aux<-data.frame(SR1=rep(0, length(CXR_B_Cd[,1])), SR2=rep(0, length(CXR_B_Cd[,1])), SR4=rep(0, length(CXR_B_Cd[,1])), SR5=rep(0, length(CXR_B_Cd[,1])))

for(i in 1:length(CXR_B_Cd[,1])){
  #coluna onde por focais
  colunaF<-which(colnames(aux)==CXR_B_Cd$Focal[i])
  #coluna onde por competidors
  colunaC<-which(colnames(aux)==CXR_B_Cd$Comp[i])
  
  #if its the same regime
  if(CXR_B_Cd$Focal[i]==CXR_B_Cd$Comp[i] & CXR_B_Cd$Dens[i]==1){
    aux[i,colunaF]<-CXR_B_Cd$Dens[i]-1
    
  }else if(CXR_B_Cd$Focal[i]==CXR_B_Cd$Comp[i]){
    aux[i,colunaF]<-CXR_B_Cd$Dens[i]-1
  }else{ #if it is heterospecific then its -1 for the competitors (because of the focal) and its one for the focal
    aux[i,colunaC]<-CXR_B_Cd$Dens[i]-1
    aux[i, colunaF]<-1
  }
  
}

CXR_B_Cd<-cbind(CXR_B_Cd, aux)

CXR_B_Cd$fitness<-sapply(c(1:length(CXR_B_Cd[,1])), function(x){
  colF<-which(colnames(CXR_B_Cd)==CXR_B_Cd$Focal[x])
  
  if(CXR_B_Cd$Focal[x]=="SR1"){
    a<-CXR_B_Cd$TuFemales[x]/CXR_B_Cd$SR1[x]
  } else if(CXR_B_Cd$Focal[x]=="SR2"){
    a<-CXR_B_Cd$TuFemales[x]/CXR_B_Cd$SR2[x]
  } else if(CXR_B_Cd$Focal[x]=="SR4"){
    a<-CXR_B_Cd$TeFemales[x]/CXR_B_Cd$SR4[x]
  } else if(CXR_B_Cd$Focal[x]=="SR5"){
    a<-CXR_B_Cd$TeFemales[x]/CXR_B_Cd$SR5[x]
  }
  
  a
})

#removing rows for which there is no data for fitness
#CXR_B_Cd<-CXR_B_Cd[-which(is.na(CXR_B_Cd$fitness)),]
#CXR_B_Cd$fitness<-CXR_B_Cd$fitness+1

CXR_B_Cd[which(CXR_B_Cd$fitness=="-Inf" | CXR_B_Cd$fitness=="Inf"),"fitness"]<-0

#0 to 1 to mainrain data
CXR_B_Cd<-CXR_B_Cd[-which(is.na(CXR_B_Cd$fitness)),]
CXR_B_Cd$fitness<-CXR_B_Cd$fitness+1



# vector that tells which are the selection regimes, the columns have to have the same name
my.reg <- c("SR1", "SR2","SR4","SR5")

# Do list per replicate and environment
R1_Cd<-list(SR1= subset(CXR_B_Cd, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(CXR_B_Cd, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(CXR_B_Cd, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(CXR_B_Cd, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R2_Cd<-list(SR1= subset(CXR_B_Cd, Rep==2 & Focal=="SR1")[,c("fitness", "SR1", "SR2","SR4", "SR5")], SR4= subset(CXR_B_Cd, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(CXR_B_Cd, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R3_Cd<-list(SR1= subset(CXR_B_Cd, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(CXR_B_Cd, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(CXR_B_Cd, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(CXR_B_Cd, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R4_Cd<-list(SR1= subset(CXR_B_Cd, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(CXR_B_Cd, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(CXR_B_Cd, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(CXR_B_Cd, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R5_Cd<-list(SR1= subset(CXR_B_Cd, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(CXR_B_Cd, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(CXR_B_Cd, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(CXR_B_Cd, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

fixed_terms_1Cd <- list(list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                        list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==2)$lambda), # focal sp 2
                        list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==4)$lambda),
                        list(lambda= subset(mean_dens1, Rep==1 & Env=="Cd" & SR==5)$lambda))

fixed_terms_2Cd <- list(list(lambda = subset(mean_dens1, Rep==2 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                        list(lambda = subset(mean_dens1, Rep==2 & Env=="Cd" & SR==4)$lambda),
                        list(lambda= subset(mean_dens1, Rep==2 & Env=="Cd" & SR==5)$lambda))

fixed_terms_3Cd <- list(list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                        list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==2)$lambda), # focal sp 2
                        list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==4)$lambda),
                        list(lambda= subset(mean_dens1, Rep==3 & Env=="Cd" & SR==5)$lambda))

fixed_terms_4Cd <- list(list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                        list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==2)$lambda), # focal sp 2
                        list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==4)$lambda),
                        list(lambda= subset(mean_dens1, Rep==4 & Env=="Cd" & SR==5)$lambda))

fixed_terms_5Cd <- list(list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                        list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==2)$lambda), # focal sp 2
                        list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==4)$lambda),
                        list(lambda= subset(mean_dens1, Rep==5 & Env=="Cd" & SR==5)$lambda))

cxr_B.R1_Cd_w0<-cxr_pm_multifit(data = R1_Cd,
                                focal_column = my.reg,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "pairwise",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_intra = 0.1,
                                                      alpha_inter = 0.1),
                                fixed_terms = fixed_terms_1Cd,
                                # no standard errors
                                bootstrap_samples = 200)

# replicate 2 below


cxr_B.R3_Cd_w0<-cxr_pm_multifit(data = R3_Cd,
                                focal_column = my.reg,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "pairwise",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_intra = 0.1,
                                                      alpha_inter = 0.1),
                                fixed_terms = fixed_terms_3Cd,
                                # no standard errors
                                bootstrap_samples =10)

cxr_B.R4_Cd_w0<-cxr_pm_multifit(data = R4_Cd,
                                focal_column = my.reg,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "pairwise",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_intra = 0.1,
                                                      alpha_inter = 0.1),
                                fixed_terms = fixed_terms_4Cd,
                                # no standard errors
                                bootstrap_samples = 200)

cxr_B.R5_Cd_w0<-cxr_pm_multifit(data = R5_Cd,
                                focal_column = my.reg,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "pairwise",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_intra = 0.1,
                                                      alpha_inter = 0.1),
                                fixed_terms = fixed_terms_5Cd,
                                # no standard errors
                                bootstrap_samples = 200)

summary(cxr_B.R1_Cd_w0)
#summary(cxr_B.R2_Cd_w0)
summary(cxr_B.R3_Cd_w0)
summary(cxr_B.R4_Cd_w0)
summary(cxr_B.R5_Cd_w0)

cxr_B.R2_Cd_w0_sr1<-cxr_pm_fit(data = R2_Cd[[1]],
                               focal_column = my.reg[1],
                               model_family = "RK",
                               covariates = NULL,
                               optimization_method = "Nelder-Mead",
                               alpha_form = "pairwise",
                               lambda_cov_form = "none",
                               alpha_cov_form = "none",
                               initial_values = list(alpha_intra = 0.1,
                                                     alpha_inter = 0.1),
                               fixed_terms = fixed_terms_2Cd[[1]],
                               # no standard errors
                               bootstrap_samples = 200)

#for replicate 2 we will do the fitting by hand because we may need to scale the parameters

cxr_B.R2_Cd_w0_sr4<-cxr_pm_fit(data = R2_Cd[[2]][which(R2_Cd[[2]][,"SR1"]==0), c("fitness", "SR4")],
                               focal_column = NULL,
                               model_family = "RK",
                               covariates = NULL,
                               optimization_method = "Nelder-Mead",
                               alpha_form = "global",
                               lambda_cov_form = "none",
                               alpha_cov_form = "none",
                               initial_values = list(alpha_inter = 0.1),
                               fixed_terms = fixed_terms_2Cd[[2]],
                               # no standard errors
                               bootstrap_samples = 200)

cxr_B.R2_Cd_w0_5<-cxr_pm_fit(data = R2_Cd[[3]][which(R2_Cd[[3]][,"SR1"]==0), c("fitness", "SR5")],
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "global",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = 0.1),
                             fixed_terms = fixed_terms_2Cd[[3]],
                             # no standard errors
                             bootstrap_samples = 200)


cxr_B.R2_Cd_w0_sr4_inter<-cxr_pm_fit(data = R2_Cd[[2]][which(R2_Cd[[2]][,"SR1"]!=0), c("fitness", "SR1")],
                                     focal_column = NULL,
                                     model_family = "RK",
                                     covariates = NULL,
                                     optimization_method = "Nelder-Mead",
                                     alpha_form = "global",
                                     lambda_cov_form = "none",
                                     alpha_cov_form = "none",
                                     initial_values = list(alpha_inter = 0.1),
                                     fixed_terms = fixed_terms_2Cd[[2]],
                                     # no standard errors
                                     bootstrap_samples = 200)

cxr_B.R2_Cd_w0_sr5_inter<-cxr_pm_fit(data = R2_Cd[[3]][which(R2_Cd[[3]][,"SR1"]!=0), c("fitness", "SR1")],
                                     focal_column = NULL,
                                     model_family = "RK",
                                     covariates = NULL,
                                     optimization_method = "Nelder-Mead",
                                     alpha_form = "global",
                                     lambda_cov_form = "none",
                                     alpha_cov_form = "none",
                                     initial_values = list(alpha_inter = 0.1),
                                     fixed_terms = fixed_terms_2Cd[[3]],
                                     # no standard errors
                                     bootstrap_samples = 200)




###### data table summary
cxr_param_BC<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_BC$Tu_lambda<-0
cxr_param_BC$Te_lambda<-0
cxr_param_BC$Tu_intra<-0
cxr_param_BC$Te_intra<-0
cxr_param_BC$Tu_inter<-0
cxr_param_BC$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_BC<-cxr_param_BC[-which(cxr_param_BC$Replicate==2 & cxr_param_BC$Tu_Regime=="SR2"),]


cxr_param_BC[which(cxr_param_BC$Replicate==1),"Tu_lambda"]<-c(cxr_B.R1_Cd_w0$fixed_terms[[1]]$lambda,cxr_B.R1_Cd_w0$fixed_terms[[2]]$lambda)
cxr_param_BC[which(cxr_param_BC$Replicate==1),"Te_lambda"]<-c(cxr_B.R1_Cd_w0$fixed_terms[[3]]$lambda,cxr_B.R1_Cd_w0$fixed_terms[[3]]$lambda, cxr_B.R1_Cd_w0$fixed_terms[[4]]$lambda,cxr_B.R1_Cd_w0$fixed_terms[[4]]$lambda)

cxr_param_BC[which(cxr_param_BC$Replicate==2),"Tu_lambda"]<-cxr_B.R2_Cd_w0_sr1$fixed_terms$lambda
cxr_param_BC[which(cxr_param_BC$Replicate==2),"Te_lambda"]<-c(cxr_B.R2_Cd_w0_sr4$fixed_terms$lambda, cxr_B.R2_Cd_w0_5$fixed_terms$lambda)

cxr_param_BC[which(cxr_param_BC$Replicate==3),"Tu_lambda"]<-c(cxr_B.R3_Cd_w0$fixed_terms[[1]]$lambda,cxr_B.R3_Cd_w0$fixed_terms[[2]]$lambda)
cxr_param_BC[which(cxr_param_BC$Replicate==3),"Te_lambda"]<-c(cxr_B.R3_Cd_w0$fixed_terms[[3]]$lambda,cxr_B.R3_Cd_w0$fixed_terms[[3]]$lambda, cxr_B.R3_Cd_w0$fixed_terms[[4]]$lambda,cxr_B.R3_Cd_w0$fixed_terms[[4]]$lambda)

cxr_param_BC[which(cxr_param_BC$Replicate==4),"Tu_lambda"]<-c(cxr_B.R4_Cd_w0$fixed_terms[[1]]$lambda,cxr_B.R4_Cd_w0$fixed_terms[[2]]$lambda)
cxr_param_BC[which(cxr_param_BC$Replicate==4),"Te_lambda"]<-c(cxr_B.R4_Cd_w0$fixed_terms[[3]]$lambda,cxr_B.R4_Cd_w0$fixed_terms[[3]]$lambda, cxr_B.R4_Cd_w0$fixed_terms[[4]]$lambda,cxr_B.R4_Cd_w0$fixed_terms[[4]]$lambda)

cxr_param_BC[which(cxr_param_BC$Replicate==5),"Tu_lambda"]<-c(cxr_B.R5_Cd_w0$fixed_terms[[1]]$lambda,cxr_B.R5_Cd_w0$fixed_terms[[2]]$lambda)
cxr_param_BC[which(cxr_param_BC$Replicate==5),"Te_lambda"]<-c(cxr_B.R5_Cd_w0$fixed_terms[[3]]$lambda,cxr_B.R5_Cd_w0$fixed_terms[[3]]$lambda, cxr_B.R5_Cd_w0$fixed_terms[[4]]$lambda,cxr_B.R5_Cd_w0$fixed_terms[[4]]$lambda)


cxr_param_BC[which(cxr_param_BC$Replicate==1),"Tu_intra"]<-rep(c(cxr_B.R1_Cd_w0$alpha_matrix[1,1], cxr_B.R1_Cd_w0$alpha_matrix[2,2]), 2)
cxr_param_BC[which(cxr_param_BC$Replicate==1),"Te_intra"]<-rep(c(cxr_B.R1_Cd_w0$alpha_matrix[3,3], cxr_B.R1_Cd_w0$alpha_matrix[4,4]), each=2)

cxr_param_BC[which(cxr_param_BC$Replicate==2),"Tu_intra"]<-cxr_B.R2_Cd_w0_sr1$alpha_intra
cxr_param_BC[which(cxr_param_BC$Replicate==2),"Te_intra"]<-c(cxr_B.R2_Cd_w0_sr4$alpha_inter,cxr_B.R2_Cd_w0_sr5_inter$alpha_inter)

cxr_param_BC[which(cxr_param_BC$Replicate==3),"Tu_intra"]<-rep(c(cxr_B.R3_Cd_w0$alpha_matrix[1,1], cxr_B.R3_Cd_w0$alpha_matrix[2,2]), 2)
cxr_param_BC[which(cxr_param_BC$Replicate==3),"Te_intra"]<-rep(c(cxr_B.R3_Cd_w0$alpha_matrix[3,3], cxr_B.R3_Cd_w0$alpha_matrix[4,4]), each=2)

cxr_param_BC[which(cxr_param_BC$Replicate==4),"Tu_intra"]<-rep(c(cxr_B.R4_Cd_w0$alpha_matrix[1,1], cxr_B.R4_Cd_w0$alpha_matrix[2,2]), 2)
cxr_param_BC[which(cxr_param_BC$Replicate==4),"Te_intra"]<-rep(c(cxr_B.R4_Cd_w0$alpha_matrix[3,3], cxr_B.R4_Cd_w0$alpha_matrix[4,4]), each=2)

cxr_param_BC[which(cxr_param_BC$Replicate==5),"Tu_intra"]<-rep(c(cxr_B.R5_Cd_w0$alpha_matrix[1,1], cxr_B.R5_Cd_w0$alpha_matrix[2,2]), 2)
cxr_param_BC[which(cxr_param_BC$Replicate==5),"Te_intra"]<-rep(c(cxr_B.R5_Cd_w0$alpha_matrix[3,3], cxr_B.R5_Cd_w0$alpha_matrix[4,4]), each=2)


cxr_param_BC[which(cxr_param_BC$Replicate==1),"Tu_inter"]<-c(cxr_B.R1_Cd_w0$alpha_matrix[1,3], cxr_B.R1_Cd_w0$alpha_matrix[2,3],cxr_B.R1_Cd_w0$alpha_matrix[1,4], cxr_B.R1_Cd_w0$alpha_matrix[2,4])
cxr_param_BC[which(cxr_param_BC$Replicate==1),"Te_inter"]<-c(cxr_B.R1_Cd_w0$alpha_matrix[3,1], cxr_B.R1_Cd_w0$alpha_matrix[3,2],cxr_B.R1_Cd_w0$alpha_matrix[4,1], cxr_B.R1_Cd_w0$alpha_matrix[4,2])

cxr_param_BC[which(cxr_param_BC$Replicate==2),"Tu_inter"]<-cxr_B.R2_Cd_w0_sr1$alpha_inter[2:3]
cxr_param_BC[which(cxr_param_BC$Replicate==2),"Te_inter"]<-c(cxr_B.R2_Cd_w0_sr4_inter$alpha_inter, cxr_B.R2_Cd_w0_sr5_inter$alpha_inter)

cxr_param_BC[which(cxr_param_BC$Replicate==3),"Tu_inter"]<-c(cxr_B.R3_Cd_w0$alpha_matrix[1,3], cxr_B.R3_Cd_w0$alpha_matrix[2,3],cxr_B.R3_Cd_w0$alpha_matrix[1,4], cxr_B.R3_Cd_w0$alpha_matrix[2,4])
cxr_param_BC[which(cxr_param_BC$Replicate==3),"Te_inter"]<-c(cxr_B.R3_Cd_w0$alpha_matrix[3,1], cxr_B.R3_Cd_w0$alpha_matrix[3,2],cxr_B.R3_Cd_w0$alpha_matrix[4,1], cxr_B.R3_Cd_w0$alpha_matrix[4,2])

cxr_param_BC[which(cxr_param_BC$Replicate==4),"Tu_inter"]<-c(cxr_B.R4_Cd_w0$alpha_matrix[1,3], cxr_B.R4_Cd_w0$alpha_matrix[2,3],cxr_B.R4_Cd_w0$alpha_matrix[1,4], cxr_B.R4_Cd_w0$alpha_matrix[2,4])
cxr_param_BC[which(cxr_param_BC$Replicate==4),"Te_inter"]<-c(cxr_B.R4_Cd_w0$alpha_matrix[3,1], cxr_B.R4_Cd_w0$alpha_matrix[3,2],cxr_B.R4_Cd_w0$alpha_matrix[4,1], cxr_B.R4_Cd_w0$alpha_matrix[4,2])

cxr_param_BC[which(cxr_param_BC$Replicate==5),"Tu_inter"]<-c(cxr_B.R5_Cd_w0$alpha_matrix[1,3], cxr_B.R5_Cd_w0$alpha_matrix[2,3],cxr_B.R5_Cd_w0$alpha_matrix[1,4], cxr_B.R5_Cd_w0$alpha_matrix[2,4])
cxr_param_BC[which(cxr_param_BC$Replicate==5),"Te_inter"]<-c(cxr_B.R5_Cd_w0$alpha_matrix[3,1], cxr_B.R5_Cd_w0$alpha_matrix[3,2],cxr_B.R5_Cd_w0$alpha_matrix[4,1], cxr_B.R5_Cd_w0$alpha_matrix[4,2])

### Lower

cxr_param_BC_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_BC_lower$Tu_lambda<-0
cxr_param_BC_lower$Te_lambda<-0
cxr_param_BC_lower$Tu_intra<-0
cxr_param_BC_lower$Te_intra<-0
cxr_param_BC_lower$Tu_inter<-0
cxr_param_BC_lower$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_BC_lower<-cxr_param_BC_lower[-which(cxr_param_BC_lower$Replicate==2 & cxr_param_BC_lower$Tu_Regime=="SR2"),]

#Since the error comes directly from the data we need to create some lists with that information
sd_1C <- list(list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==1 & Env=="Cd" & SR==5)$sd_lambda))

sd_2C <- list(list(lambda = subset(mean_dens1, Rep==2 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==2 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==2 & Env=="Cd" & SR==5)$sd_lambda))

sd_3C <- list(list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==3 & Env=="Cd" & SR==5)$sd_lambda))

sd_4C <- list(list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==4 & Env=="Cd" & SR==5)$sd_lambda))

sd_5C <- list(list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==5 & Env=="Cd" & SR==5)$sd_lambda))

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==1),"Tu_lambda"]<-c(cxr_B.R1_Cd_w0$fixed_terms[[1]]$lambda-sd_1C[[1]]$lambda,cxr_B.R1_Cd_w0$fixed_terms[[2]]$lambda-sd_1C[[2]]$lambda)
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==1),"Te_lambda"]<-c(cxr_B.R1_Cd_w0$fixed_terms[[3]]$lambda-sd_1C[[3]]$lambda,cxr_B.R1_Cd_w0$fixed_terms[[3]]$lambda-sd_1C[[3]]$lambda, cxr_B.R1_Cd_w0$fixed_terms[[4]]$lambda-sd_1C[[4]]$lambda,cxr_B.R1_Cd_w0$fixed_terms[[4]]$lambda-sd_1C[[4]]$lambda)

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==2),"Tu_lambda"]<-cxr_B.R2_Cd_w0_sr1$fixed_terms$lambda-sd_2C[[1]]$lambda
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==2),"Te_lambda"]<-c(cxr_B.R2_Cd_w0_sr4$fixed_terms$lambda-sd_2C[[2]]$lambda, cxr_B.R2_Cd_w0_5$fixed_terms$lambda-sd_2C[[3]]$lambda)

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==3),"Tu_lambda"]<-c(cxr_B.R3_Cd_w0$fixed_terms[[1]]$lambda-sd_3C[[1]]$lambda,cxr_B.R3_Cd_w0$fixed_terms[[2]]$lambda-sd_3C[[2]]$lambda)
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==3),"Te_lambda"]<-c(cxr_B.R3_Cd_w0$fixed_terms[[3]]$lambda-sd_3C[[3]]$lambda,cxr_B.R3_Cd_w0$fixed_terms[[3]]$lambda-sd_3C[[3]]$lambda, cxr_B.R3_Cd_w0$fixed_terms[[4]]$lambda-sd_3C[[4]]$lambda,cxr_B.R3_Cd_w0$fixed_terms[[4]]$lambda-sd_3C[[4]]$lambda)

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==4),"Tu_lambda"]<-c(cxr_B.R4_Cd_w0$fixed_terms[[1]]$lambda-sd_4C[[1]]$lambda,cxr_B.R4_Cd_w0$fixed_terms[[2]]$lambda-sd_4C[[2]]$lambda)
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==4),"Te_lambda"]<-c(cxr_B.R4_Cd_w0$fixed_terms[[3]]$lambda-sd_4C[[3]]$lambda,cxr_B.R4_Cd_w0$fixed_terms[[3]]$lambda-sd_4C[[3]]$lambda, cxr_B.R4_Cd_w0$fixed_terms[[4]]$lambda-sd_4C[[4]]$lambda,cxr_B.R4_Cd_w0$fixed_terms[[4]]$lambda-sd_4C[[4]]$lambda)

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==5),"Tu_lambda"]<-c(cxr_B.R5_Cd_w0$fixed_terms[[1]]$lambda-sd_5C[[1]]$lambda,cxr_B.R5_Cd_w0$fixed_terms[[2]]$lambda-sd_5C[[2]]$lambda)
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==5),"Te_lambda"]<-c(cxr_B.R5_Cd_w0$fixed_terms[[3]]$lambda-sd_5C[[3]]$lambda,cxr_B.R5_Cd_w0$fixed_terms[[3]]$lambda-sd_5C[[3]]$lambda, cxr_B.R5_Cd_w0$fixed_terms[[4]]$lambda-sd_5C[[4]]$lambda,cxr_B.R5_Cd_w0$fixed_terms[[4]]$lambda-sd_5C[[4]]$lambda)


cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==1),"Tu_intra"]<-rep(c(cxr_B.R1_Cd_w0$alpha_matrix[1,1]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[1,1], cxr_B.R1_Cd_w0$alpha_matrix[2,2]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==1),"Te_intra"]<-rep(c(cxr_B.R1_Cd_w0$alpha_matrix[3,3]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[3,3], cxr_B.R1_Cd_w0$alpha_matrix[4,4]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==2),"Tu_intra"]<-cxr_B.R2_Cd_w0_sr1$alpha_intra[1]-cxr_B.R2_Cd_w0_sr1$alpha_intra_standard_error[1]
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==2),"Te_intra"]<-c(cxr_B.R2_Cd_w0_sr4$alpha_inter[1]-cxr_B.R2_Cd_w0_sr4$alpha_inter_standard_error[1], cxr_B.R2_Cd_w0_5$alpha_inter[1]-cxr_B.R2_Cd_w0_5$alpha_inter_standard_error[1])

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==3),"Tu_intra"]<-rep(c(cxr_B.R3_Cd_w0$alpha_matrix[1,1]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[1,1], cxr_B.R3_Cd_w0$alpha_matrix[2,2]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==3),"Te_intra"]<-rep(c(cxr_B.R3_Cd_w0$alpha_matrix[3,3]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[3,3], cxr_B.R3_Cd_w0$alpha_matrix[4,4]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==4),"Tu_intra"]<-rep(c(cxr_B.R4_Cd_w0$alpha_matrix[1,1]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[1,1], cxr_B.R4_Cd_w0$alpha_matrix[2,2]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==4),"Te_intra"]<-rep(c(cxr_B.R4_Cd_w0$alpha_matrix[3,3]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[3,3], cxr_B.R4_Cd_w0$alpha_matrix[4,4]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==5),"Tu_intra"]<-rep(c(cxr_B.R5_Cd_w0$alpha_matrix[1,1]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[1,1], cxr_B.R5_Cd_w0$alpha_matrix[2,2]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==5),"Te_intra"]<-rep(c(cxr_B.R5_Cd_w0$alpha_matrix[3,3]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[3,3], cxr_B.R5_Cd_w0$alpha_matrix[4,4]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)


cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==1),"Tu_inter"]<-c(cxr_B.R1_Cd_w0$alpha_matrix[1,3]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[1,3], cxr_B.R1_Cd_w0$alpha_matrix[2,3]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[2,3],cxr_B.R1_Cd_w0$alpha_matrix[1,4]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[1,4], cxr_B.R1_Cd_w0$alpha_matrix[2,4]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==1),"Te_inter"]<-c(cxr_B.R1_Cd_w0$alpha_matrix[3,1]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[3,1], cxr_B.R1_Cd_w0$alpha_matrix[3,2]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[3,2],cxr_B.R1_Cd_w0$alpha_matrix[4,1]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[4,1], cxr_B.R1_Cd_w0$alpha_matrix[4,2]-cxr_B.R1_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==2),"Tu_inter"]<-cxr_B.R2_Cd_w0_sr1$alpha_inter[2:3]-cxr_B.R2_Cd_w0_sr1$alpha_inter_standard_error[2:3]
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==2),"Te_inter"]<-c(cxr_B.R2_Cd_w0_sr4_inter$alpha_inter[1]-cxr_B.R2_Cd_w0_sr4_inter$alpha_inter_standard_error[1], cxr_B.R2_Cd_w0_sr5_inter$alpha_inter[1]-cxr_B.R2_Cd_w0_sr5_inter$alpha_inter_standard_error[1])

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==3),"Tu_inter"]<-c(cxr_B.R3_Cd_w0$alpha_matrix[1,3]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[1,3], cxr_B.R3_Cd_w0$alpha_matrix[2,3]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[2,3],cxr_B.R3_Cd_w0$alpha_matrix[1,4]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[1,4], cxr_B.R3_Cd_w0$alpha_matrix[2,4]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==3),"Te_inter"]<-c(cxr_B.R3_Cd_w0$alpha_matrix[3,1]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[3,1], cxr_B.R3_Cd_w0$alpha_matrix[3,2]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[3,2],cxr_B.R3_Cd_w0$alpha_matrix[4,1]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[4,1], cxr_B.R3_Cd_w0$alpha_matrix[4,2]-cxr_B.R3_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==4),"Tu_inter"]<-c(cxr_B.R4_Cd_w0$alpha_matrix[1,3]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[1,3], cxr_B.R4_Cd_w0$alpha_matrix[2,3]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[2,3],cxr_B.R4_Cd_w0$alpha_matrix[1,4]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[1,4], cxr_B.R4_Cd_w0$alpha_matrix[2,4]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==4),"Te_inter"]<-c(cxr_B.R4_Cd_w0$alpha_matrix[3,1]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[3,1], cxr_B.R4_Cd_w0$alpha_matrix[3,2]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[3,2],cxr_B.R4_Cd_w0$alpha_matrix[4,1]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[4,1], cxr_B.R4_Cd_w0$alpha_matrix[4,2]-cxr_B.R4_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==5),"Tu_inter"]<-c(cxr_B.R5_Cd_w0$alpha_matrix[1,3]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[1,3], cxr_B.R5_Cd_w0$alpha_matrix[2,3]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[2,3],cxr_B.R5_Cd_w0$alpha_matrix[1,4]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[1,4], cxr_B.R5_Cd_w0$alpha_matrix[2,4]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_BC_lower[which(cxr_param_BC_lower$Replicate==5),"Te_inter"]<-c(cxr_B.R5_Cd_w0$alpha_matrix[3,1]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[3,1], cxr_B.R5_Cd_w0$alpha_matrix[3,2]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[3,2],cxr_B.R5_Cd_w0$alpha_matrix[4,1]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[4,1], cxr_B.R5_Cd_w0$alpha_matrix[4,2]-cxr_B.R5_Cd_w0$alpha_matrix_standard_error[4,2])

### upper

cxr_param_BC_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_BC_upper$Tu_lambda<-0
cxr_param_BC_upper$Te_lambda<-0
cxr_param_BC_upper$Tu_intra<-0
cxr_param_BC_upper$Te_intra<-0
cxr_param_BC_upper$Tu_inter<-0
cxr_param_BC_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_BC_upper<-cxr_param_BC_upper[-which(cxr_param_BC_upper$Replicate==2 & cxr_param_BC_upper$Tu_Regime=="SR2"),]


cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==1),"Tu_lambda"]<-c(cxr_B.R1_Cd_w0$fixed_terms[[1]]$lambda+sd_1C[[1]]$lambda,cxr_B.R1_Cd_w0$fixed_terms[[2]]$lambda+sd_1C[[2]]$lambda)
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==1),"Te_lambda"]<-c(cxr_B.R1_Cd_w0$fixed_terms[[3]]$lambda+sd_1C[[3]]$lambda,cxr_B.R1_Cd_w0$fixed_terms[[3]]$lambda+sd_1C[[3]]$lambda, cxr_B.R1_Cd_w0$fixed_terms[[4]]$lambda+sd_1C[[4]]$lambda,cxr_B.R1_Cd_w0$fixed_terms[[4]]$lambda+sd_1C[[4]]$lambda)

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==2),"Tu_lambda"]<-cxr_B.R2_Cd_w0_sr1$fixed_terms$lambda+sd_2C[[1]]$lambda
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==2),"Te_lambda"]<-c(cxr_B.R2_Cd_w0_sr4$fixed_terms$lambda+sd_2C[[2]]$lambda, cxr_B.R2_Cd_w0_5$fixed_terms$lambda+sd_2C[[3]]$lambda)

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==3),"Tu_lambda"]<-c(cxr_B.R3_Cd_w0$fixed_terms[[1]]$lambda+sd_3C[[1]]$lambda,cxr_B.R3_Cd_w0$fixed_terms[[2]]$lambda+sd_3C[[2]]$lambda)
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==3),"Te_lambda"]<-c(cxr_B.R3_Cd_w0$fixed_terms[[3]]$lambda+sd_3C[[3]]$lambda,cxr_B.R3_Cd_w0$fixed_terms[[3]]$lambda+sd_3C[[3]]$lambda, cxr_B.R3_Cd_w0$fixed_terms[[4]]$lambda+sd_3C[[4]]$lambda,cxr_B.R3_Cd_w0$fixed_terms[[4]]$lambda+sd_3C[[4]]$lambda)

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==4),"Tu_lambda"]<-c(cxr_B.R4_Cd_w0$fixed_terms[[1]]$lambda+sd_4C[[1]]$lambda,cxr_B.R4_Cd_w0$fixed_terms[[2]]$lambda+sd_4C[[2]]$lambda)
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==4),"Te_lambda"]<-c(cxr_B.R4_Cd_w0$fixed_terms[[3]]$lambda+sd_4C[[3]]$lambda,cxr_B.R4_Cd_w0$fixed_terms[[3]]$lambda+sd_4C[[3]]$lambda, cxr_B.R4_Cd_w0$fixed_terms[[4]]$lambda+sd_4C[[4]]$lambda,cxr_B.R4_Cd_w0$fixed_terms[[4]]$lambda+sd_4C[[4]]$lambda)

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==5),"Tu_lambda"]<-c(cxr_B.R5_Cd_w0$fixed_terms[[1]]$lambda+sd_5C[[1]]$lambda,cxr_B.R5_Cd_w0$fixed_terms[[2]]$lambda+sd_5C[[2]]$lambda)
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==5),"Te_lambda"]<-c(cxr_B.R5_Cd_w0$fixed_terms[[3]]$lambda+sd_5C[[3]]$lambda,cxr_B.R5_Cd_w0$fixed_terms[[3]]$lambda+sd_5C[[3]]$lambda, cxr_B.R5_Cd_w0$fixed_terms[[4]]$lambda+sd_5C[[4]]$lambda,cxr_B.R5_Cd_w0$fixed_terms[[4]]$lambda+sd_5C[[4]]$lambda)


cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==1),"Tu_intra"]<-rep(c(cxr_B.R1_Cd_w0$alpha_matrix[1,1]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[1,1], cxr_B.R1_Cd_w0$alpha_matrix[2,2]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==1),"Te_intra"]<-rep(c(cxr_B.R1_Cd_w0$alpha_matrix[3,3]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[3,3], cxr_B.R1_Cd_w0$alpha_matrix[4,4]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==2),"Tu_intra"]<-cxr_B.R2_Cd_w0_sr1$alpha_intra[1]+cxr_B.R2_Cd_w0_sr1$alpha_intra_standard_error[1]
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==2),"Te_intra"]<-c(cxr_B.R2_Cd_w0_sr4$alpha_inter[1]+cxr_B.R2_Cd_w0_sr4$alpha_inter_standard_error[1], cxr_B.R2_Cd_w0_5$alpha_inter[1]+cxr_B.R2_Cd_w0_5$alpha_inter_standard_error[1])

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==3),"Tu_intra"]<-rep(c(cxr_B.R3_Cd_w0$alpha_matrix[1,1]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[1,1], cxr_B.R3_Cd_w0$alpha_matrix[2,2]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==3),"Te_intra"]<-rep(c(cxr_B.R3_Cd_w0$alpha_matrix[3,3]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[3,3], cxr_B.R3_Cd_w0$alpha_matrix[4,4]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==4),"Tu_intra"]<-rep(c(cxr_B.R4_Cd_w0$alpha_matrix[1,1]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[1,1], cxr_B.R4_Cd_w0$alpha_matrix[2,2]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==4),"Te_intra"]<-rep(c(cxr_B.R4_Cd_w0$alpha_matrix[3,3]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[3,3], cxr_B.R4_Cd_w0$alpha_matrix[4,4]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==5),"Tu_intra"]<-rep(c(cxr_B.R5_Cd_w0$alpha_matrix[1,1]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[1,1], cxr_B.R5_Cd_w0$alpha_matrix[2,2]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==5),"Te_intra"]<-rep(c(cxr_B.R5_Cd_w0$alpha_matrix[3,3]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[3,3], cxr_B.R5_Cd_w0$alpha_matrix[4,4]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[4,4]), each=2)


cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==1),"Tu_inter"]<-c(cxr_B.R1_Cd_w0$alpha_matrix[1,3]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[1,3], cxr_B.R1_Cd_w0$alpha_matrix[2,3]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[2,3],cxr_B.R1_Cd_w0$alpha_matrix[1,4]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[1,4], cxr_B.R1_Cd_w0$alpha_matrix[2,4]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==1),"Te_inter"]<-c(cxr_B.R1_Cd_w0$alpha_matrix[3,1]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[3,1], cxr_B.R1_Cd_w0$alpha_matrix[3,2]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[3,2],cxr_B.R1_Cd_w0$alpha_matrix[4,1]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[4,1], cxr_B.R1_Cd_w0$alpha_matrix[4,2]+cxr_B.R1_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==2),"Tu_inter"]<-cxr_B.R2_Cd_w0_sr1$alpha_inter[2:3]+cxr_B.R2_Cd_w0_sr1$alpha_inter_standard_error[2:3]
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==2),"Te_inter"]<-c(cxr_B.R2_Cd_w0_sr4_inter$alpha_inter[1]+cxr_B.R2_Cd_w0_sr4_inter$alpha_inter_standard_error[1], cxr_B.R2_Cd_w0_sr5_inter$alpha_inter[1]+cxr_B.R2_Cd_w0_sr5_inter$alpha_inter_standard_error[1])

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==3),"Tu_inter"]<-c(cxr_B.R3_Cd_w0$alpha_matrix[1,3]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[1,3], cxr_B.R3_Cd_w0$alpha_matrix[2,3]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[2,3],cxr_B.R3_Cd_w0$alpha_matrix[1,4]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[1,4], cxr_B.R3_Cd_w0$alpha_matrix[2,4]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==3),"Te_inter"]<-c(cxr_B.R3_Cd_w0$alpha_matrix[3,1]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[3,1], cxr_B.R3_Cd_w0$alpha_matrix[3,2]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[3,2],cxr_B.R3_Cd_w0$alpha_matrix[4,1]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[4,1], cxr_B.R3_Cd_w0$alpha_matrix[4,2]+cxr_B.R3_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==4),"Tu_inter"]<-c(cxr_B.R4_Cd_w0$alpha_matrix[1,3]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[1,3], cxr_B.R4_Cd_w0$alpha_matrix[2,3]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[2,3],cxr_B.R4_Cd_w0$alpha_matrix[1,4]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[1,4], cxr_B.R4_Cd_w0$alpha_matrix[2,4]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==4),"Te_inter"]<-c(cxr_B.R4_Cd_w0$alpha_matrix[3,1]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[3,1], cxr_B.R4_Cd_w0$alpha_matrix[3,2]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[3,2],cxr_B.R4_Cd_w0$alpha_matrix[4,1]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[4,1], cxr_B.R4_Cd_w0$alpha_matrix[4,2]+cxr_B.R4_Cd_w0$alpha_matrix_standard_error[4,2])

cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==5),"Tu_inter"]<-c(cxr_B.R5_Cd_w0$alpha_matrix[1,3]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[1,3], cxr_B.R5_Cd_w0$alpha_matrix[2,3]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[2,3],cxr_B.R5_Cd_w0$alpha_matrix[1,4]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[1,4], cxr_B.R5_Cd_w0$alpha_matrix[2,4]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[2,4])
cxr_param_BC_upper[which(cxr_param_BC_upper$Replicate==5),"Te_inter"]<-c(cxr_B.R5_Cd_w0$alpha_matrix[3,1]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[3,1], cxr_B.R5_Cd_w0$alpha_matrix[3,2]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[3,2],cxr_B.R5_Cd_w0$alpha_matrix[4,1]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[4,1], cxr_B.R5_Cd_w0$alpha_matrix[4,2]+cxr_B.R5_Cd_w0$alpha_matrix_standard_error[4,2])




##### joining data frame
param_all_B<-as.data.frame(rbind(cxr_param_B, cxr_param_BC))

param_all_B_lower<-as.data.frame(rbind(cxr_param_B_lower, cxr_param_BC_lower))
param_all_B_upper<-as.data.frame(rbind(cxr_param_B_upper, cxr_param_BC_upper))

param_all_B_lower
param_all_B_upper

#write.csv(param_all_B, "./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed.csv")
#write.csv(param_all_B_upper, "./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed_upper.csv")
#write.csv(param_all_B_lower, "./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed_lower.csv")



##### importing data frame
param_all_B<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed.csv")
param_all_B_upper<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed_upper.csv")
param_all_B_lower<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed_lower.csv")

param_all_B<-param_all_B[,-1]
param_all_B_upper<-param_all_B_upper[,-1]
param_all_B_lower<-param_all_B_lower[,-1]



##### Plotting data
param_all_B_long<-gather(param_all_B, parameter, value,Tu_lambda:Te_inter )

param_all_B_long$category<-mapvalues(param_all_B_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

param_all_B_lower_long<-gather(param_all_B_lower, parameter, value,Tu_lambda:Te_inter )

param_all_B_lower_long$category<-mapvalues(param_all_B_lower_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

param_all_B_upper_long<-gather(param_all_B_upper, parameter, value,Tu_lambda:Te_inter )

param_all_B_upper_long$category<-mapvalues(param_all_B_upper_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

colnames(param_all_B_lower_long)[6]<-"lower"
colnames(param_all_B_upper_long)[6]<-"upper"

str(param_all_B_long)

param_all_B_long<-cbind(param_all_B_long[,1:7],param_all_B_lower_long$lower, param_all_B_upper_long$upper)

colnames(param_all_B_long)[8:9]<-c("lower","upper")


#### Predicting densities
density_aux<-seq(0, 10, by=(10/100))

pred_df_cxr_B<-as.data.frame(expand_grid(Density=density_aux, Tu_Regime=c("SR1","SR2"), Te_Regime=c("SR4","SR5"), Replicate=c(1:5), Environment=c("N", "Cd")))

pred_df_cxr_B$Tu_mean_intra<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Tu_mean_inter<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_B$Density[x])
  
  pred
})


pred_df_cxr_B$Tu_intra_L<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Tu_inter_L<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Tu_intra_U<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Tu_inter_U<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Te_mean_intra<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Te_mean_inter<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_B, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Te_intra_L<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Te_inter_L<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_B_lower, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Te_intra_U<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_B$Density[x])
  
  pred
})

pred_df_cxr_B$Te_inter_U<-sapply(c(1:length(pred_df_cxr_B[,1])), function(x){
  alpha_i<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Tu_Regime==pred_df_cxr_B$Tu_Regime[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_B_upper, Environment==pred_df_cxr_B$Environment[x] & Te_Regime==pred_df_cxr_B$Te_Regime[x] & Replicate==pred_df_cxr_B$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_B$Density[x])
  
  pred
})

# Removing Tu evolved replicate 2 because there is no data
pred_df_cxr_B<-pred_df_cxr_B[-which(pred_df_cxr_B$Tu_Regime=="SR2" & pred_df_cxr_B$Replicate==2),]



# Transforming everything bellow 0 into 0 for the lower interval

pred_df_cxr_B$Te_inter_L[which(pred_df_cxr_B$Te_inter_L<0)]<-0
pred_df_cxr_B$Te_intra_L[which(pred_df_cxr_B$Te_intra_L<0)]<-0
pred_df_cxr_B$Tu_inter_L[which(pred_df_cxr_B$Tu_inter_L<0)]<-0
pred_df_cxr_B$Tu_intra_L[which(pred_df_cxr_B$Tu_intra_L<0)]<-0

##### Predicted vs observed
str(param_all_B)

str(ca)

rk_func<- function(lambda, alpha_ii, alpha_ij, dens_i, dens_j, ...){
  gr<-lambda*exp(-alpha_ii*dens_i - alpha_ij*dens_j)
  
  return(gr)
}

red_ca_B<-ca[,c("Env", "Rep", "FocalSR", "CompSR", "Dens", "Type", "TeFemales", "TuFemales", "GrowthRateOA")]

red_ca_B

red_ca_B$Dens_Focal<-sapply(c(1:length(red_ca_B[,1])), function(x){
  if(red_ca_B$Type[x]=="INTRA"){
    a<-red_ca_B$Dens[x]-1
  }else if(red_ca_B$Type[x]=="INTER"){
    a<-1
  }
  
  a
})

red_ca_B$Dens_Comp<-sapply(c(1:length(red_ca_B[,1])), function(x){
  if(red_ca_B$Type[x]=="INTRA"){
    a<-0
  }else if(red_ca_B$Type[x]=="INTER"){
    a<-red_ca_B$Dens[x]-1
  }
  
  a
})

red_ca_B$Focal<-mapvalues(red_ca_B$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4", "SR5"))
red_ca_B$Comp<-mapvalues(red_ca_B$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4", "SR5"))

red_ca_B$pred<-sapply(c(1:length(red_ca_B[,1])), function(x){
  
  if(red_ca_B$Focal[x]=="SR1" | red_ca_B$Focal[x]=="SR2"){
    aux_data<-subset(param_all_B, Environment==red_ca_B$Env[x] & Replicate== red_ca_B$Rep[x] & as.character(Tu_Regime)==red_ca_B$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Tu_lambda[1], alpha_ii =aux_data$Tu_intra[1], alpha_ij = aux_data$Tu_inter[1], dens_i = red_ca_B$Dens_Focal[x], dens_j =  red_ca_B$Dens_Comp[x])
    
  }else if(red_ca_B$Focal[x]=="SR4" | red_ca_B$Focal[x]=="SR5"){
    aux_data<-subset(param_all_B, Environment==red_ca_B$Env[x] & Replicate== red_ca_B$Rep[x] & as.character(Te_Regime)==red_ca_B$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Te_lambda[1], alpha_ii =aux_data$Te_intra[1], alpha_ij = aux_data$Te_inter[1], dens_i = red_ca_B$Dens_Focal[x], dens_j =  red_ca_B$Dens_Comp[x])
  }
  
  aux_pred
})

red_ca_B$pred_L<-sapply(c(1:length(red_ca_B[,1])), function(x){
  
  if(red_ca_B$Focal[x]=="SR1" | red_ca_B$Focal[x]=="SR2"){
    aux_data<-subset(param_all_B_lower, Environment==red_ca_B$Env[x] & Replicate== red_ca_B$Rep[x] & as.character(Tu_Regime)==red_ca_B$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Tu_lambda[1], alpha_ii =aux_data$Tu_intra[1], alpha_ij = aux_data$Tu_inter[1], dens_i = red_ca_B$Dens_Focal[x], dens_j =  red_ca_B$Dens_Comp[x])
    
  }else if(red_ca_B$Focal[x]=="SR4" | red_ca_B$Focal[x]=="SR5"){
    aux_data<-subset(param_all_B_lower, Environment==red_ca_B$Env[x] & Replicate== red_ca_B$Rep[x] & as.character(Te_Regime)==red_ca_B$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Te_lambda[1], alpha_ii =aux_data$Te_intra[1], alpha_ij = aux_data$Te_inter[1], dens_i = red_ca_B$Dens_Focal[x], dens_j =  red_ca_B$Dens_Comp[x])
  }
  
  aux_pred
})

red_ca_B$pred_U<-sapply(c(1:length(red_ca_B[,1])), function(x){
  
  if(red_ca_B$Focal[x]=="SR1" | red_ca_B$Focal[x]=="SR2"){
    aux_data<-subset(param_all_B_upper, Environment==red_ca_B$Env[x] & Replicate== red_ca_B$Rep[x] & as.character(Tu_Regime)==red_ca_B$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Tu_lambda[1], alpha_ii =aux_data$Tu_intra[1], alpha_ij = aux_data$Tu_inter[1], dens_i = red_ca_B$Dens_Focal[x], dens_j =  red_ca_B$Dens_Comp[x])
    
  }else if(red_ca_B$Focal[x]=="SR4" | red_ca_B$Focal[x]=="SR5"){
    aux_data<-subset(param_all_B_upper, Environment==red_ca_B$Env[x] & Replicate== red_ca_B$Rep[x] & as.character(Te_Regime)==red_ca_B$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Te_lambda[1], alpha_ii =aux_data$Te_intra[1], alpha_ij = aux_data$Te_inter[1], dens_i = red_ca_B$Dens_Focal[x], dens_j =  red_ca_B$Dens_Comp[x])
  }
  
  aux_pred
})

red_ca_B$Replicate<-red_ca_B$Rep

### C - CXR nested
#To do this we have to trick the cxr (Oscar suggestion), by putting the intraspecific competitors in another column than the focal and then estimate only intra!
  
##### normal
dir.create("./Analyses/MethodComparison/cxr_lambda_fixed_nested", showWarnings = FALSE)

# modifying data frame to fit the type of setup that is need for CXR
CXR_C_N<-subset(ca, Env=="N")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "TuFemales")]

CXR_C_N$Focal<-mapvalues(CXR_C_N$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
CXR_C_N$CompSR2<-mapvalues(CXR_C_N$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

CXR_C_N$Comp<-sapply(c(1:length(CXR_C_N[,1])), function(x){
  if(is.na(CXR_C_N$CompSR2[x])){
    a<- CXR_C_N$Focal[x]
  }else{
    a<-CXR_C_N$CompSR2[x]
  }
  
  a
})

aux<-data.frame(SR1=rep(0, length(CXR_C_N[,1])), SR2=rep(0, length(CXR_C_N[,1])), SR4=rep(0, length(CXR_C_N[,1])), SR5=rep(0, length(CXR_C_N[,1])))

for(i in 1:length(CXR_C_N[,1])){
  #coluna onde por focais
  colunaF<-which(colnames(aux)==CXR_C_N$Focal[i])
  #coluna onde por competidors
  colunaC<-which(colnames(aux)==CXR_C_N$Comp[i])
  
  #if its the same regime
  if(CXR_C_N$Focal[i]==CXR_C_N$Comp[i] & CXR_C_N$Dens[i]==1){
    aux[i,colunaF]<-CXR_C_N$Dens[i]-1
    
  }else if(CXR_C_N$Focal[i]==CXR_C_N$Comp[i]){
    aux[i,colunaF]<-CXR_C_N$Dens[i]-1
  }else{ #if it is heterospecific then its -1 for the competitors (because of the focal) and its one for the focal
    aux[i,colunaC]<-CXR_C_N$Dens[i]-1
    aux[i, colunaF]<-1
  }
  
}

CXR_C_N<-cbind(CXR_C_N, aux)

CXR_C_N$fitness<-sapply(c(1:length(CXR_C_N[,1])), function(x){
  colF<-which(colnames(CXR_C_N)==CXR_C_N$Focal[x])
  
  if(CXR_C_N$Focal[x]=="SR1"){
    a<-CXR_C_N$TuFemales[x]/CXR_C_N$SR1[x]
  } else if(CXR_C_N$Focal[x]=="SR2"){
    a<-CXR_C_N$TuFemales[x]/CXR_C_N$SR2[x]
  } else if(CXR_C_N$Focal[x]=="SR4"){
    a<-CXR_C_N$TeFemales[x]/CXR_C_N$SR4[x]
  } else if(CXR_C_N$Focal[x]=="SR5"){
    a<-CXR_C_N$TeFemales[x]/CXR_C_N$SR5[x]
  }
  
  a
})

#removing rows for which there is no data for fitness
CXR_C_N<-CXR_C_N[-which(is.na(CXR_C_N$fitness)),]

# adding +1 to all data
#CXR_C_N$fitness<-CXR_C_N$fitness+1

CXR_C_N[which(CXR_C_N$fitness=="-Inf" | CXR_C_N$fitness=="Inf"),"fitness"]<-0


# all data gets +1 because of the 0 problem
CXR_C_N$fitness<-CXR_C_N$fitness+1

# vector that tells which are the selection regimes, the columns have to have the same name
my.reg <- c("SR1", "SR2","SR4","SR5")
str(CXR_C_N)

# Do list per replicate and environment
R1_intra<-list(SR1= subset(CXR_C_N, Rep==1 & Focal=="SR1" & Comp=="SR1")[,c("fitness", "SR1")], SR2= subset(CXR_C_N, Rep==1 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_N, Rep==1 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==1 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R2_intra<-list(SR1= subset(CXR_C_N, Rep==2 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR4= subset(CXR_C_N, Rep==2 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==2 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R3_intra<-list(SR1= subset(CXR_C_N, Rep==3 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_N, Rep==3 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_N, Rep==3 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==3 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R4_intra<-list(SR1= subset(CXR_C_N, Rep==4 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_N, Rep==4 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_N, Rep==4 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==4 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R5_intra<-list(SR1= subset(CXR_C_N, Rep==5 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_N, Rep==5 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_N, Rep==5 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==5 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])



####################################
######## DOING THE INTRA ESTIMATES
###################################


cxr_C.R1_intra<-cxr_pm_multifit(data = R1_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = 0.1),
                                fixed_terms = fixed_terms_1N,
                                # no standard errors
                                bootstrap_samples = 200)

cxr_C.R2_intra<-cxr_pm_multifit(data = R2_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = 0.1),
                                fixed_terms = fixed_terms_2N,
                                # no standard errors
                                bootstrap_samples = 200)

cxr_C.R3_intra<-cxr_pm_multifit(data = R3_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = 0.1),
                                fixed_terms = fixed_terms_3N,
                                # no standard errors
                                bootstrap_samples = 200)

cxr_C.R4_intra<-cxr_pm_multifit(data = R4_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = 0.1),
                                fixed_terms = fixed_terms_4N,
                                # no standard errors
                                bootstrap_samples = 200)

cxr_C.R5_intra<-cxr_pm_multifit(data = R5_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = 0.1),
                                fixed_terms = fixed_terms_5N,
                                # no standard errors
                                bootstrap_samples = 200)

summary(cxr_C.R1_intra)

#################
######Doing the inter estimates
##################

R1<-list(SR1= subset(CXR_C_N, Rep==1 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR2= subset(CXR_C_N, Rep==1 & Focal=="SR2"& Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==1 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_N, Rep==1 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R2<-list(SR1= subset(CXR_C_N, Rep==2 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==2 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1")], SR5= subset(CXR_C_N, Rep==2 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1")])

R3<-list(SR1= subset(CXR_C_N, Rep==3 & Focal=="SR1" & Comp!="SR1")[,c("fitness",  "SR4", "SR5")], SR2= subset(CXR_C_N, Rep==3 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==3 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_N, Rep==3 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R4<-list(SR1= subset(CXR_C_N, Rep==4 & Focal=="SR1" & Comp!="SR1")[,c("fitness",  "SR4", "SR5")], SR2= subset(CXR_C_N, Rep==4 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==4 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_N, Rep==4 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R5<-list(SR1= subset(CXR_C_N, Rep==5 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR2= subset(CXR_C_N, Rep==5 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==5 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_N, Rep==5 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])




cxr_C.R1<-cxr_pm_multifit(data = R1,
                          focal_column = NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(alpha_inter = 0.1),
                          fixed_terms = fixed_terms_1N,
                          # no standard errors
                          bootstrap_samples = 200)

cxr_C.R2_sr1<-cxr_pm_fit(data = R2[[1]],
                         focal_column = NULL,
                         model_family = "RK",
                         covariates = NULL,
                         optimization_method = "Nelder-Mead",
                         alpha_form = "pairwise",
                         lambda_cov_form = "none",
                         alpha_cov_form = "none",
                         initial_values = list(alpha_inter = 0.1),
                         fixed_terms = fixed_terms_2N[[1]],
                         # no standard errors
                         bootstrap_samples = 200)

cxr_C.R2_sr4<-cxr_pm_fit(data = R2[[2]],
                         focal_column = NULL,
                         model_family = "RK",
                         covariates = NULL,
                         optimization_method = "Nelder-Mead",
                         alpha_form = "global",
                         lambda_cov_form = "none",
                         alpha_cov_form = "none",
                         initial_values = list(alpha_inter = 0.1),
                         fixed_terms = fixed_terms_2N[[2]],
                         # no standard errors
                         bootstrap_samples = 200)

cxr_C.R2_sr5<-cxr_pm_fit(data = R2[[3]],
                         focal_column = NULL,
                         model_family = "RK",
                         covariates = NULL,
                         optimization_method = "Nelder-Mead",
                         alpha_form = "global",
                         lambda_cov_form = "none",
                         alpha_cov_form = "none",
                         initial_values = list(alpha_inter = 0.1),
                         fixed_terms = fixed_terms_2N[[3]],
                         # no standard errors
                         bootstrap_samples = 200)


cxr_C.R3<-cxr_pm_multifit(data = R3,
                          focal_column = NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(alpha_inter = 0.1),
                          fixed_terms = fixed_terms_3N,
                          # no standard errors
                          bootstrap_samples = 200)

cxr_C.R4<-cxr_pm_multifit(data = R4,
                          focal_column = NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(alpha_inter = 0.1),
                          fixed_terms = fixed_terms_4N,
                          # no standard errors
                          bootstrap_samples = 200)

cxr_C.R5<-cxr_pm_multifit(data = R5,
                          focal_column = NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(alpha_inter = 0.1),
                          fixed_terms = fixed_terms_5N,
                          # no standard errors
                          bootstrap_samples = 200)

cxr_C.R2_intra$alpha_matrix


###### cadmium


# modifying data frame to fit the type of setup that is need for CXR
CXR_C_Cd<-subset(ca, Env=="Cd")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "TuFemales")]

CXR_C_Cd$Focal<-mapvalues(CXR_C_Cd$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
CXR_C_Cd$CompSR2<-mapvalues(CXR_C_Cd$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

CXR_C_Cd$Comp<-sapply(c(1:length(CXR_C_Cd[,1])), function(x){
  if(is.na(CXR_C_Cd$CompSR2[x])){
    a<- CXR_C_Cd$Focal[x]
  }else{
    a<-CXR_C_Cd$CompSR2[x]
  }
  
  a
})

aux<-data.frame(SR1=rep(0, length(CXR_C_Cd[,1])), SR2=rep(0, length(CXR_C_Cd[,1])), SR4=rep(0, length(CXR_C_Cd[,1])), SR5=rep(0, length(CXR_C_Cd[,1])))

for(i in 1:length(CXR_C_Cd[,1])){
  #coluna onde por focais
  colunaF<-which(colnames(aux)==CXR_C_Cd$Focal[i])
  #coluna onde por competidors
  colunaC<-which(colnames(aux)==CXR_C_Cd$Comp[i])
  
  #if its the same regime
  if(CXR_C_Cd$Focal[i]==CXR_C_Cd$Comp[i] & CXR_C_Cd$Dens[i]==1){
    aux[i,colunaF]<-CXR_C_Cd$Dens[i]-1
    
  }else if(CXR_C_Cd$Focal[i]==CXR_C_Cd$Comp[i]){
    aux[i,colunaF]<-CXR_C_Cd$Dens[i]-1
  }else{ #if it is heterospecific then its -1 for the competitors (because of the focal) and its one for the focal
    aux[i,colunaC]<-CXR_C_Cd$Dens[i]-1
    aux[i, colunaF]<-1
  }
  
}

CXR_C_Cd<-cbind(CXR_C_Cd, aux)

CXR_C_Cd$fitness<-sapply(c(1:length(CXR_C_Cd[,1])), function(x){
  colF<-which(colnames(CXR_C_Cd)==CXR_C_Cd$Focal[x])
  
  if(CXR_C_Cd$Focal[x]=="SR1"){
    a<-CXR_C_Cd$TuFemales[x]/CXR_C_Cd$SR1[x]
  } else if(CXR_C_Cd$Focal[x]=="SR2"){
    a<-CXR_C_Cd$TuFemales[x]/CXR_C_Cd$SR2[x]
  } else if(CXR_C_Cd$Focal[x]=="SR4"){
    a<-CXR_C_Cd$TeFemales[x]/CXR_C_Cd$SR4[x]
  } else if(CXR_C_Cd$Focal[x]=="SR5"){
    a<-CXR_C_Cd$TeFemales[x]/CXR_C_Cd$SR5[x]
  }
  
  a
})

#removing rows for which there is no data for fitness
CXR_C_Cd<-CXR_C_Cd[-which(is.na(CXR_C_Cd$fitness)),]

# adding +1 to all data
#CXR_C_Cd$fitness<-CXR_C_Cd$fitness+1

CXR_C_Cd[which(CXR_C_Cd$fitness=="-Inf" | CXR_C_Cd$fitness=="Inf"),"fitness"]<-0


# all data gets +1 because of the 0 problem
CXR_C_Cd$fitness<-CXR_C_Cd$fitness+1

# vector that tells which are the selection regimes, the columns have to have the same name
my.reg <- c("SR1", "SR2","SR4","SR5")
str(CXR_C_Cd)

# Do list per replicate and environment
R1_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==1 & Focal=="SR1" & Comp=="SR1")[,c("fitness", "SR1")], SR2= subset(CXR_C_Cd, Rep==1 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_Cd, Rep==1 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==1 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R2_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==2 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR4= subset(CXR_C_Cd, Rep==2 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==2 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R3_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==3 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_Cd, Rep==3 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_Cd, Rep==3 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==3 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R4_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==4 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_Cd, Rep==4 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_Cd, Rep==4 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==4 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R5_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==5 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_Cd, Rep==5 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_Cd, Rep==5 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==5 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])


#### lambda

fixed_terms_C_1N <- list(list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                         list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==2)$lambda), # focal sp 2
                         list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==4)$lambda),
                         list(lambda= subset(mean_dens1, Rep==1 & Env=="Cd" & SR==5)$lambda))

fixed_terms_C_2N <- list(list(lambda = subset(mean_dens1, Rep==2 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                         list(lambda = subset(mean_dens1, Rep==2 & Env=="Cd" & SR==4)$lambda),
                         list(lambda= subset(mean_dens1, Rep==2 & Env=="Cd" & SR==5)$lambda))

fixed_terms_C_3N <- list(list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                         list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==2)$lambda), # focal sp 2
                         list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==4)$lambda),
                         list(lambda= subset(mean_dens1, Rep==3 & Env=="Cd" & SR==5)$lambda))

fixed_terms_C_4N <- list(list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                         list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==2)$lambda), # focal sp 2
                         list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==4)$lambda),
                         list(lambda= subset(mean_dens1, Rep==4 & Env=="Cd" & SR==5)$lambda))

fixed_terms_C_5N <- list(list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==1)$lambda ), # focal sp 1
                         list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==2)$lambda), # focal sp 2
                         list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==4)$lambda),
                         list(lambda= subset(mean_dens1, Rep==5 & Env=="Cd" & SR==5)$lambda))

####################################
######## DOING THE INTRA ESTIMATES
###################################


cxr_C.R1_cd_intra<-cxr_pm_multifit(data = R1_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = 0.1),
                                   fixed_terms = fixed_terms_C_1N,
                                   # no standard errors
                                   bootstrap_samples = 200)

cxr_C.R2_cd_intra<-cxr_pm_multifit(data = R2_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = 0.1),
                                   fixed_terms = fixed_terms_C_2N,
                                   # no standard errors
                                   bootstrap_samples = 200)

cxr_C.R3_cd_intra<-cxr_pm_multifit(data = R3_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = 0.1),
                                   fixed_terms = fixed_terms_C_3N,
                                   # no standard errors
                                   bootstrap_samples = 200)

cxr_C.R4_cd_intra<-cxr_pm_multifit(data = R4_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = 0.1),
                                   fixed_terms = fixed_terms_C_4N,
                                   # no standard errors
                                   bootstrap_samples = 200)

cxr_C.R5_cd_intra<-cxr_pm_multifit(data = R5_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = 0.1),
                                   fixed_terms = fixed_terms_C_5N,
                                   # no standard errors
                                   bootstrap_samples = 200)

summary(cxr_C.R1_cd_intra)

#################
######Doing the inter estimates
##################

R1_cd<-list(SR1= subset(CXR_C_Cd, Rep==1 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR2= subset(CXR_C_Cd, Rep==1 & Focal=="SR2"& Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==1 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_Cd, Rep==1 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R2_cd<-list(SR1= subset(CXR_C_Cd, Rep==2 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==2 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1")], SR5= subset(CXR_C_Cd, Rep==2 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1")])

R3_cd<-list(SR1= subset(CXR_C_Cd, Rep==3 & Focal=="SR1" & Comp!="SR1")[,c("fitness",  "SR4", "SR5")], SR2= subset(CXR_C_Cd, Rep==3 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==3 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_Cd, Rep==3 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R4_cd<-list(SR1= subset(CXR_C_Cd, Rep==4 & Focal=="SR1" & Comp!="SR1")[,c("fitness",  "SR4", "SR5")], SR2= subset(CXR_C_Cd, Rep==4 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==4 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_Cd, Rep==4 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R5_cd<-list(SR1= subset(CXR_C_Cd, Rep==5 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR2= subset(CXR_C_Cd, Rep==5 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==5 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_Cd, Rep==5 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])




cxr_C.R1_cd<-cxr_pm_multifit(data = R1_cd,
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = 0.1),
                             fixed_terms = fixed_terms_C_1N,
                             # no standard errors
                             bootstrap_samples = 200)

cxr_C.R2_cd_sr1<-cxr_pm_fit(data = R2_cd[[1]],
                            focal_column = NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(alpha_inter = 0.1),
                            fixed_terms = fixed_terms_C_2N[[1]],
                            # no standard errors
                            bootstrap_samples = 200)

cxr_C.R2_cd_sr4<-cxr_pm_fit(data = R2_cd[[2]],
                            focal_column = NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "global",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(alpha_inter = 0.1),
                            fixed_terms = fixed_terms_C_2N[[2]],
                            # no standard errors
                            bootstrap_samples = 200)

cxr_C.R2_cd_sr5<-cxr_pm_fit(data = R2_cd[[3]],
                            focal_column = NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "global",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(alpha_inter = 0.1),
                            fixed_terms = fixed_terms_C_2N[[3]],
                            # no standard errors
                            bootstrap_samples = 200)

cxr_C.R3_cd<-cxr_pm_multifit(data = R3_cd,
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = 0.1),
                             fixed_terms = fixed_terms_C_3N,
                             # no standard errors
                             bootstrap_samples = 200)

cxr_C.R4_cd<-cxr_pm_multifit(data = R4_cd,
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = 0.1),
                             fixed_terms = fixed_terms_C_4N,
                             # no standard errors
                             bootstrap_samples = 200)

cxr_C.R5_cd<-cxr_pm_multifit(data = R5_cd,
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = 0.1),
                             fixed_terms = fixed_terms_C_5N,
                             # no standard errors
                             bootstrap_samples = 200)

cxr_C.R1_cd$alpha_matrix



#rows in the alpha element of the returning list correspond to species i and columns to species j for each αij coefficient.

###### data table summary water



cxr_param_C<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_C$Tu_lambda<-0
cxr_param_C$Te_lambda<-0
cxr_param_C$Tu_intra<-0
cxr_param_C$Te_intra<-0
cxr_param_C$Tu_inter<-0
cxr_param_C$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_C<-cxr_param_C[-which(cxr_param_C$Replicate==2 & cxr_param_C$Tu_Regime=="SR2"),]


cxr_param_C[which(cxr_param_C$Replicate==1),"Tu_lambda"]<-c(cxr_C.R1_intra$fixed_terms[[1]]$lambda,cxr_C.R1_intra$fixed_terms[[2]]$lambda)
cxr_param_C[which(cxr_param_C$Replicate==1),"Te_lambda"]<-c(cxr_C.R1_intra$fixed_terms[[3]]$lambda,cxr_C.R1_intra$fixed_terms[[3]]$lambda, cxr_C.R1_intra$fixed_terms[[4]]$lambda,cxr_C.R1_intra$fixed_terms[[4]]$lambda)

cxr_param_C[which(cxr_param_C$Replicate==2),"Tu_lambda"]<-c(cxr_C.R2_intra$fixed_terms[[1]]$lambda,cxr_C.R2_intra$fixed_terms[[1]]$lambda)
cxr_param_C[which(cxr_param_C$Replicate==2),"Te_lambda"]<-c(cxr_C.R2_intra$fixed_terms[[2]]$lambda,cxr_C.R2_intra$fixed_terms[[3]]$lambda)

cxr_param_C[which(cxr_param_C$Replicate==3),"Tu_lambda"]<-c(cxr_C.R3_intra$fixed_terms[[1]]$lambda,cxr_C.R3_intra$fixed_terms[[2]]$lambda)
cxr_param_C[which(cxr_param_C$Replicate==3),"Te_lambda"]<-c(cxr_C.R3_intra$fixed_terms[[3]]$lambda,cxr_C.R3_intra$fixed_terms[[3]]$lambda, cxr_C.R3_intra$fixed_terms[[4]]$lambda,cxr_C.R3_intra$fixed_terms[[4]]$lambda)

cxr_param_C[which(cxr_param_C$Replicate==4),"Tu_lambda"]<-c(cxr_C.R4_intra$fixed_terms[[1]]$lambda,cxr_C.R4_intra$fixed_terms[[2]]$lambda)
cxr_param_C[which(cxr_param_C$Replicate==4),"Te_lambda"]<-c(cxr_C.R4_intra$fixed_terms[[3]]$lambda,cxr_C.R4_intra$fixed_terms[[3]]$lambda, cxr_C.R4_intra$fixed_terms[[4]]$lambda,cxr_C.R4_intra$fixed_terms[[4]]$lambda)

cxr_param_C[which(cxr_param_C$Replicate==5),"Tu_lambda"]<-c(cxr_C.R5_intra$fixed_terms[[1]]$lambda,cxr_C.R5_intra$fixed_terms[[2]]$lambda)
cxr_param_C[which(cxr_param_C$Replicate==5),"Te_lambda"]<-c(cxr_C.R5_intra$fixed_terms[[3]]$lambda,cxr_C.R5_intra$fixed_terms[[3]]$lambda, cxr_C.R5_intra$fixed_terms[[4]]$lambda,cxr_C.R5_intra$fixed_terms[[4]]$lambda)


cxr_param_C[which(cxr_param_C$Replicate==1),"Tu_intra"]<-rep(c(cxr_C.R1_intra$alpha_matrix[1,1], cxr_C.R1_intra$alpha_matrix[2,1]), 2)
cxr_param_C[which(cxr_param_C$Replicate==1),"Te_intra"]<-rep(c(cxr_C.R1_intra$alpha_matrix[3,1], cxr_C.R1_intra$alpha_matrix[4,1]), each=2)

cxr_param_C[which(cxr_param_C$Replicate==2),"Tu_intra"]<-rep(c(cxr_C.R2_intra$alpha_matrix[1,1]), 2)
cxr_param_C[which(cxr_param_C$Replicate==2),"Te_intra"]<-rep(c(cxr_C.R2_intra$alpha_matrix[2,1], cxr_C.R2_intra$alpha_matrix[3,1]))

cxr_param_C[which(cxr_param_C$Replicate==3),"Tu_intra"]<-rep(c(cxr_C.R3_intra$alpha_matrix[1,1], cxr_C.R3_intra$alpha_matrix[2,1]), 2)
cxr_param_C[which(cxr_param_C$Replicate==3),"Te_intra"]<-rep(c(cxr_C.R3_intra$alpha_matrix[3,1], cxr_C.R3_intra$alpha_matrix[4,1]), each=2)

cxr_param_C[which(cxr_param_C$Replicate==4),"Tu_intra"]<-rep(c(cxr_C.R4_intra$alpha_matrix[1,1], cxr_C.R4_intra$alpha_matrix[2,1]), 2)
cxr_param_C[which(cxr_param_C$Replicate==4),"Te_intra"]<-rep(c(cxr_C.R4_intra$alpha_matrix[3,1], cxr_C.R4_intra$alpha_matrix[4,1]), each=2)

cxr_param_C[which(cxr_param_C$Replicate==5),"Tu_intra"]<-rep(c(cxr_C.R5_intra$alpha_matrix[1,1], cxr_C.R5_intra$alpha_matrix[2,1]), 2)
cxr_param_C[which(cxr_param_C$Replicate==5),"Te_intra"]<-rep(c(cxr_C.R5_intra$alpha_matrix[3,1], cxr_C.R5_intra$alpha_matrix[4,1]), each=2)


cxr_param_C[which(cxr_param_C$Replicate==1),"Tu_inter"]<-c(cxr_C.R1$alpha_matrix[1,3], cxr_C.R1$alpha_matrix[2,3],cxr_C.R1$alpha_matrix[1,4], cxr_C.R1$alpha_matrix[2,4])
cxr_param_C[which(cxr_param_C$Replicate==1),"Te_inter"]<-c(cxr_C.R1$alpha_matrix[3,1], cxr_C.R1$alpha_matrix[3,2],cxr_C.R1$alpha_matrix[4,1], cxr_C.R1$alpha_matrix[4,2])

cxr_param_C[which(cxr_param_C$Replicate==2),"Tu_inter"]<-c(cxr_C.R2_sr1$alpha_inter[1], cxr_C.R2_sr1$alpha_inter[2])
cxr_param_C[which(cxr_param_C$Replicate==2),"Te_inter"]<-c(cxr_C.R2_sr4$alpha_inter[1],cxr_C.R2_sr1$alpha_inter[1])

cxr_param_C[which(cxr_param_C$Replicate==3),"Tu_inter"]<-c(cxr_C.R3$alpha_matrix[1,3], cxr_C.R3$alpha_matrix[2,3],cxr_C.R3$alpha_matrix[1,4], cxr_C.R3$alpha_matrix[2,4])
cxr_param_C[which(cxr_param_C$Replicate==3),"Te_inter"]<-c(cxr_C.R3$alpha_matrix[3,1], cxr_C.R3$alpha_matrix[3,2],cxr_C.R3$alpha_matrix[4,1], cxr_C.R3$alpha_matrix[4,2])

cxr_param_C[which(cxr_param_C$Replicate==4),"Tu_inter"]<-c(cxr_C.R4$alpha_matrix[1,3], cxr_C.R4$alpha_matrix[2,3],cxr_C.R4$alpha_matrix[1,4], cxr_C.R4$alpha_matrix[2,4])
cxr_param_C[which(cxr_param_C$Replicate==4),"Te_inter"]<-c(cxr_C.R4$alpha_matrix[3,1], cxr_C.R4$alpha_matrix[3,2],cxr_C.R4$alpha_matrix[4,1], cxr_C.R4$alpha_matrix[4,2])

cxr_param_C[which(cxr_param_C$Replicate==5),"Tu_inter"]<-c(cxr_C.R5$alpha_matrix[1,3], cxr_C.R5$alpha_matrix[2,3],cxr_C.R5$alpha_matrix[1,4], cxr_C.R5$alpha_matrix[2,4])
cxr_param_C[which(cxr_param_C$Replicate==5),"Te_inter"]<-c(cxr_C.R5$alpha_matrix[3,1], cxr_C.R5$alpha_matrix[3,2],cxr_C.R5$alpha_matrix[4,1], cxr_C.R5$alpha_matrix[4,2])

### Lower

cxr_param_C_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_C_lower$Tu_lambda<-0
cxr_param_C_lower$Te_lambda<-0
cxr_param_C_lower$Tu_intra<-0
cxr_param_C_lower$Te_intra<-0
cxr_param_C_lower$Tu_inter<-0
cxr_param_C_lower$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_C_lower<-cxr_param_C_lower[-which(cxr_param_C_lower$Replicate==2 & cxr_param_C_lower$Tu_Regime=="SR2"),]

#Since the error comes directly from the data we need to create some lists with that information
sd_1N <- list(list(lambda = subset(mean_dens1, Rep==1 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==1 & Env=="N" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==1 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==1 & Env=="N" & SR==5)$sd_lambda))

sd_2N <- list(list(lambda = subset(mean_dens1, Rep==2 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==2 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==2 & Env=="N" & SR==5)$sd_lambda))

sd_3N <- list(list(lambda = subset(mean_dens1, Rep==3 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==3 & Env=="N" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==3 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==3 & Env=="N" & SR==5)$sd_lambda))

sd_4N <- list(list(lambda = subset(mean_dens1, Rep==4 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==4 & Env=="N" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==4 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==4 & Env=="N" & SR==5)$sd_lambda))

sd_5N <- list(list(lambda = subset(mean_dens1, Rep==5 & Env=="N" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==5 & Env=="N" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==5 & Env=="N" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==5 & Env=="N" & SR==5)$sd_lambda))

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==1),"Tu_lambda"]<-c(cxr_C.R1$fixed_terms[[1]]$lambda-sd_1N[[1]]$lambda,cxr_C.R1$fixed_terms[[2]]$lambda-sd_1N[[2]]$lambda)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==1),"Te_lambda"]<-c(cxr_C.R1$fixed_terms[[3]]$lambda-sd_1N[[3]]$lambda,cxr_C.R1$fixed_terms[[3]]$lambda-sd_1N[[3]]$lambda, cxr_C.R1$fixed_terms[[4]]$lambda-sd_1N[[4]]$lambda,cxr_C.R1$fixed_terms[[4]]$lambda-sd_1N[[4]]$lambda)

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==2),"Tu_lambda"]<-c(cxr_C.R2_sr1$fixed_terms[[1]]-sd_2N[[1]]$lambda,cxr_C.R2_sr1$fixed_terms[[1]]-sd_2N[[1]]$lambda)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==2),"Te_lambda"]<-c(cxr_C.R2_sr4$fixed_terms[[1]]-sd_2N[[2]]$lambda,cxr_C.R2_sr5$fixed_terms[[1]]-sd_2N[[3]]$lambda)

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==3),"Tu_lambda"]<-c(cxr_C.R3$fixed_terms[[1]]$lambda-sd_3N[[1]]$lambda,cxr_C.R3$fixed_terms[[2]]$lambda-sd_3N[[2]]$lambda)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==3),"Te_lambda"]<-c(cxr_C.R3$fixed_terms[[3]]$lambda-sd_3N[[3]]$lambda,cxr_C.R3$fixed_terms[[3]]$lambda-sd_3N[[3]]$lambda, cxr_C.R3$fixed_terms[[4]]$lambda-sd_3N[[4]]$lambda,cxr_C.R3$fixed_terms[[4]]$lambda-sd_3N[[4]]$lambda)

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==4),"Tu_lambda"]<-c(cxr_C.R4$fixed_terms[[1]]$lambda-sd_4N[[1]]$lambda,cxr_C.R4$fixed_terms[[2]]$lambda-sd_4N[[2]]$lambda)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==4),"Te_lambda"]<-c(cxr_C.R4$fixed_terms[[3]]$lambda-sd_4N[[3]]$lambda,cxr_C.R4$fixed_terms[[3]]$lambda-sd_4N[[3]]$lambda, cxr_C.R4$fixed_terms[[4]]$lambda-sd_4N[[4]]$lambda,cxr_C.R4$fixed_terms[[4]]$lambda-sd_4N[[4]]$lambda)

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==5),"Tu_lambda"]<-c(cxr_C.R5$fixed_terms[[1]]$lambda-sd_5N[[1]]$lambda,cxr_C.R5$fixed_terms[[2]]$lambda-sd_5N[[2]]$lambda)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==5),"Te_lambda"]<-c(cxr_C.R5$fixed_terms[[3]]$lambda-sd_5N[[3]]$lambda,cxr_C.R5$fixed_terms[[3]]$lambda-sd_5N[[3]]$lambda, cxr_C.R5$fixed_terms[[4]]$lambda-sd_5N[[4]]$lambda,cxr_C.R5$fixed_terms[[4]]$lambda-sd_5N[[4]]$lambda)


cxr_param_C_lower[which(cxr_param_C_lower$Replicate==1),"Tu_intra"]<-rep(c(cxr_C.R1_intra$alpha_matrix[1,1]-cxr_C.R1_intra$alpha_matrix_standard_error[1,1], cxr_C.R1_intra$alpha_matrix[2,1]-cxr_C.R1_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==1),"Te_intra"]<-rep(c(cxr_C.R1_intra$alpha_matrix[3,1]-cxr_C.R1_intra$alpha_matrix_standard_error[3,1], cxr_C.R1_intra$alpha_matrix[4,1]-cxr_C.R1_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==2),"Tu_intra"]<-rep(c(cxr_C.R2_intra$alpha_matrix[1,1]-cxr_C.R2_intra$alpha_matrix_standard_error[1,1]), 2)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==2),"Te_intra"]<-c(cxr_C.R2_intra$alpha_matrix[2,1]-cxr_C.R2_intra$alpha_matrix_standard_error[2,1], cxr_C.R2_intra$alpha_matrix[3,1]-cxr_C.R2_intra$alpha_matrix_standard_error[3,1])

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==3),"Tu_intra"]<-rep(c(cxr_C.R3_intra$alpha_matrix[1,1]-cxr_C.R3_intra$alpha_matrix_standard_error[1,1], cxr_C.R3_intra$alpha_matrix[2,1]-cxr_C.R3_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==3),"Te_intra"]<-rep(c(cxr_C.R3_intra$alpha_matrix[3,1]-cxr_C.R3_intra$alpha_matrix_standard_error[3,1], cxr_C.R3_intra$alpha_matrix[4,1]-cxr_C.R3_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==4),"Tu_intra"]<-rep(c(cxr_C.R4_intra$alpha_matrix[1,1]-cxr_C.R4_intra$alpha_matrix_standard_error[1,1], cxr_C.R4_intra$alpha_matrix[2,1]-cxr_C.R4_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==4),"Te_intra"]<-rep(c(cxr_C.R4_intra$alpha_matrix[3,1]-cxr_C.R4_intra$alpha_matrix_standard_error[3,1], cxr_C.R4_intra$alpha_matrix[4,1]-cxr_C.R4_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==5),"Tu_intra"]<-rep(c(cxr_C.R5_intra$alpha_matrix[1,1]-cxr_C.R5_intra$alpha_matrix_standard_error[1,1], cxr_C.R5_intra$alpha_matrix[2,1]-cxr_C.R5_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==5),"Te_intra"]<-rep(c(cxr_C.R5_intra$alpha_matrix[3,1]-cxr_C.R5_intra$alpha_matrix_standard_error[3,1], cxr_C.R5_intra$alpha_matrix[4,1]-cxr_C.R5_intra$alpha_matrix_standard_error[4,1]), each=2)


cxr_param_C_lower[which(cxr_param_C_lower$Replicate==1),"Tu_inter"]<-c(cxr_C.R1$alpha_matrix[1,3]-cxr_C.R1$alpha_matrix_standard_error[1,3], cxr_C.R1$alpha_matrix[2,3]-cxr_C.R1$alpha_matrix_standard_error[2,3],cxr_C.R1$alpha_matrix[1,4]-cxr_C.R1$alpha_matrix_standard_error[1,4], cxr_C.R1$alpha_matrix[2,4]-cxr_C.R1$alpha_matrix_standard_error[2,4])
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==1),"Te_inter"]<-c(cxr_C.R1$alpha_matrix[3,1]-cxr_C.R1$alpha_matrix_standard_error[3,1], cxr_C.R1$alpha_matrix[3,2]-cxr_C.R1$alpha_matrix_standard_error[3,2],cxr_C.R1$alpha_matrix[4,1]-cxr_C.R1$alpha_matrix_standard_error[4,1], cxr_C.R1$alpha_matrix[4,2]-cxr_C.R1$alpha_matrix_standard_error[4,2])

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==2),"Tu_inter"]<-c(cxr_C.R2_sr1$alpha_inter[1]-cxr_C.R2_sr1$alpha_inter_standard_error[1], cxr_C.R2_sr1$alpha_inter[2]-cxr_C.R2_sr1$alpha_inter_standard_error[2])
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==2),"Te_inter"]<-c(cxr_C.R2_sr4$alpha_inter[1]-cxr_C.R2_sr4$alpha_inter_standard_error[1], cxr_C.R2_sr5$alpha_inter[1]-cxr_C.R2_sr5$alpha_inter_standard_error[1])

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==3),"Tu_inter"]<-c(cxr_C.R3$alpha_matrix[1,3]-cxr_C.R3$alpha_matrix_standard_error[1,3], cxr_C.R3$alpha_matrix[2,3]-cxr_C.R3$alpha_matrix_standard_error[2,3],cxr_C.R3$alpha_matrix[1,4]-cxr_C.R3$alpha_matrix_standard_error[1,4], cxr_C.R3$alpha_matrix[2,4]-cxr_C.R3$alpha_matrix_standard_error[2,4])
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==3),"Te_inter"]<-c(cxr_C.R3$alpha_matrix[3,1]-cxr_C.R3$alpha_matrix_standard_error[3,1], cxr_C.R3$alpha_matrix[3,2]-cxr_C.R3$alpha_matrix_standard_error[3,2],cxr_C.R3$alpha_matrix[4,1]-cxr_C.R3$alpha_matrix_standard_error[4,1], cxr_C.R3$alpha_matrix[4,2]-cxr_C.R3$alpha_matrix_standard_error[4,2])

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==4),"Tu_inter"]<-c(cxr_C.R4$alpha_matrix[1,3]-cxr_C.R4$alpha_matrix_standard_error[1,3], cxr_C.R4$alpha_matrix[2,3]-cxr_C.R4$alpha_matrix_standard_error[2,3],cxr_C.R4$alpha_matrix[1,4]-cxr_C.R4$alpha_matrix_standard_error[1,4], cxr_C.R4$alpha_matrix[2,4]-cxr_C.R4$alpha_matrix_standard_error[2,4])
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==4),"Te_inter"]<-c(cxr_C.R4$alpha_matrix[3,1]-cxr_C.R4$alpha_matrix_standard_error[3,1], cxr_C.R4$alpha_matrix[3,2]-cxr_C.R4$alpha_matrix_standard_error[3,2],cxr_C.R4$alpha_matrix[4,1]-cxr_C.R4$alpha_matrix_standard_error[4,1], cxr_C.R4$alpha_matrix[4,2]-cxr_C.R4$alpha_matrix_standard_error[4,2])

cxr_param_C_lower[which(cxr_param_C_lower$Replicate==5),"Tu_inter"]<-c(cxr_C.R5$alpha_matrix[1,3]-cxr_C.R5$alpha_matrix_standard_error[1,3], cxr_C.R5$alpha_matrix[2,3]-cxr_C.R5$alpha_matrix_standard_error[2,3],cxr_C.R5$alpha_matrix[1,4]-cxr_C.R5$alpha_matrix_standard_error[1,4], cxr_C.R5$alpha_matrix[2,4]-cxr_C.R5$alpha_matrix_standard_error[2,4])
cxr_param_C_lower[which(cxr_param_C_lower$Replicate==5),"Te_inter"]<-c(cxr_C.R5$alpha_matrix[3,1]-cxr_C.R5$alpha_matrix_standard_error[3,1], cxr_C.R5$alpha_matrix[3,2]-cxr_C.R5$alpha_matrix_standard_error[3,2],cxr_C.R5$alpha_matrix[4,1]-cxr_C.R5$alpha_matrix_standard_error[4,1], cxr_C.R5$alpha_matrix[4,2]-cxr_C.R5$alpha_matrix_standard_error[4,2])

### upper

cxr_param_C_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_C_upper$Tu_lambda<-0
cxr_param_C_upper$Te_lambda<-0
cxr_param_C_upper$Tu_intra<-0
cxr_param_C_upper$Te_intra<-0
cxr_param_C_upper$Tu_inter<-0
cxr_param_C_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_C_upper<-cxr_param_C_upper[-which(cxr_param_C_upper$Replicate==2 & cxr_param_C_upper$Tu_Regime=="SR2"),]


cxr_param_C_upper[which(cxr_param_C_upper$Replicate==1),"Tu_lambda"]<-c(cxr_C.R1$fixed_terms[[1]]$lambda+sd_1N[[1]]$lambda,cxr_C.R1$fixed_terms[[2]]$lambda+sd_1N[[2]]$lambda)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==1),"Te_lambda"]<-c(cxr_C.R1$fixed_terms[[3]]$lambda+sd_1N[[3]]$lambda,cxr_C.R1$fixed_terms[[3]]$lambda+sd_1N[[3]]$lambda, cxr_C.R1$fixed_terms[[4]]$lambda+sd_1N[[4]]$lambda,cxr_C.R1$fixed_terms[[4]]$lambda+sd_1N[[4]]$lambda)

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==2),"Tu_lambda"]<-c(cxr_C.R2_sr1$fixed_terms[[1]]+sd_2N[[1]]$lambda,cxr_C.R2_sr1$fixed_terms[[1]]+sd_2N[[1]]$lambda)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==2),"Te_lambda"]<-c(cxr_C.R2_sr4$fixed_terms[[1]]+sd_2N[[2]]$lambda,cxr_C.R2_sr5$fixed_terms[[1]]+sd_2N[[3]]$lambda)

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==3),"Tu_lambda"]<-c(cxr_C.R3$fixed_terms[[1]]$lambda+sd_3N[[1]]$lambda,cxr_C.R3$fixed_terms[[2]]$lambda+sd_3N[[2]]$lambda)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==3),"Te_lambda"]<-c(cxr_C.R3$fixed_terms[[3]]$lambda+sd_3N[[3]]$lambda,cxr_C.R3$fixed_terms[[3]]$lambda+sd_3N[[3]]$lambda, cxr_C.R3$fixed_terms[[4]]$lambda+sd_3N[[4]]$lambda,cxr_C.R3$fixed_terms[[4]]$lambda+sd_3N[[4]]$lambda)

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==4),"Tu_lambda"]<-c(cxr_C.R4$fixed_terms[[1]]$lambda+sd_4N[[1]]$lambda,cxr_C.R4$fixed_terms[[2]]$lambda+sd_4N[[2]]$lambda)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==4),"Te_lambda"]<-c(cxr_C.R4$fixed_terms[[3]]$lambda+sd_4N[[3]]$lambda,cxr_C.R4$fixed_terms[[3]]$lambda+sd_4N[[3]]$lambda, cxr_C.R4$fixed_terms[[4]]$lambda+sd_4N[[4]]$lambda,cxr_C.R4$fixed_terms[[4]]$lambda+sd_4N[[4]]$lambda)

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==5),"Tu_lambda"]<-c(cxr_C.R5$fixed_terms[[1]]$lambda+sd_5N[[1]]$lambda,cxr_C.R5$fixed_terms[[2]]$lambda+sd_5N[[2]]$lambda)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==5),"Te_lambda"]<-c(cxr_C.R5$fixed_terms[[3]]$lambda+sd_5N[[3]]$lambda,cxr_C.R5$fixed_terms[[3]]$lambda+sd_5N[[3]]$lambda, cxr_C.R5$fixed_terms[[4]]$lambda+sd_5N[[4]]$lambda,cxr_C.R5$fixed_terms[[4]]$lambda+sd_5N[[4]]$lambda)


cxr_param_C_upper[which(cxr_param_C_upper$Replicate==1),"Tu_intra"]<-rep(c(cxr_C.R1_intra$alpha_matrix[1,1]+cxr_C.R1_intra$alpha_matrix_standard_error[1,1], cxr_C.R1_intra$alpha_matrix[2,1]+cxr_C.R1_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==1),"Te_intra"]<-rep(c(cxr_C.R1_intra$alpha_matrix[3,1]+cxr_C.R1_intra$alpha_matrix_standard_error[3,1], cxr_C.R1_intra$alpha_matrix[4,1]+cxr_C.R1_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==2),"Tu_intra"]<-rep(c(cxr_C.R2_intra$alpha_matrix[1,1]+cxr_C.R2_intra$alpha_matrix_standard_error[1,1]), 2)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==2),"Te_intra"]<-c(cxr_C.R2_intra$alpha_matrix[2,1]+cxr_C.R2_intra$alpha_matrix_standard_error[2,1], cxr_C.R2_intra$alpha_matrix[3,1]+cxr_C.R2_intra$alpha_matrix_standard_error[3,1])

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==3),"Tu_intra"]<-rep(c(cxr_C.R3_intra$alpha_matrix[1,1]+cxr_C.R3_intra$alpha_matrix_standard_error[1,1], cxr_C.R3_intra$alpha_matrix[2,1]+cxr_C.R3_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==3),"Te_intra"]<-rep(c(cxr_C.R3_intra$alpha_matrix[3,1]+cxr_C.R3_intra$alpha_matrix_standard_error[3,1], cxr_C.R3_intra$alpha_matrix[4,1]+cxr_C.R3_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==4),"Tu_intra"]<-rep(c(cxr_C.R4_intra$alpha_matrix[1,1]+cxr_C.R4_intra$alpha_matrix_standard_error[1,1], cxr_C.R4_intra$alpha_matrix[2,1]+cxr_C.R4_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==4),"Te_intra"]<-rep(c(cxr_C.R4_intra$alpha_matrix[3,1]+cxr_C.R4_intra$alpha_matrix_standard_error[3,1], cxr_C.R4_intra$alpha_matrix[4,1]+cxr_C.R4_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==5),"Tu_intra"]<-rep(c(cxr_C.R5_intra$alpha_matrix[1,1]+cxr_C.R5_intra$alpha_matrix_standard_error[1,1], cxr_C.R5_intra$alpha_matrix[2,1]+cxr_C.R5_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==5),"Te_intra"]<-rep(c(cxr_C.R5_intra$alpha_matrix[3,1]+cxr_C.R5_intra$alpha_matrix_standard_error[3,1], cxr_C.R5_intra$alpha_matrix[4,1]+cxr_C.R5_intra$alpha_matrix_standard_error[4,1]), each=2)


cxr_param_C_upper[which(cxr_param_C_upper$Replicate==1),"Tu_inter"]<-c(cxr_C.R1$alpha_matrix[1,3]+cxr_C.R1$alpha_matrix_standard_error[1,3], cxr_C.R1$alpha_matrix[2,3]+cxr_C.R1$alpha_matrix_standard_error[2,3],cxr_C.R1$alpha_matrix[1,4]+cxr_C.R1$alpha_matrix_standard_error[1,4], cxr_C.R1$alpha_matrix[2,4]+cxr_C.R1$alpha_matrix_standard_error[2,4])
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==1),"Te_inter"]<-c(cxr_C.R1$alpha_matrix[3,1]+cxr_C.R1$alpha_matrix_standard_error[3,1], cxr_C.R1$alpha_matrix[3,2]+cxr_C.R1$alpha_matrix_standard_error[3,2],cxr_C.R1$alpha_matrix[4,1]+cxr_C.R1$alpha_matrix_standard_error[4,1], cxr_C.R1$alpha_matrix[4,2]+cxr_C.R1$alpha_matrix_standard_error[4,2])

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==2),"Tu_inter"]<-c(cxr_C.R2_sr1$alpha_inter[1]+cxr_C.R2_sr1$alpha_inter_standard_error[1], cxr_C.R2_sr1$alpha_inter[2]+cxr_C.R2_sr1$alpha_inter_standard_error[2])
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==2),"Te_inter"]<-c(cxr_C.R2_sr4$alpha_inter[1]+cxr_C.R2_sr4$alpha_inter_standard_error[1], cxr_C.R2_sr5$alpha_inter[1]+cxr_C.R2_sr5$alpha_inter_standard_error[1])

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==3),"Tu_inter"]<-c(cxr_C.R3$alpha_matrix[1,3]+cxr_C.R3$alpha_matrix_standard_error[1,3], cxr_C.R3$alpha_matrix[2,3]+cxr_C.R3$alpha_matrix_standard_error[2,3],cxr_C.R3$alpha_matrix[1,4]+cxr_C.R3$alpha_matrix_standard_error[1,4], cxr_C.R3$alpha_matrix[2,4]+cxr_C.R3$alpha_matrix_standard_error[2,4])
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==3),"Te_inter"]<-c(cxr_C.R3$alpha_matrix[3,1]+cxr_C.R3$alpha_matrix_standard_error[3,1], cxr_C.R3$alpha_matrix[3,2]+cxr_C.R3$alpha_matrix_standard_error[3,2],cxr_C.R3$alpha_matrix[4,1]+cxr_C.R3$alpha_matrix_standard_error[4,1], cxr_C.R3$alpha_matrix[4,2]+cxr_C.R3$alpha_matrix_standard_error[4,2])

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==4),"Tu_inter"]<-c(cxr_C.R4$alpha_matrix[1,3]+cxr_C.R4$alpha_matrix_standard_error[1,3], cxr_C.R4$alpha_matrix[2,3]+cxr_C.R4$alpha_matrix_standard_error[2,3],cxr_C.R4$alpha_matrix[1,4]+cxr_C.R4$alpha_matrix_standard_error[1,4], cxr_C.R4$alpha_matrix[2,4]+cxr_C.R4$alpha_matrix_standard_error[2,4])
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==4),"Te_inter"]<-c(cxr_C.R4$alpha_matrix[3,1]+cxr_C.R4$alpha_matrix_standard_error[3,1], cxr_C.R4$alpha_matrix[3,2]+cxr_C.R4$alpha_matrix_standard_error[3,2],cxr_C.R4$alpha_matrix[4,1]+cxr_C.R4$alpha_matrix_standard_error[4,1], cxr_C.R4$alpha_matrix[4,2]+cxr_C.R4$alpha_matrix_standard_error[4,2])

cxr_param_C_upper[which(cxr_param_C_upper$Replicate==5),"Tu_inter"]<-c(cxr_C.R5$alpha_matrix[1,3]+cxr_C.R5$alpha_matrix_standard_error[1,3], cxr_C.R5$alpha_matrix[2,3]+cxr_C.R5$alpha_matrix_standard_error[2,3],cxr_C.R5$alpha_matrix[1,4]+cxr_C.R5$alpha_matrix_standard_error[1,4], cxr_C.R5$alpha_matrix[2,4]+cxr_C.R5$alpha_matrix_standard_error[2,4])
cxr_param_C_upper[which(cxr_param_C_upper$Replicate==5),"Te_inter"]<-c(cxr_C.R5$alpha_matrix[3,1]+cxr_C.R5$alpha_matrix_standard_error[3,1], cxr_C.R5$alpha_matrix[3,2]+cxr_C.R5$alpha_matrix_standard_error[3,2],cxr_C.R5$alpha_matrix[4,1]+cxr_C.R5$alpha_matrix_standard_error[4,1], cxr_C.R5$alpha_matrix[4,2]+cxr_C.R5$alpha_matrix_standard_error[4,2])



###### data table summary cadmium
cxr_param_CC<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_CC$Tu_lambda<-0
cxr_param_CC$Te_lambda<-0
cxr_param_CC$Tu_intra<-0
cxr_param_CC$Te_intra<-0
cxr_param_CC$Tu_inter<-0
cxr_param_CC$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_CC<-cxr_param_CC[-which(cxr_param_CC$Replicate==2 & cxr_param_CC$Tu_Regime=="SR2"),]


cxr_param_CC[which(cxr_param_CC$Replicate==1),"Tu_lambda"]<-c(cxr_C.R1_cd$fixed_terms[[1]]$lambda,cxr_C.R1_cd$fixed_terms[[2]]$lambda)
cxr_param_CC[which(cxr_param_CC$Replicate==1),"Te_lambda"]<-c(cxr_C.R1_cd$fixed_terms[[3]]$lambda,cxr_C.R1_cd$fixed_terms[[3]]$lambda, cxr_C.R1_cd$fixed_terms[[4]]$lambda,cxr_C.R1_cd$fixed_terms[[4]]$lambda)

cxr_param_CC[which(cxr_param_CC$Replicate==2),"Tu_lambda"]<-c(cxr_C.R2_cd_sr1$fixed_terms[[1]],cxr_C.R2_cd_sr1$fixed_terms[[1]])
cxr_param_CC[which(cxr_param_CC$Replicate==2),"Te_lambda"]<-c(cxr_C.R2_cd_sr4$fixed_terms[[1]], cxr_C.R2_cd_sr5$fixed_terms[[1]])

cxr_param_CC[which(cxr_param_CC$Replicate==3),"Tu_lambda"]<-c(cxr_C.R3_cd$fixed_terms[[1]]$lambda,cxr_C.R3_cd$fixed_terms[[2]]$lambda)
cxr_param_CC[which(cxr_param_CC$Replicate==3),"Te_lambda"]<-c(cxr_C.R3_cd$fixed_terms[[3]]$lambda,cxr_C.R3_cd$fixed_terms[[3]]$lambda, cxr_C.R3_cd$fixed_terms[[4]]$lambda,cxr_C.R3_cd$fixed_terms[[4]]$lambda)

cxr_param_CC[which(cxr_param_CC$Replicate==4),"Tu_lambda"]<-c(cxr_C.R4_cd$fixed_terms[[1]]$lambda,cxr_C.R4_cd$fixed_terms[[2]]$lambda)
cxr_param_CC[which(cxr_param_CC$Replicate==4),"Te_lambda"]<-c(cxr_C.R4_cd$fixed_terms[[3]]$lambda,cxr_C.R4_cd$fixed_terms[[3]]$lambda, cxr_C.R4_cd$fixed_terms[[4]]$lambda,cxr_C.R4_cd$fixed_terms[[4]]$lambda)

cxr_param_CC[which(cxr_param_CC$Replicate==5),"Tu_lambda"]<-c(cxr_C.R5_cd$fixed_terms[[1]]$lambda,cxr_C.R5_cd$fixed_terms[[2]]$lambda)
cxr_param_CC[which(cxr_param_CC$Replicate==5),"Te_lambda"]<-c(cxr_C.R5_cd$fixed_terms[[3]]$lambda,cxr_C.R5_cd$fixed_terms[[3]]$lambda, cxr_C.R5_cd$fixed_terms[[4]]$lambda,cxr_C.R5_cd$fixed_terms[[4]]$lambda)


cxr_param_CC[which(cxr_param_CC$Replicate==1),"Tu_intra"]<-rep(c(cxr_C.R1_cd_intra$alpha_matrix[1,1], cxr_C.R1_cd_intra$alpha_matrix[2,1]), 2)
cxr_param_CC[which(cxr_param_CC$Replicate==1),"Te_intra"]<-rep(c(cxr_C.R1_cd_intra$alpha_matrix[3,1], cxr_C.R1_cd_intra$alpha_matrix[4,1]), each=2)

cxr_param_CC[which(cxr_param_CC$Replicate==2),"Tu_intra"]<-rep(c(cxr_C.R2_cd_intra$alpha_matrix[1,1]), 2)
cxr_param_CC[which(cxr_param_CC$Replicate==2),"Te_intra"]<-c(cxr_C.R2_cd_intra$alpha_matrix[2,1], cxr_C.R2_cd_intra$alpha_matrix[3,1])

cxr_param_CC[which(cxr_param_CC$Replicate==3),"Tu_intra"]<-rep(c(cxr_C.R3_cd_intra$alpha_matrix[1,1], cxr_C.R3_cd_intra$alpha_matrix[2,1]), 2)
cxr_param_CC[which(cxr_param_CC$Replicate==3),"Te_intra"]<-rep(c(cxr_C.R3_cd_intra$alpha_matrix[3,1], cxr_C.R3_cd_intra$alpha_matrix[4,1]), each=2)

cxr_param_CC[which(cxr_param_CC$Replicate==4),"Tu_intra"]<-rep(c(cxr_C.R4_cd_intra$alpha_matrix[1,1], cxr_C.R4_cd_intra$alpha_matrix[2,1]), 2)
cxr_param_CC[which(cxr_param_CC$Replicate==4),"Te_intra"]<-rep(c(cxr_C.R4_cd_intra$alpha_matrix[3,1], cxr_C.R4_cd_intra$alpha_matrix[4,1]), each=2)

cxr_param_CC[which(cxr_param_CC$Replicate==5),"Tu_intra"]<-rep(c(cxr_C.R5_cd_intra$alpha_matrix[1,1], cxr_C.R5_cd_intra$alpha_matrix[2,1]), 2)
cxr_param_CC[which(cxr_param_CC$Replicate==5),"Te_intra"]<-rep(c(cxr_C.R5_cd_intra$alpha_matrix[3,1], cxr_C.R5_cd_intra$alpha_matrix[4,1]), each=2)


cxr_param_CC[which(cxr_param_CC$Replicate==1),"Tu_inter"]<-c(cxr_C.R1_cd$alpha_matrix[1,3], cxr_C.R1_cd$alpha_matrix[2,3],cxr_C.R1_cd$alpha_matrix[1,4], cxr_C.R1_cd$alpha_matrix[2,4])
cxr_param_CC[which(cxr_param_CC$Replicate==1),"Te_inter"]<-c(cxr_C.R1_cd$alpha_matrix[3,1], cxr_C.R1_cd$alpha_matrix[3,2],cxr_C.R1_cd$alpha_matrix[4,1], cxr_C.R1_cd$alpha_matrix[4,2])

cxr_param_CC[which(cxr_param_CC$Replicate==2),"Tu_inter"]<-c(cxr_C.R2_cd_sr1$alpha_inter[1], cxr_C.R2_cd_sr1$alpha_inter[2])
cxr_param_CC[which(cxr_param_CC$Replicate==2),"Te_inter"]<-c(cxr_C.R2_cd_sr4$alpha_inter[1], cxr_C.R2_cd_sr5$alpha_inter[1])

cxr_param_CC[which(cxr_param_CC$Replicate==3),"Tu_inter"]<-c(cxr_C.R3_cd$alpha_matrix[1,3], cxr_C.R3_cd$alpha_matrix[2,3],cxr_C.R3_cd$alpha_matrix[1,4], cxr_C.R3_cd$alpha_matrix[2,4])
cxr_param_CC[which(cxr_param_CC$Replicate==3),"Te_inter"]<-c(cxr_C.R3_cd$alpha_matrix[3,1], cxr_C.R3_cd$alpha_matrix[3,2],cxr_C.R3_cd$alpha_matrix[4,1], cxr_C.R3_cd$alpha_matrix[4,2])

cxr_param_CC[which(cxr_param_CC$Replicate==4),"Tu_inter"]<-c(cxr_C.R4_cd$alpha_matrix[1,3], cxr_C.R4_cd$alpha_matrix[2,3],cxr_C.R4_cd$alpha_matrix[1,4], cxr_C.R4_cd$alpha_matrix[2,4])
cxr_param_CC[which(cxr_param_CC$Replicate==4),"Te_inter"]<-c(cxr_C.R4_cd$alpha_matrix[3,1], cxr_C.R4_cd$alpha_matrix[3,2],cxr_C.R4_cd$alpha_matrix[4,1], cxr_C.R4_cd$alpha_matrix[4,2])

cxr_param_CC[which(cxr_param_CC$Replicate==5),"Tu_inter"]<-c(cxr_C.R5_cd$alpha_matrix[1,3], cxr_C.R5_cd$alpha_matrix[2,3],cxr_C.R5_cd$alpha_matrix[1,4], cxr_C.R5_cd$alpha_matrix[2,4])
cxr_param_CC[which(cxr_param_CC$Replicate==5),"Te_inter"]<-c(cxr_C.R5_cd$alpha_matrix[3,1], cxr_C.R5_cd$alpha_matrix[3,2],cxr_C.R5_cd$alpha_matrix[4,1], cxr_C.R5_cd$alpha_matrix[4,2])

### Lower

cxr_param_CC_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_CC_lower$Tu_lambda<-0
cxr_param_CC_lower$Te_lambda<-0
cxr_param_CC_lower$Tu_intra<-0
cxr_param_CC_lower$Te_intra<-0
cxr_param_CC_lower$Tu_inter<-0
cxr_param_CC_lower$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_CC_lower<-cxr_param_CC_lower[-which(cxr_param_CC_lower$Replicate==2 & cxr_param_CC_lower$Tu_Regime=="SR2"),]

#Since the error comes directly from the data we need to create some lists with that information
sd_1C <- list(list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==1 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==1 & Env=="Cd" & SR==5)$sd_lambda))

sd_2C <- list(list(lambda = subset(mean_dens1, Rep==2 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==2 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==2 & Env=="Cd" & SR==5)$sd_lambda))

sd_3C <- list(list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==3 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==3 & Env=="Cd" & SR==5)$sd_lambda))

sd_4C <- list(list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==4 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==4 & Env=="Cd" & SR==5)$sd_lambda))

sd_5C <- list(list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==1)$sd_lambda ), # focal sp 1
              list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==2)$sd_lambda), # focal sp 2
              list(lambda = subset(mean_dens1, Rep==5 & Env=="Cd" & SR==4)$sd_lambda),
              list(lambda= subset(mean_dens1, Rep==5 & Env=="Cd" & SR==5)$sd_lambda))

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==1),"Tu_lambda"]<-c(cxr_C.R1_cd$fixed_terms[[1]]$lambda-sd_1C[[1]]$lambda,cxr_C.R1_cd$fixed_terms[[2]]$lambda-sd_1C[[2]]$lambda)
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==1),"Te_lambda"]<-c(cxr_C.R1_cd$fixed_terms[[3]]$lambda-sd_1C[[3]]$lambda,cxr_C.R1_cd$fixed_terms[[3]]$lambda-sd_1C[[3]]$lambda, cxr_C.R1_cd$fixed_terms[[4]]$lambda-sd_1C[[4]]$lambda,cxr_C.R1_cd$fixed_terms[[4]]$lambda-sd_1C[[4]]$lambda)

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==2),"Tu_lambda"]<-c(cxr_C.R2_cd_sr1$fixed_terms$lambda-sd_2C[[1]]$lambda,cxr_C.R2_cd_sr1$fixed_terms$lambda-sd_2C[[1]]$lambda)
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==2),"Te_lambda"]<-c(cxr_C.R2_cd_sr4$fixed_terms$lambda-sd_2C[[2]]$lambda,cxr_C.R2_cd_sr5$fixed_terms$lambda-sd_2C[[3]]$lambda)

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==3),"Tu_lambda"]<-c(cxr_C.R3_cd$fixed_terms[[1]]$lambda-sd_3C[[1]]$lambda,cxr_C.R3_cd$fixed_terms[[2]]$lambda-sd_3C[[2]]$lambda)
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==3),"Te_lambda"]<-c(cxr_C.R3_cd$fixed_terms[[3]]$lambda-sd_3C[[3]]$lambda,cxr_C.R3_cd$fixed_terms[[3]]$lambda-sd_3C[[3]]$lambda, cxr_C.R3_cd$fixed_terms[[4]]$lambda-sd_3C[[4]]$lambda,cxr_C.R3_cd$fixed_terms[[4]]$lambda-sd_3C[[4]]$lambda)

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==4),"Tu_lambda"]<-c(cxr_C.R4_cd$fixed_terms[[1]]$lambda-sd_4C[[1]]$lambda,cxr_C.R4_cd$fixed_terms[[2]]$lambda-sd_4C[[2]]$lambda)
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==4),"Te_lambda"]<-c(cxr_C.R4_cd$fixed_terms[[3]]$lambda-sd_4C[[3]]$lambda,cxr_C.R4_cd$fixed_terms[[3]]$lambda-sd_4C[[3]]$lambda, cxr_C.R4_cd$fixed_terms[[4]]$lambda-sd_4C[[4]]$lambda,cxr_C.R4_cd$fixed_terms[[4]]$lambda-sd_4C[[4]]$lambda)

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==5),"Tu_lambda"]<-c(cxr_C.R5_cd$fixed_terms[[1]]$lambda-sd_5C[[1]]$lambda,cxr_C.R5_cd$fixed_terms[[2]]$lambda-sd_5C[[2]]$lambda)
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==5),"Te_lambda"]<-c(cxr_C.R5_cd$fixed_terms[[3]]$lambda-sd_5C[[3]]$lambda,cxr_C.R5_cd$fixed_terms[[3]]$lambda-sd_5C[[3]]$lambda, cxr_C.R5_cd$fixed_terms[[4]]$lambda-sd_5C[[4]]$lambda,cxr_C.R5_cd$fixed_terms[[4]]$lambda-sd_5C[[4]]$lambda)


cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==1),"Tu_intra"]<-rep(c(cxr_C.R1_cd_intra$alpha_matrix[1,1]-cxr_C.R1_cd_intra$alpha_matrix_standard_error[1,1], cxr_C.R1_cd_intra$alpha_matrix[2,1]-cxr_C.R1_cd_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==1),"Te_intra"]<-rep(c(cxr_C.R1_cd_intra$alpha_matrix[3,1]-cxr_C.R1_cd_intra$alpha_matrix_standard_error[3,1], cxr_C.R1_cd_intra$alpha_matrix[4,1]-cxr_C.R1_cd_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==2),"Tu_intra"]<-cxr_C.R2_cd_intra$alpha_matrix[1,1]-cxr_C.R2_cd_intra$alpha_matrix_standard_error[1,1]
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==2),"Te_intra"]<-c(cxr_C.R2_cd_intra$alpha_matrix[2,1]-cxr_C.R2_cd_intra$alpha_matrix_standard_error[2,1], cxr_C.R2_cd_intra$alpha_matrix[3,1]-cxr_C.R2_cd_intra$alpha_matrix_standard_error[3,1])

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==3),"Tu_intra"]<-rep(c(cxr_C.R3_cd_intra$alpha_matrix[1,1]-cxr_C.R3_cd_intra$alpha_matrix_standard_error[1,1], cxr_C.R3_cd_intra$alpha_matrix[2,1]-cxr_C.R3_cd_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==3),"Te_intra"]<-rep(c(cxr_C.R3_cd_intra$alpha_matrix[3,1]-cxr_C.R3_cd_intra$alpha_matrix_standard_error[3,1], cxr_C.R3_cd_intra$alpha_matrix[4,1]-cxr_C.R3_cd_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==4),"Tu_intra"]<-rep(c(cxr_C.R4_cd_intra$alpha_matrix[1,1]-cxr_C.R4_cd_intra$alpha_matrix_standard_error[1,1], cxr_C.R4_cd_intra$alpha_matrix[2,1]-cxr_C.R4_cd_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==4),"Te_intra"]<-rep(c(cxr_C.R4_cd_intra$alpha_matrix[3,1]-cxr_C.R4_cd_intra$alpha_matrix_standard_error[3,1], cxr_C.R4_cd_intra$alpha_matrix[4,1]-cxr_C.R4_cd_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==5),"Tu_intra"]<-rep(c(cxr_C.R5_cd_intra$alpha_matrix[1,1]-cxr_C.R5_cd_intra$alpha_matrix_standard_error[1,1], cxr_C.R5_cd_intra$alpha_matrix[2,1]-cxr_C.R5_cd_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==5),"Te_intra"]<-rep(c(cxr_C.R5_cd_intra$alpha_matrix[3,1]-cxr_C.R5_cd_intra$alpha_matrix_standard_error[3,1], cxr_C.R5_cd_intra$alpha_matrix[4,1]-cxr_C.R5_cd_intra$alpha_matrix_standard_error[4,1]), each=2)


cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==1),"Tu_inter"]<-c(cxr_C.R1_cd$alpha_matrix[1,3]-cxr_C.R1_cd$alpha_matrix_standard_error[1,3], cxr_C.R1_cd$alpha_matrix[2,3]-cxr_C.R1_cd$alpha_matrix_standard_error[2,3],cxr_C.R1_cd$alpha_matrix[1,4]-cxr_C.R1_cd$alpha_matrix_standard_error[1,4], cxr_C.R1_cd$alpha_matrix[2,4]-cxr_C.R1_cd$alpha_matrix_standard_error[2,4])
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==1),"Te_inter"]<-c(cxr_C.R1_cd$alpha_matrix[3,1]-cxr_C.R1_cd$alpha_matrix_standard_error[3,1], cxr_C.R1_cd$alpha_matrix[3,2]-cxr_C.R1_cd$alpha_matrix_standard_error[3,2],cxr_C.R1_cd$alpha_matrix[4,1]-cxr_C.R1_cd$alpha_matrix_standard_error[4,1], cxr_C.R1_cd$alpha_matrix[4,2]-cxr_C.R1_cd$alpha_matrix_standard_error[4,2])

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==2),"Tu_inter"]<-c(cxr_C.R2_cd_sr1$alpha_inter[1]-cxr_C.R2_cd_sr1$alpha_inter_standard_error[1], cxr_C.R2_cd_sr1$alpha_inter[2]-cxr_C.R2_cd_sr1$alpha_inter_standard_error[2])
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==2),"Te_inter"]<-c(cxr_C.R2_cd_sr4$alpha_inter[1]-cxr_C.R2_cd_sr4$alpha_inter_standard_error[1], cxr_C.R2_cd_sr5$alpha_inter[1]-cxr_C.R2_cd_sr5$alpha_inter_standard_error[1])

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==3),"Tu_inter"]<-c(cxr_C.R3_cd$alpha_matrix[1,3]-cxr_C.R3_cd$alpha_matrix_standard_error[1,3], cxr_C.R3_cd$alpha_matrix[2,3]-cxr_C.R3_cd$alpha_matrix_standard_error[2,3],cxr_C.R3_cd$alpha_matrix[1,4]-cxr_C.R3_cd$alpha_matrix_standard_error[1,4], cxr_C.R3_cd$alpha_matrix[2,4]-cxr_C.R3_cd$alpha_matrix_standard_error[2,4])
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==3),"Te_inter"]<-c(cxr_C.R3_cd$alpha_matrix[3,1]-cxr_C.R3_cd$alpha_matrix_standard_error[3,1], cxr_C.R3_cd$alpha_matrix[3,2]-cxr_C.R3_cd$alpha_matrix_standard_error[3,2],cxr_C.R3_cd$alpha_matrix[4,1]-cxr_C.R3_cd$alpha_matrix_standard_error[4,1], cxr_C.R3_cd$alpha_matrix[4,2]-cxr_C.R3_cd$alpha_matrix_standard_error[4,2])

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==4),"Tu_inter"]<-c(cxr_C.R4_cd$alpha_matrix[1,3]-cxr_C.R4_cd$alpha_matrix_standard_error[1,3], cxr_C.R4_cd$alpha_matrix[2,3]-cxr_C.R4_cd$alpha_matrix_standard_error[2,3],cxr_C.R4_cd$alpha_matrix[1,4]-cxr_C.R4_cd$alpha_matrix_standard_error[1,4], cxr_C.R4_cd$alpha_matrix[2,4]-cxr_C.R4_cd$alpha_matrix_standard_error[2,4])
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==4),"Te_inter"]<-c(cxr_C.R4_cd$alpha_matrix[3,1]-cxr_C.R4_cd$alpha_matrix_standard_error[3,1], cxr_C.R4_cd$alpha_matrix[3,2]-cxr_C.R4_cd$alpha_matrix_standard_error[3,2],cxr_C.R4_cd$alpha_matrix[4,1]-cxr_C.R4_cd$alpha_matrix_standard_error[4,1], cxr_C.R4_cd$alpha_matrix[4,2]-cxr_C.R4_cd$alpha_matrix_standard_error[4,2])

cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==5),"Tu_inter"]<-c(cxr_C.R5_cd$alpha_matrix[1,3]-cxr_C.R5_cd$alpha_matrix_standard_error[1,3], cxr_C.R5_cd$alpha_matrix[2,3]-cxr_C.R5_cd$alpha_matrix_standard_error[2,3],cxr_C.R5_cd$alpha_matrix[1,4]-cxr_C.R5_cd$alpha_matrix_standard_error[1,4], cxr_C.R5_cd$alpha_matrix[2,4]-cxr_C.R5_cd$alpha_matrix_standard_error[2,4])
cxr_param_CC_lower[which(cxr_param_CC_lower$Replicate==5),"Te_inter"]<-c(cxr_C.R5_cd$alpha_matrix[3,1]-cxr_C.R5_cd$alpha_matrix_standard_error[3,1], cxr_C.R5_cd$alpha_matrix[3,2]-cxr_C.R5_cd$alpha_matrix_standard_error[3,2],cxr_C.R5_cd$alpha_matrix[4,1]-cxr_C.R5_cd$alpha_matrix_standard_error[4,1], cxr_C.R5_cd$alpha_matrix[4,2]-cxr_C.R5_cd$alpha_matrix_standard_error[4,2])

### upper

cxr_param_CC_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_CC_upper$Tu_lambda<-0
cxr_param_CC_upper$Te_lambda<-0
cxr_param_CC_upper$Tu_intra<-0
cxr_param_CC_upper$Te_intra<-0
cxr_param_CC_upper$Tu_inter<-0
cxr_param_CC_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_CC_upper<-cxr_param_CC_upper[-which(cxr_param_CC_upper$Replicate==2 & cxr_param_CC_upper$Tu_Regime=="SR2"),]


cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==1),"Tu_lambda"]<-c(cxr_C.R1_cd$fixed_terms[[1]]$lambda+sd_1C[[1]]$lambda,cxr_C.R1_cd$fixed_terms[[2]]$lambda+sd_1C[[2]]$lambda)
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==1),"Te_lambda"]<-c(cxr_C.R1_cd$fixed_terms[[3]]$lambda+sd_1C[[3]]$lambda,cxr_C.R1_cd$fixed_terms[[3]]$lambda+sd_1C[[3]]$lambda, cxr_C.R1_cd$fixed_terms[[4]]$lambda+sd_1C[[4]]$lambda,cxr_C.R1_cd$fixed_terms[[4]]$lambda+sd_1C[[4]]$lambda)

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==2),"Tu_lambda"]<-c(cxr_C.R2_cd_sr1$fixed_terms$lambda+sd_2C[[1]]$lambda,cxr_C.R2_cd_sr1$fixed_terms$lambda+sd_2C[[1]]$lambda)
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==2),"Te_lambda"]<-c(cxr_C.R2_cd_sr4$fixed_terms$lambda+sd_2C[[2]]$lambda,cxr_C.R2_cd_sr5$fixed_terms$lambda+sd_2C[[3]]$lambda)

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==3),"Tu_lambda"]<-c(cxr_C.R3_cd$fixed_terms[[1]]$lambda+sd_3C[[1]]$lambda,cxr_C.R3_cd$fixed_terms[[2]]$lambda+sd_3C[[2]]$lambda)
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==3),"Te_lambda"]<-c(cxr_C.R3_cd$fixed_terms[[3]]$lambda+sd_3C[[3]]$lambda,cxr_C.R3_cd$fixed_terms[[3]]$lambda+sd_3C[[3]]$lambda, cxr_C.R3_cd$fixed_terms[[4]]$lambda+sd_3C[[4]]$lambda,cxr_C.R3_cd$fixed_terms[[4]]$lambda+sd_3C[[4]]$lambda)

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==4),"Tu_lambda"]<-c(cxr_C.R4_cd$fixed_terms[[1]]$lambda+sd_4C[[1]]$lambda,cxr_C.R4_cd$fixed_terms[[2]]$lambda+sd_4C[[2]]$lambda)
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==4),"Te_lambda"]<-c(cxr_C.R4_cd$fixed_terms[[3]]$lambda+sd_4C[[3]]$lambda,cxr_C.R4_cd$fixed_terms[[3]]$lambda+sd_4C[[3]]$lambda, cxr_C.R4_cd$fixed_terms[[4]]$lambda+sd_4C[[4]]$lambda,cxr_C.R4_cd$fixed_terms[[4]]$lambda+sd_4C[[4]]$lambda)

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==5),"Tu_lambda"]<-c(cxr_C.R5_cd$fixed_terms[[1]]$lambda+sd_5C[[1]]$lambda,cxr_C.R5_cd$fixed_terms[[2]]$lambda+sd_5C[[2]]$lambda)
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==5),"Te_lambda"]<-c(cxr_C.R5_cd$fixed_terms[[3]]$lambda+sd_5C[[3]]$lambda,cxr_C.R5_cd$fixed_terms[[3]]$lambda+sd_5C[[3]]$lambda, cxr_C.R5_cd$fixed_terms[[4]]$lambda+sd_5C[[4]]$lambda,cxr_C.R5_cd$fixed_terms[[4]]$lambda+sd_5C[[4]]$lambda)


cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==1),"Tu_intra"]<-rep(c(cxr_C.R1_cd_intra$alpha_matrix[1,1]+cxr_C.R1_cd_intra$alpha_matrix_standard_error[1,1], cxr_C.R1_cd_intra$alpha_matrix[2,1]+cxr_C.R1_cd_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==1),"Te_intra"]<-rep(c(cxr_C.R1_cd_intra$alpha_matrix[3,1]+cxr_C.R1_cd_intra$alpha_matrix_standard_error[3,1], cxr_C.R1_cd_intra$alpha_matrix[4,1]+cxr_C.R1_cd_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==2),"Tu_intra"]<-c(cxr_C.R2_cd_intra$alpha_matrix[1,1]+cxr_C.R2_cd_intra$alpha_matrix_standard_error[1,1],cxr_C.R2_cd_intra$alpha_matrix[1,1]+cxr_C.R2_cd_intra$alpha_matrix_standard_error[1,1])
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==2),"Te_intra"]<-c(cxr_C.R2_cd_intra$alpha_matrix[2,1]+cxr_C.R2_cd_intra$alpha_matrix_standard_error[2,1],cxr_C.R2_cd_intra$alpha_matrix[3,1]+cxr_C.R2_cd_intra$alpha_matrix_standard_error[3,1])

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==3),"Tu_intra"]<-rep(c(cxr_C.R3_cd_intra$alpha_matrix[1,1]+cxr_C.R3_cd_intra$alpha_matrix_standard_error[1,1], cxr_C.R3_cd_intra$alpha_matrix[2,1]+cxr_C.R3_cd_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==3),"Te_intra"]<-rep(c(cxr_C.R3_cd_intra$alpha_matrix[3,1]+cxr_C.R3_cd_intra$alpha_matrix_standard_error[3,1], cxr_C.R3_cd_intra$alpha_matrix[4,1]+cxr_C.R3_cd_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==4),"Tu_intra"]<-rep(c(cxr_C.R4_cd_intra$alpha_matrix[1,1]+cxr_C.R4_cd_intra$alpha_matrix_standard_error[1,1], cxr_C.R4_cd_intra$alpha_matrix[2,1]+cxr_C.R4_cd_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==4),"Te_intra"]<-rep(c(cxr_C.R4_cd_intra$alpha_matrix[3,1]+cxr_C.R4_cd_intra$alpha_matrix_standard_error[3,1], cxr_C.R4_cd_intra$alpha_matrix[4,1]+cxr_C.R4_cd_intra$alpha_matrix_standard_error[4,1]), each=2)

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==5),"Tu_intra"]<-rep(c(cxr_C.R5_cd_intra$alpha_matrix[1,1]+cxr_C.R5_cd_intra$alpha_matrix_standard_error[1,1], cxr_C.R5_cd_intra$alpha_matrix[2,1]+cxr_C.R5_cd_intra$alpha_matrix_standard_error[2,1]), 2)
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==5),"Te_intra"]<-rep(c(cxr_C.R5_cd_intra$alpha_matrix[3,1]+cxr_C.R5_cd_intra$alpha_matrix_standard_error[3,1], cxr_C.R5_cd_intra$alpha_matrix[4,1]+cxr_C.R5_cd_intra$alpha_matrix_standard_error[4,1]), each=2)


cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==1),"Tu_inter"]<-c(cxr_C.R1_cd$alpha_matrix[1,3]+cxr_C.R1_cd$alpha_matrix_standard_error[1,3], cxr_C.R1_cd$alpha_matrix[2,3]+cxr_C.R1_cd$alpha_matrix_standard_error[2,3],cxr_C.R1_cd$alpha_matrix[1,4]+cxr_C.R1_cd$alpha_matrix_standard_error[1,4], cxr_C.R1_cd$alpha_matrix[2,4]+cxr_C.R1_cd$alpha_matrix_standard_error[2,4])
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==1),"Te_inter"]<-c(cxr_C.R1_cd$alpha_matrix[3,1]+cxr_C.R1_cd$alpha_matrix_standard_error[3,1], cxr_C.R1_cd$alpha_matrix[3,2]+cxr_C.R1_cd$alpha_matrix_standard_error[3,2],cxr_C.R1_cd$alpha_matrix[4,1]+cxr_C.R1_cd$alpha_matrix_standard_error[4,1], cxr_C.R1_cd$alpha_matrix[4,2]+cxr_C.R1_cd$alpha_matrix_standard_error[4,2])

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==2),"Tu_inter"]<-c(cxr_C.R2_cd_sr1$alpha_inter[1]+cxr_C.R2_cd_sr1$alpha_inter_standard_error[1], cxr_C.R2_cd_sr1$alpha_inter[2]+cxr_C.R2_cd_sr1$alpha_inter_standard_error[2])
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==2),"Te_inter"]<-c(cxr_C.R2_cd_sr4$alpha_inter[1]+cxr_C.R2_cd_sr4$alpha_inter_standard_error[1], cxr_C.R2_cd_sr5$alpha_inter[1]+cxr_C.R2_cd_sr5$alpha_inter_standard_error[1])

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==3),"Tu_inter"]<-c(cxr_C.R3_cd$alpha_matrix[1,3]+cxr_C.R3_cd$alpha_matrix_standard_error[1,3], cxr_C.R3_cd$alpha_matrix[2,3]+cxr_C.R3_cd$alpha_matrix_standard_error[2,3],cxr_C.R3_cd$alpha_matrix[1,4]+cxr_C.R3_cd$alpha_matrix_standard_error[1,4], cxr_C.R3_cd$alpha_matrix[2,4]+cxr_C.R3_cd$alpha_matrix_standard_error[2,4])
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==3),"Te_inter"]<-c(cxr_C.R3_cd$alpha_matrix[3,1]+cxr_C.R3_cd$alpha_matrix_standard_error[3,1], cxr_C.R3_cd$alpha_matrix[3,2]+cxr_C.R3_cd$alpha_matrix_standard_error[3,2],cxr_C.R3_cd$alpha_matrix[4,1]+cxr_C.R3_cd$alpha_matrix_standard_error[4,1], cxr_C.R3_cd$alpha_matrix[4,2]+cxr_C.R3_cd$alpha_matrix_standard_error[4,2])

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==4),"Tu_inter"]<-c(cxr_C.R4_cd$alpha_matrix[1,3]+cxr_C.R4_cd$alpha_matrix_standard_error[1,3], cxr_C.R4_cd$alpha_matrix[2,3]+cxr_C.R4_cd$alpha_matrix_standard_error[2,3],cxr_C.R4_cd$alpha_matrix[1,4]+cxr_C.R4_cd$alpha_matrix_standard_error[1,4], cxr_C.R4_cd$alpha_matrix[2,4]+cxr_C.R4_cd$alpha_matrix_standard_error[2,4])
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==4),"Te_inter"]<-c(cxr_C.R4_cd$alpha_matrix[3,1]+cxr_C.R4_cd$alpha_matrix_standard_error[3,1], cxr_C.R4_cd$alpha_matrix[3,2]+cxr_C.R4_cd$alpha_matrix_standard_error[3,2],cxr_C.R4_cd$alpha_matrix[4,1]+cxr_C.R4_cd$alpha_matrix_standard_error[4,1], cxr_C.R4_cd$alpha_matrix[4,2]+cxr_C.R4_cd$alpha_matrix_standard_error[4,2])

cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==5),"Tu_inter"]<-c(cxr_C.R5_cd$alpha_matrix[1,3]+cxr_C.R5_cd$alpha_matrix_standard_error[1,3], cxr_C.R5_cd$alpha_matrix[2,3]+cxr_C.R5_cd$alpha_matrix_standard_error[2,3],cxr_C.R5_cd$alpha_matrix[1,4]+cxr_C.R5_cd$alpha_matrix_standard_error[1,4], cxr_C.R5_cd$alpha_matrix[2,4]+cxr_C.R5_cd$alpha_matrix_standard_error[2,4])
cxr_param_CC_upper[which(cxr_param_CC_upper$Replicate==5),"Te_inter"]<-c(cxr_C.R5_cd$alpha_matrix[3,1]+cxr_C.R5_cd$alpha_matrix_standard_error[3,1], cxr_C.R5_cd$alpha_matrix[3,2]+cxr_C.R5_cd$alpha_matrix_standard_error[3,2],cxr_C.R5_cd$alpha_matrix[4,1]+cxr_C.R5_cd$alpha_matrix_standard_error[4,1], cxr_C.R5_cd$alpha_matrix[4,2]+cxr_C.R5_cd$alpha_matrix_standard_error[4,2])




##### joining data frame



param_all_C<-as.data.frame(rbind(cxr_param_C, cxr_param_CC))

param_all_C_lower<-as.data.frame(rbind(cxr_param_C_lower, cxr_param_CC_lower))
param_all_C_upper<-as.data.frame(rbind(cxr_param_C_upper, cxr_param_CC_upper))


param_all_C_lower
param_all_C_upper

#write.csv(param_all_C, "./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed.csv")
#write.csv(param_all_C_upper, "./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed_upper.csv")
#write.csv(param_all_C_lower, "./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed_lower.csv")




##### importing data frame
param_all_C<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed.csv")
param_all_C_upper<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed_upper.csv")
param_all_C_lower<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed_lower.csv")

param_all_C<-param_all_C[,-1]
param_all_C_upper<-param_all_C_upper[,-1]
param_all_C_lower<-param_all_C_lower[,-1]





param_all_C_long<-gather(param_all_C, parameter, value,Tu_lambda:Te_inter )

param_all_C_long$category<-mapvalues(param_all_C_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

param_all_C_lower_long<-gather(param_all_C_lower, parameter, value,Tu_lambda:Te_inter )

param_all_C_lower_long$category<-mapvalues(param_all_C_lower_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

param_all_C_upper_long<-gather(param_all_C_upper, parameter, value,Tu_lambda:Te_inter )

param_all_C_upper_long$category<-mapvalues(param_all_C_upper_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

colnames(param_all_C_lower_long)[6]<-"lower"
colnames(param_all_C_upper_long)[6]<-"upper"

str(param_all_C_long)

param_all_C_long<-cbind(param_all_C_long[,1:7],param_all_C_lower_long$lower, param_all_C_upper_long$upper)

colnames(param_all_C_long)[8:9]<-c("lower","upper")

#### Predicting densities
density_aux<-seq(0, 10, by=(10/100))

pred_df_cxr_C<-as.data.frame(expand_grid(Density=density_aux, Tu_Regime=c("SR1","SR2"), Te_Regime=c("SR4","SR5"), Replicate=c(1:5), Environment=c("N", "Cd")))

pred_df_cxr_C$Tu_mean_intra<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Tu_mean_inter<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_C$Density[x])
  
  pred
})


pred_df_cxr_C$Tu_intra_L<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Tu_inter_L<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Tu_intra_U<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Tu_inter_U<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_inter[1]
  lambda<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Te_mean_intra<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Te_mean_inter<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_C, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Te_intra_L<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Te_inter_L<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_C_lower, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Te_intra_U<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_cxr_C$Density[x])
  
  pred
})

pred_df_cxr_C$Te_inter_U<-sapply(c(1:length(pred_df_cxr_C[,1])), function(x){
  alpha_i<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Tu_Regime==pred_df_cxr_C$Tu_Regime[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_inter[1]
  lambda<-subset(param_all_C_upper, Environment==pred_df_cxr_C$Environment[x] & Te_Regime==pred_df_cxr_C$Te_Regime[x] & Replicate==pred_df_cxr_C$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_cxr_C$Density[x])
  
  pred
})

# Removing Tu evolved replicate 2 because there is no data
pred_df_cxr_C<-pred_df_cxr_C[-which(pred_df_cxr_C$Tu_Regime=="SR2" & pred_df_cxr_C$Replicate==2),]



# Transforming everything bellow 0 into 0 for the lower interval

pred_df_cxr_C$Te_inter_L[which(pred_df_cxr_C$Te_inter_L<0)]<-0
pred_df_cxr_C$Te_intra_L[which(pred_df_cxr_C$Te_intra_L<0)]<-0
pred_df_cxr_C$Tu_inter_L[which(pred_df_cxr_C$Tu_inter_L<0)]<-0
pred_df_cxr_C$Tu_intra_L[which(pred_df_cxr_C$Tu_intra_L<0)]<-0

##### Predicted vs observed
rk_func<- function(lambda, alpha_ii, alpha_ij, dens_i, dens_j, ...){
  gr<-lambda*exp(-alpha_ii*dens_i - alpha_ij*dens_j)
  
  return(gr)
}

red_ca_C<-ca[,c("Env", "Rep", "FocalSR", "CompSR", "Dens", "Type", "TeFemales", "TuFemales", "GrowthRateOA")]

red_ca_C

red_ca_C$Dens_Focal<-sapply(c(1:length(red_ca_C[,1])), function(x){
  if(red_ca_C$Type[x]=="INTRA"){
    a<-red_ca_C$Dens[x]-1
  }else if(red_ca_C$Type[x]=="INTER"){
    a<-1
  }
  
  a
})

red_ca_C$Dens_Comp<-sapply(c(1:length(red_ca_C[,1])), function(x){
  if(red_ca_C$Type[x]=="INTRA"){
    a<-0
  }else if(red_ca_C$Type[x]=="INTER"){
    a<-red_ca_C$Dens[x]-1
  }
  
  a
})

red_ca_C$Focal<-mapvalues(red_ca_C$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4", "SR5"))
red_ca_C$Comp<-mapvalues(red_ca_C$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4", "SR5"))

red_ca_C$pred<-sapply(c(1:length(red_ca_C[,1])), function(x){
  
  if(red_ca_C$Focal[x]=="SR1" | red_ca_C$Focal[x]=="SR2"){
    aux_data<-subset(param_all_C, Environment==red_ca_C$Env[x] & Replicate== red_ca_C$Rep[x] & as.character(Tu_Regime)==red_ca_C$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Tu_lambda[1], alpha_ii =aux_data$Tu_intra[1], alpha_ij = aux_data$Tu_inter[1], dens_i = red_ca_C$Dens_Focal[x], dens_j =  red_ca_C$Dens_Comp[x])
    
  }else if(red_ca_C$Focal[x]=="SR4" | red_ca_C$Focal[x]=="SR5"){
    aux_data<-subset(param_all_C, Environment==red_ca_C$Env[x] & Replicate== red_ca_C$Rep[x] & as.character(Te_Regime)==red_ca_C$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Te_lambda[1], alpha_ii =aux_data$Te_intra[1], alpha_ij = aux_data$Te_inter[1], dens_i = red_ca_C$Dens_Focal[x], dens_j =  red_ca_C$Dens_Comp[x])
  }
  
  aux_pred
})

red_ca_C$pred_L<-sapply(c(1:length(red_ca_C[,1])), function(x){
  
  if(red_ca_C$Focal[x]=="SR1" | red_ca_C$Focal[x]=="SR2"){
    aux_data<-subset(param_all_C_lower, Environment==red_ca_C$Env[x] & Replicate== red_ca_C$Rep[x] & as.character(Tu_Regime)==red_ca_C$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Tu_lambda[1], alpha_ii =aux_data$Tu_intra[1], alpha_ij = aux_data$Tu_inter[1], dens_i = red_ca_C$Dens_Focal[x], dens_j =  red_ca_C$Dens_Comp[x])
    
  }else if(red_ca_C$Focal[x]=="SR4" | red_ca_C$Focal[x]=="SR5"){
    aux_data<-subset(param_all_C_lower, Environment==red_ca_C$Env[x] & Replicate== red_ca_C$Rep[x] & as.character(Te_Regime)==red_ca_C$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Te_lambda[1], alpha_ii =aux_data$Te_intra[1], alpha_ij = aux_data$Te_inter[1], dens_i = red_ca_C$Dens_Focal[x], dens_j =  red_ca_C$Dens_Comp[x])
  }
  
  aux_pred
})

red_ca_C$pred_U<-sapply(c(1:length(red_ca_C[,1])), function(x){
  
  if(red_ca_C$Focal[x]=="SR1" | red_ca_C$Focal[x]=="SR2"){
    aux_data<-subset(param_all_C_upper, Environment==red_ca_C$Env[x] & Replicate== red_ca_C$Rep[x] & as.character(Tu_Regime)==red_ca_C$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Tu_lambda[1], alpha_ii =aux_data$Tu_intra[1], alpha_ij = aux_data$Tu_inter[1], dens_i = red_ca_C$Dens_Focal[x], dens_j =  red_ca_C$Dens_Comp[x])
    
  }else if(red_ca_C$Focal[x]=="SR4" | red_ca_C$Focal[x]=="SR5"){
    aux_data<-subset(param_all_C_upper, Environment==red_ca_C$Env[x] & Replicate== red_ca_C$Rep[x] & as.character(Te_Regime)==red_ca_C$Focal[x])
    
    aux_pred<-rk_func(lambda=aux_data$Te_lambda[1], alpha_ii =aux_data$Te_intra[1], alpha_ij = aux_data$Te_inter[1], dens_i = red_ca_C$Dens_Focal[x], dens_j =  red_ca_C$Dens_Comp[x])
  }
  
  aux_pred
})

red_ca_C$Replicate<-red_ca_C$Rep
str(red_ca_C)

### D - optim normal

##### Estimating parameters

# creating folder to put the analyses inside, this should be the same as the file path in the function
dir.create("./Analyses/MethodComparison/Optim_normal", showWarnings = FALSE)

source("./Scrips/Function_riker_27May.R")
# This matrix has all the comparisons that need to be done between regimes
comparison_mat<-matrix(nrow=4, ncol=3)
comparison_mat[1,]<-c(1,4,5)
comparison_mat[2,]<-c(2,4,5)
comparison_mat[3,]<-c(4,1,2)
comparison_mat[4,]<-c(5,1,2)

#lam2 is the data from density one corresponding to the focals populations
# data2 is the data (format) Regime (name of focal pop), background (name of competitor, the same if its intraspecific competition), focal (number of focal individuals in g0), comp (number of competitors in g0), growth rate
# Attention that for intraspecific you need to add 0 in the comp and all individuals in the focal



rep2<-mod_df(subset(ca,Rep==1 & Env=="N"))  
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/",data2=rep2, reps2=1, env="N", comparisons = comparison_mat)

rep2<-mod_df(subset(ca,Rep==1 & Env=="Cd"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=1, env="Cd", comparisons = comparison_mat)

rep2<-mod_df(subset(ca,Rep==3 & Env=="N"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=3, env="N", comparisons = comparison_mat)

rep2<-mod_df(subset(ca,Rep==3 & Env=="Cd"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=3, env="Cd", comparisons = comparison_mat)

rep2<-mod_df(subset(ca,Rep==4 & Env=="N"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=4, env="N", comparisons = comparison_mat)

rep2<-mod_df(subset(ca,Rep==4 & Env=="Cd"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=4, env="Cd", comparisons = comparison_mat)

rep2<-mod_df(subset(ca,Rep==5 & Env=="N"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=5, env="N", comparisons = comparison_mat)

rep2<-mod_df(subset(ca,Rep==5 & Env=="Cd"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=5, env="Cd", comparisons = comparison_mat)

# For two we have to change the comparison matrix
comparison_mat2<-matrix(nrow=3, ncol=3)
comparison_mat2[1,]<-c(1,4,5)
comparison_mat2[2,]<-c(4,1,NA)
comparison_mat2[3,]<-c(5,1,NA)

rep2<-mod_df(subset(ca,Rep==2 & Env=="N"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=2, env="N", comparisons = comparison_mat2)

rep2<-mod_df(subset(ca,Rep==2 & Env=="Cd"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=2, env="Cd", comparisons = comparison_mat2)



##### Importing files of alpha and lambda
alpha_file<-list.files("./Analyses/MethodComparison/Optim_normal/", pattern="alpha_estimates") #the alphas are always tu, te (row), tu, te (col)

alphaUpper_file<-list.files("./Analyses/MethodComparison/Optim_normal/", pattern="alpha_upper")

alphaLower_file<-list.files("./Analyses/MethodComparison/Optim_normal/", pattern="alpha_lower")

lambda_file<-list.files("./Analyses/MethodComparison/Optim_normal/", pattern="lambda_estimates")


alpha_list<- lapply(alpha_file, function(x) read.csv(paste("./Analyses/MethodComparison/Optim_normal/",x, sep=""), header = TRUE))
alphaUpper_list<- lapply(alphaUpper_file, function(x) read.csv(paste("./Analyses/MethodComparison/Optim_normal/",x, sep=""), header = TRUE))
alphaLower_list<- lapply(alphaLower_file, function(x) read.csv(paste("./Analyses/MethodComparison/Optim_normal/",x, sep=""), header = TRUE))
lambda_list<- lapply(lambda_file, function(x) read.csv(paste("./Analyses/MethodComparison/Optim_normal/",x, sep=""), header = TRUE))

# passing from list to data frame
# First we need to do the first iteration (to create everything)
lambda_intra_fixed<-data.frame(Regime1=rep(c(1,1,2,2),10), Regime2=rep(c(4,5,4,5), 10), Replicate=c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8)), Env=rep(c(rep("N",4), rep("Cd",4)), 5))

lambda_intra_fixed<-lambda_intra_fixed[-which(lambda_intra_fixed$Regime1==2 & lambda_intra_fixed$Replicate==2),] # to remove SR2 from replicate 2 because it does not exist

alpha_list[[1]]
lambda_list[[1]]

# passing alphas to dataframe
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

aux_alpha<-as.data.frame(alpha_list[[1]])

aux_alpha2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu=c(aux_alpha[1,2], aux_alpha[1,2], aux_alpha[2,2],aux_alpha[2,2]), intraTe=c(aux_alpha[3,2], aux_alpha[4,2], aux_alpha[3,2],aux_alpha[4,2]), interTu=c(aux_alpha[1,3], aux_alpha[1,4], aux_alpha[2,3], aux_alpha[2,4]), interTe=c(aux_alpha[3,3], aux_alpha[4,3], aux_alpha[3,4], aux_alpha[4,4]))

for(x in 2:length(lambda_list)){
  repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[x], split="_")[1])[6], split="[.]"))[1],split=""))[1]
  env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[x], split="_")[1])[6], split="[.]"))[1],split=""))[2]
  
  if(x==3 | x==4){# because there is no SR2 here
    regimeTu<-c("1","1")
    regimeTe<-c("4","5")
    Env<-rep(env, 2)
    Rep<-rep(repli,2)
    
    aux_alpha<-as.data.frame(alpha_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu=c(aux_alpha[1,2], aux_alpha[1,2]), intraTe=c(aux_alpha[2,2], aux_alpha[3,2]), interTu=c(aux_alpha[1,3], aux_alpha[1,4]), interTe=c(aux_alpha[2,3], aux_alpha[3,3]))
    
    
  }else{
    regimeTu<-c("1","1", "2","2")
    regimeTe<-c("4","5", "4","5")
    Env<-rep(env, 4)
    Rep<-rep(repli,4)
    
    aux_alpha<-as.data.frame(alpha_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu=c(aux_alpha[1,2], aux_alpha[1,2], aux_alpha[2,2],aux_alpha[2,2]), intraTe=c(aux_alpha[3,2], aux_alpha[4,2], aux_alpha[3,2],aux_alpha[4,2]), interTu=c(aux_alpha[1,3], aux_alpha[1,4], aux_alpha[2,3], aux_alpha[2,4]), interTe=c(aux_alpha[3,3], aux_alpha[4,3], aux_alpha[3,4], aux_alpha[4,4]))
  }
  
  aux_alpha2<-rbind(aux_alpha2, aux2)
}

### Alpha Lower

# passing alphas to dataframe
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

aux_alphaLower<-as.data.frame(alphaLower_list[[1]])

aux_alphaLower2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_L=c(aux_alphaLower[1,2], aux_alphaLower[1,2], aux_alphaLower[2,2],aux_alphaLower[2,2]), intraTe_L=c(aux_alphaLower[3,2], aux_alphaLower[4,2], aux_alphaLower[3,2],aux_alphaLower[4,2]), interTu_L=c(aux_alphaLower[1,3], aux_alphaLower[1,4], aux_alphaLower[2,3], aux_alphaLower[2,4]), interTe_L=c(aux_alphaLower[3,3], aux_alphaLower[4,3], aux_alphaLower[3,4], aux_alphaLower[4,4]))

for(x in 2:length(lambda_list)){
  repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[x], split="_")[1])[4], split="[.]"))[1],split=""))[1]
  env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[x], split="_")[1])[4], split="[.]"))[1],split=""))[2]
  
  if(x==3 | x==4){# because there is no SR2 here
    regimeTu<-c("1","1")
    regimeTe<-c("4","5")
    Env<-rep(env, 2)
    Rep<-rep(repli,2)
    
    aux_alphaLower<-as.data.frame(alphaLower_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_L=c(aux_alphaLower[1,2], aux_alphaLower[1,2]), intraTe_L=c(aux_alphaLower[2,2], aux_alphaLower[3,2]), interTu_L=c(aux_alphaLower[1,3], aux_alphaLower[1,4]), interTe_L=c(aux_alphaLower[2,3], aux_alphaLower[3,3]))
    
    
  }else{
    regimeTu<-c("1","1", "2","2")
    regimeTe<-c("4","5", "4","5")
    Env<-rep(env, 4)
    Rep<-rep(repli,4)
    
    aux_alphaLower<-as.data.frame(alphaLower_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_L=c(aux_alphaLower[1,2], aux_alphaLower[1,2], aux_alphaLower[2,2],aux_alphaLower[2,2]), intraTe_L=c(aux_alphaLower[3,2], aux_alphaLower[4,2], aux_alphaLower[3,2],aux_alphaLower[4,2]), interTu_L=c(aux_alphaLower[1,3], aux_alphaLower[1,4], aux_alphaLower[2,3], aux_alphaLower[2,4]), interTe_L=c(aux_alphaLower[3,3], aux_alphaLower[4,3], aux_alphaLower[3,4], aux_alphaLower[4,4]))
  }
  
  aux_alphaLower2<-rbind(aux_alphaLower2, aux2)
}

### Alpha Upper

# passing alphas to dataframe
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

aux_alphaUpper<-as.data.frame(alphaUpper_list[[1]])

aux_alphaUpper2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_U=c(aux_alphaUpper[1,2], aux_alphaUpper[1,2], aux_alphaUpper[2,2],aux_alphaUpper[2,2]), intraTe_U=c(aux_alphaUpper[3,2], aux_alphaUpper[4,2], aux_alphaUpper[3,2],aux_alphaUpper[4,2]), interTu_U=c(aux_alphaUpper[1,3], aux_alphaUpper[1,4], aux_alphaUpper[2,3], aux_alphaUpper[2,4]), interTe_U=c(aux_alphaUpper[3,3], aux_alphaUpper[4,3], aux_alphaUpper[3,4], aux_alphaUpper[4,4]))

for(x in 2:length(lambda_list)){
  repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[x], split="_")[1])[4], split="[.]"))[1],split=""))[1]
  env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[x], split="_")[1])[4], split="[.]"))[1],split=""))[2]
  
  if(x==3 | x==4){# because there is no SR2 here
    regimeTu<-c("1","1")
    regimeTe<-c("4","5")
    Env<-rep(env, 2)
    Rep<-rep(repli,2)
    
    aux_alphaUpper<-as.data.frame(alphaUpper_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_U=c(aux_alphaUpper[1,2], aux_alphaUpper[1,2]), intraTe_U=c(aux_alphaUpper[2,2], aux_alphaUpper[3,2]), interTu_U=c(aux_alphaUpper[1,3], aux_alphaUpper[1,4]), interTe_U=c(aux_alphaUpper[2,3], aux_alphaUpper[3,3]))
    
    
  }else{
    regimeTu<-c("1","1", "2","2")
    regimeTe<-c("4","5", "4","5")
    Env<-rep(env, 4)
    Rep<-rep(repli,4)
    
    aux_alphaUpper<-as.data.frame(alphaUpper_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_U=c(aux_alphaUpper[1,2], aux_alphaUpper[1,2], aux_alphaUpper[2,2],aux_alphaUpper[2,2]), intraTe_U=c(aux_alphaUpper[3,2], aux_alphaUpper[4,2], aux_alphaUpper[3,2],aux_alphaUpper[4,2]), interTu_U=c(aux_alphaUpper[1,3], aux_alphaUpper[1,4], aux_alphaUpper[2,3], aux_alphaUpper[2,4]), interTe_U=c(aux_alphaUpper[3,3], aux_alphaUpper[4,3], aux_alphaUpper[3,4], aux_alphaUpper[4,4]))
  }
  
  aux_alphaUpper2<-rbind(aux_alphaUpper2, aux2)
}

# Passing lambda to data frame
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[2]

Focal<-c("1","1","2","2","4","4","5","5")
Comp<-c("4","5","4","5","1","2","1","2")
Env<-rep(env, 8)
Rep<-rep(repli,8)

aux_lambda<-cbind(as.data.frame(lambda_list[[1]])[,c(3,4,5)],Focal,Comp, Env, Rep)

for(x in 2:length(lambda_list)){
  repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[x], split="_")[1])[6], split="[.]"))[1],split=""))[1]
  env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[x], split="_")[1])[6], split="[.]"))[1],split=""))[2]
  
  if(x==3 | x==4){# because there is no SR2 here
    Focal<-c("1","1","4","5")
    Comp<-c("4","5","1","1")
    Env<-rep(env, 4)
    Rep<-rep(repli,4)
    
    aux<-cbind(as.data.frame(lambda_list[[x]])[,c(3,4,5)],Focal,Comp, Env, Rep)
    
  }else{
    Focal<-c("1","1","2","2","4","4","5","5")
    Comp<-c("4","5","4","5","1","2","1","2")
    Env<-rep(env, 8)
    Rep<-rep(repli,8)
    
    aux<-cbind(as.data.frame(lambda_list[[x]])[,c(3,4,5)],Focal,Comp, Env, Rep)
  }
  
  aux_lambda<-rbind(aux_lambda, aux)
}


#Matching all the data

alphas_mat_D<-as.data.frame(cbind(aux_alpha2, aux_alphaLower2, aux_alphaUpper2))

str(lambda_intra_fixed)

#### adding lambda

alphas_mat_D$lambdaTu<-sapply(c(1:length(alphas_mat_D[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_D$regimeTu[x] & Comp==alphas_mat_D$regimeTe[x]) & Rep==alphas_mat_D$Rep[x] & Env==alphas_mat_D$Env[x] )
  
  auxi[1,1]
})

alphas_mat_D$lambdaTe<-sapply(c(1:length(alphas_mat_D[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_D$regimeTe[x] & Comp==alphas_mat_D$regimeTu[x]) & Rep==alphas_mat_D$Rep[x] & Env==alphas_mat_D$Env[x] )
  
  auxi[1,1]
})


alphas_mat_D$lambdaTu_L<-sapply(c(1:length(alphas_mat_D[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_D$regimeTu[x] & Comp==alphas_mat_D$regimeTe[x]) & Rep==alphas_mat_D$Rep[x] & Env==alphas_mat_D$Env[x] )
  
  auxi[1,2]
})

alphas_mat_D$lambdaTe_L<-sapply(c(1:length(alphas_mat_D[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_D$regimeTe[x] & Comp==alphas_mat_D$regimeTu[x]) & Rep==alphas_mat_D$Rep[x] & Env==alphas_mat_D$Env[x] )
  
  auxi[1,2]
})

alphas_mat_D$lambdaTu_U<-sapply(c(1:length(alphas_mat_D[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_D$regimeTu[x] & Comp==alphas_mat_D$regimeTe[x]) & Rep==alphas_mat_D$Rep[x] & Env==alphas_mat_D$Env[x] )
  
  auxi[1,3]
})

alphas_mat_D$lambdaTe_U<-sapply(c(1:length(alphas_mat_D[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_D$regimeTe[x] & Comp==alphas_mat_D$regimeTu[x]) & Rep==alphas_mat_D$Rep[x] & Env==alphas_mat_D$Env[x] )
  
  auxi[1,3]
})

alphas_mat_D$Env2<-mapvalues(alphas_mat_D$Env, c("C","N"), c("Cd","N"))

# clean up the matrix, because it has a lot of repeated columns
alphas_mat_D<-alphas_mat_D[,c(1:8, 13:16,21:30)]

alphas_mat_D



colnames(alphas_mat_D)<-c("Tu_Regime", "Te_Regime", "Environment", "Replicate", "Tu_intra", "Te_intra", "Tu_inter", "Te_inter", "Tu_intra_L", "Te_intra_L", "Tu_inter_L", "Te_inter_L", "Tu_intra_U", "Te_intra_U", "Tu_inter_U", "Te_inter_U", "Tu_lambda", "Te_lambda","Tu_lambda_L", "Te_lambda_L","Tu_lambda_U", "Te_lambda_U")


alphas_mat_D_long<-gather(alphas_mat_D, parameter, value,Tu_intra:Te_lambda_U )

alphas_mat_D_long$category<-mapvalues(alphas_mat_D_long$parameter, c("Tu_intra", "Te_intra", "Tu_inter", "Te_inter", "Tu_intra_L", "Te_intra_L", "Tu_inter_L", "Te_inter_L", "Tu_intra_U", "Te_intra_U", "Tu_inter_U", "Te_inter_U", "Tu_lambda", "Te_lambda","Tu_lambda_L", "Te_lambda_L","Tu_lambda_U", "Te_lambda_U"), c("intra", "intra", "inter", "inter", "intra_L", "intra_L", "inter_L", "inter_L","intra_U", "intra_U", "inter_U", "inter_U","lambda","lambda","lambda_L","lambda_L","lambda_U","lambda_U"))

##### Predicting data
str(alphas_mat_D)

alphas_mat_D$Env2<-mapvalues(alphas_mat_D$Environment, c("C", "N"), c("Cd","N"))

# Since the lambda is from the log data
ca$pred_D<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR[x]) & Replicate==ca$Rep[x])$Tu_intra[1]
    alpha_ij<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Te_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Tu_inter[1]
    lambda<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_lambda[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_intra[1]
    alpha_ij<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Tu_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Te_inter[1]
    lambda<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_lambda[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]
    pred<-lambda*exp(-alpha_ij*densC)
  }
  
  pred
  
})
x<-1
ca$pred_D_L<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_intra_L[1]
    alpha_ij<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Te_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Tu_inter_L[1]
    lambda<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_lambda_L[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_intra_L[1]
    alpha_ij<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Tu_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Te_inter_L[1]
    lambda<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_lambda_L[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*(densF-1))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    pred<-lambda*exp(-alpha_ij*densC)
  }
  
  pred
  
})

ca$pred_D_U<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_intra_U[1]
    alpha_ij<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Te_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Tu_inter_U[1]
    lambda<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_lambda_U[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_intra_U[1]
    alpha_ij<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Tu_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Te_inter_U[1]
    lambda<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_lambda_U[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*(densF-1))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    pred<-lambda*exp(-alpha_ij*densC)
  }
  
  pred
  
})



##### Predicting each density
density_aux<-seq(0, 10, by=(10/100))

pred_df_D<-as.data.frame(expand_grid(Density=density_aux, Tu_Regime=c(1,2), Te_Regime=c(4,5), Replicate=c(1:5), Environment=c("N", "C")))

pred_df_D$Tu_mean_intra<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_inter[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Tu_mean_inter<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_inter[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_D$Density[x])
  
  pred
})


pred_df_D$Tu_intra_L<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_intra_L[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_inter_L[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_lambda_L[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Tu_inter_L<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_intra_L[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_inter_L[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_lambda_L[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Tu_intra_U<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_intra_U[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_inter_U[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_lambda_U[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Tu_inter_U<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_intra_U[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_inter_U[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Replicate==pred_df_D$Replicate[x])$Tu_lambda_U[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Te_mean_intra<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_inter[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Te_mean_inter<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_inter[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Te_intra_L<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_intra_L[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_inter_L[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_lambda_L[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Te_inter_L<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_intra_L[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_inter_L[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_lambda_L[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Te_intra_U<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_intra_U[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_inter_U[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_lambda_U[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_D$Density[x])
  
  pred
})

pred_df_D$Te_inter_U<-sapply(c(1:length(pred_df_D[,1])), function(x){
  alpha_i<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_intra_U[1]
  alpha_ij<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Tu_Regime==pred_df_D$Tu_Regime[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_inter_U[1]
  lambda<-subset(alphas_mat_D, Environment==pred_df_D$Environment[x] & Te_Regime==pred_df_D$Te_Regime[x] & Replicate==pred_df_D$Replicate[x])$Te_lambda_U[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_D$Density[x])
  
  pred
})

# Removing Tu evolved replicate 2 because there is no data
pred_df_D<-pred_df_D[-which(pred_df_D$Tu_Regime==2 & pred_df_D$Replicate==2),]



# Transforming everything bellow 0 into 0 for the lower interval

pred_df_D$Te_inter_L[which(pred_df_D$Te_inter_L<0)]<-0
pred_df_D$Te_intra_L[which(pred_df_D$Te_intra_L<0)]<-0
pred_df_D$Tu_inter_L[which(pred_df_D$Tu_inter_L<0)]<-0
pred_df_D$Tu_intra_L[which(pred_df_D$Tu_intra_L<0)]<-0

### E - optim lambda fixed

##### Estimate parameters

# creating folder to put the analyses inside, this should be the same as the file path in the function
dir.create("./Analyses/MethodComparison/optim_lambda_fixed", showWarnings = FALSE)

source("./Scrips/Function_riker_27May.R")
# This matrix has all the comparisons that need to be done between regimes
comparison_mat<-matrix(nrow=4, ncol=3)
comparison_mat[1,]<-c(1,4,5)
comparison_mat[2,]<-c(2,4,5)
comparison_mat[3,]<-c(4,1,2)
comparison_mat[4,]<-c(5,1,2)

#lam2 is the data from density one corresponding to the focals populations
# data2 is the data (format) Regime (name of focal pop), background (name of competitor, the same if its intraspecific competition), focal (number of focal individuals in g0), comp (number of competitors in g0), growth rate
# Attention that for intraspecific you need to add 0 in the comp and all individuals in the focal

rep2<-mod_df(subset(ca,Rep==1 & Env=="N"))  
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/",data2=rep2, reps2=1, env="N", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==1 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==1 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/", data2=rep2, reps2=1, env="Cd", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==1 & Env=="Cd"))

rep2<-mod_df(subset(ca,Rep==3 & Env=="N"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/",  data2=rep2, reps2=3, env="N", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==3 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==3 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/", data2=rep2, reps2=3, env="Cd", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==3 & Env=="Cd"))

rep2<-mod_df(subset(ca,Rep==4 & Env=="N"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/", data2=rep2, reps2=4, env="N", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==4 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==4 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/", data2=rep2, reps2=4, env="Cd", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==4 & Env=="Cd"))

rep2<-mod_df(subset(ca,Rep==5 & Env=="N"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/", data2=rep2, reps2=5, env="N", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==5 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==5 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/",data2=rep2, reps2=5, env="Cd", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==5 & Env=="Cd"))

# For two we have to change the comparison matrix
comparison_mat2<-matrix(nrow=3, ncol=3)
comparison_mat2[1,]<-c(1,4,5)
comparison_mat2[2,]<-c(4,1,NA)
comparison_mat2[3,]<-c(5,1,NA)

rep2<-mod_df(subset(ca,Rep==2 & Env=="N"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/",  data2=rep2, reps2=2, env="N", comparisons = comparison_mat2, lam2=subset(mean_dens1, Rep==2 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==2 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed/",  data2=rep2, reps2=2, env="Cd", comparisons = comparison_mat2, lam2=subset(mean_dens1, Rep==2 & Env=="Cd"))






##### Importing files of alpha and lambda
alpha_file<-list.files("./Analyses/MethodComparison/optim_lambda_fixed/", pattern="alpha_estimates") #the alphas are always tu, te (row), tu, te (col)

alphaUpper_file<-list.files("./Analyses/MethodComparison/optim_lambda_fixed/", pattern="alpha_upper")

alphaLower_file<-list.files("./Analyses/MethodComparison/optim_lambda_fixed/", pattern="alpha_lower")

lambda_file<-list.files("./Analyses/MethodComparison/optim_lambda_fixed/", pattern="lambda_estimates")


alpha_list<- lapply(alpha_file, function(x) read.csv(paste("./Analyses/MethodComparison/optim_lambda_fixed/",x, sep=""), header = TRUE))
alphaUpper_list<- lapply(alphaUpper_file, function(x) read.csv(paste("./Analyses/MethodComparison/optim_lambda_fixed/",x, sep=""), header = TRUE))
alphaLower_list<- lapply(alphaLower_file, function(x) read.csv(paste("./Analyses/MethodComparison/optim_lambda_fixed/",x, sep=""), header = TRUE))
lambda_list<- lapply(lambda_file, function(x) read.csv(paste("./Analyses/MethodComparison/optim_lambda_fixed/",x, sep=""), header = TRUE))

# passing from list to data frame
# First we need to do the first iteration (to create everything)
lambda_intra_fixed<-data.frame(Regime1=rep(c(1,1,2,2),10), Regime2=rep(c(4,5,4,5), 10), Replicate=c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8)), Env=rep(c(rep("N",4), rep("Cd",4)), 5))

lambda_intra_fixed<-lambda_intra_fixed[-which(lambda_intra_fixed$Regime1==2 & lambda_intra_fixed$Replicate==2),] # to remove SR2 from replicate 2 because it does not exist

alpha_list[[1]]
lambda_list[[1]]

# passing alphas to dataframe
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

aux_alpha<-as.data.frame(alpha_list[[1]])

aux_alpha2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu=c(aux_alpha[1,2], aux_alpha[1,2], aux_alpha[2,2],aux_alpha[2,2]), intraTe=c(aux_alpha[3,2], aux_alpha[4,2], aux_alpha[3,2],aux_alpha[4,2]), interTu=c(aux_alpha[1,3], aux_alpha[1,4], aux_alpha[2,3], aux_alpha[2,4]), interTe=c(aux_alpha[3,3], aux_alpha[4,3], aux_alpha[3,4], aux_alpha[4,4]))

for(x in 2:length(lambda_list)){
  repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[x], split="_")[1])[6], split="[.]"))[1],split=""))[1]
  env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[x], split="_")[1])[6], split="[.]"))[1],split=""))[2]
  
  if(x==3 | x==4){# because there is no SR2 here
    regimeTu<-c("1","1")
    regimeTe<-c("4","5")
    Env<-rep(env, 2)
    Rep<-rep(repli,2)
    
    aux_alpha<-as.data.frame(alpha_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu=c(aux_alpha[1,2], aux_alpha[1,2]), intraTe=c(aux_alpha[2,2], aux_alpha[3,2]), interTu=c(aux_alpha[1,3], aux_alpha[1,4]), interTe=c(aux_alpha[2,3], aux_alpha[3,3]))
    
    
  }else{
    regimeTu<-c("1","1", "2","2")
    regimeTe<-c("4","5", "4","5")
    Env<-rep(env, 4)
    Rep<-rep(repli,4)
    
    aux_alpha<-as.data.frame(alpha_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu=c(aux_alpha[1,2], aux_alpha[1,2], aux_alpha[2,2],aux_alpha[2,2]), intraTe=c(aux_alpha[3,2], aux_alpha[4,2], aux_alpha[3,2],aux_alpha[4,2]), interTu=c(aux_alpha[1,3], aux_alpha[1,4], aux_alpha[2,3], aux_alpha[2,4]), interTe=c(aux_alpha[3,3], aux_alpha[4,3], aux_alpha[3,4], aux_alpha[4,4]))
  }
  
  aux_alpha2<-rbind(aux_alpha2, aux2)
}

### Alpha Lower

# passing alphas to dataframe
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

aux_alphaLower<-as.data.frame(alphaLower_list[[1]])

aux_alphaLower2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_L=c(aux_alphaLower[1,2], aux_alphaLower[1,2], aux_alphaLower[2,2],aux_alphaLower[2,2]), intraTe_L=c(aux_alphaLower[3,2], aux_alphaLower[4,2], aux_alphaLower[3,2],aux_alphaLower[4,2]), interTu_L=c(aux_alphaLower[1,3], aux_alphaLower[1,4], aux_alphaLower[2,3], aux_alphaLower[2,4]), interTe_L=c(aux_alphaLower[3,3], aux_alphaLower[4,3], aux_alphaLower[3,4], aux_alphaLower[4,4]))

for(x in 2:length(lambda_list)){
  repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[x], split="_")[1])[4], split="[.]"))[1],split=""))[1]
  env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[x], split="_")[1])[4], split="[.]"))[1],split=""))[2]
  
  if(x==3 | x==4){# because there is no SR2 here
    regimeTu<-c("1","1")
    regimeTe<-c("4","5")
    Env<-rep(env, 2)
    Rep<-rep(repli,2)
    
    aux_alphaLower<-as.data.frame(alphaLower_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_L=c(aux_alphaLower[1,2], aux_alphaLower[1,2]), intraTe_L=c(aux_alphaLower[2,2], aux_alphaLower[3,2]), interTu_L=c(aux_alphaLower[1,3], aux_alphaLower[1,4]), interTe_L=c(aux_alphaLower[2,3], aux_alphaLower[3,3]))
    
    
  }else{
    regimeTu<-c("1","1", "2","2")
    regimeTe<-c("4","5", "4","5")
    Env<-rep(env, 4)
    Rep<-rep(repli,4)
    
    aux_alphaLower<-as.data.frame(alphaLower_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_L=c(aux_alphaLower[1,2], aux_alphaLower[1,2], aux_alphaLower[2,2],aux_alphaLower[2,2]), intraTe_L=c(aux_alphaLower[3,2], aux_alphaLower[4,2], aux_alphaLower[3,2],aux_alphaLower[4,2]), interTu_L=c(aux_alphaLower[1,3], aux_alphaLower[1,4], aux_alphaLower[2,3], aux_alphaLower[2,4]), interTe_L=c(aux_alphaLower[3,3], aux_alphaLower[4,3], aux_alphaLower[3,4], aux_alphaLower[4,4]))
  }
  
  aux_alphaLower2<-rbind(aux_alphaLower2, aux2)
}

### Alpha Upper

# passing alphas to dataframe
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

aux_alphaUpper<-as.data.frame(alphaUpper_list[[1]])

aux_alphaUpper2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_U=c(aux_alphaUpper[1,2], aux_alphaUpper[1,2], aux_alphaUpper[2,2],aux_alphaUpper[2,2]), intraTe_U=c(aux_alphaUpper[3,2], aux_alphaUpper[4,2], aux_alphaUpper[3,2],aux_alphaUpper[4,2]), interTu_U=c(aux_alphaUpper[1,3], aux_alphaUpper[1,4], aux_alphaUpper[2,3], aux_alphaUpper[2,4]), interTe_U=c(aux_alphaUpper[3,3], aux_alphaUpper[4,3], aux_alphaUpper[3,4], aux_alphaUpper[4,4]))

for(x in 2:length(lambda_list)){
  repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[x], split="_")[1])[4], split="[.]"))[1],split=""))[1]
  env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[x], split="_")[1])[4], split="[.]"))[1],split=""))[2]
  
  if(x==3 | x==4){# because there is no SR2 here
    regimeTu<-c("1","1")
    regimeTe<-c("4","5")
    Env<-rep(env, 2)
    Rep<-rep(repli,2)
    
    aux_alphaUpper<-as.data.frame(alphaUpper_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_U=c(aux_alphaUpper[1,2], aux_alphaUpper[1,2]), intraTe_U=c(aux_alphaUpper[2,2], aux_alphaUpper[3,2]), interTu_U=c(aux_alphaUpper[1,3], aux_alphaUpper[1,4]), interTe_U=c(aux_alphaUpper[2,3], aux_alphaUpper[3,3]))
    
    
  }else{
    regimeTu<-c("1","1", "2","2")
    regimeTe<-c("4","5", "4","5")
    Env<-rep(env, 4)
    Rep<-rep(repli,4)
    
    aux_alphaUpper<-as.data.frame(alphaUpper_list[[x]])
    
    aux2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_U=c(aux_alphaUpper[1,2], aux_alphaUpper[1,2], aux_alphaUpper[2,2],aux_alphaUpper[2,2]), intraTe_U=c(aux_alphaUpper[3,2], aux_alphaUpper[4,2], aux_alphaUpper[3,2],aux_alphaUpper[4,2]), interTu_U=c(aux_alphaUpper[1,3], aux_alphaUpper[1,4], aux_alphaUpper[2,3], aux_alphaUpper[2,4]), interTe_U=c(aux_alphaUpper[3,3], aux_alphaUpper[4,3], aux_alphaUpper[3,4], aux_alphaUpper[4,4]))
  }
  
  aux_alphaUpper2<-rbind(aux_alphaUpper2, aux2)
}

# Passing lambda to data frame
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[2]

Focal<-c("1","1","2","2","4","4","5","5")
Comp<-c("4","5","4","5","1","2","1","2")
Env<-rep(env, 8)
Rep<-rep(repli,8)

aux_lambda<-cbind(as.data.frame(lambda_list[[1]])[,c(3,4,5)],Focal,Comp, Env, Rep)

for(x in 2:length(lambda_list)){
  repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[x], split="_")[1])[6], split="[.]"))[1],split=""))[1]
  env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[x], split="_")[1])[6], split="[.]"))[1],split=""))[2]
  
  if(x==3 | x==4){# because there is no SR2 here
    Focal<-c("1","1","4","5")
    Comp<-c("4","5","1","1")
    Env<-rep(env, 4)
    Rep<-rep(repli,4)
    
    aux<-cbind(as.data.frame(lambda_list[[x]])[,c(3,4,5)],Focal,Comp, Env, Rep)
    
  }else{
    Focal<-c("1","1","2","2","4","4","5","5")
    Comp<-c("4","5","4","5","1","2","1","2")
    Env<-rep(env, 8)
    Rep<-rep(repli,8)
    
    aux<-cbind(as.data.frame(lambda_list[[x]])[,c(3,4,5)],Focal,Comp, Env, Rep)
  }
  
  aux_lambda<-rbind(aux_lambda, aux)
}


#Matching all the data

alphas_mat_E<-as.data.frame(cbind(aux_alpha2, aux_alphaLower2, aux_alphaUpper2))

str(lambda_intra_fixed)

#### adding lambda

alphas_mat_E$lambdaTu<-sapply(c(1:length(alphas_mat_E[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_E$regimeTu[x] & Comp==alphas_mat_E$regimeTe[x]) & Rep==alphas_mat_E$Rep[x] & Env==alphas_mat_E$Env[x] )
  
  auxi[1,1]
})

alphas_mat_E$lambdaTe<-sapply(c(1:length(alphas_mat_E[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_E$regimeTe[x] & Comp==alphas_mat_E$regimeTu[x]) & Rep==alphas_mat_E$Rep[x] & Env==alphas_mat_E$Env[x] )
  
  auxi[1,1]
})


alphas_mat_E$lambdaTu_L<-sapply(c(1:length(alphas_mat_E[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_E$regimeTu[x] & Comp==alphas_mat_E$regimeTe[x]) & Rep==alphas_mat_E$Rep[x] & Env==alphas_mat_E$Env[x] )
  
  auxi[1,2]
})

alphas_mat_E$lambdaTe_L<-sapply(c(1:length(alphas_mat_E[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_E$regimeTe[x] & Comp==alphas_mat_E$regimeTu[x]) & Rep==alphas_mat_E$Rep[x] & Env==alphas_mat_E$Env[x] )
  
  auxi[1,2]
})

alphas_mat_E$lambdaTu_U<-sapply(c(1:length(alphas_mat_E[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_E$regimeTu[x] & Comp==alphas_mat_E$regimeTe[x]) & Rep==alphas_mat_E$Rep[x] & Env==alphas_mat_E$Env[x] )
  
  auxi[1,3]
})

alphas_mat_E$lambdaTe_U<-sapply(c(1:length(alphas_mat_E[,1])), function(x){
  auxi<-subset(aux_lambda, (Focal==alphas_mat_E$regimeTe[x] & Comp==alphas_mat_E$regimeTu[x]) & Rep==alphas_mat_E$Rep[x] & Env==alphas_mat_E$Env[x] )
  
  auxi[1,3]
})

alphas_mat_E$Env2<-mapvalues(alphas_mat_E$Env, c("C","N"), c("Cd","N"))

# clean up the matrix, because it has a lot of repeated columns
alphas_mat_E<-alphas_mat_E[,c(1:8, 13:16,21:30)]

alphas_mat_E


colnames(alphas_mat_E)<-c("Tu_Regime", "Te_Regime", "Environment", "Replicate", "Tu_intra", "Te_intra", "Tu_inter", "Te_inter", "Tu_intra_L", "Te_intra_L", "Tu_inter_L", "Te_inter_L", "Tu_intra_U", "Te_intra_U", "Tu_inter_U", "Te_inter_U", "Tu_lambda", "Te_lambda","Tu_lambda_L", "Te_lambda_L","Tu_lambda_U", "Te_lambda_U")


alphas_mat_E_long<-gather(alphas_mat_E, parameter, value,Tu_intra:Te_lambda_U )

alphas_mat_E_long$category<-mapvalues(alphas_mat_E_long$parameter, c("Tu_intra", "Te_intra", "Tu_inter", "Te_inter", "Tu_intra_L", "Te_intra_L", "Tu_inter_L", "Te_inter_L", "Tu_intra_U", "Te_intra_U", "Tu_inter_U", "Te_inter_U", "Tu_lambda", "Te_lambda","Tu_lambda_L", "Te_lambda_L","Tu_lambda_U", "Te_lambda_U"), c("intra", "intra", "inter", "inter", "intra_L", "intra_L", "inter_L", "inter_L","intra_U", "intra_U", "inter_U", "inter_U","lambda","lambda","lambda_L","lambda_L","lambda_U","lambda_U"))

######### Predicting data
str(alphas_mat_E)

alphas_mat_E$Env2<-mapvalues(alphas_mat_E$Environment, c("C", "N"), c("Cd","N"))
str(ca)

# Since the lambda is from the log data
ca$pred_E<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR[x]) & Replicate==ca$Rep[x])$Tu_intra[1]
    alpha_ij<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Te_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Tu_inter[1]
    lambda<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_lambda[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_intra[1]
    alpha_ij<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Tu_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Te_inter[1]
    lambda<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_lambda[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*densF)
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]
    pred<-lambda*exp(-alpha_ij*densC)
  }
  
  pred
  
})
x<-1
ca$pred_E_L<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_intra_L[1]
    alpha_ij<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Te_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Tu_inter_L[1]
    lambda<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_lambda_L[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_intra_L[1]
    alpha_ij<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Tu_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Te_inter_L[1]
    lambda<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_lambda_L[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*densF)
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    pred<-lambda*exp(-alpha_ij*densC)
  }
  
  pred
  
})

ca$pred_E_U<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_intra_U[1]
    alpha_ij<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Te_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Tu_inter_U[1]
    lambda<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Tu_lambda_U[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_intra_U[1]
    alpha_ij<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Tu_Regime==ca$CompSR[x] & Replicate==ca$Rep[x])$Te_inter_U[1]
    lambda<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR[x] & Replicate==ca$Rep[x])$Te_lambda_U[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*densF)
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    pred<-lambda*exp(-alpha_ij*densC)
  }
  
  pred
  
})



### Predicting each density
density_aux<-seq(0, 10, by=(10/100))

pred_df_E<-as.data.frame(expand_grid(Density=density_aux, Tu_Regime=c(1,2), Te_Regime=c(4,5), Replicate=c(1:5), Environment=c("N", "C")))

pred_df_E$Tu_mean_intra<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_inter[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Tu_mean_inter<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_intra[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_inter[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_E$Density[x])
  
  pred
})


pred_df_E$Tu_intra_L<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_intra_L[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_inter_L[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_lambda_L[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Tu_inter_L<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_intra_L[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_inter_L[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_lambda_L[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Tu_intra_U<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_intra_U[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_inter_U[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_lambda_U[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Tu_inter_U<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_intra_U[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_inter_U[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Replicate==pred_df_E$Replicate[x])$Tu_lambda_U[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Te_mean_intra<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_inter[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Te_mean_inter<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_intra[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_inter[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_lambda[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Te_intra_L<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_intra_L[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_inter_L[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_lambda_L[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Te_inter_L<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_intra_L[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_inter_L[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_lambda_L[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Te_intra_U<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_intra_U[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_inter_U[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_lambda_U[1]
  
  pred<-lambda*exp(-alpha_i*pred_df_E$Density[x])
  
  pred
})

pred_df_E$Te_inter_U<-sapply(c(1:length(pred_df_E[,1])), function(x){
  alpha_i<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_intra_U[1]
  alpha_ij<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Tu_Regime==pred_df_E$Tu_Regime[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_inter_U[1]
  lambda<-subset(alphas_mat_E, Environment==pred_df_E$Environment[x] & Te_Regime==pred_df_E$Te_Regime[x] & Replicate==pred_df_E$Replicate[x])$Te_lambda_U[1]
  
  pred<-lambda*exp(-alpha_ij*pred_df_E$Density[x])
  
  pred
})

# Removing Tu evolved replicate 2 because there is no data
pred_df_E<-pred_df_E[-which(pred_df_E$Tu_Regime==2 & pred_df_E$Replicate==2),]



# Transforming everything bellow 0 into 0 for the lower interval

pred_df_E$Te_inter_L[which(pred_df_E$Te_inter_L<0)]<-0
pred_df_E$Te_intra_L[which(pred_df_E$Te_intra_L<0)]<-0
pred_df_E$Tu_inter_L[which(pred_df_E$Tu_inter_L<0)]<-0
pred_df_E$Tu_intra_L[which(pred_df_E$Tu_intra_L<0)]<-0

# 4 - Testing similarities in estimation

### Compare methods visually
# Putting names of selection regimes all the same
alphas_mat_D$Tu_Regime2<-alphas_mat_D$Tu_Regime
alphas_mat_D$Te_Regime2<-alphas_mat_D$Te_Regime


alphas_mat_D$Tu_Regime<-mapvalues(alphas_mat_D$Tu_Regime2, c("1","2","4","5"), c("SR1", "SR2", "SR4","SR5"))
alphas_mat_D$Te_Regime<-mapvalues(alphas_mat_D$Te_Regime2, c("1","2","4","5"), c("SR1", "SR2", "SR4","SR5"))


alphas_mat_E$Tu_Regime2<-alphas_mat_E$Tu_Regime
alphas_mat_E$Te_Regime2<-alphas_mat_E$Te_Regime


alphas_mat_E$Tu_Regime<-mapvalues(alphas_mat_E$Tu_Regime2, c("1","2","4","5"), c("SR1", "SR2", "SR4","SR5"))
alphas_mat_E$Te_Regime<-mapvalues(alphas_mat_E$Te_Regime2, c("1","2","4","5"), c("SR1", "SR2", "SR4","SR5"))


# doing the same with the environments

alphas_mat_D$Environment2<-alphas_mat_D$Environment

alphas_mat_D$Environment<-mapvalues(alphas_mat_D$Environment2, c("N","C"), c("N", "Cd"))

alphas_mat_E$Environment2<-alphas_mat_E$Environment

alphas_mat_E$Environment<-mapvalues(alphas_mat_E$Environment2, c("N","C"), c("N", "Cd"))


### Adding variable to say who estimated what

param_all_w0$Method<-"cxr"
param_all_B$Method<-"cxr lambda fixed"
param_all_C$Method<-"cxr lambda fixed, nested"
alphas_mat_D$Method<-"optim"
alphas_mat_E$Method<-"optim lambda fixed"

cols_to_join<-c("Tu_Regime", "Te_Regime", "Environment", "Replicate","Tu_lambda","Te_lambda", "Tu_intra","Te_intra", "Tu_inter", "Te_inter", "Method" )

comparison_methods<-rbind(param_all_w0[,cols_to_join],param_all_B[,cols_to_join],param_all_C[,cols_to_join], alphas_mat_D[,cols_to_join], alphas_mat_E[,cols_to_join] )

ggplot(comparison_methods, aes(x=Method, y=Tu_lambda, colour=Environment, fill=Environment, shape=Replicate))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(group=Method, fill=Environment), alpha=0.75, outlier.colour = NA)+
  geom_point(position = position_dodge2(0.5))+
  theme_plot+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Tu lambda")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))

ggplot(comparison_methods, aes(x=Method, y=Te_lambda, colour=Environment, fill=Environment, shape=Replicate))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(group=Method, fill=Environment), alpha=0.75, outlier.colour = NA)+
  geom_point(position = position_dodge2(0.5))+
  theme_plot+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Te lambda")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))


ggplot(comparison_methods, aes(x=Method, y=Tu_intra, colour=Environment, fill=Environment, shape=Replicate))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(group=Method, fill=Environment), alpha=0.75, outlier.colour = NA)+
  geom_point(position = position_dodge2(0.5))+
  theme_plot+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Tu intra")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))

ggplot(comparison_methods, aes(x=Method, y=Te_intra, colour=Environment, fill=Environment, shape=Replicate))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(group=Method, fill=Environment), alpha=0.75, outlier.colour = NA)+
  geom_point(position = position_dodge2(0.5))+
  theme_plot+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Te intra")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))

ggplot(comparison_methods, aes(x=Method, y=Tu_inter, colour=Environment, fill=Environment, shape=Replicate))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(group=Method, fill=Environment), alpha=0.75, outlier.colour = NA)+
  geom_point(position = position_dodge2(0.5))+
  theme_plot+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Tu inter")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))

ggplot(comparison_methods, aes(x=Method, y=Te_inter, colour=Environment, fill=Environment, shape=Replicate))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(group=Method, fill=Environment), alpha=0.75, outlier.colour = NA)+
  geom_point(position = position_dodge2(0.5))+
  theme_plot+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Te inter")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))



## Estimate distance between predicted and observed

#Since I can't really know what is the best approach, I will estimate the predicted vs observed for each method and use that as metric to define which method to use in the results

##### Predict values


ca$FocalSR3<-mapvalues(ca$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

ca$CompSR3<-mapvalues(ca$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

ca$Env3<-mapvalues(ca$Env, c("N", "Cd"), c("N", "C"))

# Since the lambda is from the log data
ca$pred_A<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(param_all_w0, Environment==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Tu_intra[1]
    alpha_ij<-subset(param_all_w0, Environment==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Te_Regime==as.character(ca$CompSR3[x]) & Replicate==ca$Rep[x])$Tu_inter[1]
    lambda<-subset(param_all_w0, Environment==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Tu_lambda[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(param_all_w0, Environment==as.character(ca$Env[x]) & Te_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Te_intra[1]
    alpha_ij<-subset(param_all_w0, Environment==as.character(ca$Env[x]) & Te_Regime==as.character(ca$FocalSR3[x]) & Tu_Regime==as.character(ca$CompSR3[x]) & Replicate==ca$Rep[x])$Te_inter[1]
    lambda<-subset(param_all_w0, Environment==as.character(ca$Env[x]) & Te_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Te_lambda[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*exp(-alpha_i*(densF)-alpha_ij*densC)
  }
    
  pred
  
})


ca$pred_B<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(param_all_B, Environment==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Tu_intra[1]
    alpha_ij<-subset(param_all_B, Environment==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Te_Regime==as.character(ca$CompSR3[x]) & Replicate==ca$Rep[x])$Tu_inter[1]
    lambda<-subset(param_all_B, Environment==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Tu_lambda[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(param_all_B, Environment==as.character(ca$Env[x]) & Te_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Te_intra[1]
    alpha_ij<-subset(param_all_B, Environment==as.character(ca$Env[x]) & Te_Regime==as.character(ca$FocalSR3[x]) & Tu_Regime==as.character(ca$CompSR3[x]) & Replicate==ca$Rep[x])$Te_inter[1]
    lambda<-subset(param_all_B, Environment==as.character(ca$Env[x]) & Te_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Te_lambda[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*exp(-alpha_i*(densF)-alpha_ij*densC)
  }
    
  pred
  
})

ca$pred_C<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(param_all_C, Environment==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Tu_intra[1]
    alpha_ij<-subset(param_all_C, Environment==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Te_Regime==as.character(ca$CompSR3[x]) & Replicate==ca$Rep[x])$Tu_inter[1]
    lambda<-subset(param_all_C, Environment==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Tu_lambda[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(param_all_C, Environment==as.character(ca$Env[x]) & Te_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Te_intra[1]
    alpha_ij<-subset(param_all_C, Environment==as.character(ca$Env[x]) & Te_Regime==as.character(ca$FocalSR3[x]) & Tu_Regime==as.character(ca$CompSR3[x]) & Replicate==ca$Rep[x])$Te_inter[1]
    lambda<-subset(param_all_C, Environment==as.character(ca$Env[x]) & Te_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Te_lambda[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*exp(-alpha_i*(densF)-alpha_ij*densC)
  }
    
  pred
  
})

ca$pred_D<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Tu_intra[1]
    alpha_ij<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR3[x] & Te_Regime==ca$CompSR3[x] & Replicate==ca$Rep[x])$Tu_inter[1]
    lambda<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR3[x] & Replicate==ca$Rep[x])$Tu_lambda[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR3[x] & Replicate==ca$Rep[x])$Te_intra[1]
    alpha_ij<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR3[x] & Tu_Regime==ca$CompSR3[x] & Replicate==ca$Rep[x])$Te_inter[1]
    lambda<-subset(alphas_mat_D, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR3[x] & Replicate==ca$Rep[x])$Te_lambda[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*exp(-alpha_i*(densF)-alpha_ij*densC)
  }
    
  pred
  
})

ca$pred_E<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    alpha_i<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==as.character(ca$FocalSR3[x]) & Replicate==ca$Rep[x])$Tu_intra[1]
    alpha_ij<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR3[x] & Te_Regime==ca$CompSR3[x] & Replicate==ca$Rep[x])$Tu_inter[1]
    lambda<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Tu_Regime==ca$FocalSR3[x] & Replicate==ca$Rep[x])$Tu_lambda[1]
    
  }else if(ca$Focalfemale[x]=="Te"){
    alpha_i<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR3[x] & Replicate==ca$Rep[x])$Te_intra[1]
    alpha_ij<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR3[x] & Tu_Regime==ca$CompSR3[x] & Replicate==ca$Rep[x])$Te_inter[1]
    lambda<-subset(alphas_mat_E, Env2==as.character(ca$Env[x]) & Te_Regime==ca$FocalSR3[x] & Replicate==ca$Rep[x])$Te_lambda[1]
  }
  
  if(ca$Type[x]=="INTRA"){
    densF<-ca$Dens[x]
    pred<-lambda*exp(-alpha_i*densF)
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*exp(-alpha_i*(densF)-alpha_ij*densC)
  }
    
  pred
  
})

### Calculate distances
#Do not forget that this is the log of GR +1

euclidean <- function(a, b) sqrt(sum((a - b)^2))

ca$distA<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_A[x], ca$GrowthRateOA[x])

  euc
})

ca$distB<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_B[x], ca$GrowthRateOA[x])

  euc
})

ca$distC<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_C[x], ca$GrowthRateOA[x])

  euc
})

ca$distD<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_D[x], ca$GrowthRateOA[x])

  euc
})

ca$distE<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_E[x], ca$GrowthRateOA[x])

  euc
})

hist(ca$distA)
hist(ca$distB)
hist(ca$distC)
hist(ca$distD)
hist(ca$distE)

sum(ca$distA, na.rm = TRUE)
sum(ca$distB, na.rm = TRUE)
sum(ca$distC, na.rm = TRUE)
sum(ca$distD, na.rm = TRUE)
sum(ca$distE, na.rm = TRUE)

#The smaller sum of euclidean distance is with cxr package

##### Plotting distance
str(ca)
distance_sum<-pivot_longer(ca[, c(48:52)], cols = c(1:5),names_to = "method", values_to = "distance")

dist_sum<-distance_sum %>% group_by(method) %>% summarize(mean=mean(distance, na.rm=TRUE), se=sd(distance, na.rm=TRUE)/sqrt(n()))

ggplot(distance_sum, aes(x=method, y=distance, colour=method, fill=method))+
   geom_boxplot(colour="black", outlier.colour = NA)+
  geom_point(alpha=0.10, position=position_dodge2(0.5), colour="black", shape=21)+
  theme_bw()+
  theme_plot+
  scale_x_discrete(labels=c("cxr", "cxr lambda\nfixed", "cxr \nnested", "optim", "optim \nlambda fixed"), name="Method")+
  scale_color_brewer(palette = "Spectral")+
  scale_fill_brewer(palette = "Spectral")+
  theme(legend.position = "none")+
  geom_text(data=dist_sum, aes(x=method, label=paste(round(mean,3), round(se,3), sep="\n+/-")), y=22, colour="black")+
  scale_y_continuous(name="Estimated euclidean distance\n (predicted-observed")

save_plot("./Plots/FigS2.pdf", width = 20, height=10)
save_plot("./Plots/FigS2.png", width = 20, height=10)

