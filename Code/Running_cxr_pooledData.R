#' ---
#' title: "R Notebook"
#' output: html_notebook
#' ---
#' 
#' # Functions, general information and packages
#' 
## ---------------------------
rm(list=ls())
library(plyr)
library(tidyverse)
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
library(MASS)
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

theme_plots<-theme(axis.text = element_text(size=14), axis.title = element_text(size=14, face="bold"), legend.text = element_text(size=12), strip.text = element_text(size=14), plot.title = element_text(size=14, face="bold"), panel.grid=element_line(colour="white"), panel.background = element_rect(fill="white") , axis.line = element_line(linewidth = 0.5, linetype = "solid",colour = "black"), strip.background = element_rect(fill="white"))

save_plot<-function(dir, width=15, height=10, ...){
  ggsave(dir, width = width, height = height, units = c("cm"))
}

# Creating vectors with regime names to use for plots
regimeTu<-c("Tu \ncontrol", "Tu evolved \n in cadmium")
names(regimeTu)<-c("SR1", "SR2")

regimeTe<-c("Te \n control", "Te evolved \n in cadmium")
names(regimeTe)<-c("SR4", "SR5")

#' 
#' # Importing data
#' This chunk is used to import data, checking data and creating the columns that are necessary
#' 
## ---------------------------
ca<-read.csv(file = "../Data/CompetitiveAbility_Cd_G40_submit.csv", header=TRUE) # cdata from the competitive ability

str(ca) 
# Summary of the data to be sure that everything is ok!
summary(as.factor(ca$FocalSR))

ca$Block2<-as.factor(ca$Block)
ca$Rep2<-as.factor(ca$Rep)
ca$Disk2<-as.factor(ca$Disk)
ca$Leaf2<-as.factor(ca$Leaf)
ca$Env2<-as.factor(ca$Env)
ca$FocalSR2<-as.factor(ca$FocalSR)
ca$CompSR2<-as.factor(ca$CompSR)
ca$Type2<-as.factor(ca$Type)
ca$Focal_Female2<-as.factor(ca$Focalfemale)

#### Creating columns that are needed
# This computes the number of focal females that were alive after 3 days (Tu and Te)
ca$Nr_Focal_Females_Tu_Alive_G0<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    if(ca$Type[x]=="INTRA"){
      a<-ca$Dens[x]-ca$FocalDead[x]-ca$FocalDrowned[x]-ca$FocalMissing[x]
    }else
      a<-1-ca$FocalDead[x]-ca$FocalDrowned[x]-ca$FocalMissing[x]
    
  }else
    a<-NA
})

ca$Nr_Focal_Females_Te_Alive_G0<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Te"){
    if(ca$Type[x]=="INTRA"){
      a<-ca$Dens[x]-ca$FocalDead[x]-ca$FocalDrowned[x]-ca$FocalMissing[x]
    }else
      a<-1-ca$FocalDead[x]-ca$FocalDrowned[x]-ca$FocalMissing[x]
    
  }else
    a<-NA
})

# Computes the number of competitor females alive after 3 days (Tu amd Te)
ca$Num_Comp_Tu_Alive_G0<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Te"){
    if(ca$Type[x]=="INTER"){
      a<-ca$Dens[x]-ca$NumbDeadComp[x]-1
    }else
      a<-NA
    
  }else
    a<-NA
})


ca$Num_Comp_Te_Alive_G0<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Focalfemale[x]=="Tu"){
    if(ca$Type[x]=="INTER"){
      a<-ca$Dens[x]-ca$NumbDeadComp[x]-1
    }else
      a<-NA
    
  }else
    a<-NA
})
#This just creates the number of focal females irrespective of the species
ca$Nr_Focal_Females_G0<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Type[x]=="INTRA"){
    a<-ca$Dens[x]
  }else
    a<-1
  
})
str(ca)
summary(ca)

# which(ca$Num_Comp_Te_Alive_G0<0)
# which(ca$Num_Comp_Tu_Alive_G0<0)
# which(ca$Nr_Focal_Females_Tu_Alive_G0<0)
# which(ca$Nr_Focal_Females_Te_Alive_G0<0)

# Remove the lines where the number of competitors or focal was lower than 0
ca<-ca[-c(which(ca$Num_Comp_Te_Alive_G0<0),which(ca$Num_Comp_Tu_Alive_G0<0), which(ca$Nr_Focal_Females_Te_Alive_G0<0) ),]

#Creating the columns with the correct number of competitors. For conspecifics its always the same as the density to calculate the growth rate. For heterospecifics its always 1 female with X competitors, and DensFocal2 its also to do the same for the conspecifics 

ca$DensFocal<-sapply(c(1:dim(ca)[1]), function(x){
  if(ca$Type[x]=="INTRA"){
    a<-ca$Dens[x]-1
  }else if(ca$Type[x]=="INTER"){
    a<-1
  }
  
  a
})

ca$DensComp<-sapply(c(1:dim(ca)[1]), function(x){
  if(ca$Type[x]=="INTRA"){
    a<-ca$Dens[x]-1
  }else if(ca$Type[x]=="INTER"){
    a<-ca$Dens[x]-1
  }
  
  a
})

ca$CompSR3<-sapply(c(1:dim(ca)[1]), function(x){
  if(ca$Type[x]=="INTRA"){
    a<-ca$FocalSR[x]
  }else if(ca$Type[x]=="INTER"){
    a<-ca$CompSR[x]
  }
  
  a
})

ca$CompSR3<-as.factor(ca$CompSR3)

#### Estimate growth rate

ca[,c("Nr_Focal_Females_G0", "Dens", "Type")]


ca$GrowthRateOA<-sapply(c(1:length(ca[,1])), function(x){
  #print(x)
  if(ca$Focal_Female[x]=="Tu"){
    a<-ca$TuFemales[x]/ca$DensFocal[x]
  }else if(ca$Focal_Female[x]=="Te"){
    a<-ca$TeFemales[x]/ca$DensFocal[x]
  }else
    a<-NA
  
  a
})


#' 
#' # Running cxr
#' 
#' ### Setup cxr
#' Setting the seed number and the number of bootstrap samples.
#' 
## ---------------------------
set.seed(1809)
bootN<-5000
eval<- TRUE
dir.create("./Analyses/cxr_normal_REP_best", showWarnings = FALSE)


#' 
#' ## Set up data
#' 
#' ### Creating a data frame for the no cadmium environment
#' 
## ---------------------------
forCXR_N<-subset(ca, Env=="N")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "TuFemales")]
  
  forCXR_N$Focal<-mapvalues(forCXR_N$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
  forCXR_N$CompSR2<-mapvalues(forCXR_N$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
  
  
  # The CompSR2 column indicated the selection regime of the interspecific competitor and NA when the competitor was from the same selection regime. This creates a column indicating that has the focal selection regime for intraspecific competition and the selection regime of the other species when there is interspecific competition
  forCXR_N$Comp<-sapply(c(1:length(forCXR_N[,1])), function(x){
    if(is.na(forCXR_N$CompSR2[x])){
      a<- forCXR_N$Focal[x]
    }else{
      a<-forCXR_N$CompSR2[x]
    }
    
    a
  })
  
  # Fills in with the number of competitors present in the disk. It always correspond to -1 competitor in the case of intraspecific, for interspecific its always 0 for focal and Density -1 for the competitor.
  aux<-data.frame(SR1=rep(0, length(forCXR_N[,1])), SR2=rep(0, length(forCXR_N[,1])), SR4=rep(0, length(forCXR_N[,1])), SR5=rep(0, length(forCXR_N[,1])))
  
  for(i in 1:length(forCXR_N[,1])){
    #coluna onde por focais
    colunaF<-which(colnames(aux)==forCXR_N$Focal[i])
    #coluna onde por competidors
    colunaC<-which(colnames(aux)==forCXR_N$Comp[i])
    
    #if its the same regime
    if(forCXR_N$Focal[i]==forCXR_N$Comp[i]){
      aux[i,colunaF]<-forCXR_N$Dens[i]-1
    }else{ #if it is heterospecific then its -1 for the competitors (because of the focal)
      aux[i,colunaC]<-forCXR_N$Dens[i]-1
      aux[i, colunaF]<-0
    }
    
  }
  
  forCXR_N<-cbind(forCXR_N, aux)
  
  forCXR_N$fitness<-sapply(c(1:length(forCXR_N[,1])), function(x){
    colF<-which(colnames(forCXR_N)==forCXR_N$Focal[x])
    
    if(forCXR_N$Focal[x]=="SR1"){
      #If the density of the focal is 0, then we divide by 1
      if(forCXR_N$SR1[x]==0){
        a<-forCXR_N$TuFemales[x]/1
      }else{
        a<-forCXR_N$TuFemales[x]/forCXR_N$SR1[x]
      }
    } else if(forCXR_N$Focal[x]=="SR2"){
      if(forCXR_N$SR2[x]==0){
        a<-forCXR_N$TuFemales[x]/1
      }else{
        a<-forCXR_N$TuFemales[x]/forCXR_N$SR2[x]
      }
    } else if(forCXR_N$Focal[x]=="SR4"){
      if(forCXR_N$SR4[x]==0){
        a<-forCXR_N$TeFemales[x]/1
      }else{
        a<-forCXR_N$TeFemales[x]/forCXR_N$SR4[x]
      }
    } else if(forCXR_N$Focal[x]=="SR5"){
      if(forCXR_N$SR5[x]==0){
        a<-forCXR_N$TeFemales[x]/1
      }else{
        a<-forCXR_N$TeFemales[x]/forCXR_N$SR5[x]
      }
    }
    
    a
  })
  
  #removing rows for which there is no data for fitness
  forCXR_N<-forCXR_N[-which(is.na(forCXR_N$fitness)),]
  
  subset(forCXR_N, (Focal=="SR4" & Comp=="SR1"))
  
  # all data gets +1 because of the 0 problem
  forCXR_N$fitness<-forCXR_N$fitness+1
  
  # vector that tells which are the selection regimes, the columns have to have the same name
  my.reg <- c("SR1", "SR2","SR4","SR5")


#' 
#' ### Creating a data frame for the cadmium environment
#' 
## ---------------------------
forCXR_Cd<-subset(ca, Env=="Cd")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "TuFemales")]
  
  forCXR_Cd$Focal<-mapvalues(forCXR_Cd$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
  forCXR_Cd$CompSR2<-mapvalues(forCXR_Cd$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
  
  
  # The CompSR2 column indicated the selection regime of the interspecific competitor and NA when the competitor was from the same selection regime. This creates a column indicating that has the focal selection regime for intraspecific competition and the selection regime of the other species when there is interspecific competition
  forCXR_Cd$Comp<-sapply(c(1:nrow(forCXR_Cd)), function(x){
    if(is.na(forCXR_Cd$CompSR2[x])){
      a<- forCXR_Cd$Focal[x]
    }else{
      a<-forCXR_Cd$CompSR2[x]
    }
    
    a
  })
  
  # Fills in with the number of competitors present in the disk. It always correspond to -1 competitor in the case of intraspecific, for interspecific its always 0 for focal and Density -1 for the competitor.
  aux<-data.frame(SR1=rep(0, nrow(forCXR_Cd)), SR2=rep(0, nrow(forCXR_Cd)), SR4=rep(0, nrow(forCXR_Cd)), SR5=rep(0, nrow(forCXR_Cd)))
  
  for(i in 1:length(forCXR_Cd[,1])){
    #coluna onde por focais
    colunaF<-which(colnames(aux)==forCXR_Cd$Focal[i])
    #coluna onde por competidors
    colunaC<-which(colnames(aux)==forCXR_Cd$Comp[i])
    
    #if its the same regime
    if(forCXR_Cd$Focal[i]==forCXR_Cd$Comp[i]){
      aux[i,colunaF]<-forCXR_Cd$Dens[i]-1
    }else{ #if it is heterospecific then its -1 for the competitors (because of the focal)
      aux[i,colunaC]<-forCXR_Cd$Dens[i]-1
      aux[i, colunaF]<-0
    }
    
  }
  
  forCXR_Cd<-cbind(forCXR_Cd, aux)
  
  forCXR_Cd$fitness<-sapply(c(1:length(forCXR_Cd[,1])), function(x){
    colF<-which(colnames(forCXR_Cd)==forCXR_Cd$Focal[x])
    
    if(forCXR_Cd$Focal[x]=="SR1"){
      #If the density of the focal is 0, then we divide by 1
      if(forCXR_Cd$SR1[x]==0){
        a<-forCXR_Cd$TuFemales[x]/1
      }else{
        a<-forCXR_Cd$TuFemales[x]/forCXR_Cd$SR1[x]
      }
    } else if(forCXR_Cd$Focal[x]=="SR2"){
      if(forCXR_Cd$SR2[x]==0){
        a<-forCXR_Cd$TuFemales[x]/1
      }else{
        a<-forCXR_Cd$TuFemales[x]/forCXR_Cd$SR2[x]
      }
    } else if(forCXR_Cd$Focal[x]=="SR4"){
      if(forCXR_Cd$SR4[x]==0){
        a<-forCXR_Cd$TeFemales[x]/1
      }else{
        a<-forCXR_Cd$TeFemales[x]/forCXR_Cd$SR4[x]
      }
    } else if(forCXR_Cd$Focal[x]=="SR5"){
      if(forCXR_Cd$SR5[x]==0){
        a<-forCXR_Cd$TeFemales[x]/1
      }else{
        a<-forCXR_Cd$TeFemales[x]/forCXR_Cd$SR5[x]
      }
    }
    
    a
  })
  
  #removing rows for which there is no data for fitness
  forCXR_Cd<-forCXR_Cd[-which(is.na(forCXR_Cd$fitness)),]
  
  subset(forCXR_Cd, (Focal=="SR4" & Comp=="SR1"))
  
  # all data gets +1 because of the 0 problem
  forCXR_Cd$fitness<-forCXR_Cd$fitness+1
  
  # vector that tells which are the selection regimes, the columns have to have the same name
  my.reg <- c("SR1", "SR2","SR4","SR5")
  
 

#' 
#' ### Data for Pooled estimates for no cadmium
#' 
## ---------------------------
  # Do list per replicate and environment
  Rep<-list(SR1= subset(forCXR_N, Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N,  Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

#' 
#' ### Data for Pooled estimates for cadmium
#' 
## ---------------------------
# Do list per replicate and environment
  Rep_Cd<-list(SR1= subset(forCXR_Cd, Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd,  Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

#' 
#' ## Running cxr for pooled data
#' The initial values used were obtained from the best models identified in the code of the file: Parameter_exploration
#' #### No cadmium
#' 
## ---------------------------
obs.w0<-cxr_pm_multifit(data = Rep,
                          focal_column = my.reg,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(lambda = 1.6,
                                                alpha_intra = 0,
                                                alpha_inter = -0.1),
                          fixed_terms = NULL,
                          # no standard errors
                           bootstrap_samples = bootN)

obs.w0$lambda
obs.w0$lambda_standard_error
obs.w0$alpha_matrix
obs.w0$alpha_matrix_standard_error
obs.w0$log_likelihood


#' 
#' ###### Saving parameter and errors
#' Storing the parameter estimates and the lower and upper bounds.
#' 
## ---------------------------
cxr_param_REP<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Environment=c("N"))
  cxr_param_REP$Tu_lambda<-0
  cxr_param_REP$Te_lambda<-0
  cxr_param_REP$Tu_intra<-0
  cxr_param_REP$Te_intra<-0
  cxr_param_REP$Tu_inter<-0
  cxr_param_REP$Te_inter<-0
  
  cxr_param_REP[,"Tu_lambda"]<-obs.w0$lambda[c(1,2,1,2)]
  cxr_param_REP[,"Te_lambda"]<-obs.w0$lambda[c(3,3,4,4)]
  
  cxr_param_REP[,"Tu_intra"]<-rep(c(obs.w0$alpha_matrix[1,1], obs.w0$alpha_matrix[2,2]), 2)
  cxr_param_REP[,"Te_intra"]<-rep(c(obs.w0$alpha_matrix[3,3], obs.w0$alpha_matrix[4,4]), each=2)
  
  cxr_param_REP[,"Tu_inter"]<-c(obs.w0$alpha_matrix[1,3], obs.w0$alpha_matrix[2,3],obs.w0$alpha_matrix[1,4], obs.w0$alpha_matrix[2,4])
  cxr_param_REP[,"Te_inter"]<-c(obs.w0$alpha_matrix[3,1], obs.w0$alpha_matrix[3,2],obs.w0$alpha_matrix[4,1], obs.w0$alpha_matrix[4,2])

#' 
#' #### Data frame with lower estimates
## ---------------------------
### Here we will apply the same reasoning but now estimating the lower estimates, using the stats from cxr. Lower bounds correspond to mean-error.
  
  cxr_param_REP_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Environment=c("N"))
  cxr_param_REP_lower$Tu_lambda<-0
  cxr_param_REP_lower$Te_lambda<-0
  cxr_param_REP_lower$Tu_intra<-0
  cxr_param_REP_lower$Te_intra<-0
  cxr_param_REP_lower$Tu_inter<-0
  cxr_param_REP_lower$Te_inter<-0
  
  cxr_param_REP_lower[,"Tu_lambda"]<-rep(c(obs.w0$lambda[1]-obs.w0$lambda_standard_error[1], obs.w0$lambda[2]-obs.w0$lambda_standard_error[2]), 2)
  cxr_param_REP_lower[,"Te_lambda"]<-rep(c(obs.w0$lambda[3]-obs.w0$lambda_standard_error[3], obs.w0$lambda[4]-obs.w0$lambda_standard_error[4]), each=2)
  
  cxr_param_REP_lower[,"Tu_intra"]<-rep(c(obs.w0$alpha_matrix[1,1]-obs.w0$alpha_matrix_standard_error[1,1], obs.w0$alpha_matrix[2,2]-obs.w0$alpha_matrix_standard_error[2,2]), 2)
  cxr_param_REP_lower[,"Te_intra"]<-rep(c(obs.w0$alpha_matrix[3,3]-obs.w0$alpha_matrix_standard_error[3,3],obs.w0$alpha_matrix[4,4]-obs.w0$alpha_matrix_standard_error[4,4]), each=2)
  
  
  cxr_param_REP_lower[,"Tu_inter"]<-c(obs.w0$alpha_matrix[1,3]-obs.w0$alpha_matrix_standard_error[1,3], obs.w0$alpha_matrix[2,3]-obs.w0$alpha_matrix_standard_error[2,3], obs.w0$alpha_matrix[1,4]-obs.w0$alpha_matrix_standard_error[1,4], obs.w0$alpha_matrix[2,4]-obs.w0$alpha_matrix_standard_error[2,4])
  cxr_param_REP_lower[,"Te_inter"]<-c(obs.w0$alpha_matrix[3,1]-obs.w0$alpha_matrix_standard_error[3,1], obs.w0$alpha_matrix[3,2]-obs.w0$alpha_matrix_standard_error[3,2],obs.w0$alpha_matrix[4,1]-obs.w0$alpha_matrix_standard_error[4,1], obs.w0$alpha_matrix[4,2]-obs.w0$alpha_matrix_standard_error[4,2])

#' 
#' 
#' #### Data frame with upper estimates
## ---------------------------
### Here we will apply the same reasoning but now estimating the upper estimates, using the stats from cxr. Upper bounds correspond to mean-error.
cxr_param_REP_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Environment=c("N"))
  cxr_param_REP_upper$Tu_lambda<-0
  cxr_param_REP_upper$Te_lambda<-0
  cxr_param_REP_upper$Tu_intra<-0
  cxr_param_REP_upper$Te_intra<-0
  cxr_param_REP_upper$Tu_inter<-0
  cxr_param_REP_upper$Te_inter<-0
  
  cxr_param_REP_upper[,"Tu_lambda"]<-rep(c(obs.w0$lambda[1]+obs.w0$lambda_standard_error[1], obs.w0$lambda[2]+obs.w0$lambda_standard_error[2]), 2)
  cxr_param_REP_upper[,"Te_lambda"]<-rep(c(obs.w0$lambda[3]+obs.w0$lambda_standard_error[3], obs.w0$lambda[4]+obs.w0$lambda_standard_error[4]), each=2)
  
  cxr_param_REP_upper[,"Tu_intra"]<-rep(c(obs.w0$alpha_matrix[1,1]+obs.w0$alpha_matrix_standard_error[1,1], obs.w0$alpha_matrix[2,2]+obs.w0$alpha_matrix_standard_error[2,2]), 2)
  cxr_param_REP_upper[,"Te_intra"]<-rep(c(obs.w0$alpha_matrix[3,3]+obs.w0$alpha_matrix_standard_error[3,3], obs.w0$alpha_matrix[4,4]+obs.w0$alpha_matrix_standard_error[4,4]), each=2)
  
  
  cxr_param_REP_upper[,"Tu_inter"]<-c(obs.w0$alpha_matrix[1,3]+obs.w0$alpha_matrix_standard_error[1,3], obs.w0$alpha_matrix[2,3]+obs.w0$alpha_matrix_standard_error[2,3],obs.w0$alpha_matrix[1,4]+obs.w0$alpha_matrix_standard_error[1,4], obs.w0$alpha_matrix[2,4]+obs.w0$alpha_matrix_standard_error[2,4])
  cxr_param_REP_upper[,"Te_inter"]<-c(obs.w0$alpha_matrix[3,1]+obs.w0$alpha_matrix_standard_error[3,1], obs.w0$alpha_matrix[3,2]+obs.w0$alpha_matrix_standard_error[3,2],obs.w0$alpha_matrix[4,1]+obs.w0$alpha_matrix_standard_error[4,1], obs.w0$alpha_matrix[4,2]+obs.w0$alpha_matrix_standard_error[4,2])

#' 
#' #### Cadmium
#' Running cxr for the cadmium environment
#' 
## ---------------------------
obs.Cd_w0<-cxr_pm_multifit(data = Rep_Cd,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = 1.9,
                                                alpha_intra = 0.1,
                                                alpha_inter = 0.1),
                             fixed_terms = NULL,
                             # no standard errors
                              bootstrap_samples = bootN)

obs.Cd_w0$lambda
obs.Cd_w0$lambda_standard_error
obs.Cd_w0$alpha_matrix
obs.Cd_w0$alpha_matrix_standard_error
obs.Cd_w0$log_likelihood

#' 
#' ###### Saving parameter and errors
#' Storing data into data frames
## ---------------------------
cxr_param_REP_C<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Environment=c("Cd"))
  cxr_param_REP_C$Tu_lambda<-0
  cxr_param_REP_C$Te_lambda<-0
  cxr_param_REP_C$Tu_intra<-0
  cxr_param_REP_C$Te_intra<-0
  cxr_param_REP_C$Tu_inter<-0
  cxr_param_REP_C$Te_inter<-0
  
  cxr_param_REP_C[,"Tu_lambda"]<-obs.Cd_w0$lambda[c(1,2,1,2)]
  cxr_param_REP_C[,"Te_lambda"]<-obs.Cd_w0$lambda[c(3,3,4,4)]
  
  cxr_param_REP_C[,"Tu_intra"]<-rep(c(obs.Cd_w0$alpha_matrix[1,1], obs.Cd_w0$alpha_matrix[2,2]), 2)
  cxr_param_REP_C[,"Te_intra"]<-rep(c(obs.Cd_w0$alpha_matrix[3,3], obs.Cd_w0$alpha_matrix[4,4]), each=2)
  
  cxr_param_REP_C[,"Tu_inter"]<-c(obs.Cd_w0$alpha_matrix[1,3], obs.Cd_w0$alpha_matrix[2,3],obs.Cd_w0$alpha_matrix[1,4], obs.Cd_w0$alpha_matrix[2,4])
  cxr_param_REP_C[,"Te_inter"]<-c(obs.Cd_w0$alpha_matrix[3,1], obs.Cd_w0$alpha_matrix[3,2],obs.Cd_w0$alpha_matrix[4,1], obs.Cd_w0$alpha_matrix[4,2])

#' 
#' #### Data frame with lower estimates
## ---------------------------
### Here we will apply the same reasoning but now estimating the lower estimates, using the stats from cxr. Lower bounds correspond to mean-error.
cxr_param_REP_C_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Environment=c("Cd"))
  cxr_param_REP_C_lower$Tu_lambda<-0
  cxr_param_REP_C_lower$Te_lambda<-0
  cxr_param_REP_C_lower$Tu_intra<-0
  cxr_param_REP_C_lower$Te_intra<-0
  cxr_param_REP_C_lower$Tu_inter<-0
  cxr_param_REP_C_lower$Te_inter<-0
  
  cxr_param_REP_C_lower[,"Tu_lambda"]<-rep(c(obs.Cd_w0$lambda[1]-obs.Cd_w0$lambda_standard_error[1], obs.Cd_w0$lambda[2]-obs.Cd_w0$lambda_standard_error[2]), 2)
  cxr_param_REP_C_lower[,"Te_lambda"]<-rep(c(obs.Cd_w0$lambda[3]-obs.Cd_w0$lambda_standard_error[3], obs.Cd_w0$lambda[4]-obs.Cd_w0$lambda_standard_error[4]), each=2)
  
  
  cxr_param_REP_C_lower[,"Tu_intra"]<-rep(c(obs.Cd_w0$alpha_matrix[1,1]-obs.Cd_w0$alpha_matrix_standard_error[1,1], obs.Cd_w0$alpha_matrix[2,2]-obs.Cd_w0$alpha_matrix_standard_error[2,2]), 2)
  cxr_param_REP_C_lower[,"Te_intra"]<-rep(c(obs.Cd_w0$alpha_matrix[3,3]-obs.Cd_w0$alpha_matrix_standard_error[3,3], obs.Cd_w0$alpha_matrix[4,4]-obs.Cd_w0$alpha_matrix_standard_error[4,4]), each=2)
  
  
  cxr_param_REP_C_lower[,"Tu_inter"]<-c(obs.Cd_w0$alpha_matrix[1,3]-obs.Cd_w0$alpha_matrix_standard_error[1,3],obs.Cd_w0$alpha_matrix[2,3]-obs.Cd_w0$alpha_matrix_standard_error[2,3],obs.Cd_w0$alpha_matrix[1,4]-obs.Cd_w0$alpha_matrix_standard_error[1,4], obs.Cd_w0$alpha_matrix[2,4]-obs.Cd_w0$alpha_matrix_standard_error[2,4])
  
  cxr_param_REP_C_lower[,"Te_inter"]<-c(obs.Cd_w0$alpha_matrix[3,1]-obs.Cd_w0$alpha_matrix_standard_error[3,1], obs.Cd_w0$alpha_matrix[3,2]-obs.Cd_w0$alpha_matrix_standard_error[3,2],obs.Cd_w0$alpha_matrix[4,1]-obs.Cd_w0$alpha_matrix_standard_error[4,1], obs.Cd_w0$alpha_matrix[4,2]-obs.Cd_w0$alpha_matrix_standard_error[4,2])

#' 
#' #### Data frame with upper estimates
## ---------------------------
### Here we will apply the same reasoning but now estimating the upper estimates, using the stats from cxr. Upper bounds correspond to mean-error.
cxr_param_REP_C_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Environment=c("Cd"))
  cxr_param_REP_C_upper$Tu_lambda<-0
  cxr_param_REP_C_upper$Te_lambda<-0
  cxr_param_REP_C_upper$Tu_intra<-0
  cxr_param_REP_C_upper$Te_intra<-0
  cxr_param_REP_C_upper$Tu_inter<-0
  cxr_param_REP_C_upper$Te_inter<-0
  
  cxr_param_REP_C_upper[,"Tu_lambda"]<-rep(c(obs.Cd_w0$lambda[1]+obs.Cd_w0$lambda_standard_error[1], obs.Cd_w0$lambda[2]+obs.Cd_w0$lambda_standard_error[2]), 2)
  cxr_param_REP_C_upper[,"Te_lambda"]<-rep(c(obs.Cd_w0$lambda[3]+obs.Cd_w0$lambda_standard_error[3], obs.Cd_w0$lambda[4]+obs.Cd_w0$lambda_standard_error[4]), each=2)
  
  cxr_param_REP_C_upper[,"Tu_intra"]<-rep(c(obs.Cd_w0$alpha_matrix[1,1]+obs.Cd_w0$alpha_matrix_standard_error[1,1],obs.Cd_w0$alpha_matrix[2,2]+obs.Cd_w0$alpha_matrix_standard_error[2,2]), 2)
  cxr_param_REP_C_upper[,"Te_intra"]<-rep(c(obs.Cd_w0$alpha_matrix[3,3]+obs.Cd_w0$alpha_matrix_standard_error[3,3], obs.Cd_w0$alpha_matrix[4,4]+obs.Cd_w0$alpha_matrix_standard_error[4,4]), each=2)
  
  
  cxr_param_REP_C_upper[,"Tu_inter"]<-c(obs.Cd_w0$alpha_matrix[1,3]+obs.Cd_w0$alpha_matrix_standard_error[1,3], obs.Cd_w0$alpha_matrix[2,3]+obs.Cd_w0$alpha_matrix_standard_error[2,3],obs.Cd_w0$alpha_matrix[1,4]+obs.Cd_w0$alpha_matrix_standard_error[1,4], obs.Cd_w0$alpha_matrix[2,4]+obs.Cd_w0$alpha_matrix_standard_error[2,4])
  cxr_param_REP_C_upper[,"Te_inter"]<-c(obs.Cd_w0$alpha_matrix[3,1]+obs.Cd_w0$alpha_matrix_standard_error[3,1], obs.Cd_w0$alpha_matrix[3,2]+obs.Cd_w0$alpha_matrix_standard_error[3,2],obs.Cd_w0$alpha_matrix[4,1]+obs.Cd_w0$alpha_matrix_standard_error[4,1], obs.Cd_w0$alpha_matrix[4,2]+obs.Cd_w0$alpha_matrix_standard_error[4,2])
  
  #cxr_param_REP_C_lower
  #cxr_param_REP_C_upper

#' 
#' #### Joining the two data frames
## ---------------------------
  param_all_REP<-as.data.frame(rbind(cxr_param_REP, cxr_param_REP_C))
  
  param_all_REP_lower<-as.data.frame(rbind(cxr_param_REP_lower, cxr_param_REP_C_lower))
  param_all_REP_upper<-as.data.frame(rbind(cxr_param_REP_upper, cxr_param_REP_C_upper))
  

#' 
#' #### Saving pooled parameter estimates 
#' 
## ---------------------------
## Save the objects
save(obs.w0, file="./Analyses/cxr_N_allEqual_best.RData")
save(obs.Cd_w0, file="./Analyses/cxr_Cd_allEqual_best.RData")

#' 
#' Save the data sets to be used later
## ---------------------------
write.csv(param_all_REP, "./Analyses/cxr_normal_REP_allEqual/parameters_cxr_normal_REP_best.csv")
  write.csv(param_all_REP_upper, "./Analyses/cxr_normal_REP_allEqual/parameters_cxr_normal_REP_upper_best.csv")
  write.csv(param_all_REP_lower, "./Analyses/cxr_normal_REP_allEqual/parameters_cxr_normal_REP_lower_best.csv")


#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
