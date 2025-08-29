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

### 1 - Estimating competitive ability data
# Importing data and checking

### Importing Competitive ability

ca<-read.csv(file = "./Data/CompetitiveAbility_Cd_G40_submit.csv", header=TRUE) # cdata from the competitive ability


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

# This is useful for plots
regimeTu<-c("Tu \ncontrol", "Tu evolved \n in cadmium")
names(regimeTu)<-c("SR1", "SR2")

regimeTe<-c("Te \n control", "Te evolved \n in cadmium")
names(regimeTe)<-c("SR4", "SR5")



#### Creating columns that are needed
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

ca$Nr_Focal_Females_G0<-sapply(c(1:length(ca$Block)), function(x){
  if(ca$Type[x]=="INTRA"){
    a<-ca$Dens[x]
  }else
    a<-1
  
})
str(ca)
summary(ca)

which(ca$Num_Comp_Te_Alive_G0<0)
which(ca$Num_Comp_Tu_Alive_G0<0)
which(ca$Nr_Focal_Females_Tu_Alive_G0<0)
which(ca$Nr_Focal_Females_Te_Alive_G0<0)

ca<-ca[-c(which(ca$Num_Comp_Te_Alive_G0<0),which(ca$Num_Comp_Tu_Alive_G0<0), which(ca$Nr_Focal_Females_Te_Alive_G0<0) ),]

#Creating the columns with the correct number of competitors. For conspecifics its always the same as the density to calculate the growth rate. For heterospecifics its always 1 female with X competitors, and DensFocal2 its also to do the same for the conspecifics 

ca$DensFocal<-sapply(c(1:dim(ca)[1]), function(x){
  if(ca$Type[x]=="INTRA"){
    a<-ca$Dens[x]
  }else if(ca$Type[x]=="INTER"){
    a<-1
  }
  
  a
})

ca$DensComp<-sapply(c(1:dim(ca)[1]), function(x){
  if(ca$Type[x]=="INTRA"){
    a<-0
  }else if(ca$Type[x]=="INTER"){
    a<-ca$Dens[x]-1
  }
  
  a
})

ca$DensFocal2<-sapply(c(1:dim(ca)[1]), function(x){
  if(ca$Type[x]=="INTRA"){
    a<-ca$Dens[x]-1
  }else if(ca$Type[x]=="INTER"){
    a<-1
  }
  
  a
})

ca$DensComp2<-sapply(c(1:dim(ca)[1]), function(x){
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

#Growth rate per day
ca$GrowthRatePD<-sapply(c(1:dim(ca)[1]), function(x){
  #print(x)
  if(ca$Focal_Female[x]=="Tu"){
    a<-(ca$TuFemales[x]/ca$DensFocal[x])/3
  }else if(ca$Focal_Female[x]=="Te"){
    a<-(ca$TeFemales[x]/ca$DensFocal[x])/3
  }else{
    a<-NA}
  
  
  a
})


# Estimate parameters from cxr package
# cxr accepts a data frame with a first column called fitness with positive values and numeric columns with number of individuals.
# Each row is one individual.
# For multiple species the easier is to create a list, each with a data frame that has in the first column number of individuals produced and then the number of neighbours
# 
# this case we transformed all 0s into 1 (so that the log is 0) For that we need to add +1 to all data so that the variance is not changed.
# 
# Note that the files of the output data are available in the folder.
# To avoid having the run the code again I added eval = FALSE.
# To use the already generated ouptut files you just need to import the data.
# To generate the data again, you need to change eval=TRUE.

eval<-FALSE

if(eval){
  
  # 2.1 - Estimate All replicates together
  
  ##### normal
  dir.create("./Analyses/cxr_normal_REP", showWarnings = FALSE)
  
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
  Rep<-list(SR1= subset(forCXR_N, Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N,  Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  
  obs.w0<-cxr_pm_multifit(data = Rep,
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
                          bootstrap_samples = 10000)
  
  obs.w0$lambda_standard_error
  obs.w0$lambda
  obs.w0$alpha_matrix_standard_error
  
  ###### rows in the alpha element of the returning list correspond to species i and columns to species j for each αij coefficient.
  
  ###### data table summary
  
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
  
  ### Lower
  
  cxr_param_REP_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Environment=c("N"))
  cxr_param_REP_lower$Tu_lambda<-0
  cxr_param_REP_lower$Te_lambda<-0
  cxr_param_REP_lower$Tu_intra<-0
  cxr_param_REP_lower$Te_intra<-0
  cxr_param_REP_lower$Tu_inter<-0
  cxr_param_REP_lower$Te_inter<-0
  
  cxr_param_REP_lower[,"Tu_lambda"]<-rep(c(obs.w0$lambda[1]-obs.w0$lambda_standard_error[1], obs.w0$lambda[2]-obs.w0$lambda_standard_error[2]), 2)
  cxr_param_REP_lower[,"Te_lambda"]<-rep(c(obs.w0$lambda[3]-obs.w0$lambda_standard_error[3], obs.w0$lambda[4]-obs.w0$lambda_standard_error[4]), each=2)
  
  cxr_param_REP_lower[,"Tu_intra"]<-rep(c(sign(obs.w0$alpha_matrix[1,1])*(abs(obs.w0$alpha_matrix[1,1])-obs.w0$alpha_matrix_standard_error[1,1]), sign(obs.w0$alpha_matrix[2,2])*(abs(obs.w0$alpha_matrix[2,2])-obs.w0$alpha_matrix_standard_error[2,2])), 2)
  cxr_param_REP_lower[,"Te_intra"]<-rep(c(sign(obs.w0$alpha_matrix[3,3])*(abs(obs.w0$alpha_matrix[3,3])-obs.w0$alpha_matrix_standard_error[3,3]), sign(obs.w0$alpha_matrix[4,4])*(abs(obs.w0$alpha_matrix[4,4])-obs.w0$alpha_matrix_standard_error[4,4])), each=2)
  
  
  cxr_param_REP_lower[,"Tu_inter"]<-c(sign(obs.w0$alpha_matrix[1,3])*(abs(obs.w0$alpha_matrix[1,3])-obs.w0$alpha_matrix_standard_error[1,3]), sign(obs.w0$alpha_matrix[2,3])*(abs(obs.w0$alpha_matrix[2,3])-obs.w0$alpha_matrix_standard_error[2,3]),sign(obs.w0$alpha_matrix[1,4])*(abs(obs.w0$alpha_matrix[1,4])-obs.w0$alpha_matrix_standard_error[1,4]), sign(obs.w0$alpha_matrix[2,4])*(abs(obs.w0$alpha_matrix[2,4])-obs.w0$alpha_matrix_standard_error[2,4]))
  cxr_param_REP_lower[,"Te_inter"]<-c(sign(obs.w0$alpha_matrix[3,1])*(abs(obs.w0$alpha_matrix[3,1])-obs.w0$alpha_matrix_standard_error[3,1]), sign(obs.w0$alpha_matrix[3,2])*(abs(obs.w0$alpha_matrix[3,2])-obs.w0$alpha_matrix_standard_error[3,2]),sign(obs.w0$alpha_matrix[4,1])*(abs(obs.w0$alpha_matrix[4,1])-obs.w0$alpha_matrix_standard_error[4,1]), sign(obs.w0$alpha_matrix[4,2])*(abs(obs.w0$alpha_matrix[4,2])-obs.w0$alpha_matrix_standard_error[4,2]))
  
  
  ### upper
  
  cxr_param_REP_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Environment=c("N"))
  cxr_param_REP_upper$Tu_lambda<-0
  cxr_param_REP_upper$Te_lambda<-0
  cxr_param_REP_upper$Tu_intra<-0
  cxr_param_REP_upper$Te_intra<-0
  cxr_param_REP_upper$Tu_inter<-0
  cxr_param_REP_upper$Te_inter<-0
  
  cxr_param_REP_upper[,"Tu_lambda"]<-rep(c(obs.w0$lambda[1]+obs.w0$lambda_standard_error[1], obs.w0$lambda[2]+obs.w0$lambda_standard_error[2]), 2)
  cxr_param_REP_upper[,"Te_lambda"]<-rep(c(obs.w0$lambda[3]+obs.w0$lambda_standard_error[3], obs.w0$lambda[4]+obs.w0$lambda_standard_error[4]), each=2)
  
  cxr_param_REP_upper[,"Tu_intra"]<-rep(c(sign(obs.w0$alpha_matrix[1,1])*(abs(obs.w0$alpha_matrix[1,1])+obs.w0$alpha_matrix_standard_error[1,1]), sign(obs.w0$alpha_matrix[2,2])*(abs(obs.w0$alpha_matrix[2,2])+obs.w0$alpha_matrix_standard_error[2,2])), 2)
  cxr_param_REP_upper[,"Te_intra"]<-rep(c(sign(obs.w0$alpha_matrix[3,3])*(abs(obs.w0$alpha_matrix[3,3])+obs.w0$alpha_matrix_standard_error[3,3]), sign(obs.w0$alpha_matrix[4,4])*(abs(obs.w0$alpha_matrix[4,4])+obs.w0$alpha_matrix_standard_error[4,4])), each=2)
  
  
  cxr_param_REP_upper[,"Tu_inter"]<-c(sign(obs.w0$alpha_matrix[1,3])*(abs(obs.w0$alpha_matrix[1,3])+obs.w0$alpha_matrix_standard_error[1,3]), sign(obs.w0$alpha_matrix[2,3])*(abs(obs.w0$alpha_matrix[2,3])+obs.w0$alpha_matrix_standard_error[2,3]),sign(obs.w0$alpha_matrix[1,4])*(abs(obs.w0$alpha_matrix[1,4])+obs.w0$alpha_matrix_standard_error[1,4]), sign(obs.w0$alpha_matrix[2,4])*(abs(obs.w0$alpha_matrix[2,4])+obs.w0$alpha_matrix_standard_error[2,4]))
  cxr_param_REP_upper[,"Te_inter"]<-c(sign(obs.w0$alpha_matrix[3,1])*(abs(obs.w0$alpha_matrix[3,1])+obs.w0$alpha_matrix_standard_error[3,1]), sign(obs.w0$alpha_matrix[3,2])*(abs(obs.w0$alpha_matrix[3,2])+obs.w0$alpha_matrix_standard_error[3,2]),sign(obs.w0$alpha_matrix[4,1])*(abs(obs.w0$alpha_matrix[4,1])+obs.w0$alpha_matrix_standard_error[4,1]), sign(obs.w0$alpha_matrix[4,2])*(abs(obs.w0$alpha_matrix[4,2])+obs.w0$alpha_matrix_standard_error[4,2]))
  
  
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
  
  #removing rows for which there is no data for fitness
  #forCXR_Cd<-forCXR_Cd[-which(is.na(forCXR_Cd$fitness)),]
  #forCXR_Cd$fitness<-forCXR_Cd$fitness+1
  
  forCXR_Cd[which(forCXR_Cd$fitness=="-Inf" | forCXR_Cd$fitness=="Inf"),"fitness"]<-0
  
  # 0 to 1 to maintain data hat is missing
  forCXR_Cd<-forCXR_Cd[-which(is.na(forCXR_Cd$fitness)),]
  forCXR_Cd$fitness<-forCXR_Cd$fitness+1
  
  
  
  # vector that tells which are the selection regimes, the columns have to have the same name
  my.reg <- c("SR1", "SR2","SR4","SR5")
  
  # Do list per replicate and environment
  Rep_Cd<-list(SR1= subset(forCXR_Cd, Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  
  obs.Cd_w0<-cxr_pm_multifit(data = Rep_Cd,
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
                             bootstrap_samples = 10000)
  
  
  ###### data table summary
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
  
  ### Lower
  
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
  
  
  ### upper
  
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
  
  cxr_param_REP_C_lower
  cxr_param_REP_C_upper
  
  
  ###### joining data frame
  param_all_REP<-as.data.frame(rbind(cxr_param_REP, cxr_param_REP_C))
  
  param_all_REP_lower<-as.data.frame(rbind(cxr_param_REP_lower, cxr_param_REP_C_lower))
  param_all_REP_upper<-as.data.frame(rbind(cxr_param_REP_upper, cxr_param_REP_C_upper))
  
  
  param_all_REP_lower
  param_all_REP_upper
  
  write.csv(param_all_REP, "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP.csv")
  write.csv(param_all_REP_upper, "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_upper.csv")
  write.csv(param_all_REP_lower, "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_lower.csv")
  
  
  # 2.2 - Estimate each replicate separately
  
  ##### normal
  
  
  dir.create("./Analyses/cxr_normal", showWarnings = FALSE)
  
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
                             bootstrap_samples = 10000)
  
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
                             bootstrap_samples = 10000)
  
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
                             bootstrap_samples = 10000)
  
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
                             bootstrap_samples = 10000)
  
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
                            bootstrap_samples = 10000)
  
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
                            bootstrap_samples = 10000)
  
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
                                  bootstrap_samples = 10000)
  
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
                            bootstrap_samples = 10000)
  
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
                                  bootstrap_samples = 10000)
  
  
  #ab<-abundance_projection(obs.R1_w0, timesteps = 1, initial_abundances = c(3,3,3,3))
  
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
                                bootstrap_samples = 10000)
  
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
                                bootstrap_samples = 10000)
  
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
                                bootstrap_samples = 10000)
  
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
                               bootstrap_samples = 10000)
  
  #for replicate 2 we will do the fitting by hand because we may need to scale the parameters
  
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
                               bootstrap_samples = 10000)
  
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
                               bootstrap_samples = 10000)
  
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
                                     bootstrap_samples = 10000)
  
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
                                     bootstrap_samples = 10000)
  
  
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
  
  
  cxr_param_w0C
  
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
  
  param_all_w0<-as.data.frame(rbind(cxr_param_w0, cxr_param_w0C))
  
  param_all_w0_lower<-as.data.frame(rbind(cxr_param_w0_lower, cxr_param_w0C_lower))
  param_all_w0_upper<-as.data.frame(rbind(cxr_param_w0_upper, cxr_param_w0C_upper))
  
  
  param_all_w0_lower
  param_all_w0_upper
  
  write.csv(param_all_w0, "./Analyses/cxr_normal/parameters_cxr_normal.csv")
  write.csv(param_all_w0_upper, "./Analyses/cxr_normal/parameters_cxr_normal_upper.csv")
  write.csv(param_all_w0_lower, "./Analyses/cxr_normal/parameters_cxr_normal_lower.csv")
  
}

# 2 -  Testing evolution

###### Importing data files for replicate estimation, importing files in case the cxr code was not run before, it will not import if eval=TRUE
if(!eval){
## Importing
param_all_w0<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal.csv")
param_all_w0_upper<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal_upper.csv")
param_all_w0_lower<-read.csv( "./Analyses/cxr_normal/parameters_cxr_normal_lower.csv")

param_all_w0<-param_all_w0[,-1]
param_all_w0_upper<-param_all_w0_upper[,-1]
param_all_w0_lower<-param_all_w0_lower[,-1]

}
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
# This takes a lot of time!
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


# 3 - Estimating parameteres from structural stability
# Importing data in case the previous code was not run

if(!eval){
  ## Importing pooled data
  param_all_REP<-read.csv("./Analyses/cxr_normal_REP/parameters_cxr_normal_REP.csv")
  param_all_REP_upper<-read.csv("./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_upper.csv")
  param_all_REP_lower<-read.csv( "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_lower.csv")
  param_all_REP<-param_all_REP[,-1]
  param_all_REP_upper<-param_all_REP_upper[,-1]
  param_all_REP_lower<-param_all_REP_lower[,-1]
  
  ## Importing per replicate
  param_all_w0<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal.csv")
  param_all_w0_upper<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal_upper.csv")
  param_all_w0_lower<-read.csv( "./Analyses/cxr_normal/parameters_cxr_normal_lower.csv")
  
  param_all_w0<-param_all_w0[,-1]
  param_all_w0_upper<-param_all_w0_upper[,-1]
  param_all_w0_lower<-param_all_w0_lower[,-1]
  
}

# Predicting coexistence (structural stability)

### Structural stability

### Defining functions

#input parameters:
#alpha = competition strength matrix 
#r = vector of intrinsic growth rates

#structural niche difference (output on a log scale)
Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha, tol = 1e-40)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}

#vector defining the centroid of the feasibility domain
r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}


#structural fitness difference (in degree)
theta <- function(alpha,r){
  r_c <- r_centroid(alpha)
  out <- acos(sum(r_c*r, na.rm = TRUE)/(sqrt(sum(r^2, na.rm = TRUE))*sqrt(sum(r_c^2, na.rm = TRUE))))*90/pi
  return(out)
}


#test if a system (alpha and r) is feasible (output 1 = feasible, 0 = not feasible)
test_feasibility <- function(alpha,r){
  out <- prod(solve(alpha,r)>0)
  return(out)
}



### 3.1.1 - estimating structural stability for pooled replicates
#x<-2
struct_mat_REP<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_REP[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_REP$Te_intra[x],param_all_REP$Te_inter[x],param_all_REP$Tu_inter[x], param_all_REP$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_REP$Te_lambda[x],param_all_REP$Tu_lambda[x])
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))


colnames(struct_mat_REP)<-c("ND", "FD", "Feasibility")

#For the lower we use the higher alphas with lower lambda, and for upper the other way around
# Since we have facilitation we have to actually test what is the lowest value

struct_mat_REP_U<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_REP_lower[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_REP_upper$Te_intra[x],param_all_REP_upper$Te_inter[x], param_all_REP_upper$Tu_inter[x],param_all_REP_upper$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_REP$Te_lambda[x],param_all_REP$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))


struct_mat_REP_L<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_REP_upper[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_REP_lower$Te_intra[x],param_all_REP_lower$Te_inter[x],param_all_REP_lower$Tu_inter[x],param_all_REP_lower$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_REP$Te_lambda[x],param_all_REP$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

colnames(struct_mat_REP_U)<-c("ND_U", "FD_U", "Feasibility_U")
colnames(struct_mat_REP_L)<-c("ND_L", "FD_L", "Feasibility_L")

struct_mat_REP<-cbind(param_all_REP, struct_mat_REP,struct_mat_REP_L,struct_mat_REP_U)

# To create the boundaries
bound_struct_rk_w0<-data.frame(ND=seq(0,1, 0.01))
bound_struct_rk_w0$FD<-45*bound_struct_rk_w0$ND

struct_mat_REP3<-struct_mat_REP

### 3.1.2 - estimating structural stability per replicate
struct_mat_w0<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_w0[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_w0$Te_intra[x],param_all_w0$Te_inter[x],param_all_w0$Tu_inter[x], param_all_w0$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_w0$Te_lambda[x],param_all_w0$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

colnames(struct_mat_w0)<-c("ND", "FD", "Feasibility")

struct_mat_w0_L<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_w0_lower[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_w0_upper$Te_intra[x],param_all_w0_upper$Te_inter[x],param_all_w0_upper$Tu_inter[x], param_all_w0_upper$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_w0$Te_lambda[x],param_all_w0$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

colnames(struct_mat_w0_L)<-c("ND_L", "FD_L", "Feasibility_L")

struct_mat_w0_U<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_w0_upper[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_w0_lower$Te_intra[x],param_all_w0_lower$Te_inter[x],param_all_w0_lower$Tu_inter[x], param_all_w0_lower$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_w0$Te_lambda[x],param_all_w0$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

colnames(struct_mat_w0_U)<-c("ND_U", "FD_U", "Feasibility_U")

struct_mat_w0<-cbind(param_all_w0, struct_mat_w0,struct_mat_w0_L,struct_mat_w0_U)

bound_struct_rk_w0<-data.frame(ND=seq(0,1, 0.01))
bound_struct_rk_w0$FD<-45*bound_struct_rk_w0$ND


### Distance to the edge and who wins

struct_mat_REP$a21_a11<-struct_mat_REP$Te_inter/struct_mat_REP$Tu_intra
struct_mat_REP$a22_a12<-struct_mat_REP$Te_intra/struct_mat_REP$Tu_inter

## Now to calculate the upper and lower bounds we have to see which is the lowest value of competition (and those create the upper bounds) or the highest values of competition (those create the lower boundaries)

struct_mat_REP$a21_a11_upper<-param_all_REP_lower$Te_inter/param_all_REP_lower$Tu_intra
struct_mat_REP$a22_a12_upper<-param_all_REP_lower$Te_intra/param_all_REP_lower$Tu_inter

struct_mat_REP$a21_a11_lower<-param_all_REP_upper$Te_inter/param_all_REP_upper$Tu_intra
struct_mat_REP$a22_a12_lower<-param_all_REP_upper$Te_intra/param_all_REP_upper$Tu_inter

struct_mat_REP$Tu_lambda_lower<-param_all_REP_lower$Tu_lambda
struct_mat_REP$Te_lambda_lower<-param_all_REP_lower$Te_lambda
struct_mat_REP$Tu_lambda_upper<-param_all_REP_upper$Tu_lambda
struct_mat_REP$Te_lambda_upper<-param_all_REP_upper$Te_lambda

#### order, because of facilitation 
struct_mat_REP$min_a21_a11<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  min(c(abs(struct_mat_REP$a21_a11_lower[x]),abs(struct_mat_REP$a21_a11_upper[x])))})

struct_mat_REP$min_a22_a12<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  min(c(abs(struct_mat_REP$a22_a12_lower[x]),abs(struct_mat_REP$a22_a12_upper[x])))})

struct_mat_REP$max_a21_a11<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  max(c(abs(struct_mat_REP$a21_a11_lower[x]),abs(struct_mat_REP$a21_a11_upper[x])))})

struct_mat_REP$max_a22_a12<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  max(c(abs(struct_mat_REP$a22_a12_lower[x]),abs(struct_mat_REP$a22_a12_upper[x])))})

### per replicate
struct_mat_w0$a21_a11<-struct_mat_w0$Te_inter/struct_mat_w0$Tu_intra
struct_mat_w0$a22_a12<-struct_mat_w0$Te_intra/struct_mat_w0$Tu_inter

struct_mat_w0$a21_a11_lower<-param_all_w0_lower$Te_inter/param_all_w0_lower$Tu_intra
struct_mat_w0$a22_a12_lower<-param_all_w0_lower$Te_intra/param_all_w0_lower$Tu_inter

struct_mat_w0$a21_a11_upper<-param_all_w0_upper$Te_inter/param_all_w0_upper$Tu_intra
struct_mat_w0$a22_a12_upper<-param_all_w0_upper$Te_intra/param_all_w0_upper$Tu_inter

if(eval){
  write.csv(struct_mat_REP, "./Analyses/structural_REP.csv")
  write.csv(struct_mat_w0, "./Analyses/structural_REP_w0.csv")
}

# 4 - Coexistence

if(!eval){ #Importing data if the code above has not been run yet
# Importing data
param_all_REP<-read.csv("./Analyses/cxr_normal_REP/parameters_cxr_normal_REP.csv")
param_all_REP_upper<-read.csv("./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_upper.csv")
param_all_REP_lower<-read.csv( "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_lower.csv")
param_all_REP<-param_all_REP[,-1]
param_all_REP_upper<-param_all_REP_upper[,-1]
param_all_REP_lower<-param_all_REP_lower[,-1]

param_all_w0<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal.csv")
param_all_w0_upper<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal_upper.csv")
param_all_w0_lower<-read.csv( "./Analyses/cxr_normal/parameters_cxr_normal_lower.csv")

param_all_w0<-param_all_w0[,-1]
param_all_w0_upper<-param_all_w0_upper[,-1]
param_all_w0_lower<-param_all_w0_lower[,-1]
}

coex_g42<-read.csv("./Data/Coexistence_Cd_G42_submit.csv", header=TRUE) # Data from the coexistence experiment

coex_g42$Rep2<-as.factor(coex_g42$Rep)
coex_g42$X1st.pair<-as.factor(coex_g42$X1st_pair)
coex_g42$X2nd.pair<-as.factor(coex_g42$X2nd_pair)
coex_g42$SRTu<-as.factor(coex_g42$SRTu)
coex_g42$SRTe<-as.factor(coex_g42$SRTe)
coex_g42$Box2<-as.factor(coex_g42$Box)

### summary data per leaf (because the leaflets are not attributable)
coex_g42_res<-gather(coex_g42, leaf, females, Leaf_2_Up_Tu:Leaf_5_Down_Te, factor_key=TRUE)
str(coex_g42_res)

coex_g42_res$char<-as.character(coex_g42_res$leaf)

aux4<-as.data.frame(t(as.data.frame(sapply(c(1:length(coex_g42_res$Rep)), function(x){
  a<-strsplit(coex_g42_res$char[x], split="_")[[1]]
  
  c(a[2:4])
}))))

colnames(aux4)<-c("Leaf2", "Direction", "Species")

coex_g42_res<-cbind(coex_g42_res, aux4)

str(coex_g42_res)

sum_coex_g42<-coex_g42_res %>%
  group_by(Rep2, SRTu, SRTe, Box2, Leaf2,X1st.pair,X2nd.pair, Direction, Species, Env) %>%
  summarize(av_females=sum(females, na.rm=TRUE))

sum_coex_g42
sum_coex_g42$Direction<-as.factor(sum_coex_g42$Direction)

sum_coex_g42_res<-as.data.frame(spread(sum_coex_g42, key=Species, value=av_females))

sum_coex_g42_res2<-sum_coex_g42_res %>%
  group_by(Rep2, SRTu, SRTe, Box2, Leaf2,X1st.pair,X2nd.pair, Env) %>%
  summarize(av_Te=sum(Te, na.rm=TRUE), av_Tu=sum(Tu, na.rm=TRUE)) %>% as.data.frame()

str(sum_coex_g42_res2)

coex_g42_rep<-sum_coex_g42_res2 %>%
  group_by(Rep2, Leaf2, SRTu, SRTe, Env, Box2) %>%
  summarize( sum_Te=sum(av_Te, na.rm=TRUE), sum_Tu=mean(av_Tu, na.rm=TRUE)) %>% as.data.frame()


coex_g42_rep$Env[which(coex_g42_rep$Env=="Cd")]<-"Cd" 


#Can we predict the outcome of species interactions?
str(sum_coex_g42_res2)

coex_no_het2<-sum_coex_g42_res2 %>%
  group_by(SRTu, SRTe, Rep2, Box2, Env) %>%
  summarize(sumTe=sum(av_Te, na.rm=TRUE), sumTu=sum(av_Tu, na.rm=TRUE)) %>% as.data.frame()

coex_no_het2$Te_ratio<-coex_no_het2$sumTe/(coex_no_het2$sumTe+coex_no_het2$sumTu)

coex_no_het2$Te_ratio[which(coex_no_het2$Te_ratio=="NaN")]<-0

coex_no_het2$Tu_Regime<-mapvalues(coex_no_het2$SRTu, c("Tu2", "Tu1"), c("SR2", "SR1"))

coex_no_het2$Te_Regime<-mapvalues(coex_no_het2$SRTe, c("Te4", "Te5"), c("SR4", "SR5"))

#### Predicting once

pred_coex_RK_REP<-expand_grid(Te=c("SR4","SR5"), Tu=c("SR1", "SR2"), Environment= c("N", "Cd"))

pred_coex_RK_REP$predTu1<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_alphas$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6 - aux_alphas$Tu_inter[1]*6)
  
  bl
})

pred_coex_RK_REP$predTe1<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})


pred_coex_RK_REP$predTu2<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_alphas$Tu_lambda[1]*pred_coex_RK_REP$predTu1[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_REP$predTu1[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_REP$predTe1[x])
  
  bl
})

pred_coex_RK_REP$predTe2<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*pred_coex_RK_REP$predTe1[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_REP$predTe1[x] - aux_alphas$Te_inter[1]*pred_coex_RK_REP$predTu1[x])
  
  bl
})

#x<-1
# lower - stronger alpha and lower lambda
pred_coex_RK_REP$predTu1_L<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6- aux_alphas$Tu_inter[1]*6)
  
  bl
})

pred_coex_RK_REP$predTe1_L<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})


pred_coex_RK_REP$predTu2_L<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Tu_lambda[1]*pred_coex_RK_REP$predTu1_L[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_REP$predTu1_L[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_REP$predTe1_L[x])
  
  bl
})

pred_coex_RK_REP$predTe2_L<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Te_lambda[1]*pred_coex_RK_REP$predTe1_L[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_REP$predTe1_L[x] - aux_alphas$Te_inter[1]*pred_coex_RK_REP$predTu1_L[x])
  
  bl
})

# upper
pred_coex_RK_REP$predTu1_U<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6- aux_alphas$Tu_inter[1]*6)
  
  bl
})

pred_coex_RK_REP$predTe1_U<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})


pred_coex_RK_REP$predTu2_U<-sapply(c(1:length(pred_coex_RK_REP$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Tu_lambda[1]* pred_coex_RK_REP$predTu1_U[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_REP$predTu1_U[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_REP$predTe1_U[x])
  
  bl
})

pred_coex_RK_REP$predTe2_U<-sapply(c(1:length(pred_coex_RK_REP$Te)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  aux_lambda<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex_RK_REP$Tu[x]) & Te_Regime==as.character(pred_coex_RK_REP$Te[x]) & Environment==as.character(pred_coex_RK_REP$Environment[x]))
  
  bl<-aux_lambda$Te_lambda[1]*pred_coex_RK_REP$predTe1_U[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_REP$predTe1_U[x] - aux_alphas$Te_inter[1]*pred_coex_RK_REP$predTu1_U[x])
  
  bl
})


names(pred_coex_RK_REP)[1:3]<-c("SRTe", "SRTu", "Env")

pred_coex_RK_REP<-as.data.frame(pred_coex_RK_REP)

## Predicting per replicate
pred_coex_RK_w0<-as.data.frame(expand_grid(Te=c("SR4","SR5"), Tu=c("SR1", "SR2"), Environment= c("N", "Cd"), Replicate=c(1,2,3,4,5)))

pred_coex_RK_w0<- pred_coex_RK_w0[- which(pred_coex_RK_w0$Replicate==2 & pred_coex_RK_w0$Tu=="SR2" & pred_coex_RK_w0$Environment=="Cd"),]


pred_coex_RK_w0$predTu1<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_alphas$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6 - aux_alphas$Tu_inter[1]*6)
  
  bl
})

pred_coex_RK_w0$predTe1<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_alphas$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})


pred_coex_RK_w0$predTu2<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_alphas$Tu_lambda[1]*pred_coex_RK_w0$predTu1[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_w0$predTu1[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_w0$predTe1[x])
  
  bl
})

pred_coex_RK_w0$predTe2<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_alphas$Te_lambda[1]*pred_coex_RK_w0$predTe1[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_w0$predTe1[x] - aux_alphas$Te_inter[1]*pred_coex_RK_w0$predTu1[x])
  
  bl
})

#x<-1
# lower - stronger alpha and lower lambda
pred_coex_RK_w0$predTu1_L<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6- aux_alphas$Tu_inter[1]*6)
  
  bl
})

pred_coex_RK_w0$predTe1_L<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})


pred_coex_RK_w0$predTu2_L<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Tu_lambda[1]*pred_coex_RK_w0$predTu1_L[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_w0$predTu1_L[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_w0$predTe1_L[x])
  
  bl
})

pred_coex_RK_w0$predTe2_L<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Te_lambda[1]*pred_coex_RK_w0$predTe1_L[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_w0$predTe1_L[x] - aux_alphas$Te_inter[1]*pred_coex_RK_w0$predTu1_L[x])
  
  bl
})

# upper
pred_coex_RK_w0$predTu1_U<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Tu_lambda[1]*6* exp(-aux_alphas$Tu_intra[1]*6- aux_alphas$Tu_inter[1]*6)
  
  bl
})

pred_coex_RK_w0$predTe1_U<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Te_lambda[1]*6* exp(-aux_alphas$Te_intra[1]*6 - aux_alphas$Te_inter[1]*6)
  
  bl
})


pred_coex_RK_w0$predTu2_U<-sapply(c(1:length(pred_coex_RK_w0$Tu)), function(x){
  aux_alphas<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Tu_lambda[1]* pred_coex_RK_w0$predTu1_U[x]* exp(-aux_alphas$Tu_intra[1]*pred_coex_RK_w0$predTu1_U[x]- aux_alphas$Tu_inter[1]*pred_coex_RK_w0$predTe1_U[x])
  
  bl
})

pred_coex_RK_w0$predTe2_U<-sapply(c(1:length(pred_coex_RK_w0$Te)), function(x){
  aux_alphas<-subset(param_all_w0_lower, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  aux_lambda<-subset(param_all_w0_upper, Tu_Regime==as.character(pred_coex_RK_w0$Tu[x]) & Te_Regime==as.character(pred_coex_RK_w0$Te[x]) & Environment==as.character(pred_coex_RK_w0$Environment[x]) & Replicate==pred_coex_RK_w0$Replicate[x])
  
  bl<-aux_lambda$Te_lambda[1]*pred_coex_RK_w0$predTe1_U[x]* exp(-aux_alphas$Te_intra[1]*pred_coex_RK_w0$predTe1_U[x] - aux_alphas$Te_inter[1]*pred_coex_RK_w0$predTu1_U[x])
  
  bl
})


names(pred_coex_RK_w0)[1:3]<-c("SRTe", "SRTu", "Env")

pred_coex_RK_w0<-as.data.frame(pred_coex_RK_w0)

#### Empirical

coex_g42_rep$Te_ratio<-sapply(c(1:dim(coex_g42_rep)[1]), function(x) coex_g42_rep$sum_Te[x]/sum(coex_g42_rep$sum_Tu[x],coex_g42_rep$sum_Te[x]))

#All the NAN were caused by division by 0, so we put it as 0
coex_g42_rep$Te_ratio[which(coex_g42_rep$Te_ratio=="NaN")]<-0

coex_g42_rep2<-coex_g42_rep %>%
  group_by(SRTu, SRTe, Rep2, Box2, Env) %>%
  summarize(sumTe=sum(sum_Te, na.rm=TRUE), sumTu=sum(sum_Tu, na.rm=TRUE))

coex_g42_rep2$Te_ratio<-sapply(c(1:dim(coex_g42_rep2)[1]), function(x) coex_g42_rep2$sumTe[x]/sum(coex_g42_rep2$sumTu[x],coex_g42_rep2$sumTe[x]))

coex_g42_rep3<-coex_g42_rep2 %>%
  group_by(SRTu, SRTe, Rep2, Env) %>%
  summarize(sum_Te=sum(sumTe, na.rm=TRUE), sum_Tu=sum(sumTu, na.rm=TRUE), sdTe=sd(sumTe, na.rm=TRUE), sdTu=sd(sumTu, na.rm=TRUE), meanTeRatio=mean(Te_ratio, na.rm=TRUE))

coex_g42_rep3$Te_ratio<-sapply(c(1:dim(coex_g42_rep3)[1]), function(x) coex_g42_rep3$sum_Te[x]/sum(coex_g42_rep3$sum_Tu[x],coex_g42_rep3$sum_Te[x]))

### summarizing data


str(coex_g42_rep3)

coex_g42_rep3$Env<-plyr::mapvalues(coex_g42_rep3$Env, c("Cd","water"), c("Cd", "N"))
coex_g42_rep3$SRTu2<-plyr::mapvalues(coex_g42_rep3$SRTu, c("Tu1","Tu2"), c("SR1", "SR2"))
coex_g42_rep3$SRTe2<-plyr::mapvalues(coex_g42_rep3$SRTe, c("Te4","Te5"), c("SR4", "SR5"))

sum_observed_coex<-coex_g42_rep3
str(sum_observed_coex)

#write.csv(sum_observed_coex, "./TableS3.csv")


### Comparing to data

sum_observed_coex2<-sum_observed_coex %>%
  group_by(SRTu2, SRTe2, Env)%>%
  summarise(sumTe=sum(sum_Te, na.rm=TRUE),sumTu=mean(sum_Tu, na.rm=TRUE), sdTe2=sd(sum_Te, na.rm=TRUE)/sqrt(5), sdTu2=sd(sum_Tu, na.rm=TRUE)/sqrt(5), TeRatio=mean(meanTeRatio, na.rm=TRUE), sdTeRatio2=sd(meanTeRatio, na.rm=TRUE)/sqrt(5)) %>% as.data.frame()

sum_observed_coex2$TeRatio_L<-sum_observed_coex2$TeRatio-sum_observed_coex2$sdTeRatio2
sum_observed_coex2$TeRatio_U<-sum_observed_coex2$TeRatio+sum_observed_coex2$sdTeRatio2


## Testing predictions proportions

str(pred_coex_RK_REP)

pred_coex_RK_REP$TeRatio<-sapply(c(1:dim(pred_coex_RK_REP)[1]), function(x){
  pred_coex_RK_REP$predTe2[x]/(pred_coex_RK_REP$predTe2[x]+pred_coex_RK_REP$predTu2[x])
})

pred_coex_RK_w0$TeRatio<-sapply(c(1:dim(pred_coex_RK_w0)[1]), function(x){
  pred_coex_RK_w0$predTe2[x]/(pred_coex_RK_w0$predTe2[x]+pred_coex_RK_w0$predTu2[x])
})

pred_coex_RK_REP$TeRatio_L<-sapply(c(1:dim(pred_coex_RK_REP)[1]), function(x){
  pred_coex_RK_REP$predTe2_L[x]/(pred_coex_RK_REP$predTe2_L[x]+pred_coex_RK_REP$predTu2_L[x])
})

pred_coex_RK_w0$TeRatio_L<-sapply(c(1:dim(pred_coex_RK_w0)[1]), function(x){
  pred_coex_RK_w0$predTe2_L[x]/(pred_coex_RK_w0$predTe2_L[x]+pred_coex_RK_w0$predTu2_L[x])
})

pred_coex_RK_REP$TeRatio_U<-sapply(c(1:dim(pred_coex_RK_REP)[1]), function(x){
  pred_coex_RK_REP$predTe2_U[x]/(pred_coex_RK_REP$predTe2_U[x]+pred_coex_RK_REP$predTu2_U[x])
})

pred_coex_RK_w0$TeRatio_U<-sapply(c(1:dim(pred_coex_RK_w0)[1]), function(x){
  pred_coex_RK_w0$predTe2_U[x]/(pred_coex_RK_w0$predTe2_U[x]+pred_coex_RK_w0$predTu2_U[x])
})


str(pred_coex_RK_REP)

sum_pred_coex_RK_REP<-pred_coex_RK_REP %>%
  group_by(SRTu, SRTe, Env)%>%
  summarise(predTe=mean(predTe2, na.rm=TRUE),predTu=mean(predTu2, na.rm=TRUE), sumTeRatio=(sum(predTe2, na.rm=TRUE)/(sum(predTe2, na.rm=TRUE)+sum(predTu2, na.rm=TRUE)))) %>% as.data.frame()

## Comparing to proportions

str(sum_observed_coex)
sum_observed_coex_rep2<-sum_observed_coex %>%
  group_by(SRTe, SRTu, Env) %>%
  summarize(obs_TeRatio=mean(meanTeRatio), SE_obs=sd(meanTeRatio)/sqrt(5), meanTe=mean(sum_Te, na.rm=TRUE), meanTu=mean(sum_Tu, na.rm=TRUE)) %>% as.data.frame()


sum_observed_coex_rep2$SRTe2<-(plyr::mapvalues(as.character(sum_observed_coex_rep2$SRTe), c("Te4","Te5"), c("SR4", "SR5")))
sum_observed_coex_rep2$SRTu2<-(plyr::mapvalues(as.character(sum_observed_coex_rep2$SRTu), c("Tu1","Tu2"), c("SR1", "SR2")))
colnames(sum_observed_coex_rep2)[c(1:2, 8,9)]<-c("SRTe2", "SRTu2","SRTe", "SRTu" )

sum_observed_coex_rep2<-sum_observed_coex_rep2[,c(8,9,3:7)]
colnames(sum_observed_coex)[c(1,2,3,11:12)]<-c("SRTu2","SRTe2","Replicate","SRTu","SRTe")

sum_observed_coex_rep<-inner_join(sum_observed_coex_rep2, pred_coex_RK_REP, by=c("SRTe", "SRTu", "Env"))

sum_observed_coex_rep<-as.data.frame(sum_observed_coex_rep[,c("SRTe", "SRTu", "Env", "obs_TeRatio","SE_obs", "TeRatio", "TeRatio_L", "TeRatio_U", "meanTe","meanTu", "predTu2","predTe2")])

str(sum_observed_coex_rep)
colnames(sum_observed_coex_rep)[6:8]<-c("pred_T1", "T1_L","T1_U")

pred_coex_RK_w0$Replicate<-as.factor(pred_coex_RK_w0$Replicate)

sum_observed_coex_ALL<-inner_join(sum_observed_coex, pred_coex_RK_w0, by=c("SRTe", "SRTu", "Env", "Replicate"))

sum2_observed_coex_ALL<-as.data.frame(sum_observed_coex_ALL[,c(11,12,3,4,5,6,9,15,16)])

str(sum2_observed_coex_ALL)

sum2_observed_coex_ALL$predTeRatio<-sapply(c(1:dim(sum2_observed_coex_ALL)[1]), function(x){sum2_observed_coex_ALL$predTe2[x]/sum(sum2_observed_coex_ALL$predTe2[x],sum2_observed_coex_ALL$predTu2[x])})

#### Testing pooled replicates

m3<- glm(cbind(meanTe, meanTu)~pred_T1, data=sum_observed_coex_rep, family="binomial")

m4<- glm(obs_TeRatio~pred_T1, data=sum_observed_coex_rep, family="binomial")

summary(m3)
summary(m4)
#small differences in the estimates, but the p-value changes dramatically 

#This forces the line to pass by the 0,0, as the reviewer suggested
m5<-glm(cbind(meanTe, meanTu)~0+pred_T1, data=sum_observed_coex_rep, family="binomial")
summary(m5)
emtrends(m5, var="pred_T1", type="response")





