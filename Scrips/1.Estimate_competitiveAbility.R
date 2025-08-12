rm(list=ls())
library(plyr)
library(tidyverse)
library(car)
library(fitdistrplus)
library(tidyr)
library(cxr)
library(MASS)
library(mvtnorm)


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


# 2 - Estimate parameters from cxr package
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
                        bootstrap_samples = 1000)

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
                           bootstrap_samples = 1000)


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

# write.csv(param_all_REP, "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP.csv")
# write.csv(param_all_REP_upper, "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_upper.csv")
# write.csv(param_all_REP_lower, "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_lower.csv")


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
                           bootstrap_samples = 1000)

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
                           bootstrap_samples = 1000)

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
                           bootstrap_samples = 1000)

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
                           bootstrap_samples = 1000)

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
                          bootstrap_samples = 1000)

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
                          bootstrap_samples = 1000)

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
                                bootstrap_samples = 1000)

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
                          bootstrap_samples = 1000)

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
                                bootstrap_samples = 1000)


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
                              bootstrap_samples = 1000)

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
                              bootstrap_samples = 1000)

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
                              bootstrap_samples = 1000)

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

# write.csv(param_all_w0, "./Analyses/cxr_normal/parameters_cxr_normal.csv")
# write.csv(param_all_w0_upper, "./Analyses/cxr_normal/parameters_cxr_normal_upper.csv")
# write.csv(param_all_w0_lower, "./Analyses/cxr_normal/parameters_cxr_normal_lower.csv")
