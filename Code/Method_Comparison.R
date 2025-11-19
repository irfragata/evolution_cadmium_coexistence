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


# Packages and functions----


if(!require("ggplot2")){
  install.packages("ggplot2")
}
if(!require("plyr")){
  install.packages("plyr")
}
if(!require("dplyr")){
  install.packages("dplyr")
}
if(!require("car")){
  install.packages("car")
}
if(!require("tidyr")){
  install.packages("tidyr")
}
if(!require("cxr")){
  install.packages("cxr")
}

library(ggplot2)
library(plyr)
library(dplyr)
library(car)
library(tidyr)
library(cxr)


theme_plots<-theme(axis.text = element_text(size=14), axis.title = element_text(size=14, face="bold"), legend.text = element_text(size=12), strip.text = element_text(size=14), plot.title = element_text(size=14, face="bold"), panel.grid=element_line(colour="white"), panel.background = element_rect(fill="white") , axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"), strip.background = element_rect(fill="white"))

# functions
save_plot<-function(dir, width=15, height=10, ...){
  ggsave(dir, width = width, height = height, units = c("cm"))
}

rk_func<- function(lambda, alpha_ii, alpha_ij, dens_i, dens_j, ...){
  gr<-lambda*dens_i*exp(-alpha_ii*dens_i - alpha_ij*dens_j)
  
  return(gr)
}

# Importing functions for running the methods D and E
source("./Code/Function_riker.R")


## Naming vectors to use for facets
regimeTu<-c("Tu \ncontrol", "Tu evolved \n in cadmium")
names(regimeTu)<-c("SR1", "SR2")

regimeTe<-c("Te \n control", "Te evolved \n in cadmium")
names(regimeTe)<-c("SR4", "SR5")

Env<-c("Water", "Cadmium")
names(Env)<-c("N", "Cd")

# setting the seed and the number of bootstrap samples
#set.seed(42)
set.seed(123)
Nboot <- 2


#' # Evaluation
#' Some pieces of the code take a lot of time to run. We provide the intermediate results to speed up the process. To run the whole code change eval from FALSE to TRUE.
#' 
# Evaluation------------------------------------
eval<-FALSE
if(!dir.exists("./Analyses")){
  dir.create("./Analyses/")
}

if(!dir.exists("./Plots")){
  dir.create("./Plots/")
}

if(!dir.exists("./Analyses/MethodComparison")){
  dir.create("./Analyses/MethodComparison/")
}

# Importing data ---- 
# Importing data and checking it 
ca_raw<-read.csv(file = "./Data/CompetitiveAbility_Cd_G40_submit.csv", header=TRUE) # cdata from the competitive ability

str(ca_raw) 
# Summary of the data to be sure that everything is ok!
summary(as.factor(ca_raw$FocalSR))

# Transforming character values into factors
ca_raw$Block2<-as.factor(ca_raw$Block)
ca_raw$Rep2<-as.factor(ca_raw$Rep)
ca_raw$Disk2<-as.factor(ca_raw$Disk)
ca_raw$Leaf2<-as.factor(ca_raw$Leaf)
ca_raw$Env2<-as.factor(ca_raw$Env)
ca_raw$FocalSR2<-as.factor(ca_raw$FocalSR)
ca_raw$CompSR2<-as.factor(ca_raw$CompSR)
ca_raw$Type2<-as.factor(ca_raw$Type)
ca_raw$Focal_Female2<-as.factor(ca_raw$Focalfemale)

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

# 2 - Estimate growth rate per generation

# loop to estimate the growth rate
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

str(ca)

# Estimate competitive ability & predict data -----

#Here we have two different methods, using the cxr package or with the optim. We will also vary if we estimate lambda from the data or from the model and if using cxr with the nested approach is better or not. So the different hypothesis are

#A - CXR normal: using cxr with the normal approach
#B - CXR lambda fixed: using cxr but lambda comes from the data
#C - CXR nested: lambda comes the data, and we use the same nested approach as the optim --> for that we can put intra as another species (column)
#D - optim normal: the same approach as used in Fragata 2022
#E - optim lambda fixed: using optim, but lambda is fixed 

#In all the models we will use density -1 for the intra, which basically corresponds to the number of competitors.

## A - CXR normal ------------------------------------------------------------

#cxr accepts a data frame with a first column called fitness with positive values and numeric columns with number of individuals. Each row is one individual. For multiple species the easier is to create a list, each with a data frame that has in the first column number of individuals produced and then the number of neighbours
#this case we transformed all 0s into 1 (so that the log is 0) For that we need to add +1 to all data so that the variance is not changed

print("Running Method A: CXR")
if(!dir.exists("./Analyses/MethodComparison/cxr_normal")){
  dir.create("./Analyses/MethodComparison/cxr_normal", showWarnings = FALSE)
}

## Create data frames ----

### No cadmium --------------------
# modifying data frame to fit the type of setup that is need for CXR
forCXR_N<-subset(ca, Env=="N")[,c("Rep", "FocalSR", "CompSR", "Dens", "TeFemales", "TuFemales")]

forCXR_N$Focal<-mapvalues(forCXR_N$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))
forCXR_N$CompSR2<-mapvalues(forCXR_N$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

# creating column to store the number of competitors
forCXR_N$Comp<-sapply(c(1:length(forCXR_N[,1])), function(x){
  if(is.na(forCXR_N$CompSR2[x])){
    a<- forCXR_N$Focal[x]
  }else{
    a<-forCXR_N$CompSR2[x]
  }
  
  a
})

# data frame to store the number of conspecific and heterospecific 
aux<-data.frame(SR1=rep(0, length(forCXR_N[,1])), SR2=rep(0, length(forCXR_N[,1])), SR4=rep(0, length(forCXR_N[,1])), SR5=rep(0, length(forCXR_N[,1])))

for(i in 1:length(forCXR_N[,1])){
  #column where to put focals
  colunaF<-which(colnames(aux)==forCXR_N$Focal[i])
  #column where to put competitors
  colunaC<-which(colnames(aux)==forCXR_N$Comp[i])
  
  #if its the same regime
  if(forCXR_N$Focal[i]==forCXR_N$Comp[i] & forCXR_N$Dens[i]==1){
    aux[i,colunaF]<-forCXR_N$Dens[i]-1
    
  }else if(forCXR_N$Focal[i]==forCXR_N$Comp[i]){
    aux[i,colunaF]<-forCXR_N$Dens[i]-1
  }else{ #if it is heterospecific then its -1 for the competitors (because of the focal) and its 0 for the focal
    aux[i,colunaC]<-forCXR_N$Dens[i]-1
    aux[i, colunaF]<-0
  }
  
}

#joining data frames
forCXR_N<-cbind(forCXR_N, aux)

# storing the fitness data 
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

# removing rows without data
forCXR_N[which(forCXR_N$fitness=="-Inf" | forCXR_N$fitness=="Inf"),"fitness"]<-0


# Adding +1 to all growth rate data to include 0s in the data
forCXR_N$fitness<-forCXR_N$fitness+1

# vector that tells which are the selection regimes, the columns have to have the same name
my.reg <- c("SR1", "SR2","SR4","SR5")

# Do list per replicate and environment
R1<-list(SR1= subset(forCXR_N, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R2<-list(SR1= subset(forCXR_N, Rep==2 & Focal=="SR1")[,c("fitness", "SR1", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR4", "SR5")])

R3<-list(SR1= subset(forCXR_N, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R4<-list(SR1= subset(forCXR_N, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R5<-list(SR1= subset(forCXR_N, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])


### Cadmium --------------------
##### Running cxr for the Cadmium environment
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

# data frame to store the number of conspecific and heterospecific 
aux<-data.frame(SR1=rep(0, length(forCXR_Cd[,1])), SR2=rep(0, length(forCXR_Cd[,1])), SR4=rep(0, length(forCXR_Cd[,1])), SR5=rep(0, length(forCXR_Cd[,1])))

for(i in 1:length(forCXR_Cd[,1])){
  #column to put the focals
  colunaF<-which(colnames(aux)==forCXR_Cd$Focal[i])
  #column to put the competitors
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

#joining data frames
forCXR_Cd<-cbind(forCXR_Cd, aux)

# Create fitness column
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
forCXR_Cd<-forCXR_Cd[-which(is.na(forCXR_Cd$fitness)),]

# removing rows without data
forCXR_Cd[which(forCXR_Cd$fitness=="-Inf" | forCXR_Cd$fitness=="Inf"),"fitness"]<-0


# Adding +1 to all growth rate data to include 0s in the data
forCXR_Cd$fitness<-forCXR_Cd$fitness+1

# vector that tells which are the selection regimes, the columns have to have the same name
my.reg <- c("SR1", "SR2","SR4","SR5")

# Do list per replicate and environment
R1_Cd<-list(SR1= subset(forCXR_Cd, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R2_Cd<-list(SR1= subset(forCXR_Cd, Rep==2 & Focal=="SR1")[,c("fitness", "SR1", "SR2","SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R3_Cd<-list(SR1= subset(forCXR_Cd, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R4_Cd<-list(SR1= subset(forCXR_Cd, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

R5_Cd<-list(SR1= subset(forCXR_Cd, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

if(!eval){
#### Importing parameters
 
 param_all_w0<-read.csv("./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal.csv")
 param_all_w0_upper<-read.csv("./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal_upper.csv")
 param_all_w0_lower<-read.csv("./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal_lower.csv")
#  
 param_all_w0<-param_all_w0[,-1]
 param_all_w0_upper<-param_all_w0_upper[,-1]
 param_all_w0_lower<-param_all_w0_lower[,-1]
}else{

## Run cxr no cadmium ------
# running the cxr to perform parameter estimates
lambda_start<- 1.8
alpha_intra_start<- -0.1
alpha_inter_start<- -0.1

print("R1 no cadmium")
obs.R1_w0<-cxr_pm_multifit(data = R1,
                           focal_column = my.reg,
                           model_family = "RK",
                           covariates = NULL,
                           optimization_method = "Nelder-Mead",
                           alpha_form = "pairwise",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none",
                           initial_values = list(lambda = lambda_start,
                                                 alpha_intra = alpha_intra_start,
                                                 alpha_inter = alpha_inter_start),
                           fixed_terms = NULL,
                           # no standard errors
                           bootstrap_samples = Nboot)

print("R3 no cadmium")
obs.R3_w0<-cxr_pm_multifit(data = R3,
                           focal_column = my.reg,
                           model_family = "RK",
                           covariates = NULL,
                           optimization_method = "Nelder-Mead",
                           alpha_form = "pairwise",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none",
                           initial_values = list(lambda = lambda_start,
                                                 alpha_intra = alpha_intra_start,
                                                 alpha_inter = alpha_inter_start),
                           fixed_terms = NULL,
                           # no standard errors
                           bootstrap_samples = Nboot)
print("R4 no cadmium")
obs.R4_w0<-cxr_pm_multifit(data = R4,
                           focal_column = my.reg,
                           model_family = "RK",
                           covariates = NULL,
                           optimization_method = "Nelder-Mead",
                           alpha_form = "pairwise",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none",
                           initial_values = list(lambda = lambda_start,
                                                 alpha_intra = alpha_intra_start,
                                                 alpha_inter = alpha_inter_start),
                           fixed_terms = NULL,
                           # no standard errors
                           bootstrap_samples = Nboot)
print("R5 no cadmium")
obs.R5_w0<-cxr_pm_multifit(data = R5,
                           focal_column = my.reg,
                           model_family = "RK",
                           covariates = NULL,
                           optimization_method = "Nelder-Mead",
                           alpha_form = "pairwise",
                           lambda_cov_form = "none",
                           alpha_cov_form = "none",
                           initial_values = list(lambda = lambda_start,
                                                 alpha_intra = alpha_intra_start,
                                                 alpha_inter = alpha_inter_start),
                           fixed_terms = NULL,
                           # no standard errors
                           bootstrap_samples = Nboot)
print("R2 no cadmium")
# For replicate 2 we need to do it differently
obs.R2_w0_sr1<-cxr_pm_fit(data = R2[[1]],
                          focal_column = my.reg[1],
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(lambda = lambda_start,
                                                alpha_intra = alpha_intra_start,
                                                alpha_inter = alpha_inter_start),
                          fixed_terms = NULL,
                          # no standard errors
                          bootstrap_samples = Nboot)

obs.R2_w0_sr4<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]==0),c("fitness", "SR4")],
                          focal_column =NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "global",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(lambda = lambda_start,
                                                alpha_inter = alpha_intra_start),
                          fixed_terms = NULL,
                          # no standard errors
                          bootstrap_samples = Nboot)

obs.R2_w0_sr4_inter<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]!=0),c("fitness", "SR1")],
                                focal_column =NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = alpha_inter_start),
                                fixed_terms = list(lambda=obs.R2_w0_sr4$lambda),
                                # no standard errors
                                bootstrap_samples = Nboot)

obs.R2_w0_sr5<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]==0),c("fitness", "SR5")],
                          focal_column =NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "global",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(lambda = lambda_start,
                                                alpha_inter = alpha_intra_start),
                          fixed_terms = NULL,
                          # no standard errors
                          bootstrap_samples = Nboot)

obs.R2_w0_sr5_inter<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]!=0),c("fitness", "SR1")],
                                focal_column =NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = alpha_inter_start),
                                fixed_terms = list(lambda=obs.R2_w0_sr5$lambda),
                                # no standard errors
                                bootstrap_samples = Nboot)


#rows in the alpha element of the returning list correspond to species i and columns to species j for each αij coefficient.

###### Storing mean parameter estimates
cxr_param_w0<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_w0$Tu_lambda<-0
cxr_param_w0$Te_lambda<-0
cxr_param_w0$Tu_intra<-0
cxr_param_w0$Te_intra<-0
cxr_param_w0$Tu_inter<-0
cxr_param_w0$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0<-cxr_param_w0[-which(cxr_param_w0$Replicate==2 & cxr_param_w0$Tu_Regime=="SR2"),]

# Growth rate
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

# Intraspecific competition
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

# Interspecific competition
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

### Storing lower estimates (mean-error)

cxr_param_w0_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_w0_lower$Tu_lambda<-0
cxr_param_w0_lower$Te_lambda<-0
cxr_param_w0_lower$Tu_intra<-0
cxr_param_w0_lower$Te_intra<-0
cxr_param_w0_lower$Tu_inter<-0
cxr_param_w0_lower$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0_lower<-cxr_param_w0_lower[-which(cxr_param_w0_lower$Replicate==2 & cxr_param_w0_lower$Tu_Regime=="SR2"),]

# Growth rate
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

# Intraspecific competition
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

# Interspecific competition
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

### Storing upper estimates (mean+error)

cxr_param_w0_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_w0_upper$Tu_lambda<-0
cxr_param_w0_upper$Te_lambda<-0
cxr_param_w0_upper$Tu_intra<-0
cxr_param_w0_upper$Te_intra<-0
cxr_param_w0_upper$Tu_inter<-0
cxr_param_w0_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0_upper<-cxr_param_w0_upper[-which(cxr_param_w0_upper$Replicate==2 & cxr_param_w0_upper$Tu_Regime=="SR2"),]

# Growth rate
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

# Intraspecific competition
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

# Interspecific competition
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

## Run cxr Cadmium------

print("R1 cadmium")
obs.R1_Cd_w0<-cxr_pm_multifit(data = R1_Cd,
                              focal_column = my.reg,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = lambda_start,
                                                    alpha_intra = alpha_intra_start,
                                                    alpha_inter = alpha_inter_start),
                              fixed_terms = NULL,
                              # no standard errors
                              bootstrap_samples = Nboot)

# replicate 2 below

print("R3 cadmium")
obs.R3_Cd_w0<-cxr_pm_multifit(data = R3_Cd,
                              focal_column = my.reg,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = lambda_start,
                                                    alpha_intra = alpha_intra_start,
                                                    alpha_inter = alpha_inter_start),
                              fixed_terms = NULL,
                              # no standard errors
                              bootstrap_samples = Nboot)
print("R4 cadmium")
obs.R4_Cd_w0<-cxr_pm_multifit(data = R4_Cd,
                              focal_column = my.reg,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = lambda_start,
                                                    alpha_intra = alpha_intra_start,
                                                    alpha_inter = alpha_inter_start),
                              fixed_terms = NULL,
                              # no standard errors
                              bootstrap_samples = Nboot)
print("R5 cadmium")
obs.R5_Cd_w0<-cxr_pm_multifit(data = R5_Cd,
                              focal_column = my.reg,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = lambda_start,
                                                    alpha_intra = alpha_intra_start,
                                                    alpha_inter = alpha_inter_start),
                              fixed_terms = NULL,
                              # no standard errors
                              bootstrap_samples = Nboot)


# Replicate 2
print("R2 cadmium")
obs.R2_Cd_w0_sr1<-cxr_pm_fit(data = R2_Cd[[1]],
                             focal_column = my.reg[1],
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = lambda_start,
                                                   alpha_intra = alpha_intra_start,
                                                   alpha_inter = alpha_inter_start),
                             fixed_terms = NULL,
                             # no standard errors
                             bootstrap_samples = Nboot)

obs.R2_Cd_w0_sr4<-cxr_pm_fit(data = R2_Cd[[2]][which(R2_Cd[[2]][,"SR1"]==0),c("fitness", "SR4")],
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "global",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = lambda_start,
                                                   alpha_inter = alpha_intra_start),
                             fixed_terms = NULL,
                             # no standard errors
                             bootstrap_samples = Nboot)

obs.R2_Cd_w0_sr5<-cxr_pm_fit(data = R2_Cd[[3]][which(R2_Cd[[3]][,"SR1"]==0),c("fitness", "SR5")],
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "global",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = lambda_start,
                                                   alpha_inter = alpha_intra_start),
                             fixed_terms = NULL,
                             # no standard errors
                             bootstrap_samples = Nboot)

obs.R2_Cd_w0_sr4_inter<-cxr_pm_fit(data = R2_Cd[[2]][which(R2_Cd[[2]][,"SR1"]!=0),c("fitness", "SR1")],
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list( alpha_inter = alpha_inter_start),
                                   fixed_terms = list(lambda=obs.R2_Cd_w0_sr4$lambda),
                                   # no standard errors
                                   bootstrap_samples = Nboot)

obs.R2_Cd_w0_sr5_inter<-cxr_pm_fit(data = R2_Cd[[3]][which(R2_Cd[[3]][,"SR1"]!=0),c("fitness", "SR1")],
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list( alpha_inter = alpha_inter_start),
                                   fixed_terms = list(lambda=obs.R2_Cd_w0_sr5$lambda),
                                   # no standard errors
                                   bootstrap_samples = Nboot)

###### Storing mean parameter estimates
cxr_param_w0C<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_w0C$Tu_lambda<-0
cxr_param_w0C$Te_lambda<-0
cxr_param_w0C$Tu_intra<-0
cxr_param_w0C$Te_intra<-0
cxr_param_w0C$Tu_inter<-0
cxr_param_w0C$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0C<-cxr_param_w0C[-which(cxr_param_w0C$Replicate==2 & cxr_param_w0C$Tu_Regime=="SR2"),]

# Growth rate
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

# Intraspecific competition
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

# Interspecific competition
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


### Storing lower estimates (mean-error)

cxr_param_w0C_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_w0C_lower$Tu_lambda<-0
cxr_param_w0C_lower$Te_lambda<-0
cxr_param_w0C_lower$Tu_intra<-0
cxr_param_w0C_lower$Te_intra<-0
cxr_param_w0C_lower$Tu_inter<-0
cxr_param_w0C_lower$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0C_lower<-cxr_param_w0C_lower[-which(cxr_param_w0C_lower$Replicate==2 & cxr_param_w0C_lower$Tu_Regime=="SR2"),]

# Growth rate
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

# Intraspecific competition
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

# Interspecific competition
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

### Storing upper estimates (mean+error)
cxr_param_w0C_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_w0C_upper$Tu_lambda<-0
cxr_param_w0C_upper$Te_lambda<-0
cxr_param_w0C_upper$Tu_intra<-0
cxr_param_w0C_upper$Te_intra<-0
cxr_param_w0C_upper$Tu_inter<-0
cxr_param_w0C_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_w0C_upper<-cxr_param_w0C_upper[-which(cxr_param_w0C_upper$Replicate==2 & cxr_param_w0C_upper$Tu_Regime=="SR2"),]

# Growth rate
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

# Intraspecific competition
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

# Interspecific competition
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


# Save parameters
write.csv(param_all_w0, "./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal.csv")
write.csv(param_all_w0_upper, "./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal_upper.csv")
write.csv(param_all_w0_lower, "./Analyses/MethodComparison/cxr_normal/parameters_cxr_normal_lower.csv")

}

# B - CXR lambda fixed --------------------
### 

#cxr accepts a data frame with a first column called fitness with positive values and numeric columns with number of individuals. Each row is one individual. For multiple species the easier is to create a list, each with a data frame that has in the first column number of individuals produced and then the number of neighbours
#this case we transformed all 0s into 1 (so that the log is 0) For that we need to add +1 to all data so that the variance is not changed
print("Running Method B: CXR with lambda fixed")
if(!dir.exists("./Analyses/MethodComparison/cxr_lambda_fixed_log")){
  dir.create("./Analyses/MethodComparison/cxr_lambda_fixed_log", showWarnings = FALSE)
}

## Create data frame -------

# estimating the mean growth rate for the cxr
mean_dens1<-data.frame(SR=c(rep(1,10), rep(2,8), rep(4,10),rep(5,10)), Env=c(rep("N", 5),rep("Cd", 5), rep("N", 4),rep("Cd", 4),rep("N", 5),rep("Cd", 5),rep("N", 5),rep("Cd", 5)), Rep=c(rep(c(1,2,3,4,5),2),rep(c(1,3,4,5),2),rep(c(1,2,3,4,5),2),rep(c(1,2,3,4,5),2)))


#since in the model we use the log of data +1, here we also have to use the +1 to estimate the lambda
# Estimating mean growth rate
mean_dens1$lambda<-sapply(c(1:length(mean_dens1[,1])), function(x){
  mean(subset(ca, FocalSR==mean_dens1$SR[x] & Dens==1 & Env==mean_dens1$Env[x] & Rep==mean_dens1$Rep[x] & GrowthRateOA!=0)$GrowthRateOA, na.rm=TRUE)
})

#Estimating standard deviation of the growth rate
mean_dens1$sd_lambda<-sapply(c(1:length(mean_dens1[,1])), function(x){
  sd(subset(ca, FocalSR==mean_dens1$SR[x] & Dens==1 & Env==mean_dens1$Env[x] & Rep==mean_dens1$Rep[x] & GrowthRateOA!=0)$GrowthRateOA, na.rm=TRUE)
})

# If the sd is NA or 0 we replace it with 0.01
mean_dens1$sd_lambda[which(is.na(mean_dens1$sd_lambda))]<-0.01
mean_dens1$sd_lambda[which(mean_dens1$sd_lambda==0)]<- 0.01

### No Cadmium ----------
#### Creating a list to store initial values for lambda

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

### Cadmium-------

# Estimate the fixed terms 
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

if(!eval){
##### importing data frame
param_all_B<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed.csv")
param_all_B_upper<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed_upper.csv")
param_all_B_lower<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed_lower.csv")

param_all_B<-param_all_B[,-1]
param_all_B_upper<-param_all_B_upper[,-1]
param_all_B_lower<-param_all_B_lower[,-1]
}else{
## Run cxr no cadmium --------


alpha_start<- -1
  
print("R1 no cadmium")
cxr_B.R1_w0<-cxr_pm_multifit(data = R1,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_intra = alpha_start,
                                                   alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_1N,
                             # no standard errors
                             bootstrap_samples = Nboot)


print("R3 no cadmium")
cxr_B.R3_w0<-cxr_pm_multifit(data = R3,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_intra = alpha_start,
                                                   alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_3N,
                             # no standard errors
                             bootstrap_samples = Nboot)
print("R4 no cadmium")
cxr_B.R4_w0<-cxr_pm_multifit(data = R4,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_intra = alpha_start,
                                                   alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_3N,
                             # no standard errors
                             bootstrap_samples = Nboot)
print("R5 no cadmium")
cxr_B.R5_w0<-cxr_pm_multifit(data = R5,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_intra = alpha_start,
                                                   alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_5N,
                             # no standard errors
                             bootstrap_samples = Nboot)

print("R2 no cadmium")
#for replicate 2 we will do the fitting by hand because we may need to scale the parameters

cxr_B.R2_w0_sr1<-cxr_pm_fit(data = R2[[1]],
                               focal_column = my.reg[1],
                               model_family = "RK",
                               covariates = NULL,
                               optimization_method = "Nelder-Mead",
                               alpha_form = "pairwise",
                               lambda_cov_form = "none",
                               alpha_cov_form = "none",
                               initial_values = list(alpha_intra = alpha_start,
                                                     alpha_inter = alpha_start),
                               fixed_terms = fixed_terms_2N[[1]],
                               # no standard errors
                               bootstrap_samples = Nboot)

cxr_B.R2_w0_sr4<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]==0), c("fitness", "SR4")],
                               focal_column = NULL,
                               model_family = "RK",
                               covariates = NULL,
                               optimization_method = "Nelder-Mead",
                               alpha_form = "global",
                               lambda_cov_form = "none",
                               alpha_cov_form = "none",
                               initial_values = list(alpha_inter = alpha_start),
                               fixed_terms = fixed_terms_2N[[2]],
                               # no standard errors
                               bootstrap_samples = Nboot)

cxr_B.R2_w0_5<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]==0), c("fitness", "SR5")],
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "global",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_2N[[3]],
                             # no standard errors
                             bootstrap_samples = Nboot)


cxr_B.R2_w0_sr4_inter<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]!=0), c("fitness", "SR1")],
                                     focal_column = NULL,
                                     model_family = "RK",
                                     covariates = NULL,
                                     optimization_method = "Nelder-Mead",
                                     alpha_form = "global",
                                     lambda_cov_form = "none",
                                     alpha_cov_form = "none",
                                     initial_values = list(alpha_inter = alpha_start),
                                     fixed_terms = fixed_terms_2N[[2]],
                                     # no standard errors
                                     bootstrap_samples = Nboot)

cxr_B.R2_w0_sr5_inter<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]!=0), c("fitness", "SR1")],
                                     focal_column = NULL,
                                     model_family = "RK",
                                     covariates = NULL,
                                     optimization_method = "Nelder-Mead",
                                     alpha_form = "global",
                                     lambda_cov_form = "none",
                                     alpha_cov_form = "none",
                                     initial_values = list(alpha_inter = alpha_start),
                                     fixed_terms = fixed_terms_2N[[3]],
                                     # no standard errors
                                     bootstrap_samples = Nboot)

# rows in the alpha element of the returning list correspond to species i and columns to species j for each αij coefficient.

###### Create data table summary
cxr_param_B<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_B$Tu_lambda<-0
cxr_param_B$Te_lambda<-0
cxr_param_B$Tu_intra<-0
cxr_param_B$Te_intra<-0
cxr_param_B$Tu_inter<-0
cxr_param_B$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_B<-cxr_param_B[-which(cxr_param_B$Replicate==2 & cxr_param_B$Tu_Regime=="SR2"),]

# Store growth rate
cxr_param_B[which(cxr_param_B$Replicate==1),"Tu_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[1]]$lambda,cxr_B.R1_w0$fixed_terms[[2]]$lambda)
cxr_param_B[which(cxr_param_B$Replicate==1),"Te_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[3]]$lambda,cxr_B.R1_w0$fixed_terms[[3]]$lambda, cxr_B.R1_w0$fixed_terms[[4]]$lambda,cxr_B.R1_w0$fixed_terms[[4]]$lambda)

cxr_param_B[which(cxr_param_B$Replicate==2),"Tu_lambda"]<-cxr_B.R2_w0_sr1$fixed_terms$lambda
cxr_param_B[which(cxr_param_B$Replicate==2),"Te_lambda"]<-c(cxr_B.R2_w0_sr4$fixed_terms$lambda, cxr_B.R2_w0_5$fixed_terms$lambda)

cxr_param_B[which(cxr_param_B$Replicate==3),"Tu_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[1]]$lambda,cxr_B.R3_w0$fixed_terms[[2]]$lambda)
cxr_param_B[which(cxr_param_B$Replicate==3),"Te_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[3]]$lambda,cxr_B.R3_w0$fixed_terms[[3]]$lambda, cxr_B.R3_w0$fixed_terms[[4]]$lambda,cxr_B.R3_w0$fixed_terms[[4]]$lambda)

cxr_param_B[which(cxr_param_B$Replicate==4),"Tu_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[1]]$lambda,cxr_B.R4_w0$fixed_terms[[2]]$lambda)
cxr_param_B[which(cxr_param_B$Replicate==4),"Te_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[3]]$lambda,cxr_B.R4_w0$fixed_terms[[3]]$lambda, cxr_B.R4_w0$fixed_terms[[4]]$lambda,cxr_B.R4_w0$fixed_terms[[4]]$lambda)

cxr_param_B[which(cxr_param_B$Replicate==5),"Tu_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[1]]$lambda,cxr_B.R5_w0$fixed_terms[[2]]$lambda)
cxr_param_B[which(cxr_param_B$Replicate==5),"Te_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[3]]$lambda,cxr_B.R5_w0$fixed_terms[[3]]$lambda, cxr_B.R5_w0$fixed_terms[[4]]$lambda,cxr_B.R5_w0$fixed_terms[[4]]$lambda)

# Store intraspecific competition
cxr_param_B[which(cxr_param_B$Replicate==1),"Tu_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[1,1], cxr_B.R1_w0$alpha_matrix[2,2]), 2)
cxr_param_B[which(cxr_param_B$Replicate==1),"Te_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[3,3], cxr_B.R1_w0$alpha_matrix[4,4]), each=2)

cxr_param_B[which(cxr_param_B$Replicate==2),"Tu_intra"]<-cxr_B.R2_w0_sr1$alpha_intra
cxr_param_B[which(cxr_param_B$Replicate==2),"Te_intra"]<-c(cxr_B.R2_w0_sr4$alpha_inter,cxr_B.R2_w0_sr5_inter$alpha_inter)

cxr_param_B[which(cxr_param_B$Replicate==3),"Tu_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[1,1], cxr_B.R3_w0$alpha_matrix[2,2]), 2)
cxr_param_B[which(cxr_param_B$Replicate==3),"Te_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[3,3], cxr_B.R3_w0$alpha_matrix[4,4]), each=2)

cxr_param_B[which(cxr_param_B$Replicate==4),"Tu_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[1,1], cxr_B.R4_w0$alpha_matrix[2,2]), 2)
cxr_param_B[which(cxr_param_B$Replicate==4),"Te_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[3,3], cxr_B.R4_w0$alpha_matrix[4,4]), each=2)

cxr_param_B[which(cxr_param_B$Replicate==5),"Tu_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[1,1], cxr_B.R5_w0$alpha_matrix[2,2]), 2)
cxr_param_B[which(cxr_param_B$Replicate==5),"Te_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[3,3], cxr_B.R5_w0$alpha_matrix[4,4]), each=2)

# Store interspecific competition
cxr_param_B[which(cxr_param_B$Replicate==1),"Tu_inter"]<-c(cxr_B.R1_w0$alpha_matrix[1,3], cxr_B.R1_w0$alpha_matrix[2,3],cxr_B.R1_w0$alpha_matrix[1,4], cxr_B.R1_w0$alpha_matrix[2,4])
cxr_param_B[which(cxr_param_B$Replicate==1),"Te_inter"]<-c(cxr_B.R1_w0$alpha_matrix[3,1], cxr_B.R1_w0$alpha_matrix[3,2],cxr_B.R1_w0$alpha_matrix[4,1], cxr_B.R1_w0$alpha_matrix[4,2])

cxr_param_B[which(cxr_param_B$Replicate==2),"Tu_inter"]<-cxr_B.R2_w0_sr1$alpha_inter[2:3]
cxr_param_B[which(cxr_param_B$Replicate==2),"Te_inter"]<-c(cxr_B.R2_w0_sr4_inter$alpha_inter, cxr_B.R2_w0_sr5_inter$alpha_inter)

cxr_param_B[which(cxr_param_B$Replicate==3),"Tu_inter"]<-c(cxr_B.R3_w0$alpha_matrix[1,3], cxr_B.R3_w0$alpha_matrix[2,3],cxr_B.R3_w0$alpha_matrix[1,4], cxr_B.R3_w0$alpha_matrix[2,4])
cxr_param_B[which(cxr_param_B$Replicate==3),"Te_inter"]<-c(cxr_B.R3_w0$alpha_matrix[3,1], cxr_B.R3_w0$alpha_matrix[3,2],cxr_B.R3_w0$alpha_matrix[4,1], cxr_B.R3_w0$alpha_matrix[4,2])

cxr_param_B[which(cxr_param_B$Replicate==4),"Tu_inter"]<-c(cxr_B.R4_w0$alpha_matrix[1,3], cxr_B.R4_w0$alpha_matrix[2,3],cxr_B.R4_w0$alpha_matrix[1,4], cxr_B.R4_w0$alpha_matrix[2,4])
cxr_param_B[which(cxr_param_B$Replicate==4),"Te_inter"]<-c(cxr_B.R4_w0$alpha_matrix[3,1], cxr_B.R4_w0$alpha_matrix[3,2],cxr_B.R4_w0$alpha_matrix[4,1], cxr_B.R4_w0$alpha_matrix[4,2])

cxr_param_B[which(cxr_param_B$Replicate==5),"Tu_inter"]<-c(cxr_B.R5_w0$alpha_matrix[1,3], cxr_B.R5_w0$alpha_matrix[2,3],cxr_B.R5_w0$alpha_matrix[1,4], cxr_B.R5_w0$alpha_matrix[2,4])
cxr_param_B[which(cxr_param_B$Replicate==5),"Te_inter"]<-c(cxr_B.R5_w0$alpha_matrix[3,1], cxr_B.R5_w0$alpha_matrix[3,2],cxr_B.R5_w0$alpha_matrix[4,1], cxr_B.R5_w0$alpha_matrix[4,2])

### Create data summary for lower boundaries
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

# Storing lower boundaries for the growth rate 
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Tu_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[1]]$lambda-sd_1N[[1]]$lambda,cxr_B.R1_w0$fixed_terms[[2]]$lambda-sd_1N[[2]]$lambda)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Te_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[3]]$lambda-sd_1N[[3]]$lambda,cxr_B.R1_w0$fixed_terms[[3]]$lambda-sd_1N[[3]]$lambda, cxr_B.R1_w0$fixed_terms[[4]]$lambda-sd_1N[[4]]$lambda,cxr_B.R1_w0$fixed_terms[[4]]$lambda-sd_1N[[4]]$lambda)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Tu_lambda"]<-cxr_B.R2_w0_sr1$fixed_terms$lambda-sd_2N[[1]]$lambda
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Te_lambda"]<-c(cxr_B.R2_w0_sr4$fixed_terms$lambda-sd_2N[[2]]$lambda, cxr_B.R2_w0_5$fixed_terms$lambda-sd_2N[[3]]$lambda)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Tu_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[1]]$lambda-sd_3N[[1]]$lambda,cxr_B.R3_w0$fixed_terms[[2]]$lambda-sd_3N[[2]]$lambda)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Te_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[3]]$lambda-sd_3N[[3]]$lambda,cxr_B.R3_w0$fixed_terms[[3]]$lambda-sd_3N[[3]]$lambda, cxr_B.R3_w0$fixed_terms[[4]]$lambda-sd_3N[[4]]$lambda,cxr_B.R3_w0$fixed_terms[[4]]$lambda-sd_3N[[4]]$lambda)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Tu_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[1]]$lambda-sd_4N[[1]]$lambda,cxr_B.R4_w0$fixed_terms[[2]]$lambda-sd_4N[[2]]$lambda)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Te_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[3]]$lambda-sd_4N[[3]]$lambda,cxr_B.R4_w0$fixed_terms[[3]]$lambda-sd_4N[[3]]$lambda, cxr_B.R4_w0$fixed_terms[[4]]$lambda-sd_4N[[4]]$lambda,cxr_B.R4_w0$fixed_terms[[4]]$lambda-sd_4N[[4]]$lambda)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Tu_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[1]]$lambda-sd_5N[[1]]$lambda,cxr_B.R5_w0$fixed_terms[[2]]$lambda-sd_5N[[2]]$lambda)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Te_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[3]]$lambda-sd_5N[[3]]$lambda,cxr_B.R5_w0$fixed_terms[[3]]$lambda-sd_5N[[3]]$lambda, cxr_B.R5_w0$fixed_terms[[4]]$lambda-sd_5N[[4]]$lambda,cxr_B.R5_w0$fixed_terms[[4]]$lambda-sd_5N[[4]]$lambda)

# Storing lower boundaries for the intraspecific competition 
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Tu_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[1,1]-cxr_B.R1_w0$alpha_matrix_standard_error[1,1], cxr_B.R1_w0$alpha_matrix[2,2]-cxr_B.R1_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Te_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[3,3]-cxr_B.R1_w0$alpha_matrix_standard_error[3,3], cxr_B.R1_w0$alpha_matrix[4,4]-cxr_B.R1_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Tu_intra"]<-cxr_B.R2_w0_sr1$alpha_intra[1]-cxr_B.R2_w0_sr1$alpha_intra_standard_error[1]
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Te_intra"]<-c(cxr_B.R2_w0_sr4$alpha_inter[1]-cxr_B.R2_w0_sr4$alpha_inter_standard_error[1], cxr_B.R2_w0_5$alpha_inter[1]-cxr_B.R2_w0_5$alpha_inter_standard_error[1])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Tu_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[1,1]-cxr_B.R3_w0$alpha_matrix_standard_error[1,1], cxr_B.R3_w0$alpha_matrix[2,2]-cxr_B.R3_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Te_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[3,3]-cxr_B.R3_w0$alpha_matrix_standard_error[3,3], cxr_B.R3_w0$alpha_matrix[4,4]-cxr_B.R3_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Tu_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[1,1]-cxr_B.R4_w0$alpha_matrix_standard_error[1,1], cxr_B.R4_w0$alpha_matrix[2,2]-cxr_B.R4_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Te_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[3,3]-cxr_B.R4_w0$alpha_matrix_standard_error[3,3], cxr_B.R4_w0$alpha_matrix[4,4]-cxr_B.R4_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Tu_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[1,1]-cxr_B.R5_w0$alpha_matrix_standard_error[1,1], cxr_B.R5_w0$alpha_matrix[2,2]-cxr_B.R5_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Te_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[3,3]-cxr_B.R5_w0$alpha_matrix_standard_error[3,3], cxr_B.R5_w0$alpha_matrix[4,4]-cxr_B.R5_w0$alpha_matrix_standard_error[4,4]), each=2)

# Storing lower boundaries for the interspecific competition 
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Tu_inter"]<-c(cxr_B.R1_w0$alpha_matrix[1,3]-cxr_B.R1_w0$alpha_matrix_standard_error[1,3], cxr_B.R1_w0$alpha_matrix[2,3]-cxr_B.R1_w0$alpha_matrix_standard_error[2,3],cxr_B.R1_w0$alpha_matrix[1,4]-cxr_B.R1_w0$alpha_matrix_standard_error[1,4], cxr_B.R1_w0$alpha_matrix[2,4]-cxr_B.R1_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==1),"Te_inter"]<-c(cxr_B.R1_w0$alpha_matrix[3,1]-cxr_B.R1_w0$alpha_matrix_standard_error[3,1], cxr_B.R1_w0$alpha_matrix[3,2]-cxr_B.R1_w0$alpha_matrix_standard_error[3,2],cxr_B.R1_w0$alpha_matrix[4,1]-cxr_B.R1_w0$alpha_matrix_standard_error[4,1], cxr_B.R1_w0$alpha_matrix[4,2]-cxr_B.R1_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Tu_inter"]<-cxr_B.R2_w0_sr1$alpha_inter[2:3]-cxr_B.R2_w0_sr1$alpha_inter_standard_error[2:3]
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==2),"Te_inter"]<-c(cxr_B.R2_w0_sr4_inter$alpha_inter[1]-cxr_B.R2_w0_sr4_inter$alpha_inter_standard_error[1], cxr_B.R2_w0_sr5_inter$alpha_inter[1]-cxr_B.R2_w0_sr5_inter$alpha_inter_standard_error[1])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Tu_inter"]<-c(cxr_B.R3_w0$alpha_matrix[1,3]-cxr_B.R3_w0$alpha_matrix_standard_error[1,3], cxr_B.R3_w0$alpha_matrix[2,3]-cxr_B.R3_w0$alpha_matrix_standard_error[2,3],cxr_B.R3_w0$alpha_matrix[1,4]-cxr_B.R3_w0$alpha_matrix_standard_error[1,4], cxr_B.R3_w0$alpha_matrix[2,4]-cxr_B.R3_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==3),"Te_inter"]<-c(cxr_B.R3_w0$alpha_matrix[3,1]-cxr_B.R3_w0$alpha_matrix_standard_error[3,1], cxr_B.R3_w0$alpha_matrix[3,2]-cxr_B.R3_w0$alpha_matrix_standard_error[3,2],cxr_B.R3_w0$alpha_matrix[4,1]-cxr_B.R3_w0$alpha_matrix_standard_error[4,1], cxr_B.R3_w0$alpha_matrix[4,2]-cxr_B.R3_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Tu_inter"]<-c(cxr_B.R4_w0$alpha_matrix[1,3]-cxr_B.R4_w0$alpha_matrix_standard_error[1,3], cxr_B.R4_w0$alpha_matrix[2,3]-cxr_B.R4_w0$alpha_matrix_standard_error[2,3],cxr_B.R4_w0$alpha_matrix[1,4]-cxr_B.R4_w0$alpha_matrix_standard_error[1,4], cxr_B.R4_w0$alpha_matrix[2,4]-cxr_B.R4_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==4),"Te_inter"]<-c(cxr_B.R4_w0$alpha_matrix[3,1]-cxr_B.R4_w0$alpha_matrix_standard_error[3,1], cxr_B.R4_w0$alpha_matrix[3,2]-cxr_B.R4_w0$alpha_matrix_standard_error[3,2],cxr_B.R4_w0$alpha_matrix[4,1]-cxr_B.R4_w0$alpha_matrix_standard_error[4,1], cxr_B.R4_w0$alpha_matrix[4,2]-cxr_B.R4_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Tu_inter"]<-c(cxr_B.R5_w0$alpha_matrix[1,3]-cxr_B.R5_w0$alpha_matrix_standard_error[1,3], cxr_B.R5_w0$alpha_matrix[2,3]-cxr_B.R5_w0$alpha_matrix_standard_error[2,3],cxr_B.R5_w0$alpha_matrix[1,4]-cxr_B.R5_w0$alpha_matrix_standard_error[1,4], cxr_B.R5_w0$alpha_matrix[2,4]-cxr_B.R5_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_lower[which(cxr_param_B_lower$Replicate==5),"Te_inter"]<-c(cxr_B.R5_w0$alpha_matrix[3,1]-cxr_B.R5_w0$alpha_matrix_standard_error[3,1], cxr_B.R5_w0$alpha_matrix[3,2]-cxr_B.R5_w0$alpha_matrix_standard_error[3,2],cxr_B.R5_w0$alpha_matrix[4,1]-cxr_B.R5_w0$alpha_matrix_standard_error[4,1], cxr_B.R5_w0$alpha_matrix[4,2]-cxr_B.R5_w0$alpha_matrix_standard_error[4,2])

### Data frame to store upper boundaries
cxr_param_B_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_B_upper$Tu_lambda<-0
cxr_param_B_upper$Te_lambda<-0
cxr_param_B_upper$Tu_intra<-0
cxr_param_B_upper$Te_intra<-0
cxr_param_B_upper$Tu_inter<-0
cxr_param_B_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_B_upper<-cxr_param_B_upper[-which(cxr_param_B_upper$Replicate==2 & cxr_param_B_upper$Tu_Regime=="SR2"),]

# Storing upper boundaries for the growth rate
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Tu_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[1]]$lambda+sd_1N[[1]]$lambda,cxr_B.R1_w0$fixed_terms[[2]]$lambda+sd_1N[[2]]$lambda)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Te_lambda"]<-c(cxr_B.R1_w0$fixed_terms[[3]]$lambda+sd_1N[[3]]$lambda,cxr_B.R1_w0$fixed_terms[[3]]$lambda+sd_1N[[3]]$lambda, cxr_B.R1_w0$fixed_terms[[4]]$lambda+sd_1N[[4]]$lambda,cxr_B.R1_w0$fixed_terms[[4]]$lambda+sd_1N[[4]]$lambda)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Tu_lambda"]<-cxr_B.R2_w0_sr1$fixed_terms$lambda+sd_2N[[1]]$lambda
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Te_lambda"]<-c(cxr_B.R2_w0_sr4$fixed_terms$lambda+sd_2N[[2]]$lambda, cxr_B.R2_w0_5$fixed_terms$lambda+sd_2N[[3]]$lambda)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Tu_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[1]]$lambda+sd_3N[[1]]$lambda,cxr_B.R3_w0$fixed_terms[[2]]$lambda+sd_3N[[2]]$lambda)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Te_lambda"]<-c(cxr_B.R3_w0$fixed_terms[[3]]$lambda+sd_3N[[3]]$lambda,cxr_B.R3_w0$fixed_terms[[3]]$lambda+sd_3N[[3]]$lambda, cxr_B.R3_w0$fixed_terms[[4]]$lambda+sd_3N[[4]]$lambda,cxr_B.R3_w0$fixed_terms[[4]]$lambda+sd_3N[[4]]$lambda)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Tu_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[1]]$lambda+sd_4N[[1]]$lambda,cxr_B.R4_w0$fixed_terms[[2]]$lambda+sd_4N[[2]]$lambda)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Te_lambda"]<-c(cxr_B.R4_w0$fixed_terms[[3]]$lambda+sd_4N[[3]]$lambda,cxr_B.R4_w0$fixed_terms[[3]]$lambda+sd_4N[[3]]$lambda, cxr_B.R4_w0$fixed_terms[[4]]$lambda+sd_4N[[4]]$lambda,cxr_B.R4_w0$fixed_terms[[4]]$lambda+sd_4N[[4]]$lambda)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Tu_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[1]]$lambda+sd_5N[[1]]$lambda,cxr_B.R5_w0$fixed_terms[[2]]$lambda+sd_5N[[2]]$lambda)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Te_lambda"]<-c(cxr_B.R5_w0$fixed_terms[[3]]$lambda+sd_5N[[3]]$lambda,cxr_B.R5_w0$fixed_terms[[3]]$lambda+sd_5N[[3]]$lambda, cxr_B.R5_w0$fixed_terms[[4]]$lambda+sd_5N[[4]]$lambda,cxr_B.R5_w0$fixed_terms[[4]]$lambda+sd_5N[[4]]$lambda)

# Storing lower boundaries for the intraspecific competition
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Tu_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[1,1]+cxr_B.R1_w0$alpha_matrix_standard_error[1,1], cxr_B.R1_w0$alpha_matrix[2,2]+cxr_B.R1_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Te_intra"]<-rep(c(cxr_B.R1_w0$alpha_matrix[3,3]+cxr_B.R1_w0$alpha_matrix_standard_error[3,3], cxr_B.R1_w0$alpha_matrix[4,4]+cxr_B.R1_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Tu_intra"]<-cxr_B.R2_w0_sr1$alpha_intra[1]+cxr_B.R2_w0_sr1$alpha_intra_standard_error[1]
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Te_intra"]<-c(cxr_B.R2_w0_sr4$alpha_inter[1]+cxr_B.R2_w0_sr4$alpha_inter_standard_error[1], cxr_B.R2_w0_5$alpha_inter[1]+cxr_B.R2_w0_5$alpha_inter_standard_error[1])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Tu_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[1,1]+cxr_B.R3_w0$alpha_matrix_standard_error[1,1], cxr_B.R3_w0$alpha_matrix[2,2]+cxr_B.R3_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Te_intra"]<-rep(c(cxr_B.R3_w0$alpha_matrix[3,3]+cxr_B.R3_w0$alpha_matrix_standard_error[3,3], cxr_B.R3_w0$alpha_matrix[4,4]+cxr_B.R3_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Tu_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[1,1]+cxr_B.R4_w0$alpha_matrix_standard_error[1,1], cxr_B.R4_w0$alpha_matrix[2,2]+cxr_B.R4_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Te_intra"]<-rep(c(cxr_B.R4_w0$alpha_matrix[3,3]+cxr_B.R4_w0$alpha_matrix_standard_error[3,3], cxr_B.R4_w0$alpha_matrix[4,4]+cxr_B.R4_w0$alpha_matrix_standard_error[4,4]), each=2)

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Tu_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[1,1]+cxr_B.R5_w0$alpha_matrix_standard_error[1,1], cxr_B.R5_w0$alpha_matrix[2,2]+cxr_B.R5_w0$alpha_matrix_standard_error[2,2]), 2)
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Te_intra"]<-rep(c(cxr_B.R5_w0$alpha_matrix[3,3]+cxr_B.R5_w0$alpha_matrix_standard_error[3,3], cxr_B.R5_w0$alpha_matrix[4,4]+cxr_B.R5_w0$alpha_matrix_standard_error[4,4]), each=2)

# Storing lower boundaries for the interspecific competition
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Tu_inter"]<-c(cxr_B.R1_w0$alpha_matrix[1,3]+cxr_B.R1_w0$alpha_matrix_standard_error[1,3], cxr_B.R1_w0$alpha_matrix[2,3]+cxr_B.R1_w0$alpha_matrix_standard_error[2,3],cxr_B.R1_w0$alpha_matrix[1,4]+cxr_B.R1_w0$alpha_matrix_standard_error[1,4], cxr_B.R1_w0$alpha_matrix[2,4]+cxr_B.R1_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==1),"Te_inter"]<-c(cxr_B.R1_w0$alpha_matrix[3,1]+cxr_B.R1_w0$alpha_matrix_standard_error[3,1], cxr_B.R1_w0$alpha_matrix[3,2]+cxr_B.R1_w0$alpha_matrix_standard_error[3,2],cxr_B.R1_w0$alpha_matrix[4,1]+cxr_B.R1_w0$alpha_matrix_standard_error[4,1], cxr_B.R1_w0$alpha_matrix[4,2]+cxr_B.R1_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Tu_inter"]<-cxr_B.R2_w0_sr1$alpha_inter[2:3]+cxr_B.R2_w0_sr1$alpha_inter_standard_error[2:3]
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==2),"Te_inter"]<-c(cxr_B.R2_w0_sr4_inter$alpha_inter[1]+cxr_B.R2_w0_sr4_inter$alpha_inter_standard_error[1], cxr_B.R2_w0_sr5_inter$alpha_inter[1]+cxr_B.R2_w0_sr5_inter$alpha_inter_standard_error[1])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Tu_inter"]<-c(cxr_B.R3_w0$alpha_matrix[1,3]+cxr_B.R3_w0$alpha_matrix_standard_error[1,3], cxr_B.R3_w0$alpha_matrix[2,3]+cxr_B.R3_w0$alpha_matrix_standard_error[2,3],cxr_B.R3_w0$alpha_matrix[1,4]+cxr_B.R3_w0$alpha_matrix_standard_error[1,4], cxr_B.R3_w0$alpha_matrix[2,4]+cxr_B.R3_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==3),"Te_inter"]<-c(cxr_B.R3_w0$alpha_matrix[3,1]+cxr_B.R3_w0$alpha_matrix_standard_error[3,1], cxr_B.R3_w0$alpha_matrix[3,2]+cxr_B.R3_w0$alpha_matrix_standard_error[3,2],cxr_B.R3_w0$alpha_matrix[4,1]+cxr_B.R3_w0$alpha_matrix_standard_error[4,1], cxr_B.R3_w0$alpha_matrix[4,2]+cxr_B.R3_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Tu_inter"]<-c(cxr_B.R4_w0$alpha_matrix[1,3]+cxr_B.R4_w0$alpha_matrix_standard_error[1,3], cxr_B.R4_w0$alpha_matrix[2,3]+cxr_B.R4_w0$alpha_matrix_standard_error[2,3],cxr_B.R4_w0$alpha_matrix[1,4]+cxr_B.R4_w0$alpha_matrix_standard_error[1,4], cxr_B.R4_w0$alpha_matrix[2,4]+cxr_B.R4_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==4),"Te_inter"]<-c(cxr_B.R4_w0$alpha_matrix[3,1]+cxr_B.R4_w0$alpha_matrix_standard_error[3,1], cxr_B.R4_w0$alpha_matrix[3,2]+cxr_B.R4_w0$alpha_matrix_standard_error[3,2],cxr_B.R4_w0$alpha_matrix[4,1]+cxr_B.R4_w0$alpha_matrix_standard_error[4,1], cxr_B.R4_w0$alpha_matrix[4,2]+cxr_B.R4_w0$alpha_matrix_standard_error[4,2])

cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Tu_inter"]<-c(cxr_B.R5_w0$alpha_matrix[1,3]+cxr_B.R5_w0$alpha_matrix_standard_error[1,3], cxr_B.R5_w0$alpha_matrix[2,3]+cxr_B.R5_w0$alpha_matrix_standard_error[2,3],cxr_B.R5_w0$alpha_matrix[1,4]+cxr_B.R5_w0$alpha_matrix_standard_error[1,4], cxr_B.R5_w0$alpha_matrix[2,4]+cxr_B.R5_w0$alpha_matrix_standard_error[2,4])
cxr_param_B_upper[which(cxr_param_B_upper$Replicate==5),"Te_inter"]<-c(cxr_B.R5_w0$alpha_matrix[3,1]+cxr_B.R5_w0$alpha_matrix_standard_error[3,1], cxr_B.R5_w0$alpha_matrix[3,2]+cxr_B.R5_w0$alpha_matrix_standard_error[3,2],cxr_B.R5_w0$alpha_matrix[4,1]+cxr_B.R5_w0$alpha_matrix_standard_error[4,1], cxr_B.R5_w0$alpha_matrix[4,2]+cxr_B.R5_w0$alpha_matrix_standard_error[4,2])


## Run cxr for Cadmium  --------------------
print("R1 cadmium")
cxr_B.R1_Cd_w0<-cxr_pm_multifit(data = R1_Cd,
                                focal_column = my.reg,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "pairwise",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_intra = alpha_start,
                                                      alpha_inter = alpha_start),
                                fixed_terms = fixed_terms_1Cd,
                                # no standard errors
                                bootstrap_samples = Nboot)

# replicate 2 below

print("R3 cadmium")
cxr_B.R3_Cd_w0<-cxr_pm_multifit(data = R3_Cd,
                                focal_column = my.reg,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "pairwise",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_intra = alpha_start,
                                                      alpha_inter = alpha_start),
                                fixed_terms = fixed_terms_3Cd,
                                # no standard errors
                                bootstrap_samples = Nboot)
print("R4 cadmium")
cxr_B.R4_Cd_w0<-cxr_pm_multifit(data = R4_Cd,
                                focal_column = my.reg,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "pairwise",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_intra = alpha_start,
                                                      alpha_inter = alpha_start),
                                fixed_terms = fixed_terms_4Cd,
                                # no standard errors
                                bootstrap_samples = Nboot)

print("R5 cadmium")
cxr_B.R5_Cd_w0<-cxr_pm_multifit(data = R5_Cd,
                                focal_column = my.reg,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "pairwise",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_intra = alpha_start,
                                                      alpha_inter = alpha_start),
                                fixed_terms = fixed_terms_5Cd,
                                # no standard errors
                                bootstrap_samples = Nboot)

print("R2 cadmium")
cxr_B.R2_Cd_w0_sr1<-cxr_pm_fit(data = R2_Cd[[1]],
                               focal_column = my.reg[1],
                               model_family = "RK",
                               covariates = NULL,
                               optimization_method = "Nelder-Mead",
                               alpha_form = "pairwise",
                               lambda_cov_form = "none",
                               alpha_cov_form = "none",
                               initial_values = list(alpha_intra = alpha_start,
                                                     alpha_inter = alpha_start),
                               fixed_terms = fixed_terms_2Cd[[1]],
                               # no standard errors
                               bootstrap_samples = Nboot)

#for replicate 2 we will do the fitting by hand because we may need to scale the parameters

cxr_B.R2_Cd_w0_sr4<-cxr_pm_fit(data = R2_Cd[[2]][which(R2_Cd[[2]][,"SR1"]==0), c("fitness", "SR4")],
                               focal_column = NULL,
                               model_family = "RK",
                               covariates = NULL,
                               optimization_method = "Nelder-Mead",
                               alpha_form = "global",
                               lambda_cov_form = "none",
                               alpha_cov_form = "none",
                               initial_values = list(alpha_inter = alpha_start),
                               fixed_terms = fixed_terms_2Cd[[2]],
                               # no standard errors
                               bootstrap_samples = Nboot)

cxr_B.R2_Cd_w0_5<-cxr_pm_fit(data = R2_Cd[[3]][which(R2_Cd[[3]][,"SR1"]==0), c("fitness", "SR5")],
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "global",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_2Cd[[3]],
                             # no standard errors
                             bootstrap_samples = Nboot)


cxr_B.R2_Cd_w0_sr4_inter<-cxr_pm_fit(data = R2_Cd[[2]][which(R2_Cd[[2]][,"SR1"]!=0), c("fitness", "SR1")],
                                     focal_column = NULL,
                                     model_family = "RK",
                                     covariates = NULL,
                                     optimization_method = "Nelder-Mead",
                                     alpha_form = "global",
                                     lambda_cov_form = "none",
                                     alpha_cov_form = "none",
                                     initial_values = list(alpha_inter = alpha_start),
                                     fixed_terms = fixed_terms_2Cd[[2]],
                                     # no standard errors
                                     bootstrap_samples = Nboot)

cxr_B.R2_Cd_w0_sr5_inter<-cxr_pm_fit(data = R2_Cd[[3]][which(R2_Cd[[3]][,"SR1"]!=0), c("fitness", "SR1")],
                                     focal_column = NULL,
                                     model_family = "RK",
                                     covariates = NULL,
                                     optimization_method = "Nelder-Mead",
                                     alpha_form = "global",
                                     lambda_cov_form = "none",
                                     alpha_cov_form = "none",
                                     initial_values = list(alpha_inter = alpha_start),
                                     fixed_terms = fixed_terms_2Cd[[3]],
                                     # no standard errors
                                     bootstrap_samples = Nboot)




###### Create data table summary to store parameter estimates
cxr_param_BC<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_BC$Tu_lambda<-0
cxr_param_BC$Te_lambda<-0
cxr_param_BC$Tu_intra<-0
cxr_param_BC$Te_intra<-0
cxr_param_BC$Tu_inter<-0
cxr_param_BC$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_BC<-cxr_param_BC[-which(cxr_param_BC$Replicate==2 & cxr_param_BC$Tu_Regime=="SR2"),]

# Storing growth rate estimates
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

# Storing intraspecific competition estimates
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

# Storing interspecific competition estimates
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

### Data frame to store lower boundary estimates
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

# Storing lower boundaries for growth rate
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

# Storing lower boundaries for intraspecific competition
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

# Storing lower boundaries for interspecific competition
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

### Create a data frame to store upper boundaries 
cxr_param_BC_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_BC_upper$Tu_lambda<-0
cxr_param_BC_upper$Te_lambda<-0
cxr_param_BC_upper$Tu_intra<-0
cxr_param_BC_upper$Te_intra<-0
cxr_param_BC_upper$Tu_inter<-0
cxr_param_BC_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_BC_upper<-cxr_param_BC_upper[-which(cxr_param_BC_upper$Replicate==2 & cxr_param_BC_upper$Tu_Regime=="SR2"),]

# Storing upper boundaries for growth rate
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

# Storing upper boundaries for intraspecific competition
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

# Storing upper boundaries for interspecific competition
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


##### Joining the two data frames
param_all_B<-as.data.frame(rbind(cxr_param_B, cxr_param_BC))

param_all_B_lower<-as.data.frame(rbind(cxr_param_B_lower, cxr_param_BC_lower))
param_all_B_upper<-as.data.frame(rbind(cxr_param_B_upper, cxr_param_BC_upper))

# Save parameters
write.csv(param_all_B, "./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed.csv")
write.csv(param_all_B_upper, "./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed_upper.csv")
write.csv(param_all_B_lower, "./Analyses/MethodComparison/cxr_lambda_fixed_log/parameters_cxr_lambda_fixed_lower.csv")
}

# C - CXR nested --------------------

#To do this we have to "trick" cxr, by putting the intraspecific competitors in another column than the focal and then estimate only intra!

print("Running Method C: CXR with lambda fixed and nested approach")

##### normal
if(!dir.exists("./Analyses/MethodComparison/cxr_lambda_fixed_nested")){
  dir.create("./Analyses/MethodComparison/cxr_lambda_fixed_nested", showWarnings = FALSE)
}

## Create data frame------

### No cadmium --------------------
CXR_C_N<-forCXR_N

# Do list per replicate and environment
R1_intra<-list(SR1= subset(CXR_C_N, Rep==1 & Focal=="SR1" & Comp=="SR1")[,c("fitness", "SR1")], SR2= subset(CXR_C_N, Rep==1 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_N, Rep==1 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==1 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R2_intra<-list(SR1= subset(CXR_C_N, Rep==2 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR4= subset(CXR_C_N, Rep==2 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==2 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R3_intra<-list(SR1= subset(CXR_C_N, Rep==3 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_N, Rep==3 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_N, Rep==3 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==3 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R4_intra<-list(SR1= subset(CXR_C_N, Rep==4 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_N, Rep==4 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_N, Rep==4 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==4 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R5_intra<-list(SR1= subset(CXR_C_N, Rep==5 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_N, Rep==5 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_N, Rep==5 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_N, Rep==5 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

# Interspecific data frame
R1<-list(SR1= subset(CXR_C_N, Rep==1 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR2= subset(CXR_C_N, Rep==1 & Focal=="SR2"& Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==1 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_N, Rep==1 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R2<-list(SR1= subset(CXR_C_N, Rep==2 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==2 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1")], SR5= subset(CXR_C_N, Rep==2 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1")])

R3<-list(SR1= subset(CXR_C_N, Rep==3 & Focal=="SR1" & Comp!="SR1")[,c("fitness",  "SR4", "SR5")], SR2= subset(CXR_C_N, Rep==3 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==3 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_N, Rep==3 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R4<-list(SR1= subset(CXR_C_N, Rep==4 & Focal=="SR1" & Comp!="SR1")[,c("fitness",  "SR4", "SR5")], SR2= subset(CXR_C_N, Rep==4 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==4 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_N, Rep==4 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R5<-list(SR1= subset(CXR_C_N, Rep==5 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR2= subset(CXR_C_N, Rep==5 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_N, Rep==5 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_N, Rep==5 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

### Cadmium -----

CXR_C_Cd<-forCXR_Cd

# Do list per replicate and environment
R1_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==1 & Focal=="SR1" & Comp=="SR1")[,c("fitness", "SR1")], SR2= subset(CXR_C_Cd, Rep==1 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_Cd, Rep==1 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==1 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R2_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==2 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR4= subset(CXR_C_Cd, Rep==2 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==2 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R3_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==3 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_Cd, Rep==3 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_Cd, Rep==3 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==3 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R4_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==4 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_Cd, Rep==4 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_Cd, Rep==4 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==4 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])

R5_cd_intra<-list(SR1= subset(CXR_C_Cd, Rep==5 & Focal=="SR1" & Comp=="SR1" )[,c("fitness", "SR1")], SR2= subset(CXR_C_Cd, Rep==5 & Focal=="SR2" & Comp=="SR2")[,c("fitness",  "SR2")], SR4= subset(CXR_C_Cd, Rep==5 & Focal=="SR4" & Comp=="SR4")[,c("fitness",  "SR4")], SR5= subset(CXR_C_Cd, Rep==5 & Focal=="SR5" & Comp=="SR5")[,c("fitness", "SR5")])


R1_cd<-list(SR1= subset(CXR_C_Cd, Rep==1 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR2= subset(CXR_C_Cd, Rep==1 & Focal=="SR2"& Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==1 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_Cd, Rep==1 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R2_cd<-list(SR1= subset(CXR_C_Cd, Rep==2 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==2 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1")], SR5= subset(CXR_C_Cd, Rep==2 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1")])

R3_cd<-list(SR1= subset(CXR_C_Cd, Rep==3 & Focal=="SR1" & Comp!="SR1")[,c("fitness",  "SR4", "SR5")], SR2= subset(CXR_C_Cd, Rep==3 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==3 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_Cd, Rep==3 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R4_cd<-list(SR1= subset(CXR_C_Cd, Rep==4 & Focal=="SR1" & Comp!="SR1")[,c("fitness",  "SR4", "SR5")], SR2= subset(CXR_C_Cd, Rep==4 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==4 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_Cd, Rep==4 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])

R5_cd<-list(SR1= subset(CXR_C_Cd, Rep==5 & Focal=="SR1" & Comp!="SR1")[,c("fitness", "SR4", "SR5")], SR2= subset(CXR_C_Cd, Rep==5 & Focal=="SR2" & Comp!="SR2")[,c("fitness", "SR4", "SR5")], SR4= subset(CXR_C_Cd, Rep==5 & Focal=="SR4" & Comp!="SR4")[,c("fitness", "SR1", "SR2")], SR5= subset(CXR_C_Cd, Rep==5 & Focal=="SR5" & Comp!="SR5")[,c("fitness", "SR1", "SR2")])


if(!eval){
  ##### importing data frame
  param_all_C<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed.csv")
  param_all_C_upper<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed_upper.csv")
  param_all_C_lower<-read.csv("./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed_lower.csv")
  
  param_all_C<-param_all_C[,-1]
  param_all_C_upper<-param_all_C_upper[,-1]
  param_all_C_lower<-param_all_C_lower[,-1]
  
}else{

## Run cxr for no cadmium ------

alpha_start<- -0.1
  
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

##### Intraspecific------------ 
print("R1 no cadmium")
cxr_C.R1_intra<-cxr_pm_multifit(data = R1_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = alpha_start),
                                fixed_terms = fixed_terms_1N,
                                # no standard errors
                                bootstrap_samples = Nboot)

print("R2 no cadmium")
cxr_C.R2_intra<-cxr_pm_multifit(data = R2_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = alpha_start),
                                fixed_terms = fixed_terms_2N,
                                # no standard errors
                                bootstrap_samples = Nboot)

print("R3 no cadmium")
cxr_C.R3_intra<-cxr_pm_multifit(data = R3_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = alpha_start),
                                fixed_terms = fixed_terms_3N,
                                # no standard errors
                                bootstrap_samples = Nboot)

print("R4 no cadmium")
cxr_C.R4_intra<-cxr_pm_multifit(data = R4_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = alpha_start),
                                fixed_terms = fixed_terms_4N,
                                # no standard errors
                                bootstrap_samples = Nboot)

print("R5 no cadmium")
cxr_C.R5_intra<-cxr_pm_multifit(data = R5_intra,
                                focal_column = NULL,
                                model_family = "RK",
                                covariates = NULL,
                                optimization_method = "Nelder-Mead",
                                alpha_form = "global",
                                lambda_cov_form = "none",
                                alpha_cov_form = "none",
                                initial_values = list(alpha_inter = alpha_start),
                                fixed_terms = fixed_terms_5N,
                                # no standard errors
                                bootstrap_samples = Nboot)


##### Interspecific ------

print("R1 no cadmium")
cxr_C.R1<-cxr_pm_multifit(data = R1,
                          focal_column = NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(alpha_inter = alpha_start),
                          fixed_terms = fixed_terms_1N,
                          # no standard errors
                          bootstrap_samples = Nboot)
print("R2 no cadmium")
cxr_C.R2_sr1<-cxr_pm_fit(data = R2[[1]],
                         focal_column = NULL,
                         model_family = "RK",
                         covariates = NULL,
                         optimization_method = "Nelder-Mead",
                         alpha_form = "pairwise",
                         lambda_cov_form = "none",
                         alpha_cov_form = "none",
                         initial_values = list(alpha_inter = alpha_start),
                         fixed_terms = fixed_terms_2N[[1]],
                         # no standard errors
                         bootstrap_samples = Nboot)

cxr_C.R2_sr4<-cxr_pm_fit(data = R2[[2]],
                         focal_column = NULL,
                         model_family = "RK",
                         covariates = NULL,
                         optimization_method = "Nelder-Mead",
                         alpha_form = "global",
                         lambda_cov_form = "none",
                         alpha_cov_form = "none",
                         initial_values = list(alpha_inter = alpha_start),
                         fixed_terms = fixed_terms_2N[[2]],
                         # no standard errors
                         bootstrap_samples = Nboot)

cxr_C.R2_sr5<-cxr_pm_fit(data = R2[[3]],
                         focal_column = NULL,
                         model_family = "RK",
                         covariates = NULL,
                         optimization_method = "Nelder-Mead",
                         alpha_form = "global",
                         lambda_cov_form = "none",
                         alpha_cov_form = "none",
                         initial_values = list(alpha_inter = alpha_start),
                         fixed_terms = fixed_terms_2N[[3]],
                         # no standard errors
                         bootstrap_samples = Nboot)

print("R3 no cadmium")
cxr_C.R3<-cxr_pm_multifit(data = R3,
                          focal_column = NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(alpha_inter = alpha_start),
                          fixed_terms = fixed_terms_3N,
                          # no standard errors
                          bootstrap_samples = Nboot)
print("R4 no cadmium")
cxr_C.R4<-cxr_pm_multifit(data = R4,
                          focal_column = NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(alpha_inter = alpha_start),
                          fixed_terms = fixed_terms_4N,
                          # no standard errors
                          bootstrap_samples = Nboot)
print("R5 no cadmium")
cxr_C.R5<-cxr_pm_multifit(data = R5,
                          focal_column = NULL,
                          model_family = "RK",
                          covariates = NULL,
                          optimization_method = "Nelder-Mead",
                          alpha_form = "pairwise",
                          lambda_cov_form = "none",
                          alpha_cov_form = "none",
                          initial_values = list(alpha_inter = alpha_start),
                          fixed_terms = fixed_terms_5N,
                          # no standard errors
                          bootstrap_samples = Nboot)


## Run cxr for Cadmium ------------

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


### Intraspecific------ 
print("R1 cadmium")
cxr_C.R1_cd_intra<-cxr_pm_multifit(data = R1_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = alpha_start),
                                   fixed_terms = fixed_terms_C_1N,
                                   # no standard errors
                                   bootstrap_samples = Nboot)

print("R2 cadmium")
cxr_C.R2_cd_intra<-cxr_pm_multifit(data = R2_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = alpha_start),
                                   fixed_terms = fixed_terms_C_2N,
                                   # no standard errors
                                   bootstrap_samples = Nboot)

print("R3 cadmium")
cxr_C.R3_cd_intra<-cxr_pm_multifit(data = R3_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = alpha_start),
                                   fixed_terms = fixed_terms_C_3N,
                                   # no standard errors
                                   bootstrap_samples = Nboot)

print("R4 cadmium")
cxr_C.R4_cd_intra<-cxr_pm_multifit(data = R4_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = alpha_start),
                                   fixed_terms = fixed_terms_C_4N,
                                   # no standard errors
                                   bootstrap_samples = Nboot)

print("R5 cadmium")
cxr_C.R5_cd_intra<-cxr_pm_multifit(data = R5_cd_intra,
                                   focal_column = NULL,
                                   model_family = "RK",
                                   covariates = NULL,
                                   optimization_method = "Nelder-Mead",
                                   alpha_form = "global",
                                   lambda_cov_form = "none",
                                   alpha_cov_form = "none",
                                   initial_values = list(alpha_inter = alpha_start),
                                   fixed_terms = fixed_terms_C_5N,
                                   # no standard errors
                                   bootstrap_samples = Nboot)


### Interspecific  --------

print("R1 cadmium")
cxr_C.R1_cd<-cxr_pm_multifit(data = R1_cd,
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_C_1N,
                             # no standard errors
                             bootstrap_samples = Nboot)

print("R2 cadmium")
cxr_C.R2_cd_sr1<-cxr_pm_fit(data = R2_cd[[1]],
                            focal_column = NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(alpha_inter = alpha_start),
                            fixed_terms = fixed_terms_C_2N[[1]],
                            # no standard errors
                            bootstrap_samples = Nboot)

cxr_C.R2_cd_sr4<-cxr_pm_fit(data = R2_cd[[2]],
                            focal_column = NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "global",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(alpha_inter = alpha_start),
                            fixed_terms = fixed_terms_C_2N[[2]],
                            # no standard errors
                            bootstrap_samples = Nboot)

cxr_C.R2_cd_sr5<-cxr_pm_fit(data = R2_cd[[3]],
                            focal_column = NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "global",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(alpha_inter = alpha_start),
                            fixed_terms = fixed_terms_C_2N[[3]],
                            # no standard errors
                            bootstrap_samples = Nboot)

print("R3 cadmium")
cxr_C.R3_cd<-cxr_pm_multifit(data = R3_cd,
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_C_3N,
                             # no standard errors
                             bootstrap_samples = Nboot)

print("R4 cadmium")
cxr_C.R4_cd<-cxr_pm_multifit(data = R4_cd,
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_C_4N,
                             # no standard errors
                             bootstrap_samples = Nboot)

print("R5 cadmium")
cxr_C.R5_cd<-cxr_pm_multifit(data = R5_cd,
                             focal_column = NULL,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(alpha_inter = alpha_start),
                             fixed_terms = fixed_terms_C_5N,
                             # no standard errors
                             bootstrap_samples = Nboot)


#rows in the alpha element of the returning list correspond to species i and columns to species j for each αij coefficient.

###### Create data table to store parameter estimates
# No cadmium environment
cxr_param_C<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_C$Tu_lambda<-0
cxr_param_C$Te_lambda<-0
cxr_param_C$Tu_intra<-0
cxr_param_C$Te_intra<-0
cxr_param_C$Tu_inter<-0
cxr_param_C$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_C<-cxr_param_C[-which(cxr_param_C$Replicate==2 & cxr_param_C$Tu_Regime=="SR2"),]

# Storing growth rates 
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

# Storing intraspecific competition 
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

# Storing interspecific competition
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

### Store estimates for lower boundaries
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

# Storing lower boundaries for growth rates 
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

# Storing lower boundaries for intraspecific competition 
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

# Storing lower boundaries for interspecific competition
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

### Storing upper estimates
cxr_param_C_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
cxr_param_C_upper$Tu_lambda<-0
cxr_param_C_upper$Te_lambda<-0
cxr_param_C_upper$Tu_intra<-0
cxr_param_C_upper$Te_intra<-0
cxr_param_C_upper$Tu_inter<-0
cxr_param_C_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_C_upper<-cxr_param_C_upper[-which(cxr_param_C_upper$Replicate==2 & cxr_param_C_upper$Tu_Regime=="SR2"),]

# Storing upper boundaries for growth rates 
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

# Storing upper boundaries for intraspecific competition
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

# Storing upper boundaries for interspecific competition
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


###### Storing data table summary for the cadmium environment
cxr_param_CC<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_CC$Tu_lambda<-0
cxr_param_CC$Te_lambda<-0
cxr_param_CC$Tu_intra<-0
cxr_param_CC$Te_intra<-0
cxr_param_CC$Tu_inter<-0
cxr_param_CC$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_CC<-cxr_param_CC[-which(cxr_param_CC$Replicate==2 & cxr_param_CC$Tu_Regime=="SR2"),]

# Storing the growth rate estimates
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

# Storing the intraspecific competition estimates
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

# Storing the interspecific competition estimates
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

### Storing Lower boundaries estimates for the cadmium environment

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

# Storing the lower boundaries for the growth rate estimates
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

# Storing the lower boundaries for the intraspecific competition estimates
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

# Storing the lower boundaries for the interspecific competition estimates
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

### Storing the upper boundary estimates
cxr_param_CC_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
cxr_param_CC_upper$Tu_lambda<-0
cxr_param_CC_upper$Te_lambda<-0
cxr_param_CC_upper$Tu_intra<-0
cxr_param_CC_upper$Te_intra<-0
cxr_param_CC_upper$Tu_inter<-0
cxr_param_CC_upper$Te_inter<-0

#removing SR2 for replicate 2
cxr_param_CC_upper<-cxr_param_CC_upper[-which(cxr_param_CC_upper$Replicate==2 & cxr_param_CC_upper$Tu_Regime=="SR2"),]

# Storing the upper boundaries for the growth rate estimates
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

# Storing the upper boundaries for the intraspecific competition estimates
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

# Storing the upper boundaries for the interspecific competition estimates
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


write.csv(param_all_C, "./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed.csv")
write.csv(param_all_C_upper, "./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed_upper.csv")
write.csv(param_all_C_lower, "./Analyses/MethodComparison/cxr_lambda_fixed_nested/parameters_cxr_lambda_fixed_lower.csv")
}
# D - optim normal ------

##### Estimating parameters
print("Running Method D: Optim nested")

# creating folder to store the analyses, this should be the same as the file path in the function

if(!dir.exists("./Analyses/MethodComparison/Optim_normal")){
  dir.create("./Analyses/MethodComparison/Optim_normal", showWarnings = FALSE)
}

## Run optim ----
# This matrix has all the comparisons that need to be done between regimes
comparison_mat<-matrix(nrow=4, ncol=3)
comparison_mat[1,]<-c(1,4,5)
comparison_mat[2,]<-c(2,4,5)
comparison_mat[3,]<-c(4,1,2)
comparison_mat[4,]<-c(5,1,2)

#lam2 is the data from density one corresponding to the focals populations
# data2 is the data (format) Regime (name of focal pop), background (name of competitor, the same if its intraspecific competition), focal (number of focal individuals in g0), comp (number of competitors in g0), growth rate
# Attention that for intraspecific you need to add 0 in the comp and all individuals -1 in the focal

# Estimating the different parameters for all replicates and the two environments.
if(eval){

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

# For replicate two we have to change the comparison matrix because there is no replicate 2 for SR2
comparison_mat2<-matrix(nrow=3, ncol=3)
comparison_mat2[1,]<-c(1,4,5)
comparison_mat2[2,]<-c(4,1,NA)
comparison_mat2[3,]<-c(5,1,NA)

rep2<-mod_df(subset(ca,Rep==2 & Env=="N"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=2, env="N", comparisons = comparison_mat2)

rep2<-mod_df(subset(ca,Rep==2 & Env=="Cd"))
magic_rk(filepath2 = "./Analyses/MethodComparison/Optim_normal/", lam2=dataForLambda, data2=rep2, reps2=2, env="Cd", comparisons = comparison_mat2)
}

## Importing file parameters ---- 
##### Importing files of alpha and lambda
#the alpha matrices are always tu, te (row), tu, te (col)
# In this case we need to import the files and parse them

# First step get a list of files for the different parameters
alpha_file<-list.files("./Analyses/MethodComparison/Optim_normal/", pattern="alpha_estimates") 

alphaUpper_file<-list.files("./Analyses/MethodComparison/Optim_normal/", pattern="alpha_upper")

alphaLower_file<-list.files("./Analyses/MethodComparison/Optim_normal/", pattern="alpha_lower")

lambda_file<-list.files("./Analyses/MethodComparison/Optim_normal/", pattern="lambda_estimates")

# Parsing the list of files and importing them
alpha_list<- lapply(alpha_file, function(x) read.csv(paste("./Analyses/MethodComparison/Optim_normal/",x, sep=""), header = TRUE))
alphaUpper_list<- lapply(alphaUpper_file, function(x) read.csv(paste("./Analyses/MethodComparison/Optim_normal/",x, sep=""), header = TRUE))
alphaLower_list<- lapply(alphaLower_file, function(x) read.csv(paste("./Analyses/MethodComparison/Optim_normal/",x, sep=""), header = TRUE))
lambda_list<- lapply(lambda_file, function(x) read.csv(paste("./Analyses/MethodComparison/Optim_normal/",x, sep=""), header = TRUE))

# passing from list to data frame
# First we need to do the first iteration (to the structure to store the parameters)
lambda_intra_fixed<-data.frame(Regime1=rep(c(1,1,2,2),10), Regime2=rep(c(4,5,4,5), 10), Replicate=c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8)), Env=rep(c(rep("N",4), rep("Cd",4)), 5))

lambda_intra_fixed<-lambda_intra_fixed[-which(lambda_intra_fixed$Regime1==2 & lambda_intra_fixed$Replicate==2),] # to remove SR2 from replicate 2 because it does not exist

# Getting the names of replicates and environments
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[2]

# Creating vectors to store the names of the selection regimes, environments and replicates
regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

# Creating the first iteration of the data frame
aux_alpha<-as.data.frame(alpha_list[[1]])

aux_alpha2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu=c(aux_alpha[1,2], aux_alpha[1,2], aux_alpha[2,2],aux_alpha[2,2]), intraTe=c(aux_alpha[3,2], aux_alpha[4,2], aux_alpha[3,2],aux_alpha[4,2]), interTu=c(aux_alpha[1,3], aux_alpha[1,4], aux_alpha[2,3], aux_alpha[2,4]), interTe=c(aux_alpha[3,3], aux_alpha[4,3], aux_alpha[3,4], aux_alpha[4,4]))

# Passing information from list to data frame
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


### Doing the same for the lower boundaries for competition (Alpha Lower)

# Setting up the items needed to do the structure
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

# Creating the first iteration of the data frame to store the lower boundary estimates
aux_alphaLower<-as.data.frame(alphaLower_list[[1]])

aux_alphaLower2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_L=c(aux_alphaLower[1,2], aux_alphaLower[1,2], aux_alphaLower[2,2],aux_alphaLower[2,2]), intraTe_L=c(aux_alphaLower[3,2], aux_alphaLower[4,2], aux_alphaLower[3,2],aux_alphaLower[4,2]), interTu_L=c(aux_alphaLower[1,3], aux_alphaLower[1,4], aux_alphaLower[2,3], aux_alphaLower[2,4]), interTe_L=c(aux_alphaLower[3,3], aux_alphaLower[4,3], aux_alphaLower[3,4], aux_alphaLower[4,4]))

# Passing information from list to data frame
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

### Doing the same for the upper boundaries for competition (Alpha Lower)

# Setting up the items needed to do the structure
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

# Creating the first iteration of the data frame to store the upper boundary estimates

aux_alphaUpper<-as.data.frame(alphaUpper_list[[1]])

aux_alphaUpper2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_U=c(aux_alphaUpper[1,2], aux_alphaUpper[1,2], aux_alphaUpper[2,2],aux_alphaUpper[2,2]), intraTe_U=c(aux_alphaUpper[3,2], aux_alphaUpper[4,2], aux_alphaUpper[3,2],aux_alphaUpper[4,2]), interTu_U=c(aux_alphaUpper[1,3], aux_alphaUpper[1,4], aux_alphaUpper[2,3], aux_alphaUpper[2,4]), interTe_U=c(aux_alphaUpper[3,3], aux_alphaUpper[4,3], aux_alphaUpper[3,4], aux_alphaUpper[4,4]))

# Passing information from list to data frame

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

# Passing lambda list to data frame
## Setting up the items needed to do the structure

repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[2]

Focal<-c("1","1","2","2","4","4","5","5")
Comp<-c("4","5","4","5","1","2","1","2")
Env<-rep(env, 8)
Rep<-rep(repli,8)

# Creating the first iteration of the data frame 
aux_lambda<-cbind(as.data.frame(lambda_list[[1]])[,c(3,4,5)],Focal,Comp, Env, Rep)

# Passing from list to data frame
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


#Joining the different data frames

alphas_mat_D<-as.data.frame(cbind(aux_alpha2, aux_alphaLower2, aux_alphaUpper2))

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

str(alphas_mat_D)

colnames(alphas_mat_D)<-c("Tu_Regime", "Te_Regime", "Environment", "Replicate", "Tu_intra", "Te_intra", "Tu_inter", "Te_inter", "Tu_intra_L", "Te_intra_L", "Tu_inter_L", "Te_inter_L", "Tu_intra_U", "Te_intra_U", "Tu_inter_U", "Te_inter_U", "Tu_lambda", "Te_lambda","Tu_lambda_L", "Te_lambda_L","Tu_lambda_U", "Te_lambda_U")


# E - optim lambda fixed ---------
print("Running Method E: Optim nested with lambda fixed")

##### Estimate parameters

# creating folder to put the analyses inside, this should be the same as the file path in the function
if(!dir.exists("./Analyses/MethodComparison/optim_lambda_fixed_nested")){
  dir.create("./Analyses/MethodComparison/optim_lambda_fixed_nested", showWarnings = FALSE)
}


# This matrix has all the comparisons that need to be done between regimes
comparison_mat<-matrix(nrow=4, ncol=3)
comparison_mat[1,]<-c(1,4,5)
comparison_mat[2,]<-c(2,4,5)
comparison_mat[3,]<-c(4,1,2)
comparison_mat[4,]<-c(5,1,2)

## Run optim ----

#lam2 is the data from density one corresponding to the focals populations
# data2 is the data (format) Regime (name of focal pop), background (name of competitor, the same if its intraspecific competition), focal (number of focal individuals in g0), comp (number of competitors in g0), growth rate
# Attention that for intraspecific you need to add 0 in the comp and all individuals in the focal
if(eval){
rep2<-mod_df(subset(ca,Rep==1 & Env=="N"))  
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/",data2=rep2, reps2=1, env="N", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==1 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==1 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/", data2=rep2, reps2=1, env="Cd", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==1 & Env=="Cd"))

rep2<-mod_df(subset(ca,Rep==3 & Env=="N"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/",  data2=rep2, reps2=3, env="N", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==3 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==3 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/", data2=rep2, reps2=3, env="Cd", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==3 & Env=="Cd"))

rep2<-mod_df(subset(ca,Rep==4 & Env=="N"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/", data2=rep2, reps2=4, env="N", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==4 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==4 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/", data2=rep2, reps2=4, env="Cd", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==4 & Env=="Cd"))

rep2<-mod_df(subset(ca,Rep==5 & Env=="N"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/", data2=rep2, reps2=5, env="N", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==5 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==5 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/",data2=rep2, reps2=5, env="Cd", comparisons = comparison_mat, lam2=subset(mean_dens1, Rep==5 & Env=="Cd"))

# For two we have to change the comparison matrix
comparison_mat2<-matrix(nrow=3, ncol=3)
comparison_mat2[1,]<-c(1,4,5)
comparison_mat2[2,]<-c(4,1,NA)
comparison_mat2[3,]<-c(5,1,NA)

rep2<-mod_df(subset(ca,Rep==2 & Env=="N"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/",  data2=rep2, reps2=2, env="N", comparisons = comparison_mat2, lam2=subset(mean_dens1, Rep==2 & Env=="N"))

rep2<-mod_df(subset(ca,Rep==2 & Env=="Cd"))
magic_rk_lambda(filepath2 = "./Analyses/MethodComparison/optim_lambda_fixed_nested/",  data2=rep2, reps2=2, env="Cd", comparisons = comparison_mat2, lam2=subset(mean_dens1, Rep==2 & Env=="Cd"))
}

## Importing file parameters ----
##### Importing files of alpha and lambda
#the alpha matrices are always tu, te (row), tu, te (col)
# In this case we need to import the files and parse them

# First step get a list of files for the different parameters

alpha_file<-list.files("./Analyses/MethodComparison/optim_lambda_fixed_nested/", pattern="alpha_estimates") #the alphas are always tu, te (row), tu, te (col)

alphaUpper_file<-list.files("./Analyses/MethodComparison/optim_lambda_fixed_nested/", pattern="alpha_upper")

alphaLower_file<-list.files("./Analyses/MethodComparison/optim_lambda_fixed_nested/", pattern="alpha_lower")

lambda_file<-list.files("./Analyses/MethodComparison/optim_lambda_fixed_nested/", pattern="lambda_estimates")

# Parsing the list of files and importing them

alpha_list<- lapply(alpha_file, function(x) read.csv(paste("./Analyses/MethodComparison/optim_lambda_fixed_nested/",x, sep=""), header = TRUE))
alphaUpper_list<- lapply(alphaUpper_file, function(x) read.csv(paste("./Analyses/MethodComparison/optim_lambda_fixed_nested/",x, sep=""), header = TRUE))
alphaLower_list<- lapply(alphaLower_file, function(x) read.csv(paste("./Analyses/MethodComparison/optim_lambda_fixed_nested/",x, sep=""), header = TRUE))
lambda_list<- lapply(lambda_file, function(x) read.csv(paste("./Analyses/MethodComparison/optim_lambda_fixed_nested/",x, sep=""), header = TRUE))

# passing from list to data frame
# First we need to do the first iteration (to create everything)
lambda_intra_fixed<-data.frame(Regime1=rep(c(1,1,2,2),10), Regime2=rep(c(4,5,4,5), 10), Replicate=c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8)), Env=rep(c(rep("N",4), rep("Cd",4)), 5))

lambda_intra_fixed<-lambda_intra_fixed[-which(lambda_intra_fixed$Regime1==2 & lambda_intra_fixed$Replicate==2),] # to remove SR2 from replicate 2 because it does not exist

# Setting up variables to create the needed data structure
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

# Transforming the first position of the list into data frame to set up the structure
aux_alpha<-as.data.frame(alpha_list[[1]])

aux_alpha2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu=c(aux_alpha[1,2], aux_alpha[1,2], aux_alpha[2,2],aux_alpha[2,2]), intraTe=c(aux_alpha[3,2], aux_alpha[4,2], aux_alpha[3,2],aux_alpha[4,2]), interTu=c(aux_alpha[1,3], aux_alpha[1,4], aux_alpha[2,3], aux_alpha[2,4]), interTe=c(aux_alpha[3,3], aux_alpha[4,3], aux_alpha[3,4], aux_alpha[4,4]))

# Doing a loop to pass the list to data frame
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

### Alpha Lower - Same as above but now creating the data frame with the lower boundaries

# Setting up variables to create the needed data structure
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaLower_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

# Transforming the first position of the list into data frame to set up the structure
aux_alphaLower<-as.data.frame(alphaLower_list[[1]])

aux_alphaLower2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_L=c(aux_alphaLower[1,2], aux_alphaLower[1,2], aux_alphaLower[2,2],aux_alphaLower[2,2]), intraTe_L=c(aux_alphaLower[3,2], aux_alphaLower[4,2], aux_alphaLower[3,2],aux_alphaLower[4,2]), interTu_L=c(aux_alphaLower[1,3], aux_alphaLower[1,4], aux_alphaLower[2,3], aux_alphaLower[2,4]), interTe_L=c(aux_alphaLower[3,3], aux_alphaLower[4,3], aux_alphaLower[3,4], aux_alphaLower[4,4]))

# Doing a loop to pass the list to data frame

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

# Setting up variables to create the needed data structure
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alphaUpper_file[1], split="_")[1])[4], split="[.]"))[1],split=""))[2]

regimeTu<-c("1","1", "2","2")
regimeTe<-c("4","5", "4","5")
Env<-rep(env, 4)
Rep<-rep(repli,4)

# Transforming the first position of the list into data frame to set up the structure
aux_alphaUpper<-as.data.frame(alphaUpper_list[[1]])

aux_alphaUpper2<-data.frame(regimeTu, regimeTe, Env, Rep, intraTu_U=c(aux_alphaUpper[1,2], aux_alphaUpper[1,2], aux_alphaUpper[2,2],aux_alphaUpper[2,2]), intraTe_U=c(aux_alphaUpper[3,2], aux_alphaUpper[4,2], aux_alphaUpper[3,2],aux_alphaUpper[4,2]), interTu_U=c(aux_alphaUpper[1,3], aux_alphaUpper[1,4], aux_alphaUpper[2,3], aux_alphaUpper[2,4]), interTe_U=c(aux_alphaUpper[3,3], aux_alphaUpper[4,3], aux_alphaUpper[3,4], aux_alphaUpper[4,4]))

# Doing a loop to pass the list to data frame
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

# Creating the data frame to store the lambda (growth rate) parameters

# Setting up variables to create the needed data structure
repli<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[1]
env<-unlist(strsplit(unlist(strsplit(unlist(strsplit(alpha_file[1], split="_")[1])[6], split="[.]"))[1],split=""))[2]

Focal<-c("1","1","2","2","4","4","5","5")
Comp<-c("4","5","4","5","1","2","1","2")
Env<-rep(env, 8)
Rep<-rep(repli,8)

# Transforming the first position of the list into data frame to set up the structure
aux_lambda<-cbind(as.data.frame(lambda_list[[1]])[,c(3,4,5)],Focal,Comp, Env, Rep)

# Doing a loop to pass the list to data frame
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


# Joining all data frames
alphas_mat_E<-as.data.frame(cbind(aux_alpha2, aux_alphaLower2, aux_alphaUpper2))


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

colnames(alphas_mat_E)<-c("Tu_Regime", "Te_Regime", "Environment", "Replicate", "Tu_intra", "Te_intra", "Tu_inter", "Te_inter", "Tu_intra_L", "Te_intra_L", "Tu_inter_L", "Te_inter_L", "Tu_intra_U", "Te_intra_U", "Tu_inter_U", "Te_inter_U", "Tu_lambda", "Te_lambda","Tu_lambda_L", "Te_lambda_L","Tu_lambda_U", "Te_lambda_U")

#################################################'
# Compare methods----
#################################################'


# Creating columns to ensure that th names of selection regimes match across methods
alphas_mat_D$Tu_Regime2<-alphas_mat_D$Tu_Regime
alphas_mat_D$Te_Regime2<-alphas_mat_D$Te_Regime


alphas_mat_D$Tu_Regime<-mapvalues(alphas_mat_D$Tu_Regime2, c("1","2","4","5"), c("SR1", "SR2", "SR4","SR5"))
alphas_mat_D$Te_Regime<-mapvalues(alphas_mat_D$Te_Regime2, c("1","2","4","5"), c("SR1", "SR2", "SR4","SR5"))


alphas_mat_E$Tu_Regime2<-alphas_mat_E$Tu_Regime
alphas_mat_E$Te_Regime2<-alphas_mat_E$Te_Regime


alphas_mat_E$Tu_Regime<-mapvalues(alphas_mat_E$Tu_Regime2, c("1","2","4","5"), c("SR1", "SR2", "SR4","SR5"))
alphas_mat_E$Te_Regime<-mapvalues(alphas_mat_E$Te_Regime2, c("1","2","4","5"), c("SR1", "SR2", "SR4","SR5"))


# Creating columns to be sure that the names of the environments match betweeb methods

alphas_mat_D$Environment2<-alphas_mat_D$Environment

## Predict observed data ----
alphas_mat_D$Env2<-mapvalues(alphas_mat_D$Environment, c("C", "N"), c("Cd","N"))

alphas_mat_D$Environment<-mapvalues(alphas_mat_D$Environment2, c("N","C"), c("N", "Cd"))

alphas_mat_E$Environment2<-alphas_mat_E$Environment

alphas_mat_E$Environment<-mapvalues(alphas_mat_E$Environment2, c("N","C"), c("N", "Cd"))


### Adding column with the information about the medthos used to estimate

param_all_w0$Method<-"cxr normal"
param_all_B$Method<-"cxr lambda fixed"
param_all_C$Method<-"cxr lambda fixed, nested"
alphas_mat_D$Method<-"optim"
alphas_mat_E$Method<-"optim lambda fixed"

# Joining the data frames
cols_to_join<-c("Tu_Regime", "Te_Regime", "Environment", "Replicate","Tu_lambda","Te_lambda", "Tu_intra","Te_intra", "Tu_inter", "Te_inter", "Method" )

comparison_methods<-rbind(param_all_w0[,cols_to_join],param_all_B[,cols_to_join],param_all_C[,cols_to_join], alphas_mat_D[,cols_to_join], alphas_mat_E[,cols_to_join] )

### Plots -----------
# Plots comparing the parameter estimates for the different methods
ggplot(comparison_methods, aes(x=Method, y=Tu_lambda))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(fill=Environment), alpha=0.35, outlier.colour = NA)+
  geom_point(aes(colour=Environment, fill=Environment, shape=as.factor(Replicate)), position = position_dodge2(0.5))+
  theme_plots+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Tu lambda")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))+
  scale_fill_manual(values=c("darkblue", "darkred"))

ggplot(comparison_methods, aes(x=Method, y=Te_lambda))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(fill=Environment), alpha=0.75, outlier.colour = NA)+
  geom_point(aes(colour=Environment, fill=Environment, shape=Replicate), position = position_dodge2(0.5))+
  theme_plots+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Te lambda")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))+
  scale_fill_manual(values=c("darkblue", "darkred"))


ggplot(comparison_methods, aes(x=Method, y=Tu_intra))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(fill=Environment), alpha=0.35, outlier.colour = NA)+
  geom_point(aes(colour=Environment, fill=Environment, shape=Replicate), position = position_dodge2(0.5))+
  theme_plots+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Tu intra")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))+
  scale_fill_manual(values=c("darkblue", "darkred"))

ggplot(comparison_methods, aes(x=Method, y=Te_intra))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(fill=Environment), alpha=0.35, outlier.colour = NA)+
  geom_point(aes(colour=Environment, fill=Environment, shape=Replicate), position = position_dodge2(0.5))+  theme_plots+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Te intra")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))+
  scale_fill_manual(values=c("darkblue", "darkred"))

ggplot(comparison_methods, aes(x=Method, y=Tu_inter))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(fill=Environment), alpha=0.35, outlier.colour = NA)+
  geom_point(aes(colour=Environment, fill=Environment, shape=Replicate), position = position_dodge2(0.5))+  theme_plots+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Tu inter")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))+
  scale_fill_manual(values=c("darkblue", "darkred"))

ggplot(comparison_methods, aes(x=Method, y=Te_inter))+
  facet_grid(Tu_Regime~Te_Regime)+
  geom_boxplot(aes(fill=Environment), alpha=0.35, outlier.colour = NA)+
  geom_point(aes(colour=Environment, fill=Environment, shape=Replicate), position = position_dodge2(0.5))+  theme_plots+
  theme_bw()+
  xlab("Methods used to estimate data")+
  ylab("Te inter")+
  scale_x_discrete(labels=c("cxr", "cxr\nlambda","cxr\nnested", "optim", "optim\nlambda"))+
  scale_colour_manual(values=c("darkblue", "darkred"))+
  scale_fill_manual(values=c("darkblue", "darkred"))


# Comparing performance between methods ----
print("Comparing predicted values")

## Estimate distance between predicted and observed-----

#We will estimate the predicted vs observed for each method and use that as metric to define which method to use going forward

## Predict observed data ----

# Setting up the columns
ca$FocalSR3<-mapvalues(ca$FocalSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

ca$CompSR3<-mapvalues(ca$CompSR, c(1,2,4,5), c("SR1", "SR2","SR4","SR5"))

ca$Env3<-mapvalues(ca$Env, c("N", "Cd"), c("N", "C"))

# For each method we will predict the number of females that would be observed according to the estimates obtained
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
    densF<-ca$Dens[x]-1
    pred<-lambda*densF*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*densF*exp(-alpha_ij*densC)
  }
  # if the predicted number would be negative it becomes 0
  #print(pred)
  if(pred<0){
    pred<-0
  }
    
  log(pred+1)
  
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
    densF<-ca$Dens[x]-1
    pred<-lambda*densF*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*densF*exp(-alpha_ij*densC)
  }
  #print(pred)
  # if the predicted number would be negative it becomes 0
  if(is.na(pred) | pred<0){
    pred<-0
  }
    
  log(pred+1)
  
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
    densF<-ca$Dens[x]-1
    pred<-lambda*densF*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*densF*exp(-alpha_ij*densC)
  }
  # if the predicted number would be negative it becomes 0
  if(pred<0){
    pred<-0
  }
    
  log(pred+1)
  
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
    densF<-ca$Dens[x]-1
    pred<-lambda*densF*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*densF*exp(-alpha_ij*densC)
  }
  # if the predicted number would be negative it becomes 0
  if(pred<0){
    pred<-0
  }
    
  log(pred+1)
  
})


alphas_mat_E$Env2<-mapvalues(alphas_mat_E$Environment, c("C", "N"), c("Cd","N"))

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
    densF<-ca$Dens[x]-1
    pred<-lambda*densF*exp(-alpha_i*(densF))
    
  }else if(ca$Type[x]=="INTER"){
    densC<-ca$Dens[x]-1
    densF<-1
    pred<-lambda*densF*exp(-alpha_ij*densC)
  }
  # if the predicted number would be negative it becomes 0
  if(pred<0){
    pred<-0
  }
    
  log(pred+1)
  
})



### Calculate distances ---- 
#Do not forget that this is the log of GR +1

euclidean <- function(a, b) sqrt(sum((a - b)^2))

ca$distA<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_A[x], log(ca$GrowthRate[x]+1))

  euc
})

ca$distB<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_B[x], log(ca$GrowthRate[x]+1))

  euc
})

ca$distC<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_C[x], log(ca$GrowthRate[x]+1))

  euc
})

ca$distD<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_D[x], log(ca$GrowthRate[x]+1))

  euc
})

ca$distE<-sapply(c(1:length(ca$Block)), function(x){
  euc<-euclidean(ca$pred_E[x], log(ca$GrowthRate[x]+1))

  euc
})

sum(ca$distA, na.rm = TRUE)
sum(ca$distB, na.rm = TRUE)
sum(ca$distC, na.rm = TRUE)
sum(ca$distD, na.rm = TRUE)
sum(ca$distE, na.rm = TRUE)
# 
 length(which(!is.na(ca$distA)))
 length(which(!is.na(ca$distB)))
 length(which(!is.na(ca$distC)))
 length(which(!is.na(ca$distD)))
 length(which(!is.na(ca$distE)))


distance_sum<-pivot_longer(ca[, c(43:47)], cols = c(1:5),names_to = "method", values_to = "distance")

dist_sum<-distance_sum %>% group_by(method) %>% summarize(mean=mean(distance, na.rm=T), sd=sd(distance, na.rm=TRUE))

dist_sum$se<-c(dist_sum$sd[1]/sqrt(3563),dist_sum$sd[2]/sqrt(3563), dist_sum$sd[3]/sqrt(3563), dist_sum$sd[4]/sqrt(3563), dist_sum$sd[5]/sqrt(3563) )

### Plotting distance----

ggplot(distance_sum, aes(x=method, y=distance, colour=method, fill=method))+
   geom_boxplot(colour="black", outlier.colour = NA)+
  geom_point(alpha=0.10, position=position_dodge2(0.35), colour="black", shape=21)+
  theme_bw()+
  theme_plots+
  scale_x_discrete(labels=c("cxr", "cxr lambda\nfixed", "cxr \nnested", "optim", "optim \nlambda fixed"), name="Method")+
  scale_color_brewer(palette = "Spectral")+
  scale_fill_brewer(palette = "Spectral")+
  theme(legend.position = "none")+
  geom_text(data=dist_sum, aes(x=method, label=paste(round(mean,3), round(se,3), sep="\n+/-")), y=3, colour="black")+
  scale_y_continuous(name="Estimated euclidean distance\n (predicted-observed")


if(eval){
save_plot("./Plots/FigS2.pdf", width = 20, height=10)
save_plot("./Plots/FigS2.png", width = 20, height=10)
}

