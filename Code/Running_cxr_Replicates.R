#' ---
#' title: "R Notebook"
#' output: html_notebook
#' ---
#' 

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

# Functions, general information and packages---------------------------
#' 
## 
if(!require("plyr")){
  install.packages("plyr")
}
if(!require("dplyr")){
  install.packages("dplyr")
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
if(!require("DescTools")){
  install.packages("DescTools")
}
if(!require("LSAfun")){
  install.packages("LSAfun")
}

library(plyr)
library(dplyr)
library(tidyr)
library(cxr)
library(MASS)
library(mvtnorm)
library(DescTools)
library(LSAfun)

# Creating vectors with regime names to use for plots
regimeTu<-c("Tu \ncontrol", "Tu evolved \n in cadmium")
names(regimeTu)<-c("SR1", "SR2")

regimeTe<-c("Te \n control", "Te evolved \n in cadmium")
names(regimeTe)<-c("SR4", "SR5")

#' 
# Import data---------------------------
#' This chunk is used to import data, checking data and creating the columns that are necessary
#' 
## 
ca<-read.csv(file = "./Data/CompetitiveAbility_Cd_G40_submit.csv", header=TRUE) # data from the competitive ability

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
# Running cxr---------------------------
#' 
#' ### Setup cxr
#' Setting the seed number and the number of bootstrap samples.
#' 
## 
set.seed(1809)
bootN<-2000
eval<- TRUE

# Creating directory to store the results
# If the directory does not exist create it 
if(!dir.exists("./Analyses/")){
  dir.create("./Analyses/")
}
dir.create("./Analyses/cxr_normal_Replicates", showWarnings = FALSE)

#' 
#' ## Set up data
#' 
#' ### Creating a data frame for the no cadmium environment
#' 
## Data setup no cadmium ---------------------------
print("Creating data files for the no cadmium environment")

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
    #column where to put the focals
    colunaF<-which(colnames(aux)==forCXR_N$Focal[i])
    #column where to put the competitors
    colunaC<-which(colnames(aux)==forCXR_N$Comp[i])
    
    #if its the same regime
    if(forCXR_N$Focal[i]==forCXR_N$Comp[i]){
      aux[i,colunaF]<-forCXR_N$Dens[i]-1
    }else{ #if it is heterospecific then its -1 for the competitors (because of the focal)
      aux[i,colunaC]<-forCXR_N$Dens[i]-1
      aux[i, colunaF]<-0
    }
    
  }

  #joining the two data frames
  forCXR_N<-cbind(forCXR_N, aux)
  
  # Loop that creates the column corresponding to the fitness of each observation
  forCXR_N$fitness<-sapply(c(1:length(forCXR_N[,1])), function(x){
    colF<-which(colnames(forCXR_N)==forCXR_N$Focal[x])
    # Depending on the focal selection regime we need to divide by different number of individuals
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

  # Adding +1 to all data to allow for the 0s to be included in the cxr estimates
  forCXR_N$fitness<-forCXR_N$fitness+1
  
  # vector that tells which are the selection regimes, the columns have to have the same name
  my.reg <- c("SR1", "SR2","SR4","SR5")


#' 
#' ### Creating a data frame for the cadmium environment
#' This is the same procedure as above
## Data setup cadmium environment ---------------------------

print("Creating data files for the cadmium environment")
# subsetting data to only contain cadmium

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

  # Adding +1 to all data to allow for the 0s to be included in the cxr estimates
  forCXR_Cd$fitness<-forCXR_Cd$fitness+1
  
  # vector that tells which are the selection regimes, the columns have to have the same name
  my.reg <- c("SR1", "SR2","SR4","SR5")

#' 
#' ### Data for replicates no cadmium
#' 
## Replicate list no cadmium---------------------------
# Do list per replicate and environment
  R1<-list(SR1= subset(forCXR_N, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  R2<-list(SR1= subset(forCXR_N, Rep==2 & Focal=="SR1")[,c("fitness", "SR1", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR4")], SR5= subset(forCXR_N, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR5")])
  
  R3<-list(SR1= subset(forCXR_N, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  R4<-list(SR1= subset(forCXR_N, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  R5<-list(SR1= subset(forCXR_N, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  

#' 
#' ### Data for replicates cadmium
#' 
## Replicate list cadmium---------------------------
Rep_R1_Cd<-list(SR1= subset(forCXR_Cd, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  Rep_R2_Cd<-list(SR1= subset(forCXR_Cd, Rep==2 & Focal=="SR1")[,c("fitness", "SR1","SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR4")], SR5= subset(forCXR_Cd, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR5")])
  
  Rep_R3_Cd<-list(SR1= subset(forCXR_Cd, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  Rep_R4_Cd<-list(SR1= subset(forCXR_Cd, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  Rep_R5_Cd<-list(SR1= subset(forCXR_Cd, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])

#' 
#' ## Data for the initial values
#' Creating a data frame with the initial values. The values to use have been tested in the parameter exploration code file.
replicates_initial<-as.data.frame(expand.grid(Replicates=c(1,2,3,4,5), Environment=c("N", "Cd"), lambda=0, alpha_intra=0, alpha_inter=0))

replicates_initial$lambda<-1.1001
replicates_initial$alpha_intra<- 0
replicates_initial$alpha_inter<- -0.1

#' 
#' # Running cxr for the no cadmium environment
#' For each replicate we use cxr multifit function and then store the results.
#' 
## Run cxr for no cadmium env---------------------------
print("Running cxr to get parameter estimates for the no cadmium environment")
### Replicate 1-----
obs.R1_w0<-cxr_pm_multifit(data = R1,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==1 & replicates_initial$Environment=="N"),"lambda"],
                                                   alpha_intra = replicates_initial[which(replicates_initial$Replicates==1 & replicates_initial$Environment=="N"),"alpha_intra"],
                                                   alpha_inter = replicates_initial[which(replicates_initial$Replicates==1 & replicates_initial$Environment=="N"),"alpha_inter"]),
                             fixed_terms = NULL,
                             # no standard errors
                              bootstrap_samples = bootN)

# save observation
save(obs.R1_w0,file= "./Analyses/cxr_normal_Replicates/cxr_R1_N.RData")

### Replicate 3-----
obs.R3_w0<-cxr_pm_multifit(data = R3,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==3 & replicates_initial$Environment=="N"),"lambda"],
                                                   alpha_intra = replicates_initial[which(replicates_initial$Replicates==3 & replicates_initial$Environment=="N"),"alpha_intra"],
                                                   alpha_inter = replicates_initial[which(replicates_initial$Replicates==3 & replicates_initial$Environment=="N"),"alpha_inter"]),
                             fixed_terms = NULL,
                             # no standard errors
                              bootstrap_samples = bootN)

#save data
save(obs.R3_w0,file= "./Analyses/cxr_normal_Replicates/cxr_R3_N.RData")

### Replicate 4-----
obs.R4_w0<-cxr_pm_multifit(data = R4,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==4 & replicates_initial$Environment=="N"),"lambda"],
                                                   alpha_intra = replicates_initial[which(replicates_initial$Replicates==4 & replicates_initial$Environment=="N"),"alpha_intra"],
                                                   alpha_inter = replicates_initial[which(replicates_initial$Replicates==4 & replicates_initial$Environment=="N"),"alpha_inter"]),
                             fixed_terms = NULL,
                             # no standard errors
                              bootstrap_samples = bootN)
# save data
save(obs.R4_w0,file= "./Analyses/cxr_normal_Replicates/cxr_R4_N.RData")

### Replicate 5-----
obs.R5_w0<-cxr_pm_multifit(data = R5,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==5 & replicates_initial$Environment=="N"),"lambda"],
                                                   alpha_intra = replicates_initial[which(replicates_initial$Replicates==5 & replicates_initial$Environment=="N"),"alpha_intra"],
                                                   alpha_inter = replicates_initial[which(replicates_initial$Replicates==5 & replicates_initial$Environment=="N"),"alpha_inter"]),
                             fixed_terms = NULL,
                             # no standard errors
                              bootstrap_samples = bootN)
# save data
save(obs.R5_w0,file= "./Analyses/cxr_normal_Replicates/cxr_R5_N.RData")

# Because there is no SR2 for replicate 2 we need to estimate SR1 and SR4 and SR5 in a different way
### Replicate 2-----

# fitting for SR1
obs.R2_w0_sr1<-cxr_pm_fit(data = R2[[1]],
                            focal_column = my.reg[1],
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="N"),"lambda"],
                                                  alpha_intra = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="N"),"alpha_intra"],
                                                  alpha_inter = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="N"),"alpha_inter"]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = bootN)
    # Fitting for sr4 lambda and intra
    obs.R2_w0_sr4<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]==0),c("fitness", "SR4")],
                            focal_column =NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "global",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="N"),"lambda"],
                                                  alpha_inter = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="N"),"alpha_intra"]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = bootN)
    
    # fitting for sr4 inter
    obs.R2_w0_sr4_inter<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]!=0),c("fitness", "SR1")],
                                  focal_column =NULL,
                                  model_family = "RK",
                                  covariates = NULL,
                                  optimization_method = "Nelder-Mead",
                                  alpha_form = "global",
                                  lambda_cov_form = "none",
                                  alpha_cov_form = "none",
                                  initial_values = list(alpha_inter=replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="N"),"alpha_inter"]),
                                  fixed_terms = list(lambda=obs.R2_w0_sr4$lambda),
                                  # no standard errors
                                   bootstrap_samples = bootN)
    # fitting for SR5 lambda and intra
    obs.R2_w0_sr5<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]==0),c("fitness", "SR5")],
                            focal_column =NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "global",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="N"),"lambda"],
                                                  alpha_inter = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="N"),"alpha_intra"]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = bootN)
    
    #fitting for SR5 inter
    obs.R2_w0_sr5_inter<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]!=0),c("fitness", "SR1")],
                                  focal_column =NULL,
                                  model_family = "RK",
                                  covariates = NULL,
                                  optimization_method = "Nelder-Mead",
                                  alpha_form = "global",
                                  lambda_cov_form = "none",
                                  alpha_cov_form = "none",
                                  initial_values = list(alpha_inter=replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="N"),"alpha_inter"]),
                                  fixed_terms = list(lambda=obs.R2_w0_sr5$lambda),
                                  # no standard errors
                                   bootstrap_samples = bootN)
# saving data
save(obs.R2_w0_sr1, file="./Analyses/cxr_normal_Replicates/cxr_R2_SR1_N.RData")
save(obs.R2_w0_sr4,file= "./Analyses/cxr_normal_Replicates/cxr_R2_SR4_N.RData")
save(obs.R2_w0_sr5, file="./Analyses/cxr_normal_Replicates/cxr_R2_SR5_N.RData")
save(obs.R2_w0_sr4_inter, file="./Analyses/cxr_normal_Replicates/cxr_R2_SR4_inter_N.RData")
save(obs.R2_w0_sr5_inter, file="./Analyses/cxr_normal_Replicates/cxr_R2_SR5_inter_N.RData")


#' 
### Storing parameter estimates (mean, lower and upper)---------------------------
#' Storing the parameter estimates and the lower and upper bounds in a data frame so we can use it later.

print("Creating data frame to store mean parameter estimates for the no cadmium environment")
#### Mean estimates ------
# creating data frame where to store the data
cxr_param_w0<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
  cxr_param_w0$Tu_lambda<-0
  cxr_param_w0$Te_lambda<-0
  cxr_param_w0$Tu_intra<-0
  cxr_param_w0$Te_intra<-0
  cxr_param_w0$Tu_inter<-0
  cxr_param_w0$Te_inter<-0
  
  #removing SR2 for replicate 2
  cxr_param_w0<-cxr_param_w0[-which(cxr_param_w0$Replicate==2 & cxr_param_w0$Tu_Regime=="SR2"),]
  
# Storing parameters: intrinsic growth rate
  cxr_param_w0[which(cxr_param_w0$Replicate==1),"Tu_lambda"]<-obs.R1_w0$lambda[c(1,2,1,2)]
  cxr_param_w0[which(cxr_param_w0$Replicate==1),"Te_lambda"]<-obs.R1_w0$lambda[c(3,3,4,4)]
  
  cxr_param_w0[which(cxr_param_w0$Replicate==2),"Tu_lambda"]<-obs.R2_w0_sr1$lambda
  cxr_param_w0[which(cxr_param_w0$Replicate==2),"Te_lambda"]<-c(obs.R2_w0_sr4$lambda,obs.R2_w0_sr5$lambda)
  
  cxr_param_w0[which(cxr_param_w0$Replicate==3),"Tu_lambda"]<-obs.R3_w0$lambda[1:2]
  cxr_param_w0[which(cxr_param_w0$Replicate==3),"Te_lambda"]<-obs.R3_w0$lambda[c(3,3,4,4)]
  
  cxr_param_w0[which(cxr_param_w0$Replicate==4),"Tu_lambda"]<-obs.R4_w0$lambda[c(1,2,1,2)]
  cxr_param_w0[which(cxr_param_w0$Replicate==4),"Te_lambda"]<-obs.R4_w0$lambda[c(3,3,4,4)]
  
  cxr_param_w0[which(cxr_param_w0$Replicate==5),"Tu_lambda"]<-obs.R5_w0$lambda[c(1,2,1,2)]
  cxr_param_w0[which(cxr_param_w0$Replicate==5),"Te_lambda"]<-obs.R5_w0$lambda[c(3,3,4,4)]
  
 # Storing parameters: intraspecific competition
  cxr_param_w0[which(cxr_param_w0$Replicate==1),"Tu_intra"]<-rep(c(obs.R1_w0$alpha_matrix[1,1], obs.R1_w0$alpha_matrix[2,2]), 2)
  cxr_param_w0[which(cxr_param_w0$Replicate==1),"Te_intra"]<-rep(c(obs.R1_w0$alpha_matrix[3,3], obs.R1_w0$alpha_matrix[4,4]), each=2)
  
  cxr_param_w0[which(cxr_param_w0$Replicate==2),"Tu_intra"]<-obs.R2_w0_sr1$alpha_intra
  #Although the parameter is called alpha_inter, it was estimated using only intraspecific competitors.
  cxr_param_w0[which(cxr_param_w0$Replicate==2),"Te_intra"]<-c(obs.R2_w0_sr4$alpha_inter, obs.R2_w0_sr5$alpha_inter)
  
  cxr_param_w0[which(cxr_param_w0$Replicate==3),"Tu_intra"]<-rep(c(obs.R3_w0$alpha_matrix[1,1], obs.R3_w0$alpha_matrix[2,2]), 2)
  cxr_param_w0[which(cxr_param_w0$Replicate==3),"Te_intra"]<-rep(c(obs.R3_w0$alpha_matrix[3,3], obs.R3_w0$alpha_matrix[4,4]), each=2)
  
  cxr_param_w0[which(cxr_param_w0$Replicate==4),"Tu_intra"]<-rep(c(obs.R4_w0$alpha_matrix[1,1], obs.R4_w0$alpha_matrix[2,2]), 2)
  cxr_param_w0[which(cxr_param_w0$Replicate==4),"Te_intra"]<-rep(c(obs.R4_w0$alpha_matrix[3,3], obs.R4_w0$alpha_matrix[4,4]), each=2)
  
  cxr_param_w0[which(cxr_param_w0$Replicate==5),"Tu_intra"]<-rep(c(obs.R5_w0$alpha_matrix[1,1], obs.R5_w0$alpha_matrix[2,2]), 2)
  cxr_param_w0[which(cxr_param_w0$Replicate==5),"Te_intra"]<-rep(c(obs.R5_w0$alpha_matrix[3,3], obs.R5_w0$alpha_matrix[4,4]), each=2)
  
  # Storing parameters: interspecific competition
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

#' 
#### Lower estimates ------
### Here we will apply the same reasoning but now estimating the lower estimates, using the stats from cxr. Lower bounds correspond to mean-error.
  print("Creating data frame to store lower parameter estimates for the no cadmium environment")
  
  cxr_param_w0_lower<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
  cxr_param_w0_lower$Tu_lambda<-0
  cxr_param_w0_lower$Te_lambda<-0
  cxr_param_w0_lower$Tu_intra<-0
  cxr_param_w0_lower$Te_intra<-0
  cxr_param_w0_lower$Tu_inter<-0
  cxr_param_w0_lower$Te_inter<-0
  
  #removing SR2 for replicate 2
  cxr_param_w0_lower<-cxr_param_w0_lower[-which(cxr_param_w0_lower$Replicate==2 & cxr_param_w0_lower$Tu_Regime=="SR2"),]
   # Storing parameters: intrinsic growth rate
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
  
   # Storing parameters: intraspecific competition
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
  
   # Storing parameters: interspecific competition
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

#' 
#' #### Data frame with upper estimates
#' 
#### Upper estimates ---------------------------
### Here we will apply the same reasoning but now estimating the upper estimates, using the stats from cxr. Upper bounds correspond to mean-error.
print("Creating data frame to store upper parameter estimates for the no cadmium environment")
    
  cxr_param_w0_upper<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("N"))
  cxr_param_w0_upper$Tu_lambda<-0
  cxr_param_w0_upper$Te_lambda<-0
  cxr_param_w0_upper$Tu_intra<-0
  cxr_param_w0_upper$Te_intra<-0
  cxr_param_w0_upper$Tu_inter<-0
  cxr_param_w0_upper$Te_inter<-0
  
  #removing SR2 for replicate 2
  cxr_param_w0_upper<-cxr_param_w0_upper[-which(cxr_param_w0_upper$Replicate==2 & cxr_param_w0_upper$Tu_Regime=="SR2"),]
  
   # Storing parameters: intrinsic growth rate
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
  
   # Storing parameters: intraspecific competition
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
  
   # Storing parameters: interspecific competition
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
  

#' 
### Saving parameter estimates--------------------------- 
#' saving intermediate files in case something goes wrong
## 

# mean estimates  
write.csv(cxr_param_w0, "./Analyses/cxr_normal_Replicates/parameters_cxr_normal_N.csv")

# upper estimates
  write.csv(cxr_param_w0_upper, "./Analyses/cxr_normal_Replicates/parameters_cxr_normal_upper_N.csv")
  
# Lower estimates
  write.csv(cxr_param_w0_lower, "./Analyses/cxr_normal_Replicates/parameters_cxr_normal_lower_N.csv")


#' 
## Running cxr for cadmium environment---------------------------
#' Here we apply exactly the same reasoning as above but for the data from cadmium environment
## 
print("Running cxr to get parameter estimates for the cadmium environment")  
### Replicate 1----- 
obs.R1_Cd_w0<-cxr_pm_multifit(data = Rep_R1_Cd,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==1 & replicates_initial$Environment=="Cd"),"lambda"],
                                                   alpha_intra = replicates_initial[which(replicates_initial$Replicates==1 & replicates_initial$Environment=="Cd"),"alpha_intra"],
                                                   alpha_inter = replicates_initial[which(replicates_initial$Replicates==1 & replicates_initial$Environment=="Cd"),"alpha_inter"]),
                             fixed_terms = NULL,
                             # no standard errors
                              bootstrap_samples = bootN)
### Replicate 3-----
obs.R3_Cd_w0<-cxr_pm_multifit(data = Rep_R3_Cd,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==3 & replicates_initial$Environment=="Cd"),"lambda"],
                                                   alpha_intra = replicates_initial[which(replicates_initial$Replicates==3 & replicates_initial$Environment=="Cd"),"alpha_intra"],
                                                   alpha_inter = replicates_initial[which(replicates_initial$Replicates==3 & replicates_initial$Environment=="Cd"),"alpha_inter"]),
                             fixed_terms = NULL,
                             # no standard errors
                              bootstrap_samples = bootN)
### Replicate 4-----
obs.R4_Cd_w0<-cxr_pm_multifit(data = Rep_R4_Cd,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==4 & replicates_initial$Environment=="Cd"),"lambda"],
                                                   alpha_intra = replicates_initial[which(replicates_initial$Replicates==4 & replicates_initial$Environment=="Cd"),"alpha_intra"],
                                                   alpha_inter = replicates_initial[which(replicates_initial$Replicates==4 & replicates_initial$Environment=="Cd"),"alpha_inter"]),
                             fixed_terms = NULL,
                             # no standard errors
                              bootstrap_samples = bootN)

### Replicate 5-----
obs.R5_Cd_w0<-cxr_pm_multifit(data = Rep_R5_Cd,
                             focal_column = my.reg,
                             model_family = "RK",
                             covariates = NULL,
                             optimization_method = "Nelder-Mead",
                             alpha_form = "pairwise",
                             lambda_cov_form = "none",
                             alpha_cov_form = "none",
                             initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==5 & replicates_initial$Environment=="Cd"),"lambda"],
                                                   alpha_intra = replicates_initial[which(replicates_initial$Replicates==5 & replicates_initial$Environment=="Cd"),"alpha_intra"],
                                                   alpha_inter = replicates_initial[which(replicates_initial$Replicates==5 & replicates_initial$Environment=="Cd"),"alpha_inter"]),
                             fixed_terms = NULL,
                             # no standard errors
                              bootstrap_samples = bootN)

### Replicate 2-----
# Because there is no SR2 for replicate 2 we need to estimate SR1 and SR4 and SR5 in a different way

obs.R2_Cd_w0_sr1<-cxr_pm_fit(data = Rep_R2_Cd[[1]],
                            focal_column = my.reg[1],
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="Cd"),"lambda"],
                                                  alpha_intra = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="Cd"),"alpha_intra"],
                                                  alpha_inter = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="Cd"),"alpha_inter"]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = bootN)
    
    obs.R2_Cd_w0_sr4<-cxr_pm_fit(data = Rep_R2_Cd[[2]][which(Rep_R2_Cd[[2]][,"SR1"]==0),c("fitness", "SR4")],
                            focal_column =NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "global",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="Cd"),"lambda"],
                                                  alpha_inter = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="Cd"),"alpha_intra"]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = bootN)
    
    obs.R2_Cd_w0_sr4_inter<-cxr_pm_fit(data = Rep_R2_Cd[[2]][which(Rep_R2_Cd[[2]][,"SR1"]!=0),c("fitness", "SR1")],
                                  focal_column =NULL,
                                  model_family = "RK",
                                  covariates = NULL,
                                  optimization_method = "Nelder-Mead",
                                  alpha_form = "global",
                                  lambda_cov_form = "none",
                                  alpha_cov_form = "none",
                                  initial_values = list(alpha_inter=replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="Cd"),"alpha_inter"]),
                                  fixed_terms = list(lambda=obs.R2_Cd_w0_sr4$lambda),
                                  # no standard errors
                                   bootstrap_samples = bootN)
    
    obs.R2_Cd_w0_sr5<-cxr_pm_fit(data = Rep_R2_Cd[[3]][which(Rep_R2_Cd[[3]][,"SR1"]==0),c("fitness", "SR5")],
                            focal_column =NULL,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "global",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="Cd"),"lambda"],
                                                  alpha_inter = replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="Cd"),"alpha_intra"]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = bootN)
    
    obs.R2_Cd_w0_sr5_inter<-cxr_pm_fit(data = Rep_R2_Cd[[3]][which(Rep_R2_Cd[[3]][,"SR1"]!=0),c("fitness", "SR1")],
                                  focal_column =NULL,
                                  model_family = "RK",
                                  covariates = NULL,
                                  optimization_method = "Nelder-Mead",
                                  alpha_form = "global",
                                  lambda_cov_form = "none",
                                  alpha_cov_form = "none",
                                  initial_values = list(alpha_inter=replicates_initial[which(replicates_initial$Replicates==2 & replicates_initial$Environment=="Cd"),"alpha_inter"]),
                                  fixed_terms = list(lambda=obs.R2_Cd_w0_sr5$lambda),
                                  # no standard errors
                                   bootstrap_samples = bootN)


#' 
### Storing parameter estimates (mean, lower and upper)---------------------------
#' Storing the parameter estimates and the lower and upper bounds in a data frame so we can use it later.
#### Mean estimates ------
# creating data frame where to store the data
    
print("Creating data frame to store mean parameter estimates for the cadmium environment")
  cxr_param_w0C<-expand.grid(Tu_Regime=c("SR1", "SR2"), Te_Regime=c("SR4", "SR5"), Replicate=c(1,2,3,4,5), Environment=c("Cd"))
  cxr_param_w0C$Tu_lambda<-0
  cxr_param_w0C$Te_lambda<-0
  cxr_param_w0C$Tu_intra<-0
  cxr_param_w0C$Te_intra<-0
  cxr_param_w0C$Tu_inter<-0
  cxr_param_w0C$Te_inter<-0
  
  #removing SR2 for replicate 2
  cxr_param_w0C<-cxr_param_w0C[-which(cxr_param_w0C$Replicate==2 & cxr_param_w0C$Tu_Regime=="SR2"),]
  
  # Storing parameters: intrinsic growth rate
  cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Tu_lambda"]<-obs.R1_Cd_w0$lambda[c(1,2,1,2)]
  cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Te_lambda"]<-obs.R1_Cd_w0$lambda[c(3,3,4,4)]
  
  cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Tu_lambda"]<-obs.R2_Cd_w0_sr1$lambda
  cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Te_lambda"]<-c(obs.R2_Cd_w0_sr4$lambda, obs.R2_Cd_w0_sr5$lambda)
  
  cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Tu_lambda"]<-obs.R3_Cd_w0$lambda[c(1,2,1,2)]
  cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Te_lambda"]<-obs.R3_Cd_w0$lambda[c(3,3,4,4)]
  
  cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Tu_lambda"]<-obs.R4_Cd_w0$lambda[c(1,2,1,2)]
  cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Te_lambda"]<-obs.R4_Cd_w0$lambda[c(3,3,4,4)]
  
  cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Tu_lambda"]<-obs.R5_Cd_w0$lambda[c(1,2,1,2)]
  cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Te_lambda"]<-obs.R5_Cd_w0$lambda[c(3,3,4,4)]
  
  # Storing parameters: intraspecific competition
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
  
  # Storing parameters: interspecific competition
  cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Tu_inter"]<-c(obs.R1_Cd_w0$alpha_matrix[1,3], obs.R1_Cd_w0$alpha_matrix[2,3],obs.R1_Cd_w0$alpha_matrix[1,4], obs.R1_Cd_w0$alpha_matrix[2,4])
  cxr_param_w0C[which(cxr_param_w0C$Replicate==1),"Te_inter"]<-c(obs.R1_Cd_w0$alpha_matrix[3,1], obs.R1_Cd_w0$alpha_matrix[3,2],obs.R1_Cd_w0$alpha_matrix[4,1], obs.R1_Cd_w0$alpha_matrix[4,2])
  
  cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Tu_inter"]<-obs.R2_Cd_w0_sr1$alpha_inter[1:2]
  cxr_param_w0C[which(cxr_param_w0C$Replicate==2),"Te_inter"]<-c(obs.R2_Cd_w0_sr4_inter$alpha_inter,obs.R2_Cd_w0_sr5_inter$alpha_inter)
  
  cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Tu_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[1,3], obs.R3_Cd_w0$alpha_matrix[2,3],obs.R3_Cd_w0$alpha_matrix[1,4], obs.R3_Cd_w0$alpha_matrix[2,4])
  cxr_param_w0C[which(cxr_param_w0C$Replicate==3),"Te_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[3,1], obs.R3_Cd_w0$alpha_matrix[3,2],obs.R3_Cd_w0$alpha_matrix[4,1], obs.R3_Cd_w0$alpha_matrix[4,2])
  
  cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Tu_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[1,3], obs.R4_Cd_w0$alpha_matrix[2,3],obs.R4_Cd_w0$alpha_matrix[1,4], obs.R4_Cd_w0$alpha_matrix[2,4])
  cxr_param_w0C[which(cxr_param_w0C$Replicate==4),"Te_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[3,1], obs.R4_Cd_w0$alpha_matrix[3,2],obs.R4_Cd_w0$alpha_matrix[4,1], obs.R4_Cd_w0$alpha_matrix[4,2])
  
  cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Tu_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[1,3], obs.R5_Cd_w0$alpha_matrix[2,3],obs.R5_Cd_w0$alpha_matrix[1,4], obs.R5_Cd_w0$alpha_matrix[2,4])
  cxr_param_w0C[which(cxr_param_w0C$Replicate==5),"Te_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[3,1], obs.R5_Cd_w0$alpha_matrix[3,2],obs.R5_Cd_w0$alpha_matrix[4,1], obs.R5_Cd_w0$alpha_matrix[4,2])
  
  
  cxr_param_w0C

#' 
#' #### Data frame with lower estimates
#' 
#### Lower estimates ---------------------------
### Here we will apply the same reasoning but now estimating the lower estimates, using the stats from cxr. Lower bounds correspond to mean-error.
  
  print("Creating data frame to store lower parameter estimates for the cadmium environment")
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
  
  cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==2),"Tu_inter"]<-obs.R2_Cd_w0_sr1$alpha_inter[1:2]-obs.R2_Cd_w0_sr1$alpha_inter_standard_error[1:2]
  cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==2),"Te_inter"]<-c(obs.R2_Cd_w0_sr4_inter$alpha_inter-obs.R2_Cd_w0_sr4_inter$alpha_inter_standard_error, obs.R2_Cd_w0_sr5_inter$alpha_inter-obs.R2_Cd_w0_sr5_inter$alpha_inter_standard_error)
  
  cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==3),"Tu_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[1,3]-obs.R3_Cd_w0$alpha_matrix_standard_error[1,3], obs.R3_Cd_w0$alpha_matrix[2,3]-obs.R3_Cd_w0$alpha_matrix_standard_error[2,3],obs.R3_Cd_w0$alpha_matrix[1,4]-obs.R3_Cd_w0$alpha_matrix_standard_error[1,4], obs.R3_Cd_w0$alpha_matrix[2,4]-obs.R3_Cd_w0$alpha_matrix_standard_error[2,4])
  cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==3),"Te_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[3,1]-obs.R3_Cd_w0$alpha_matrix_standard_error[3,1], obs.R3_Cd_w0$alpha_matrix[3,2]-obs.R3_Cd_w0$alpha_matrix_standard_error[3,2],obs.R3_Cd_w0$alpha_matrix[4,1]-obs.R3_Cd_w0$alpha_matrix_standard_error[4,1], obs.R3_Cd_w0$alpha_matrix[4,2]-obs.R3_Cd_w0$alpha_matrix_standard_error[4,2])
  
  cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==4),"Tu_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[1,3]-obs.R4_Cd_w0$alpha_matrix_standard_error[1,3], obs.R4_Cd_w0$alpha_matrix[2,3]-obs.R4_Cd_w0$alpha_matrix_standard_error[2,3],obs.R4_Cd_w0$alpha_matrix[1,4]-obs.R4_Cd_w0$alpha_matrix_standard_error[1,4], obs.R4_Cd_w0$alpha_matrix[2,4]-obs.R4_Cd_w0$alpha_matrix_standard_error[2,4])
  cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==4),"Te_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[3,1]-obs.R4_Cd_w0$alpha_matrix_standard_error[3,1], obs.R4_Cd_w0$alpha_matrix[3,2]-obs.R4_Cd_w0$alpha_matrix_standard_error[3,2],obs.R4_Cd_w0$alpha_matrix[4,1]-obs.R4_Cd_w0$alpha_matrix_standard_error[4,1], obs.R4_Cd_w0$alpha_matrix[4,2]-obs.R4_Cd_w0$alpha_matrix_standard_error[4,2])
  
  cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==5),"Tu_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[1,3]-obs.R5_Cd_w0$alpha_matrix_standard_error[1,3], obs.R5_Cd_w0$alpha_matrix[2,3]-obs.R5_Cd_w0$alpha_matrix_standard_error[2,3],obs.R5_Cd_w0$alpha_matrix[1,4]-obs.R5_Cd_w0$alpha_matrix_standard_error[1,4], obs.R5_Cd_w0$alpha_matrix[2,4]-obs.R5_Cd_w0$alpha_matrix_standard_error[2,4])
  cxr_param_w0C_lower[which(cxr_param_w0C_lower$Replicate==5),"Te_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[3,1]-obs.R5_Cd_w0$alpha_matrix_standard_error[3,1], obs.R5_Cd_w0$alpha_matrix[3,2]-obs.R5_Cd_w0$alpha_matrix_standard_error[3,2],obs.R5_Cd_w0$alpha_matrix[4,1]-obs.R5_Cd_w0$alpha_matrix_standard_error[4,1], obs.R5_Cd_w0$alpha_matrix[4,2]-obs.R5_Cd_w0$alpha_matrix_standard_error[4,2])

#' 
#' #### Data frame with upper estimates
#' 
#### Upper estimates ---------------------------
### Here we will apply the same reasoning but now estimating the upper estimates, using the stats from cxr. Upper bounds correspond to mean-error.
print("Creating data frame to store upper parameter estimates for the cadmium environment")
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
  
  cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==2),"Tu_inter"]<-obs.R2_Cd_w0_sr1$alpha_inter[1:2]+obs.R2_Cd_w0_sr1$alpha_inter_standard_error[1:2]
  cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==2),"Te_inter"]<-c(obs.R2_Cd_w0_sr4_inter$alpha_inter+obs.R2_Cd_w0_sr4_inter$alpha_inter_standard_error, obs.R2_Cd_w0_sr5_inter$alpha_inter+obs.R2_Cd_w0_sr5_inter$alpha_inter_standard_error)
  
  cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==3),"Tu_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[1,3]+obs.R3_Cd_w0$alpha_matrix_standard_error[1,3], obs.R3_Cd_w0$alpha_matrix[2,3]+obs.R3_Cd_w0$alpha_matrix_standard_error[2,3],obs.R3_Cd_w0$alpha_matrix[1,4]+obs.R3_Cd_w0$alpha_matrix_standard_error[1,4], obs.R3_Cd_w0$alpha_matrix[2,4]+obs.R3_Cd_w0$alpha_matrix_standard_error[2,4])
  cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==3),"Te_inter"]<-c(obs.R3_Cd_w0$alpha_matrix[3,1]+obs.R3_Cd_w0$alpha_matrix_standard_error[3,1], obs.R3_Cd_w0$alpha_matrix[3,2]+obs.R3_Cd_w0$alpha_matrix_standard_error[3,2],obs.R3_Cd_w0$alpha_matrix[4,1]+obs.R3_Cd_w0$alpha_matrix_standard_error[4,1], obs.R3_Cd_w0$alpha_matrix[4,2]+obs.R3_Cd_w0$alpha_matrix_standard_error[4,2])
  
  cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==4),"Tu_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[1,3]+obs.R4_Cd_w0$alpha_matrix_standard_error[1,3], obs.R4_Cd_w0$alpha_matrix[2,3]+obs.R4_Cd_w0$alpha_matrix_standard_error[2,3],obs.R4_Cd_w0$alpha_matrix[1,4]+obs.R4_Cd_w0$alpha_matrix_standard_error[1,4], obs.R4_Cd_w0$alpha_matrix[2,4]+obs.R4_Cd_w0$alpha_matrix_standard_error[2,4])
  cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==4),"Te_inter"]<-c(obs.R4_Cd_w0$alpha_matrix[3,1]+obs.R4_Cd_w0$alpha_matrix_standard_error[3,1], obs.R4_Cd_w0$alpha_matrix[3,2]+obs.R4_Cd_w0$alpha_matrix_standard_error[3,2],obs.R4_Cd_w0$alpha_matrix[4,1]+obs.R4_Cd_w0$alpha_matrix_standard_error[4,1], obs.R4_Cd_w0$alpha_matrix[4,2]+obs.R4_Cd_w0$alpha_matrix_standard_error[4,2])
  
  cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==5),"Tu_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[1,3]+obs.R5_Cd_w0$alpha_matrix_standard_error[1,3], obs.R5_Cd_w0$alpha_matrix[2,3]+obs.R5_Cd_w0$alpha_matrix_standard_error[2,3],obs.R5_Cd_w0$alpha_matrix[1,4]+obs.R5_Cd_w0$alpha_matrix_standard_error[1,4], obs.R5_Cd_w0$alpha_matrix[2,4]+obs.R5_Cd_w0$alpha_matrix_standard_error[2,4])
  cxr_param_w0C_upper[which(cxr_param_w0C_upper$Replicate==5),"Te_inter"]<-c(obs.R5_Cd_w0$alpha_matrix[3,1]+obs.R5_Cd_w0$alpha_matrix_standard_error[3,1], obs.R5_Cd_w0$alpha_matrix[3,2]+obs.R5_Cd_w0$alpha_matrix_standard_error[3,2],obs.R5_Cd_w0$alpha_matrix[4,1]+obs.R5_Cd_w0$alpha_matrix_standard_error[4,1], obs.R5_Cd_w0$alpha_matrix[4,2]+obs.R5_Cd_w0$alpha_matrix_standard_error[4,2])
  

#' 
#' 
### Saving pooled parameter estimates ---------------------------
#' intermediate files in case something goes wrong
## 
## Save the objects
save(obs.R1_Cd_w0,file= "./Analyses/cxr_normal_Replicates/cxr_R1_Cd.RData")
save(obs.R3_Cd_w0, file="./Analyses/cxr_normal_Replicates/cxr_R3_Cd.RData")
save(obs.R4_Cd_w0, file="./Analyses/cxr_normal_Replicates/cxr_R4_Cd.RData")
save(obs.R5_Cd_w0, file="./Analyses/cxr_normal_Replicates/cxr_R5_Cd.RData")
save(obs.R2_Cd_w0_sr1, file="./Analyses/cxr_normal_Replicates/cxr_R2_SR1_Cd.RData")
save(obs.R2_Cd_w0_sr4,file= "./Analyses/cxr_normal_Replicates/cxr_R2_SR4_Cd.RData")
save(obs.R2_Cd_w0_sr5, file="./Analyses/cxr_normal_Replicates/cxr_R2_SR5_Cd.RData")
save(obs.R2_Cd_w0_sr4_inter, file="./Analyses/cxr_normal_Replicates/cxr_R2_SR4_inter_Cd.RData")
save(obs.R2_Cd_w0_sr5_inter, file="./Analyses/cxr_normal_Replicates/cxr_R2_SR5_inter_Cd.RData")

#' 
#' 
##Save parameter estimates ---------------------------

# Mean estimates
write.csv(cxr_param_w0C, file = "./Analyses/cxr_normal_Replicates/parameters_cxr_normal_Cd.csv")
# Upper estimates
  write.csv(cxr_param_w0C_upper, file = "./Analyses/cxr_normal_Replicates/parameters_cxr_normal_upper_Cd.csv")

# Lower estimates
  write.csv(cxr_param_w0C_lower, file = "./Analyses/cxr_normal_Replicates/parameters_cxr_normal_lower_Cd.csv")

print("Joining data frames")
#' 
# Joining data frames---------------------------

## 
param_all_w0<-as.data.frame(rbind(cxr_param_w0, cxr_param_w0C))

param_all_w0_lower<-as.data.frame(rbind(cxr_param_w0_lower, cxr_param_w0C_lower))
param_all_w0_upper<-as.data.frame(rbind(cxr_param_w0_upper, cxr_param_w0C_upper))

print("Saving files")
# Save files--------------------------- 
#' These are the final files that are then used for the rest of the analyses
  write.csv(param_all_w0, "./Analyses/cxr_normal_Replicates/parameters_cxr_normal.csv")
  write.csv(param_all_w0_upper, "./Analyses/cxr_normal_Replicates/parameters_cxr_normal_upper.csv")
  write.csv(param_all_w0_lower, "./Analyses/cxr_normal_Replicates/parameters_cxr_normal_lower.csv")

#' 
