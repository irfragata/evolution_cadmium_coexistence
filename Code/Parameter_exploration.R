  #' ---
  #' title: "R Notebook"
  #' output: html_notebook
  #' ---
  #' 
 
  #' 
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

# Functions, general information and packages------
  if(!require(dplyr)){
    install.packages("dplyr")
  }
  if(!require(tidyr)){
    install.packages("tidyr")
  }
  if(!require(cxr)){
    install.packages("cxr")
  }
  if(!require(MASS)){
    install.packages("MASS")
  }
  if(!require(plyr)){
    install.packages("plyr")
  }
  
  library(dplyr)
  library(tidyr)
  library(cxr)
  library(MASS)
  library(plyr)
  
  rssq_function<-function(df, lambda, alpha_intra, alpha_inter,... ){
    pred<-sapply(c(1:nrow(df)), function(a){
      m<-lambda*exp(-alpha_intra*df[a,2]-alpha_inter*df[a,3])
      # Estimating the diff between predicted and observed
      re<-(m-df[a,1])^2
      
      re
    })
    return(sum(pred))
  }
  #' 
  #' WARNING: This code takes a large amount of time to run. It is provided not as part of the main analyses, but to show how the initial parameter values were selected.
  #' 
  # Importing data------
  #' This chunk is used to import data, checking data and creating the columns that are necessary
  #' 
  
  
  ca<-read.csv(file = "./Data/CompetitiveAbility_Cd_G40_submit.csv", header=TRUE) # cdata from the competitive ability
  
  str(ca) 
  # Summary of the data to be sure that everything is ok!
  summary(as.factor(ca$Foca_rawlSR))
  
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
  #' # Evaluation
  #' Some pieces of the code take a lot of time to run. We provide the intermediate results to speed up the process. To run the whole code change eval from FALSE to TRUE.
  #' 
  # Evaluation------------------------------------
  
  evaluation <- FALSE
  evaluation2<- FALSE
  if(!file.exists("./Analyses")){
    dir.create("./Analyses/")
  }
  
  #' 
  # Finding the best optim values------
  #' 
  #' When using the cxr package we detected at least two likelihood peaks in some selection regime.
  #' We will perform a grid search with variation of the initial parameter values.
  #' For each combination of values we will:
  #' - Estimate the loglikelihood
  #' - Estimate RMSE between predicted and observed values
  #' 
  #' To make it faster no bootstrap samples will be run
  #' 
  
  set.seed(42)
  bootN<-0
  lambda_range<-seq(0.0001,2,0.1)
  alpha_range<-seq(-1,1, 0.1)
  
  
  mat_ini_val<-as.data.frame(expand_grid(lambda=lambda_range, alpha_intra=alpha_range, alpha_inter=alpha_range))
  

  #' 
  #' # Set up data
  #' 
  ## Creating a data frame for the no cadmium environment ------------------------------------
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
  ## Creating a data frame for the cadmium environment------------------------------------
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
  
    # Do list per replicate and environment
    Rep<-list(SR1= subset(forCXR_N, Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N,  Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  #' 
  #'  ### Data for Pooled estimates for cadmium
  #' 
  
  # Do list per replicate and environment
    Rep_Cd<-list(SR1= subset(forCXR_Cd, Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd,  Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  #' 
  #' ### Data for replicates no cadmium
  
  
  # Do list per replicate and environment
    R1<-list(SR1= subset(forCXR_N, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
    
    R2<-list(SR1= subset(forCXR_N, Rep==2 & Focal=="SR1")[,c("fitness", "SR1", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR4")], SR5= subset(forCXR_N, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR5")])
    
    R3<-list(SR1= subset(forCXR_N, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
    
    R4<-list(SR1= subset(forCXR_N, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
    
    R5<-list(SR1= subset(forCXR_N, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_N, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_N, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_N, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
    
  
  #' 
  #' ### Data for replicates cadmium
  #' 
  
  Rep_R1_Cd<-list(SR1= subset(forCXR_Cd, Rep==1 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==1 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==1 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==1 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
    
    Rep_R2_Cd<-list(SR1= subset(forCXR_Cd, Rep==2 & Focal=="SR1")[,c("fitness", "SR1","SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==2 & Focal=="SR4")[,c("fitness", "SR1", "SR4")], SR5= subset(forCXR_Cd, Rep==2 & Focal=="SR5")[,c("fitness", "SR1", "SR5")])
    
    Rep_R3_Cd<-list(SR1= subset(forCXR_Cd, Rep==3 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==3 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==3 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==3 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
    
    Rep_R4_Cd<-list(SR1= subset(forCXR_Cd, Rep==4 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==4 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==4 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==4 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
    
    Rep_R5_Cd<-list(SR1= subset(forCXR_Cd, Rep==5 & Focal=="SR1")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR2= subset(forCXR_Cd, Rep==5 & Focal=="SR2")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR4= subset(forCXR_Cd, Rep==5 & Focal=="SR4")[,c("fitness", "SR1", "SR2", "SR4", "SR5")], SR5= subset(forCXR_Cd, Rep==5 & Focal=="SR5")[,c("fitness", "SR1", "SR2", "SR4", "SR5")])
  
  #' 
  # Running parameter grid exploration------
  #' 
  ### Pooled no cadmium ------------------------------------
if(evaluation){
  print("Running cxr for pooled no cadmium data")
    
  aux_N<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = Rep,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_N)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  
  pooled_N<-as.data.frame(cbind(mat_ini_val, aux_N))
  
  #save(pooled_N, file="./Pooled_N.RData")
  
  print("Running cxr for pooled cadmium data")
  #' 
  ### Pooled Cadmium------------------------------------
  aux_Cd<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = Rep_Cd,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_Cd)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  str(aux_Cd)
  
  pooled_Cd<-as.data.frame(cbind(mat_ini_val, aux_Cd))
  #save(pooled_Cd, file="./Pooled_Cd.RData")
  
  
  #' 
  ### Replicates no cadmium------
  #' 
  ##### R1 N------------------------------------
  
  print("Running cxr for replicate no cadmium data - R1")
  
  aux_R1_N<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = R1,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_R1_N)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  str(aux_R1_N)
  
  R1_N<-as.data.frame(cbind(mat_ini_val, aux_R1_N))
  
  #save(R1_N, file="./R1_N.RData")
  
  
  #' 
  ##### R3 N------------------------------------
  print("Running cxr for replicate no cadmium data - R3")
  
  aux_R3_N<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = R3,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_R3_N)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  str(aux_R3_N)
  
  R3_N<-as.data.frame(cbind(mat_ini_val, aux_R3_N))
  
  #save(R3_N, file="./R3_N.RData")
  
  
  #' 
  ##### R4 N ------------------------------------
  print("Running cxr for replicate no cadmium data - R4")
  
  aux_R4_N<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = R4,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_R4_N)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  str(aux_R4_N)
  
  R4_N<-as.data.frame(cbind(mat_ini_val, aux_R4_N))
  
  #save(R4_N, file="./R4_N.RData")
  
  
  #' 
  ##### R5 N------------------------------------
  print("Running cxr for replicate no cadmium data - R5")
  
  aux_R5_N<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = R5,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_R5_N)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  str(aux_R5_N)
  
  R5_N<-as.data.frame(cbind(mat_ini_val, aux_R5_N))
  
  #save(R5_N, file="./R5_N.RData")
  
  
  #' 
  ##### R2 N-------------------
  print("Running cxr for replicate no cadmium data - R2")
  
  aux_R2_N<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.R2_w0_sr1<-cxr_pm_fit(data = R2[[1]],
                              focal_column = my.reg[1],
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = mat_ini_val$lambda[x],
                                                    alpha_intra = mat_ini_val$alpha_intra[x],
                                                    alpha_inter = mat_ini_val$alpha_inter[x]),
                              fixed_terms = NULL,
                              # no standard errors
                               bootstrap_samples = 0)
      
      obs.R2_w0_sr4<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]==0),c("fitness", "SR4")],
                              focal_column =NULL,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "global",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = mat_ini_val$lambda[x],
                                                    alpha_inter = mat_ini_val$alpha_intra[x]),
                              fixed_terms = NULL,
                              # no standard errors
                               bootstrap_samples = 0)
      
      obs.R2_w0_sr4_inter<-cxr_pm_fit(data = R2[[2]][which(R2[[2]][,"SR1"]!=0),c("fitness", "SR1")],
                                    focal_column =NULL,
                                    model_family = "RK",
                                    covariates = NULL,
                                    optimization_method = "Nelder-Mead",
                                    alpha_form = "global",
                                    lambda_cov_form = "none",
                                    alpha_cov_form = "none",
                                    initial_values = list(alpha_inter=mat_ini_val$alpha_inter[x]),
                                    fixed_terms = list(lambda=obs.R2_w0_sr4$lambda),
                                    # no standard errors
                                     bootstrap_samples = 0)
      
      obs.R2_w0_sr5<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]==0),c("fitness", "SR5")],
                              focal_column =NULL,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "global",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = mat_ini_val$lambda[x],
                                                    alpha_inter = mat_ini_val$alpha_intra[x]),
                              fixed_terms = NULL,
                              # no standard errors
                               bootstrap_samples = 0)
      
      obs.R2_w0_sr5_inter<-cxr_pm_fit(data = R2[[3]][which(R2[[3]][,"SR1"]!=0),c("fitness", "SR1")],
                                    focal_column =NULL,
                                    model_family = "RK",
                                    covariates = NULL,
                                    optimization_method = "Nelder-Mead",
                                    alpha_form = "global",
                                    lambda_cov_form = "none",
                                    alpha_cov_form = "none",
                                    initial_values = list(alpha_inter=mat_ini_val$alpha_inter[x]),
                                    fixed_terms = list(lambda=obs.R2_w0_sr5$lambda),
                                    # no standard errors
                                     bootstrap_samples = 0)
      
      auxi<-c(obs.R2_w0_sr1$lambda[1],NA,obs.R2_w0_sr4$lambda[1],obs.R2_w0_sr5$lambda[1], obs.R2_w0_sr1$alpha_intra[1],NA,obs.R2_w0_sr4$alpha_inter[1],obs.R2_w0_sr5$alpha_inter[1], obs.R2_w0_sr1$alpha_inter[1], obs.R2_w0_sr1$alpha_inter[2], NA, NA, obs.R2_w0_sr4_inter$alpha_inter[1], NA, obs.R2_w0_sr5_inter$alpha_inter[1], NA, obs.R2_w0_sr1$log_likelihood[1], NA, obs.R2_w0_sr4$log_likelihood[1], obs.R2_w0_sr4_inter$log_likelihood[1], obs.R2_w0_sr5$log_likelihood[1], obs.R2_w0_sr5_inter$log_likelihood[1])
  
      auxi
      })))
  
  colnames(aux_R2_N)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR4inter", "LL_SR5", "LL_SR5inter")
  
  str(aux_R2_N)
  
  R2_N<-as.data.frame(cbind(mat_ini_val, aux_R2_N))
  #save(R2_N, file="./R2_N.RData")
  
  #' 
  #' 
  ### Replicates no cadmium------
  #' 
  ##### R1 Cd------------------------------------
  print("Running cxr for replicate cadmium data - R1")
  
  aux_R1_Cd<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = Rep_R1_Cd,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_R1_Cd)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  str(aux_R1_Cd)
  
  R1_Cd<-as.data.frame(cbind(mat_ini_val, aux_R1_Cd))
  #save(R1_Cd, file="./R1_Cd.RData")
  
  #' 
  ##### R3 Cd------------------------------------
  print("Running cxr for replicate cadmium data - R3")
  
  aux_R3_Cd<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = Rep_R3_Cd,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_R3_Cd)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  str(aux_R3_Cd)
  
  R3_Cd<-as.data.frame(cbind(mat_ini_val, aux_R3_Cd))
  #save(R3_Cd, file="./R3_Cd.RData")
  
  #' 
  ##### R4 Cd ------------------------------------
  print("Running cxr for replicate cadmium data - R4")
  
  aux_R4_Cd<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = Rep_R4_Cd,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_R4_Cd)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  str(aux_R4_Cd)
  
  R4_Cd<-as.data.frame(cbind(mat_ini_val, aux_R4_Cd))
  #save(R4_Cd, file="./R4_Cd.RData")
  
  #' 
  ##### R5 Cd------------------------------------
  print("Running cxr for replicate cadmium data - R5")
  
  aux_R5_Cd<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.w0_0<-cxr_pm_multifit(data = Rep_R5_Cd,
                            focal_column = my.reg,
                            model_family = "RK",
                            covariates = NULL,
                            optimization_method = "Nelder-Mead",
                            alpha_form = "pairwise",
                            lambda_cov_form = "none",
                            alpha_cov_form = "none",
                            initial_values = list(lambda = mat_ini_val$lambda[x],
                                                  alpha_intra = mat_ini_val$alpha_intra[x],
                                                  alpha_inter = mat_ini_val$alpha_inter[x]),
                            fixed_terms = NULL,
                            # no standard errors
                             bootstrap_samples = 0)
      
      auxi<-c(obs.w0_0$lambda[1],obs.w0_0$lambda[2],obs.w0_0$lambda[3],obs.w0_0$lambda[4], obs.w0_0$alpha_matrix[1,1],obs.w0_0$alpha_matrix[2,2],obs.w0_0$alpha_matrix[3,3],obs.w0_0$alpha_matrix[4,4],obs.w0_0$alpha_matrix[1,3], obs.w0_0$alpha_matrix[1,4], obs.w0_0$alpha_matrix[2,3], obs.w0_0$alpha_matrix[2,4], obs.w0_0$alpha_matrix[3,1], obs.w0_0$alpha_matrix[3,2], obs.w0_0$alpha_matrix[4,1], obs.w0_0$alpha_matrix[4,2], obs.w0_0$log_likelihood[1], obs.w0_0$log_likelihood[2], obs.w0_0$log_likelihood[3], obs.w0_0$log_likelihood[4])
  
      auxi
      })))
  
  colnames(aux_R5_Cd)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR5")
  
  str(aux_R5_Cd)
  
  R5_Cd<-as.data.frame(cbind(mat_ini_val, aux_R5_Cd))
  #save(R5_Cd, file="./R5_Cd.RData")
  
  #' 
  ##### R2 Cd-------------------
  print("Running cxr for replicate cadmium data - R2")
  
  aux_R2_Cd<-as.data.frame(t(sapply(c(1:nrow(mat_ini_val)), function(x){
    
      obs.R2_w0_sr1<-cxr_pm_fit(data = Rep_R2_Cd[[1]],
                              focal_column = my.reg[1],
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "pairwise",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = mat_ini_val$lambda[x],
                                                    alpha_intra = mat_ini_val$alpha_intra[x],
                                                    alpha_inter = mat_ini_val$alpha_inter[x]),
                              fixed_terms = NULL,
                              # no standard errors
                               bootstrap_samples = 0)
      
      obs.R2_w0_sr4<-cxr_pm_fit(data = Rep_R2_Cd[[2]][which(Rep_R2_Cd[[2]][,"SR1"]==0),c("fitness", "SR4")],
                              focal_column =NULL,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "global",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = mat_ini_val$lambda[x],
                                                    alpha_inter = mat_ini_val$alpha_intra[x]),
                              fixed_terms = NULL,
                              # no standard errors
                               bootstrap_samples = 0)
      
      obs.R2_w0_sr4_inter<-cxr_pm_fit(data = Rep_R2_Cd[[2]][which(Rep_R2_Cd[[2]][,"SR1"]!=0),c("fitness", "SR1")],
                                    focal_column =NULL,
                                    model_family = "RK",
                                    covariates = NULL,
                                    optimization_method = "Nelder-Mead",
                                    alpha_form = "global",
                                    lambda_cov_form = "none",
                                    alpha_cov_form = "none",
                                    initial_values = list(alpha_inter=mat_ini_val$alpha_inter[x]),
                                    fixed_terms = list(lambda=obs.R2_w0_sr4$lambda),
                                    # no standard errors
                                     bootstrap_samples = 0)
      
      obs.R2_w0_sr5<-cxr_pm_fit(data = Rep_R2_Cd[[3]][which(Rep_R2_Cd[[3]][,"SR1"]==0),c("fitness", "SR5")],
                              focal_column =NULL,
                              model_family = "RK",
                              covariates = NULL,
                              optimization_method = "Nelder-Mead",
                              alpha_form = "global",
                              lambda_cov_form = "none",
                              alpha_cov_form = "none",
                              initial_values = list(lambda = mat_ini_val$lambda[x],
                                                    alpha_inter = mat_ini_val$alpha_intra[x]),
                              fixed_terms = NULL,
                              # no standard errors
                               bootstrap_samples = 0)
      
      obs.R2_w0_sr5_inter<-cxr_pm_fit(data = Rep_R2_Cd[[3]][which(Rep_R2_Cd[[3]][,"SR1"]!=0),c("fitness", "SR1")],
                                    focal_column =NULL,
                                    model_family = "RK",
                                    covariates = NULL,
                                    optimization_method = "Nelder-Mead",
                                    alpha_form = "global",
                                    lambda_cov_form = "none",
                                    alpha_cov_form = "none",
                                    initial_values = list(alpha_inter=mat_ini_val$alpha_inter[x]),
                                    fixed_terms = list(lambda=obs.R2_w0_sr5$lambda),
                                    # no standard errors
                                     bootstrap_samples = 0)
      
      auxi<-c(obs.R2_w0_sr1$lambda[1],NA,obs.R2_w0_sr4$lambda[1],obs.R2_w0_sr5$lambda[1], obs.R2_w0_sr1$alpha_intra[1],NA,obs.R2_w0_sr4$alpha_inter[1],obs.R2_w0_sr5$alpha_inter[1], obs.R2_w0_sr1$alpha_inter[1], obs.R2_w0_sr1$alpha_inter[2], NA, NA, obs.R2_w0_sr4_inter$alpha_inter[1], NA, obs.R2_w0_sr5_inter$alpha_inter[1], NA, obs.R2_w0_sr1$log_likelihood[1], NA, obs.R2_w0_sr4$log_likelihood[1], obs.R2_w0_sr4_inter$log_likelihood[1], obs.R2_w0_sr5$log_likelihood[1], obs.R2_w0_sr5_inter$log_likelihood[1])
  
      auxi
      })))
  
  colnames(aux_R2_Cd)<-c("lambda_SR1","lambda_SR2", "lambda_SR4", "lambda_SR5", "intra_SR1", "intra_SR2", "intra_SR4", "intra_SR5", "inter_SR1.SR4", "inter_SR1.SR5", "inter_SR2.SR4", "inter_SR2.SR5", "inter_SR4.SR1", "inter_SR4.SR2", "inter_SR5.SR1", "inter_SR5.SR2", "LL_SR1", "LL_SR2", "LL_SR4", "LL_SR4inter", "LL_SR5", "LL_SR5inter")
  
  str(aux_R2_Cd)
  
  R2_Cd<-as.data.frame(cbind(mat_ini_val, aux_R2_Cd))
  #save(R2_Cd, file="./R2_Cd.RData")
  
  ## Saving the information ------------------------------------
  save.image(file='./Analyses/Parameter_exploration.RData')
  
    }else{
    load(file='./Analyses/Parameter_exploration.RData') 
  }
  
  
  #' 
  #' 
  # Model selection------
  #' 
  ## Pooled data------
  #' 
  ### Analysing likelihood information ------------------------------------
  print("Estimate mean likelihood for pooled data")
  
  
    pooled_N$meanLL<-sapply(c(1:nrow(pooled_N)), function(x) mean(c(pooled_N$LL_SR1[x], pooled_N$LL_SR2[x],pooled_N$LL_SR4[x], pooled_N$LL_SR5[x]), na.rm=TRUE))
  
  pooled_N$varLL<-sapply(c(1:nrow(pooled_N)), function(x) sd(c(pooled_N$LL_SR1[x], pooled_N$LL_SR2[x],pooled_N$LL_SR4[x], pooled_N$LL_SR5[x]), na.rm=TRUE))
  
  pooled_Cd$meanLL<-sapply(c(1:nrow(pooled_Cd)), function(x) mean(c(pooled_Cd$LL_SR1[x], pooled_Cd$LL_SR2[x],pooled_Cd$LL_SR4[x], pooled_Cd$LL_SR5[x]), na.rm=TRUE))
  
  pooled_Cd$varLL<-sapply(c(1:nrow(pooled_Cd)), function(x) sd(c(pooled_Cd$LL_SR1[x], pooled_Cd$LL_SR2[x],pooled_Cd$LL_SR4[x], pooled_Cd$LL_SR5[x]), na.rm=TRUE))
  
  # ordering
  pooled_N[order(pooled_N$meanLL, decreasing=FALSE)[1:20],]
  
  # another option
  subset(pooled_N, meanLL<=(min(pooled_N$meanLL)+2))
  
  colnames(pooled_N)[6:7]<-c("lambda_SR4","lambda_SR5")
  colnames(pooled_Cd)[6:7]<-c("lambda_SR4","lambda_SR5")
  #' 
  # Estimating root mean squared error (RMSE) ------------------------------------
  # Creates a function that estimates the residual mean sum of squares between the predicted data and the observed
print(evaluation)
 
if(!evaluation2){
  load(file='./Analyses/Parameter_exploration_RMSE.RData')
}else{
  print("Estimate root mean square error for pooled data")
  
  pooled_N$RMSE_SR1<-sapply(c(1:nrow(pooled_N)), function(x){
    rssq_function(Rep[[1]][,c("fitness","SR1","SR4")], pooled_N$lambda_SR1[x], pooled_N$intra_SR1[x],pooled_N$inter_SR1.SR4[x])+rssq_function(Rep[[1]][,c("fitness","SR1","SR5")], pooled_N$lambda_SR1[x], pooled_N$intra_SR1[x],pooled_N$inter_SR1.SR5[x])
  })
  
  pooled_N$RMSE_SR2<-sapply(c(1:nrow(pooled_N)), function(x){
    rssq_function(Rep[[2]][,c("fitness","SR2","SR4")], pooled_N$lambda_SR2[x], pooled_N$intra_SR2[x],pooled_N$inter_SR2.SR4[x])+rssq_function(Rep[[2]][,c("fitness","SR2","SR5")], pooled_N$lambda_SR2[x], pooled_N$intra_SR2[x],pooled_N$inter_SR2.SR5[x])
  })
  
  pooled_N$RMSE_SR4<-sapply(c(1:nrow(pooled_N)), function(x){
    rssq_function(Rep[[3]][,c("fitness","SR4","SR1")], pooled_N$lambda_SR4[x], pooled_N$intra_SR4[x],pooled_N$inter_SR4.SR1[x])+rssq_function(Rep[[3]][,c("fitness","SR4","SR1")], pooled_N$lambda_SR4[x], pooled_N$intra_SR4[x],pooled_N$inter_SR4.SR2[x])
  })
  
  pooled_N$RMSE_SR5<-sapply(c(1:nrow(pooled_N)), function(x){
    rssq_function(Rep[[4]][,c("fitness","SR5","SR1")], pooled_N$lambda_SR5[x], pooled_N$intra_SR5[x],pooled_N$inter_SR5.SR1[x])+rssq_function(Rep[[4]][,c("fitness","SR5","SR2")], pooled_N$lambda_SR5[x], pooled_N$intra_SR5[x],pooled_N$inter_SR5.SR2[x])
  })
  
  #' ## Pooled Cd
  
  pooled_Cd$RMSE_SR1<-sapply(c(1:nrow(pooled_Cd)), function(x){
    rssq_function(Rep_Cd[[1]][,c("fitness","SR1","SR4")], pooled_Cd$lambda_SR1[x], pooled_Cd$intra_SR1[x],pooled_Cd$inter_SR1.SR4[x])+rssq_function(Rep_Cd[[1]][,c("fitness","SR1","SR5")], pooled_Cd$lambda_SR1[x], pooled_Cd$intra_SR1[x],pooled_Cd$inter_SR1.SR5[x])
  })
  
  pooled_Cd$RMSE_SR2<-sapply(c(1:nrow(pooled_Cd)), function(x){
    rssq_function(Rep_Cd[[2]][,c("fitness","SR2","SR4")], pooled_Cd$lambda_SR2[x], pooled_Cd$intra_SR2[x],pooled_Cd$inter_SR2.SR4[x])+rssq_function(Rep_Cd[[2]][,c("fitness","SR2","SR5")], pooled_Cd$lambda_SR2[x], pooled_Cd$intra_SR2[x],pooled_Cd$inter_SR2.SR5[x])
  })
  
  pooled_Cd$RMSE_SR4<-sapply(c(1:nrow(pooled_Cd)), function(x){
    rssq_function(Rep_Cd[[3]][,c("fitness","SR4","SR1")], pooled_Cd$lambda_SR4[x], pooled_Cd$intra_SR4[x],pooled_Cd$inter_SR4.SR1[x])+rssq_function(Rep_Cd[[3]][,c("fitness","SR4","SR1")], pooled_Cd$lambda_SR4[x], pooled_Cd$intra_SR4[x],pooled_Cd$inter_SR4.SR2[x])
  })
  
  pooled_Cd$RMSE_SR5<-sapply(c(1:nrow(pooled_Cd)), function(x){
    rssq_function(Rep_Cd[[4]][,c("fitness","SR5","SR1")], pooled_Cd$lambda_SR5[x], pooled_Cd$intra_SR5[x],pooled_Cd$inter_SR5.SR1[x])+rssq_function(Rep_Cd[[4]][,c("fitness","SR5","SR2")], pooled_Cd$lambda_SR5[x], pooled_Cd$intra_SR5[x],pooled_Cd$inter_SR5.SR2[x])
  })
  
  #' Estimate the sum or mean of RMSE
  
  pooled_N$Sum_RMSE<-sapply(c(1:nrow(pooled_N)), function(x) sum(c(pooled_N$RMSE_SR1[x],pooled_N$RMSE_SR2[x],pooled_N$RMSE_SR4[x],pooled_N$RMSE_SR5[x])))
  
  pooled_N$Mean_RMSE<-sapply(c(1:nrow(pooled_N)), function(x) mean(c(pooled_N$RMSE_SR1[x],pooled_N$RMSE_SR2[x],pooled_N$RMSE_SR4[x],pooled_N$RMSE_SR5[x])))
  
  pooled_Cd$Sum_RMSE<-sapply(c(1:nrow(pooled_Cd)), function(x) sum(c(pooled_Cd$RMSE_SR1[x],pooled_Cd$RMSE_SR2[x],pooled_Cd$RMSE_SR4[x],pooled_Cd$RMSE_SR5[x])))
  
  pooled_Cd$Mean_RMSE<-sapply(c(1:nrow(pooled_Cd)), function(x) mean(c(pooled_Cd$RMSE_SR1[x],pooled_Cd$RMSE_SR2[x],pooled_Cd$RMSE_SR4[x],pooled_Cd$RMSE_SR5[x])))
  
  #' Using the mean or the sum renders the same best
  #' However the estimates for alpha and lambda show a large variation, specially in the cadmium environment.
  head(pooled_N[order(pooled_N$Sum_RMSE, decreasing=FALSE),], n=10)
  head(pooled_N[order(pooled_N$Mean_RMSE, decreasing=FALSE),], n=10)
  
  head(pooled_Cd[order(pooled_Cd$Sum_RMSE, decreasing=FALSE),], n=10)
  head(pooled_Cd[order(pooled_Cd$Mean_RMSE, decreasing=FALSE),], n=10)
  
  
  #' 
  # Estimating mean likelihood and RMSE for each Replicate ------------------------------------
  
  # colnames(R1_N)[6:7]<-c("lambda_SR4","lambda_SR5")
  # colnames(R1_Cd)[6:7]<-c("lambda_SR4","lambda_SR5")
  # colnames(R2_N)[6:7]<-c("lambda_SR4","lambda_SR5")
  # colnames(R2_Cd)[6:7]<-c("lambda_SR4","lambda_SR5")
  # colnames(R3_N)[6:7]<-c("lambda_SR4","lambda_SR5")
  # colnames(R3_Cd)[6:7]<-c("lambda_SR4","lambda_SR5")
  # colnames(R4_N)[6:7]<-c("lambda_SR4","lambda_SR5")
  # colnames(R4_Cd)[6:7]<-c("lambda_SR4","lambda_SR5")
  # colnames(R5_N)[6:7]<-c("lambda_SR4","lambda_SR5")
  # colnames(R5_Cd)[6:7]<-c("lambda_SR4","lambda_SR5")
  
  print("Estimate mean likelihood for replicate data")
  
  # Estimating mean log likelihood for the different replicates
  R1_N$meanLL<-sapply(c(1:nrow(R1_N)), function(x) mean(c(R1_N$LL_SR1[x], R1_N$LL_SR2[x],R1_N$LL_SR4[x], R1_N$LL_SR5[x]), na.rm=TRUE))
  
  R1_Cd$meanLL<-sapply(c(1:nrow(R1_Cd)), function(x) mean(c(R1_Cd$LL_SR1[x], R1_Cd$LL_SR2[x],R1_Cd$LL_SR4[x], R1_Cd$LL_SR5[x]), na.rm=TRUE))
  
  R2_N$meanLL<-sapply(c(1:nrow(R2_N)), function(x) mean(c(R2_N$LL_SR1[x], R2_N$LL_SR4[x], R2_N$LL_SR5[x]), na.rm=TRUE))
  
  R2_Cd$meanLL<-sapply(c(1:nrow(R2_Cd)), function(x) mean(c(R2_Cd$LL_SR1[x],R2_Cd$LL_SR4[x], R2_Cd$LL_SR5[x]), na.rm=TRUE))
  
  R3_N$meanLL<-sapply(c(1:nrow(R3_N)), function(x) mean(c(R3_N$LL_SR1[x], R3_N$LL_SR2[x],R3_N$LL_SR4[x], R3_N$LL_SR5[x]), na.rm=TRUE))
  
  R3_Cd$meanLL<-sapply(c(1:nrow(R3_Cd)), function(x) mean(c(R3_Cd$LL_SR1[x], R3_Cd$LL_SR2[x],R3_Cd$LL_SR4[x], R3_Cd$LL_SR5[x]), na.rm=TRUE))
  
  R4_N$meanLL<-sapply(c(1:nrow(R4_N)), function(x) mean(c(R4_N$LL_SR1[x], R4_N$LL_SR2[x],R4_N$LL_SR4[x], R4_N$LL_SR5[x]), na.rm=TRUE))
  
  R4_Cd$meanLL<-sapply(c(1:nrow(R4_Cd)), function(x) mean(c(R4_Cd$LL_SR1[x], R4_Cd$LL_SR2[x],R4_Cd$LL_SR4[x], R4_Cd$LL_SR5[x]), na.rm=TRUE))
  
  R5_N$meanLL<-sapply(c(1:nrow(R5_N)), function(x) mean(c(R5_N$LL_SR1[x], R5_N$LL_SR2[x],R5_N$LL_SR4[x], R5_N$LL_SR5[x]), na.rm=TRUE))
  
  R5_Cd$meanLL<-sapply(c(1:nrow(R5_Cd)), function(x) mean(c(R5_Cd$LL_SR1[x], R5_Cd$LL_SR2[x],R5_Cd$LL_SR4[x], R5_Cd$LL_SR5[x]), na.rm=TRUE))
  
  
  #' 
  # Estimating root mean squared error (RMSE)------
  #' 
  ##### Replicate 1 ------------------------------------
  
  print("Estimate root mean square error for replicate data")
  
  R1_N$RMSE_SR1<-sapply(c(1:nrow(R1_N)), function(x){
    rssq_function(R1[[1]][,c("fitness","SR1","SR4")], R1_N$lambda_SR1[x], R1_N$intra_SR1[x],R1_N$inter_SR1.SR4[x])+rssq_function(R1[[1]][,c("fitness","SR1","SR5")], R1_N$lambda_SR1[x], R1_N$intra_SR1[x],R1_N$inter_SR1.SR5[x])
  })
  
  R1_N$RMSE_SR2<-sapply(c(1:nrow(R1_N)), function(x){
    rssq_function(R1[[2]][,c("fitness","SR2","SR4")], R1_N$lambda_SR2[x], R1_N$intra_SR2[x],R1_N$inter_SR2.SR4[x])+rssq_function(R1[[2]][,c("fitness","SR2","SR5")], R1_N$lambda_SR2[x], R1_N$intra_SR2[x],R1_N$inter_SR2.SR5[x])
  })
  
  R1_N$RMSE_SR4<-sapply(c(1:nrow(R1_N)), function(x){
    rssq_function(R1[[3]][,c("fitness","SR4","SR1")], R1_N$lambda_SR4[x], R1_N$intra_SR4[x],R1_N$inter_SR4.SR1[x])+rssq_function(R1[[3]][,c("fitness","SR4","SR1")], R1_N$lambda_SR4[x], R1_N$intra_SR4[x],R1_N$inter_SR4.SR2[x])
  })
  
  R1_N$RMSE_SR5<-sapply(c(1:nrow(R1_N)), function(x){
    rssq_function(R1[[4]][,c("fitness","SR5","SR1")], R1_N$lambda_SR5[x], R1_N$intra_SR5[x],R1_N$inter_SR5.SR1[x])+rssq_function(R1[[4]][,c("fitness","SR5","SR2")], R1_N$lambda_SR5[x], R1_N$intra_SR5[x],R1_N$inter_SR5.SR2[x])
  })
  
  ### Replicates Cd
  
  R1_Cd$RMSE_SR1<-sapply(c(1:nrow(R1_Cd)), function(x){
    rssq_function(Rep_R1_Cd[[1]][,c("fitness","SR1","SR4")], R1_Cd$lambda_SR1[x], R1_Cd$intra_SR1[x],R1_Cd$inter_SR1.SR4[x])+rssq_function(Rep_R1_Cd[[1]][,c("fitness","SR1","SR5")], R1_Cd$lambda_SR1[x], R1_Cd$intra_SR1[x],R1_Cd$inter_SR1.SR5[x])
  })
  
  R1_Cd$RMSE_SR2<-sapply(c(1:nrow(R1_Cd)), function(x){
    rssq_function(Rep_R1_Cd[[2]][,c("fitness","SR2","SR4")], R1_Cd$lambda_SR2[x], R1_Cd$intra_SR2[x],R1_Cd$inter_SR2.SR4[x])+rssq_function(Rep_R1_Cd[[2]][,c("fitness","SR2","SR5")], R1_Cd$lambda_SR2[x], R1_Cd$intra_SR2[x],R1_Cd$inter_SR2.SR5[x])
  })
  
  R1_Cd$RMSE_SR4<-sapply(c(1:nrow(R1_Cd)), function(x){
    rssq_function(Rep_R1_Cd[[3]][,c("fitness","SR4","SR1")], R1_Cd$lambda_SR4[x], R1_Cd$intra_SR4[x],R1_Cd$inter_SR4.SR1[x])+rssq_function(Rep_R1_Cd[[3]][,c("fitness","SR4","SR1")], R1_Cd$lambda_SR4[x], R1_Cd$intra_SR4[x],R1_Cd$inter_SR4.SR2[x])
  })
  
  R1_Cd$RMSE_SR5<-sapply(c(1:nrow(R1_Cd)), function(x){
    rssq_function(Rep_R1_Cd[[4]][,c("fitness","SR5","SR1")], R1_Cd$lambda_SR5[x], R1_Cd$intra_SR5[x],R1_Cd$inter_SR5.SR1[x])+rssq_function(Rep_R1_Cd[[4]][,c("fitness","SR5","SR2")], R1_Cd$lambda_SR5[x], R1_Cd$intra_SR5[x],R1_Cd$inter_SR5.SR2[x])
  })
  
  # Estimate the sum or mean of RMSE
  
  R1_N$Sum_RMSE<-sapply(c(1:nrow(R1_N)), function(x) sum(c(R1_N$RMSE_SR1[x],R1_N$RMSE_SR2[x],R1_N$RMSE_SR4[x],R1_N$RMSE_SR5[x])))
  
  R1_N$Mean_RMSE<-sapply(c(1:nrow(R1_N)), function(x) mean(c(R1_N$RMSE_SR1[x],R1_N$RMSE_SR2[x],R1_N$RMSE_SR4[x],R1_N$RMSE_SR5[x])))
  
  R1_Cd$Sum_RMSE<-sapply(c(1:nrow(R1_Cd)), function(x) sum(c(R1_Cd$RMSE_SR1[x],R1_Cd$RMSE_SR2[x],R1_Cd$RMSE_SR4[x],R1_Cd$RMSE_SR5[x])))
  
  R1_Cd$Mean_RMSE<-sapply(c(1:nrow(R1_Cd)), function(x) mean(c(R1_Cd$RMSE_SR1[x],R1_Cd$RMSE_SR2[x],R1_Cd$RMSE_SR4[x],R1_Cd$RMSE_SR5[x])))
  
  # Using the mean or the sum renders the same best
  # However the estimates for alpha and lambda show a large variation, specially in the cadmium environment. So we will fit the data for the 10 best candidates from RMSE and LL and for each selection regime
  head(R1_N[order(R1_N$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R1_N[order(R1_N$Mean_RMSE, decreasing=FALSE),], n=10)
  
  head(R1_Cd[order(R1_Cd$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R1_Cd[order(R1_Cd$Mean_RMSE, decreasing=FALSE),], n=10)
  
  
  #' 
  ##### Replicate 2------------------------------------
  
  R2_N$RMSE_SR1<-sapply(c(1:nrow(R2_N)), function(x){
    rssq_function(R2[[1]][,c("fitness","SR1","SR4")], R2_N$lambda_SR1[x], R2_N$intra_SR1[x],R2_N$inter_SR1.SR4[x])+rssq_function(R2[[1]][,c("fitness","SR1","SR5")], R2_N$lambda_SR1[x], R2_N$intra_SR1[x],R2_N$inter_SR1.SR5[x])
  })
  
  
  R2_N$RMSE_SR4<-sapply(c(1:nrow(R2_N)), function(x){
    rssq_function(R2[[2]][,c("fitness","SR4","SR1")], R2_N$lambda_SR4[x], R2_N$intra_SR4[x],R2_N$inter_SR4.SR1[x])
  })
  
  R2_N$RMSE_SR5<-sapply(c(1:nrow(R2_N)), function(x){
    rssq_function(R2[[3]][,c("fitness","SR5","SR1")], R2_N$lambda_SR5[x], R2_N$intra_SR5[x],R2_N$inter_SR5.SR1[x])
  })
  
  ### Replicates Cd
  
  R2_Cd$RMSE_SR1<-sapply(c(1:nrow(R2_Cd)), function(x){
    rssq_function(Rep_R2_Cd[[1]][,c("fitness","SR1","SR4")], R2_Cd$lambda_SR1[x], R2_Cd$intra_SR1[x],R2_Cd$inter_SR1.SR4[x])+rssq_function(Rep_R2_Cd[[1]][,c("fitness","SR1","SR5")], R2_Cd$lambda_SR1[x], R2_Cd$intra_SR1[x],R2_Cd$inter_SR1.SR5[x])
  })
  
  
  R2_Cd$RMSE_SR4<-sapply(c(1:nrow(R2_Cd)), function(x){
    rssq_function(Rep_R2_Cd[[2]][,c("fitness","SR4","SR1")], R2_Cd$lambda_SR4[x], R2_Cd$intra_SR4[x],R2_Cd$inter_SR4.SR1[x])
  })
  
  R2_Cd$RMSE_SR5<-sapply(c(1:nrow(R2_Cd)), function(x){
    rssq_function(Rep_R2_Cd[[3]][,c("fitness","SR5","SR1")], R2_Cd$lambda_SR5[x], R2_Cd$intra_SR5[x],R2_Cd$inter_SR5.SR1[x])
  })
  
  # Estimate the sum or mean of RMSE
  
  R2_N$Sum_RMSE<-sapply(c(1:nrow(R2_N)), function(x) sum(c(R2_N$RMSE_SR1[x],R2_N$RMSE_SR4[x],R2_N$RMSE_SR5[x])))
  
  R2_N$Mean_RMSE<-sapply(c(1:nrow(R2_N)), function(x) mean(c(R2_N$RMSE_SR1[x],R2_N$RMSE_SR4[x],R2_N$RMSE_SR5[x])))
  
  R2_Cd$Sum_RMSE<-sapply(c(1:nrow(R2_Cd)), function(x) sum(c(R2_Cd$RMSE_SR1[x],R2_Cd$RMSE_SR4[x],R2_Cd$RMSE_SR5[x])))
  
  R2_Cd$Mean_RMSE<-sapply(c(1:nrow(R2_Cd)), function(x) mean(c(R2_Cd$RMSE_SR1[x],R2_Cd$RMSE_SR4[x],R2_Cd$RMSE_SR5[x])))
  
  # Using the mean or the sum renders the same best
  head(R2_N[order(R2_N$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R2_N[order(R2_N$Mean_RMSE, decreasing=FALSE),], n=10)
  
  head(R2_Cd[order(R2_Cd$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R2_Cd[order(R2_Cd$Mean_RMSE, decreasing=FALSE),], n=10)
  
  
  #' 
  ##### Replicate 3------------------------------------
  
  R3_N$RMSE_SR1<-sapply(c(1:nrow(R3_N)), function(x){
    rssq_function(R3[[1]][,c("fitness","SR1","SR4")], R3_N$lambda_SR1[x], R3_N$intra_SR1[x],R3_N$inter_SR1.SR4[x])+rssq_function(R3[[1]][,c("fitness","SR1","SR5")], R3_N$lambda_SR1[x], R3_N$intra_SR1[x],R3_N$inter_SR1.SR5[x])
  })
  
  R3_N$RMSE_SR2<-sapply(c(1:nrow(R3_N)), function(x){
    rssq_function(R3[[2]][,c("fitness","SR2","SR4")], R3_N$lambda_SR2[x], R3_N$intra_SR2[x],R3_N$inter_SR2.SR4[x])+rssq_function(R3[[2]][,c("fitness","SR2","SR5")], R3_N$lambda_SR2[x], R3_N$intra_SR2[x],R3_N$inter_SR2.SR5[x])
  })
  
  R3_N$RMSE_SR4<-sapply(c(1:nrow(R3_N)), function(x){
    rssq_function(R3[[3]][,c("fitness","SR4","SR1")], R3_N$lambda_SR4[x], R3_N$intra_SR4[x],R3_N$inter_SR4.SR1[x])+rssq_function(R3[[3]][,c("fitness","SR4","SR1")], R3_N$lambda_SR4[x], R3_N$intra_SR4[x],R3_N$inter_SR4.SR2[x])
  })
  
  R3_N$RMSE_SR5<-sapply(c(1:nrow(R3_N)), function(x){
    rssq_function(R3[[4]][,c("fitness","SR5","SR1")], R3_N$lambda_SR5[x], R3_N$intra_SR5[x],R3_N$inter_SR5.SR1[x])+rssq_function(R3[[4]][,c("fitness","SR5","SR2")], R3_N$lambda_SR5[x], R3_N$intra_SR5[x],R3_N$inter_SR5.SR2[x])
  })
  
  ### Replicates Cd
  
  R3_Cd$RMSE_SR1<-sapply(c(1:nrow(R3_Cd)), function(x){
    rssq_function(Rep_R3_Cd[[1]][,c("fitness","SR1","SR4")], R3_Cd$lambda_SR1[x], R3_Cd$intra_SR1[x],R3_Cd$inter_SR1.SR4[x])+rssq_function(Rep_R3_Cd[[1]][,c("fitness","SR1","SR5")], R3_Cd$lambda_SR1[x], R3_Cd$intra_SR1[x],R3_Cd$inter_SR1.SR5[x])
  })
  
  R3_Cd$RMSE_SR2<-sapply(c(1:nrow(R3_Cd)), function(x){
    rssq_function(Rep_R3_Cd[[2]][,c("fitness","SR2","SR4")], R3_Cd$lambda_SR2[x], R3_Cd$intra_SR2[x],R3_Cd$inter_SR2.SR4[x])+rssq_function(Rep_R3_Cd[[2]][,c("fitness","SR2","SR5")], R3_Cd$lambda_SR2[x], R3_Cd$intra_SR2[x],R3_Cd$inter_SR2.SR5[x])
  })
  
  R3_Cd$RMSE_SR4<-sapply(c(1:nrow(R3_Cd)), function(x){
    rssq_function(Rep_R3_Cd[[3]][,c("fitness","SR4","SR1")], R3_Cd$lambda_SR4[x], R3_Cd$intra_SR4[x],R3_Cd$inter_SR4.SR1[x])+rssq_function(Rep_R3_Cd[[3]][,c("fitness","SR4","SR1")], R3_Cd$lambda_SR4[x], R3_Cd$intra_SR4[x],R3_Cd$inter_SR4.SR2[x])
  })
  
  R3_Cd$RMSE_SR5<-sapply(c(1:nrow(R3_Cd)), function(x){
    rssq_function(Rep_R3_Cd[[4]][,c("fitness","SR5","SR1")], R3_Cd$lambda_SR5[x], R3_Cd$intra_SR5[x],R3_Cd$inter_SR5.SR1[x])+rssq_function(Rep_R3_Cd[[4]][,c("fitness","SR5","SR2")], R3_Cd$lambda_SR5[x], R3_Cd$intra_SR5[x],R3_Cd$inter_SR5.SR2[x])
  })
  
  # Estimate the sum or mean of RMSE
  
  R3_N$Sum_RMSE<-sapply(c(1:nrow(R3_N)), function(x) sum(c(R3_N$RMSE_SR1[x],R3_N$RMSE_SR2[x],R3_N$RMSE_SR4[x],R3_N$RMSE_SR5[x])))
  
  R3_N$Mean_RMSE<-sapply(c(1:nrow(R3_N)), function(x) mean(c(R3_N$RMSE_SR1[x],R3_N$RMSE_SR2[x],R3_N$RMSE_SR4[x],R3_N$RMSE_SR5[x])))
  
  R3_Cd$Sum_RMSE<-sapply(c(1:nrow(R3_Cd)), function(x) sum(c(R3_Cd$RMSE_SR1[x],R3_Cd$RMSE_SR2[x],R3_Cd$RMSE_SR4[x],R3_Cd$RMSE_SR5[x])))
  
  R3_Cd$Mean_RMSE<-sapply(c(1:nrow(R3_Cd)), function(x) mean(c(R3_Cd$RMSE_SR1[x],R3_Cd$RMSE_SR2[x],R3_Cd$RMSE_SR4[x],R3_Cd$RMSE_SR5[x])))
  
  # Using the mean or the sum renders the same best
  head(R3_N[order(R3_N$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R3_N[order(R3_N$Mean_RMSE, decreasing=FALSE),], n=10)
  
  head(R3_Cd[order(R3_Cd$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R3_Cd[order(R3_Cd$Mean_RMSE, decreasing=FALSE),], n=10)
  
  
  #' 
  ##### Replicate 4------------------------------------
  
  R4_N$RMSE_SR1<-sapply(c(1:nrow(R4_N)), function(x){
    rssq_function(R4[[1]][,c("fitness","SR1","SR4")], R4_N$lambda_SR1[x], R4_N$intra_SR1[x],R4_N$inter_SR1.SR4[x])+rssq_function(R4[[1]][,c("fitness","SR1","SR5")], R4_N$lambda_SR1[x], R4_N$intra_SR1[x],R4_N$inter_SR1.SR5[x])
  })
  
  R4_N$RMSE_SR2<-sapply(c(1:nrow(R4_N)), function(x){
    rssq_function(R4[[2]][,c("fitness","SR2","SR4")], R4_N$lambda_SR2[x], R4_N$intra_SR2[x],R4_N$inter_SR2.SR4[x])+rssq_function(R4[[2]][,c("fitness","SR2","SR5")], R4_N$lambda_SR2[x], R4_N$intra_SR2[x],R4_N$inter_SR2.SR5[x])
  })
  
  R4_N$RMSE_SR4<-sapply(c(1:nrow(R4_N)), function(x){
    rssq_function(R4[[3]][,c("fitness","SR4","SR1")], R4_N$lambda_SR4[x], R4_N$intra_SR4[x],R4_N$inter_SR4.SR1[x])+rssq_function(R4[[3]][,c("fitness","SR4","SR1")], R4_N$lambda_SR4[x], R4_N$intra_SR4[x],R4_N$inter_SR4.SR2[x])
  })
  
  R4_N$RMSE_SR5<-sapply(c(1:nrow(R4_N)), function(x){
    rssq_function(R4[[4]][,c("fitness","SR5","SR1")], R4_N$lambda_SR5[x], R4_N$intra_SR5[x],R4_N$inter_SR5.SR1[x])+rssq_function(R4[[4]][,c("fitness","SR5","SR2")], R4_N$lambda_SR5[x], R4_N$intra_SR5[x],R4_N$inter_SR5.SR2[x])
  })
  
  ### Replicates Cd
  
  R4_Cd$RMSE_SR1<-sapply(c(1:nrow(R4_Cd)), function(x){
    rssq_function(Rep_R4_Cd[[1]][,c("fitness","SR1","SR4")], R4_Cd$lambda_SR1[x], R4_Cd$intra_SR1[x],R4_Cd$inter_SR1.SR4[x])+rssq_function(Rep_R4_Cd[[1]][,c("fitness","SR1","SR5")], R4_Cd$lambda_SR1[x], R4_Cd$intra_SR1[x],R4_Cd$inter_SR1.SR5[x])
  })
  
  R4_Cd$RMSE_SR2<-sapply(c(1:nrow(R4_Cd)), function(x){
    rssq_function(Rep_R4_Cd[[2]][,c("fitness","SR2","SR4")], R4_Cd$lambda_SR2[x], R4_Cd$intra_SR2[x],R4_Cd$inter_SR2.SR4[x])+rssq_function(Rep_R4_Cd[[2]][,c("fitness","SR2","SR5")], R4_Cd$lambda_SR2[x], R4_Cd$intra_SR2[x],R4_Cd$inter_SR2.SR5[x])
  })
  
  R4_Cd$RMSE_SR4<-sapply(c(1:nrow(R4_Cd)), function(x){
    rssq_function(Rep_R4_Cd[[3]][,c("fitness","SR4","SR1")], R4_Cd$lambda_SR4[x], R4_Cd$intra_SR4[x],R4_Cd$inter_SR4.SR1[x])+rssq_function(Rep_R4_Cd[[3]][,c("fitness","SR4","SR1")], R4_Cd$lambda_SR4[x], R4_Cd$intra_SR4[x],R4_Cd$inter_SR4.SR2[x])
  })
  
  R4_Cd$RMSE_SR5<-sapply(c(1:nrow(R4_Cd)), function(x){
    rssq_function(Rep_R4_Cd[[4]][,c("fitness","SR5","SR1")], R4_Cd$lambda_SR5[x], R4_Cd$intra_SR5[x],R4_Cd$inter_SR5.SR1[x])+rssq_function(Rep_R4_Cd[[4]][,c("fitness","SR5","SR2")], R4_Cd$lambda_SR5[x], R4_Cd$intra_SR5[x],R4_Cd$inter_SR5.SR2[x])
  })
  
  # Estimate the sum or mean of RMSE
  
  R4_N$Sum_RMSE<-sapply(c(1:nrow(R4_N)), function(x) sum(c(R4_N$RMSE_SR1[x],R4_N$RMSE_SR2[x],R4_N$RMSE_SR4[x],R4_N$RMSE_SR5[x])))
  
  R4_N$Mean_RMSE<-sapply(c(1:nrow(R4_N)), function(x) mean(c(R4_N$RMSE_SR1[x],R4_N$RMSE_SR2[x],R4_N$RMSE_SR4[x],R4_N$RMSE_SR5[x])))
  
  R4_Cd$Sum_RMSE<-sapply(c(1:nrow(R4_Cd)), function(x) sum(c(R4_Cd$RMSE_SR1[x],R4_Cd$RMSE_SR2[x],R4_Cd$RMSE_SR4[x],R4_Cd$RMSE_SR5[x])))
  
  R4_Cd$Mean_RMSE<-sapply(c(1:nrow(R4_Cd)), function(x) mean(c(R4_Cd$RMSE_SR1[x],R4_Cd$RMSE_SR2[x],R4_Cd$RMSE_SR4[x],R4_Cd$RMSE_SR5[x])))
  
  # Using the mean or the sum renders the same best
  # However the estimates for alpha and lambda show a large variation, specially in the cadmium environment. So we will fit the data for the 10 best candidates from RMSE and LL and for each selection regime
  head(R4_N[order(R4_N$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R4_N[order(R4_N$Mean_RMSE, decreasing=FALSE),], n=10)
  
  head(R4_Cd[order(R4_Cd$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R4_Cd[order(R4_Cd$Mean_RMSE, decreasing=FALSE),], n=10)
  
  
  #' 
  ##### Replicate 5------------------------------------
  
  R5_N$RMSE_SR1<-sapply(c(1:nrow(R5_N)), function(x){
    rssq_function(R5[[1]][,c("fitness","SR1","SR4")], R5_N$lambda_SR1[x], R5_N$intra_SR1[x],R5_N$inter_SR1.SR4[x])+rssq_function(R5[[1]][,c("fitness","SR1","SR5")], R5_N$lambda_SR1[x], R5_N$intra_SR1[x],R5_N$inter_SR1.SR5[x])
  })
  
  R5_N$RMSE_SR2<-sapply(c(1:nrow(R5_N)), function(x){
    rssq_function(R5[[2]][,c("fitness","SR2","SR4")], R5_N$lambda_SR2[x], R5_N$intra_SR2[x],R5_N$inter_SR2.SR4[x])+rssq_function(R5[[2]][,c("fitness","SR2","SR5")], R5_N$lambda_SR2[x], R5_N$intra_SR2[x],R5_N$inter_SR2.SR5[x])
  })
  
  R5_N$RMSE_SR4<-sapply(c(1:nrow(R5_N)), function(x){
    rssq_function(R5[[3]][,c("fitness","SR4","SR1")], R5_N$lambda_SR4[x], R5_N$intra_SR4[x],R5_N$inter_SR4.SR1[x])+rssq_function(R5[[3]][,c("fitness","SR4","SR1")], R5_N$lambda_SR4[x], R5_N$intra_SR4[x],R5_N$inter_SR4.SR2[x])
  })
  
  R5_N$RMSE_SR5<-sapply(c(1:nrow(R5_N)), function(x){
    rssq_function(R5[[4]][,c("fitness","SR5","SR1")], R5_N$lambda_SR5[x], R5_N$intra_SR5[x],R5_N$inter_SR5.SR1[x])+rssq_function(R5[[4]][,c("fitness","SR5","SR2")], R5_N$lambda_SR5[x], R5_N$intra_SR5[x],R5_N$inter_SR5.SR2[x])
  })
  
  ### Replicates Cd
  
  R5_Cd$RMSE_SR1<-sapply(c(1:nrow(R5_Cd)), function(x){
    rssq_function(Rep_R5_Cd[[1]][,c("fitness","SR1","SR4")], R5_Cd$lambda_SR1[x], R5_Cd$intra_SR1[x],R5_Cd$inter_SR1.SR4[x])+rssq_function(Rep_R5_Cd[[1]][,c("fitness","SR1","SR5")], R5_Cd$lambda_SR1[x], R5_Cd$intra_SR1[x],R5_Cd$inter_SR1.SR5[x])
  })
  
  R5_Cd$RMSE_SR2<-sapply(c(1:nrow(R5_Cd)), function(x){
    rssq_function(Rep_R5_Cd[[2]][,c("fitness","SR2","SR4")], R5_Cd$lambda_SR2[x], R5_Cd$intra_SR2[x],R5_Cd$inter_SR2.SR4[x])+rssq_function(Rep_R5_Cd[[2]][,c("fitness","SR2","SR5")], R5_Cd$lambda_SR2[x], R5_Cd$intra_SR2[x],R5_Cd$inter_SR2.SR5[x])
  })
  
  R5_Cd$RMSE_SR4<-sapply(c(1:nrow(R5_Cd)), function(x){
    rssq_function(Rep_R5_Cd[[3]][,c("fitness","SR4","SR1")], R5_Cd$lambda_SR4[x], R5_Cd$intra_SR4[x],R5_Cd$inter_SR4.SR1[x])+rssq_function(Rep_R5_Cd[[3]][,c("fitness","SR4","SR1")], R5_Cd$lambda_SR4[x], R5_Cd$intra_SR4[x],R5_Cd$inter_SR4.SR2[x])
  })
  
  R5_Cd$RMSE_SR5<-sapply(c(1:nrow(R5_Cd)), function(x){
    rssq_function(Rep_R5_Cd[[4]][,c("fitness","SR5","SR1")], R5_Cd$lambda_SR5[x], R5_Cd$intra_SR5[x],R5_Cd$inter_SR5.SR1[x])+rssq_function(Rep_R5_Cd[[4]][,c("fitness","SR5","SR2")], R5_Cd$lambda_SR5[x], R5_Cd$intra_SR5[x],R5_Cd$inter_SR5.SR2[x])
  })
  
  # Estimate the sum or mean of RMSE
  
  R5_N$Sum_RMSE<-sapply(c(1:nrow(R5_N)), function(x) sum(c(R5_N$RMSE_SR1[x],R5_N$RMSE_SR2[x],R5_N$RMSE_SR4[x],R5_N$RMSE_SR5[x])))
  
  R5_N$Mean_RMSE<-sapply(c(1:nrow(R5_N)), function(x) mean(c(R5_N$RMSE_SR1[x],R5_N$RMSE_SR2[x],R5_N$RMSE_SR4[x],R5_N$RMSE_SR5[x])))
  
  R5_Cd$Sum_RMSE<-sapply(c(1:nrow(R5_Cd)), function(x) sum(c(R5_Cd$RMSE_SR1[x],R5_Cd$RMSE_SR2[x],R5_Cd$RMSE_SR4[x],R5_Cd$RMSE_SR5[x])))
  
  R5_Cd$Mean_RMSE<-sapply(c(1:nrow(R5_Cd)), function(x) mean(c(R5_Cd$RMSE_SR1[x],R5_Cd$RMSE_SR2[x],R5_Cd$RMSE_SR4[x],R5_Cd$RMSE_SR5[x])))
  
  # Using the mean or the sum renders the same best
  
  head(R5_N[order(R5_N$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R5_N[order(R5_N$Mean_RMSE, decreasing=FALSE),], n=10)
  
  head(R5_Cd[order(R5_Cd$Sum_RMSE, decreasing=FALSE),], n=10)
  head(R5_Cd[order(R5_Cd$Mean_RMSE, decreasing=FALSE),], n=10)
  
  str(R5_Cd)
  
  #' 
  ### Save or load------
    save.image(file='./Analyses/Parameter_exploration_RMSE.RData')
  }

  #' 
  #' 

  # Final list of parameters to explore------
  #' 
  #### No cadmium environment------
  #' Checking the best 10 values across the different 
  #' 
  
  
  head(pooled_N[order(pooled_N$meanLL, decreasing=FALSE),], n=10)
  
  head(pooled_Cd[order(pooled_Cd$meanLL, decreasing=FALSE),], n=10)
  
  #' 
  ## Estimate loglikelihood across the 5 replicates------
  #' From here we estimate the initial points for the replicates and pooled data
  
  print("Estimate mean likelihood across the five replicates")
  
  mean_replicates<-as.data.frame(R1_N[,c("lambda", "alpha_intra", "alpha_inter")])
  
  mean_replicates$mean_LL_N<-sapply(c(1:nrow(mean_replicates)), function(x){
    aux1<-R1_N$meanLL[which(R1_N$lambda==mean_replicates$lambda[x] & R1_N$alpha_intra==mean_replicates$alpha_intra[x] & R1_N$alpha_inter==mean_replicates$alpha_inter[x])]
    
    aux2<-R2_N$meanLL[which(R2_N$lambda==mean_replicates$lambda[x] & R2_N$alpha_intra==mean_replicates$alpha_intra[x] & R2_N$alpha_inter==mean_replicates$alpha_inter[x])]
    
    aux3<-R3_N$meanLL[which(R3_N$lambda==mean_replicates$lambda[x] & R3_N$alpha_intra==mean_replicates$alpha_intra[x] & R3_N$alpha_inter==mean_replicates$alpha_inter[x])]
    
    aux4<-R4_N$meanLL[which(R4_N$lambda==mean_replicates$lambda[x] & R4_N$alpha_intra==mean_replicates$alpha_intra[x] & R4_N$alpha_inter==mean_replicates$alpha_inter[x])]
    
    aux5<-R5_N$meanLL[which(R5_N$lambda==mean_replicates$lambda[x] & R5_N$alpha_intra==mean_replicates$alpha_intra[x] & R5_N$alpha_inter==mean_replicates$alpha_inter[x])]
    
    mean(c(aux1, aux2,aux3, aux4, aux5), na.rm=TRUE)
    
  })
  
  mean_replicates$mean_LL_Cd<-sapply(c(1:nrow(mean_replicates)), function(x){
    
    aux1<-R1_Cd$meanLL[which(R1_Cd$lambda==mean_replicates$lambda[x] & R1_Cd$alpha_intra==mean_replicates$alpha_intra[x] & R1_Cd$alpha_inter==mean_replicates$alpha_inter[x])]
    
    aux2<-R2_Cd$meanLL[which(R2_Cd$lambda==mean_replicates$lambda[x] & R2_Cd$alpha_intra==mean_replicates$alpha_intra[x] & R2_Cd$alpha_inter==mean_replicates$alpha_inter[x])]
    
    aux3<-R3_Cd$meanLL[which(R3_Cd$lambda==mean_replicates$lambda[x] & R3_Cd$alpha_intra==mean_replicates$alpha_intra[x] & R3_Cd$alpha_inter==mean_replicates$alpha_inter[x])]
    
    aux4<-R4_Cd$meanLL[which(R4_Cd$lambda==mean_replicates$lambda[x] & R4_Cd$alpha_intra==mean_replicates$alpha_intra[x] & R4_Cd$alpha_inter==mean_replicates$alpha_inter[x])]
    
    aux5<-R5_Cd$meanLL[which(R5_Cd$lambda==mean_replicates$lambda[x] & R5_Cd$alpha_intra==mean_replicates$alpha_intra[x] & R5_Cd$alpha_inter==mean_replicates$alpha_inter[x])]
    
    mean(c(aux1, aux2,aux3, aux4, aux5), na.rm=TRUE)
    
  })
  
  mean_replicates[which(mean_replicates$mean_LL_N==min(mean_replicates$mean_LL_N)),]
  
  mean_replicates[which(mean_replicates$mean_LL_Cd==min(mean_replicates$mean_LL_Cd)),]
  
  best_pooledN<-pooled_N[which(pooled_N$meanLL==min(pooled_N$meanLL)),c("lambda", "alpha_intra", "alpha_inter")]
  best_pooledCd<-pooled_Cd[which(pooled_Cd$meanLL==min(pooled_Cd$meanLL)),c("lambda", "alpha_intra", "alpha_inter")]
  
  pooled_N[which(pooled_N$lambda==mean_replicates$lambda[which(mean_replicates$mean_LL_N==min(mean_replicates$mean_LL_N))] & pooled_N$alpha_intra==mean_replicates$alpha_intra[which(mean_replicates$mean_LL_N==min(mean_replicates$mean_LL_N))] & pooled_N$alpha_inter==mean_replicates$alpha_inter[which(mean_replicates$mean_LL_N==min(mean_replicates$mean_LL_N))]),]
  
  pooled_Cd[which(pooled_Cd$lambda==mean_replicates$lambda[which(mean_replicates$mean_LL_Cd==min(mean_replicates$mean_LL_Cd))] & pooled_Cd$alpha_intra==mean_replicates$alpha_intra[which(mean_replicates$mean_LL_Cd==min(mean_replicates$mean_LL_Cd))] & pooled_Cd$alpha_inter==mean_replicates$alpha_inter[which(mean_replicates$mean_LL_Cd==min(mean_replicates$mean_LL_Cd))]),]
  
  # The best parameters are not the same for both environments. But we can see which combination is best across both
  
  #For pooled estimates the two values are very similar. We will just use the pooled N initial values
  subset(pooled_Cd, (lambda==best_pooledCd$lambda[1] & alpha_inter==best_pooledCd$alpha_inter[1] & alpha_intra==best_pooledCd$alpha_intra[1]) | (lambda==best_pooledN$lambda[1] & alpha_inter==best_pooledN$alpha_inter[1] & alpha_intra==best_pooledN$alpha_intra[1]))
  
  subset(pooled_N, (lambda==best_pooledCd$lambda[1] & alpha_inter==best_pooledCd$alpha_inter[1] & alpha_intra==best_pooledCd$alpha_intra[1]) | (lambda==best_pooledN$lambda[1] & alpha_inter==best_pooledN$alpha_inter[1] & alpha_intra==best_pooledN$alpha_intra[1]))
  
  # For replicates
  mean_replicates[which(mean_replicates$mean_LL_N==min(mean_replicates$mean_LL_N) | mean_replicates$mean_LL_Cd==min(mean_replicates$mean_LL_Cd)) ,] 
  
  
  best_replicates<-mean_replicates[which(mean_replicates$mean_LL_N==min(mean_replicates$mean_LL_N) | mean_replicates$mean_LL_Cd==min(mean_replicates$mean_LL_Cd)) ,c("lambda", "alpha_intra", "alpha_inter")] 
  # Again very similar values, similar to the pooled estimates we will use the N estimates
  
  
   best_pooledCd
   best_pooledN
   best_replicates
  
  # See the parameter values that correspond to the best estimate to see if there are large differences between the two
  
  # For the non-cadmium environment
  subset(R1_N, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  subset(R2_N, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  subset(R3_N, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  subset(R4_N, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  subset(R5_N, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  # For the cadmium environment
  subset(R1_Cd, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  subset(R2_Cd, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  subset(R3_Cd, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  subset(R4_Cd, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  subset(R5_Cd, lambda==best_replicates$lambda & alpha_intra==best_replicates$alpha_intra & alpha_inter==best_replicates$alpha_inter)

  
  
  #' 
  #' 
  #' 
  #' 

 
