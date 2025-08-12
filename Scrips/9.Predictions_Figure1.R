rm(list=ls())

library(tidyverse)
library(plyr)
library(dplyr)
library(car)
library(fitdistrplus)
library(tidyr)
library(ggtext)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(ggbreak)
library(effects)
library(cowplot)
library(ggeffects)
library(marginaleffects)
library(ggtext)
library(gridExtra)
library(RColorBrewer)

theme_plot<-theme(axis.text = element_text(size=14), axis.title = element_text(size=14, face="bold"), legend.text = element_text(size=12), strip.text = element_text(size=14), plot.title = element_text(size=14, face="bold"), panel.grid=element_line(colour="white"), panel.background = element_rect(fill="white") , axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"), strip.background = element_rect(fill="white"))

save_plot<-function(dir, width=15, height=10, ...){
  ggsave(dir, width = width, height = height, units = c("cm"))
}


Env<-c("No cadmium", "Cadmium")
names(Env)<-c("N", "Cd")

regimeTu<-c("Tu no cadmium", "Tu cadmium")
names(regimeTu)<-c("SR1", "SR2")

regimeTe<-c("Te no cadmium", "Te cadmium")
names(regimeTe)<-c("SR4", "SR5")

colors_comb<-brewer.pal(name = "Spectral", 4)

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



# Simulations for figure 1

pred_coex1Gen<-as.data.frame(expand_grid(Te=c("SR4","SR5"), Tu=c("SR1", "SR2"), Environment= c("N", "Cd")))


pred_coex1Gen$predTu_onlyLambda<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Tu_lambda[1]*10
  
  bl
})

pred_coex1Gen$predTu_Lambda_INTRA<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10)
  
  bl
})

pred_coex1Gen$predTu_ALL<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10- aux_alphas$Tu_inter[1]*10)
  
  bl
})


pred_coex1Gen$predTe_onlyLambda<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10
  
  bl
})

pred_coex1Gen$predTe_Lambda_INTRA<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10)
  
  bl
})

pred_coex1Gen$predTe_ALL<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10- aux_alphas$Te_inter[1]*10)
  
  bl
})

pred_coex1Gen$Control_lambdaTu<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Tu=="SR1")
  
  pred_coex1Gen$predTu_onlyLambda[x]/cont$predTu_onlyLambda[1]
  
})

pred_coex1Gen$Control_lambdaTe<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4")
  
  pred_coex1Gen$predTe_onlyLambda[x]/cont$predTe_onlyLambda[1]
  
})

pred_coex1Gen$Control_lambdaIntraTu<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTu_Lambda_INTRA[x]/cont$predTu_Lambda_INTRA[1]
  
})

pred_coex1Gen$Control_lambdaIntraTe<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTe_Lambda_INTRA[x]/cont$predTe_Lambda_INTRA[1]
  
})

pred_coex1Gen$Control_ALLTu<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTu_ALL[x]/cont$predTu_ALL[1]
  
})

pred_coex1Gen$Control_ALLTe<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTe_ALL[x]/cont$predTe_ALL[1]
  
})

### Lower

pred_coex1Gen$predTu_onlyLambda_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Tu_lambda[1]*10
  
  bl
})

pred_coex1Gen$predTu_Lambda_INTRA_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_lambdas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_lambdas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10)
  
  bl
})

pred_coex1Gen$predTu_ALL_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_lambdas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_lambdas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10- aux_alphas$Tu_inter[1]*10)
  
  bl
})


pred_coex1Gen$predTe_onlyLambda_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10
  
  bl
})

pred_coex1Gen$predTe_Lambda_INTRA_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_lambdas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_lambdas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10)
  
  bl
})

pred_coex1Gen$predTe_ALL_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_lambdas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_lambdas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10- aux_alphas$Te_inter[1]*10)
  
  bl
})

pred_coex1Gen$Control_lambdaTu_L<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Tu=="SR1")
  
  pred_coex1Gen$predTu_onlyLambda_L[x]/cont$predTu_onlyLambda_L[1]
  
})

pred_coex1Gen$Control_lambdaTe_L<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4")
  
  pred_coex1Gen$predTe_onlyLambda_L[x]/cont$predTe_onlyLambda_L[1]
  
})

pred_coex1Gen$Control_lambdaIntraTu_L<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTu_Lambda_INTRA_L[x]/cont$predTu_Lambda_INTRA_L[1]
  
})

pred_coex1Gen$Control_lambdaIntraTe_L<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTe_Lambda_INTRA_L[x]/cont$predTe_Lambda_INTRA_L[1]
  
})

pred_coex1Gen$Control_ALLTu_L<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTu_ALL_L[x]/cont$predTu_ALL_L[1]
  
})

pred_coex1Gen$Control_ALLTe_L<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTe_ALL_L[x]/cont$predTe_ALL_L[1]
  
})


### Upper

pred_coex1Gen$predTu_onlyLambda_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Tu_lambda[1]*10
  
  bl
})

pred_coex1Gen$predTu_lambda_INTRA_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_Lambdas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_Lambdas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10)
  
  bl
})

pred_coex1Gen$predTu_ALL_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_Lambdas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_Lambdas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10- aux_alphas$Tu_inter[1]*10)
  
  bl
})


pred_coex1Gen$predTe_onlyLambda_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10
  
  bl
})

pred_coex1Gen$predTe_lambda_INTRA_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_Lambdas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_Lambdas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10)
  
  bl
})

pred_coex1Gen$predTe_ALL_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_Lambdas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_Lambdas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10- aux_alphas$Te_inter[1]*10)
  
  bl
})

pred_coex1Gen$Control_LambdaTu_U<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Tu=="SR1")
  
  pred_coex1Gen$predTu_onlyLambda_U[x]/cont$predTu_onlyLambda_U[1]
  
})

pred_coex1Gen$Control_LambdaTe_U<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4")
  
  pred_coex1Gen$predTe_onlyLambda_U[x]/cont$predTe_onlyLambda_U[1]
  
})

pred_coex1Gen$Control_LambdaIntraTu_U<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTu_lambda_INTRA_U[x]/cont$predTu_lambda_INTRA_U[1]
  
})

pred_coex1Gen$Control_LambdaIntraTe_U<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTe_lambda_INTRA_U[x]/cont$predTe_lambda_INTRA_U[1]
  
})

pred_coex1Gen$Control_ALLTu_U<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTu_ALL_U[x]/cont$predTu_ALL_U[1]
  
})

pred_coex1Gen$Control_ALLTe_U<-sapply(c(1:dim(pred_coex1Gen)[1]), function(x){
  cont<-subset(pred_coex1Gen, Environment==pred_coex1Gen$Environment[x] & Te=="SR4" & Tu=="SR1")
  
  pred_coex1Gen$predTe_ALL_U[x]/cont$predTe_ALL_U[1]
  
})

pred_coex1Gen[,c("predTu_ALL", "predTu_ALL_U", "predTu_onlyLambda", "predTu_onlyLambda_U", "predTu_Lambda_INTRA", "predTu_lambda_INTRA_U")]


#### reshaping

pred_coex1Gen_long<-gather(pred_coex1Gen[,c("Te","Tu","Environment","Control_lambdaTe" , "Control_lambdaTu","Control_lambdaIntraTe", "Control_lambdaIntraTu","Control_ALLTe" ,"Control_ALLTu",  "predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL" )], parameter, value, c("Control_lambdaTe" , "Control_lambdaTu","Control_lambdaIntraTe", "Control_lambdaIntraTu","Control_ALLTe" ,"Control_ALLTu",  "predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL" ))

pred_coex1Gen_long_L<-gather(pred_coex1Gen[,c("Te","Tu","Environment","Control_lambdaTe_L" , "Control_lambdaTu_L","Control_lambdaIntraTe_L", "Control_lambdaIntraTu_L","Control_ALLTe_L" ,"Control_ALLTu_L",  "predTe_onlyLambda_L","predTe_Lambda_INTRA_L" ,"predTe_ALL_L" ,"predTu_onlyLambda_L","predTu_Lambda_INTRA_L","predTu_ALL_L" )], parameter, value_L, c("Control_lambdaTe_L" , "Control_lambdaTu_L","Control_lambdaIntraTe_L", "Control_lambdaIntraTu_L","Control_ALLTe_L" ,"Control_ALLTu_L",  "predTe_onlyLambda_L","predTe_Lambda_INTRA_L" ,"predTe_ALL_L" ,"predTu_onlyLambda_L","predTu_Lambda_INTRA_L","predTu_ALL_L"))

pred_coex1Gen_long_L$parameter2<-mapvalues(pred_coex1Gen_long_L$parameter,c("Control_lambdaTe_L" , "Control_lambdaTu_L","Control_lambdaIntraTe_L", "Control_lambdaIntraTu_L","Control_ALLTe_L" ,"Control_ALLTu_L",  "predTe_onlyLambda_L","predTe_Lambda_INTRA_L" ,"predTe_ALL_L" ,"predTu_onlyLambda_L","predTu_Lambda_INTRA_L","predTu_ALL_L"), c("Control_lambdaTe" , "Control_lambdaTu","Control_lambdaIntraTe", "Control_lambdaIntraTu","Control_ALLTe" ,"Control_ALLTu",  "predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL") )


pred_coex1Gen_long_U<-gather(pred_coex1Gen[,c("Te","Tu","Environment","Control_LambdaTe_U" , "Control_LambdaTu_U","Control_LambdaIntraTe_U", "Control_LambdaIntraTu_U","Control_ALLTe_U" ,"Control_ALLTu_U",  "predTe_onlyLambda_U","predTe_lambda_INTRA_U" ,"predTe_ALL_U" ,"predTu_onlyLambda_U","predTu_lambda_INTRA_U","predTu_ALL_U" )], parameter, value_U, c("Control_LambdaTe_U" , "Control_LambdaTu_U","Control_LambdaIntraTe_U", "Control_LambdaIntraTu_U","Control_ALLTe_U" ,"Control_ALLTu_U",  "predTe_onlyLambda_U","predTe_lambda_INTRA_U" ,"predTe_ALL_U" ,"predTu_onlyLambda_U","predTu_lambda_INTRA_U","predTu_ALL_U") )

pred_coex1Gen_long_U$parameter2<-mapvalues(pred_coex1Gen_long_U$parameter,c("Control_LambdaTe_U" , "Control_LambdaTu_U","Control_LambdaIntraTe_U", "Control_LambdaIntraTu_U","Control_ALLTe_U" ,"Control_ALLTu_U",  "predTe_onlyLambda_U","predTe_lambda_INTRA_U" ,"predTe_ALL_U" ,"predTu_onlyLambda_U","predTu_lambda_INTRA_U","predTu_ALL_U"), c("Control_lambdaTe" , "Control_lambdaTu","Control_lambdaIntraTe", "Control_lambdaIntraTu","Control_ALLTe" ,"Control_ALLTu",  "predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL") )

pred_coex1Gen_long$parameter2<-pred_coex1Gen_long$parameter

pred_coex1Gen_long<-left_join(pred_coex1Gen_long, pred_coex1Gen_long_L, by=c("Te","Tu", "parameter2","Environment"))

pred_coex1Gen_long<-left_join(pred_coex1Gen_long, pred_coex1Gen_long_U, by=c("Te","Tu", "parameter2","Environment"))

colnames(pred_coex1Gen_long)<-c("Te", "Tu", "Environment", "parameter", "value", "parameter2","parameter_L", "value_L", "parameter_U", "value_U" )

pred_coex1Gen_long$parameter3<-factor(pred_coex1Gen_long$parameter, c("Control_lambdaTe" , "Control_lambdaTu","Control_lambdaIntraTe", "Control_lambdaIntraTu","Control_ALLTe" ,"Control_ALLTu",  "predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL"))
str(pred_coex1Gen)


### Figure1

ggplot(subset(pred_coex1Gen_long, parameter=="predTe_onlyLambda" |  parameter=="predTe_Lambda_INTRA" |  parameter=="predTe_ALL"), aes(y=parameter3, x=value))+
  facet_grid(.~Environment, labeller=labeller(Environment=Env))+
  geom_errorbarh(aes(xmin=value_L, xmax=value_U, group=interaction(Te, Tu)), colour="black", height=0.5, position=position_dodge2(0.5))+
  geom_point(aes(fill=interaction(Te, Tu)),size=2.5, position=position_dodge2(0.5), stat="identity", shape=21)+
  geom_vline(xintercept = 1, colour="lightgray", linetype="dashed")+
  theme_bw()+
  theme_plot+
  scale_fill_brewer(palette = "Spectral", labels=c("Te no cadmium:Tu no cadmium", "Te cadmium:Tu no cadmium", "Te no cadmium:Tu cadmium", "Te cadmium:Tu cadmium"), name="")+
  guides(fill=guide_legend(nrow=2))+
  xlab(expression(paste("Predicted offspring production for ", italic("T. evansi"))))+
  scale_y_discrete(labels=c(expression(lambda+ alpha[ii] + alpha [ij]),expression(lambda+ alpha[ii]), expression(lambda)), limits=rev(levels(droplevels(subset(pred_coex1Gen_long, parameter=="predTe_onlyLambda" |  parameter=="predTe_Lambda_INTRA" |  parameter=="predTe_ALL"))$parameter3)))+
  theme(legend.position = "bottom", axis.text = element_text(size=12), axis.title = element_text(face="plain", size=12))+
  ylab("")
save_plot("./Plots/Fig1A.pdf", width=17.5, height=10)



ggplot(subset(pred_coex1Gen_long, parameter=="predTu_onlyLambda" |  parameter=="predTu_Lambda_INTRA" |  parameter=="predTu_ALL"), aes(y=parameter3, x=value))+
  facet_grid(.~Environment, labeller=labeller(Environment=Env))+
  geom_errorbarh(aes(xmin=value_L, xmax=value_U, group=interaction(Te, Tu)), colour="black", height=0.5, position=position_dodge2(0.5))+
  geom_point(aes(fill=interaction(Te, Tu)),size=2.5, position=position_dodge2(0.5), stat="identity", shape=21)+
  geom_vline(xintercept = 1, colour="lightgray", linetype="dashed")+
  theme_bw()+
  theme_plot+
  scale_fill_brewer(palette = "Spectral", labels=c("Te no cadmium:Tu no cadmium", "Te cadmium:Tu no cadmium", "Te no cadmium:Tu cadmium", "Te cadmium:Tu cadmium"), name="")+
  xlab(expression(paste("Predicted offspring production for ", italic("T. urticae"))))+
  guides(fill=guide_legend(nrow=2))+
  scale_y_discrete(labels=c(expression(lambda+ alpha[ii] + alpha [ij]),expression(lambda+ alpha[ii]), expression(lambda)), limits=rev(levels(droplevels(subset(pred_coex1Gen_long, parameter=="predTu_onlyLambda" |  parameter=="predTu_Lambda_INTRA" |  parameter=="predTu_ALL"))$parameter3)))+
  theme(legend.position = "bottom", axis.text = element_text(size=12), axis.title = element_text(face="plain", size=12))+
  ylab("")
save_plot("./Plots/Fig1B.pdf", width=17.5, height=10)
