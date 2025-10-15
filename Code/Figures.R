#' ---
#' title: "R Notebook"
#' output:
#'   html_document:
#'     df_print: paged
#' ---
## ---------------------------
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

theme_plots<-theme(axis.text = element_text(size=14), axis.title = element_text(size=14, face="bold"), legend.text = element_text(size=12), strip.text = element_text(size=14), plot.title = element_text(size=14, face="bold"), panel.grid=element_line(colour="white"), panel.background = element_rect(fill="white") , axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "black"), strip.background = element_rect(fill="white"))

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

dir.create("./Plots")

#' 
#' # Figures
#' 
#' This script has all the necessary code to redo the figures from the manuscript "Evolution in response to an abiotic stress shapes species coexistence".
#' The only figure that is missing is figure S2, which is created through the file: Method_comparison.R
#' 
#' 
#' # Figure 1
#' 
#' ### Importing data
#' Importing the data to get the parameter estimates for pooled data
## ---------------------------
param_all_REP<-read.csv("./Analyses/cxr_normal_REP_allEqual/parameters_cxr_normal_REP_best.csv")
param_all_REP_upper<-read.csv("./Analyses/cxr_normal_REP_allEqual/parameters_cxr_normal_REP_upper_best.csv")
param_all_REP_lower<-read.csv( "./Analyses/cxr_normal_REP_allEqual/parameters_cxr_normal_REP_lower_best.csv")
param_all_REP<-param_all_REP[,-1]
param_all_REP_upper<-param_all_REP_upper[,-1]
param_all_REP_lower<-param_all_REP_lower[,-1]

#' 
#' Importing the data to get the parameter estimates for each replicate
## ---------------------------
param_all_w0<-read.csv("./Analyses/cxr_normal_allEqual/parameters_cxr_normal_new.csv")
param_all_w0_upper<-read.csv("./Analyses/cxr_normal_allEqual/parameters_cxr_normal_upper_new.csv")
param_all_w0_lower<-read.csv( "./Analyses/cxr_normal_allEqual/parameters_cxr_normal_lower_new.csv")

param_all_w0<-param_all_w0[,-1]
param_all_w0_upper<-param_all_w0_upper[,-1]
param_all_w0_lower<-param_all_w0_lower[,-1]

#' 
#' ### Simulations for figure 1
#' Using the parameter estimates we will predict the number offspring produced by 10 females when considering only the effect of growth rate, growth rate + intraspecific competition, growth rate + intra and interspecific competition.
#' We will use the same approach for each of the two species and for upper and lower bounds.
#' 
## ---------------------------
# create the data frame to store the data
pred_coex1Gen<-as.data.frame(expand_grid(Te=c("SR4","SR5"), Tu=c("SR1", "SR2"), Environment= c("N", "Cd")))

# Loop to create the prediction for the number of offspring for TU considering only growth rate
pred_coex1Gen$predTu_onlyLambda<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  # subset the data frame corresponding to the regimes 
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  #Predicting number of offspring with only growth rate
  bl<-aux_alphas$Tu_lambda[1]*10
  
  bl
})
# Loop to create the prediction for the number of offspring for TU considering growth rate + intraspecific competition
pred_coex1Gen$predTu_Lambda_INTRA<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  # subset the data frame corresponding to the regimes 
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  #Predicting number of offspring with growth rate and the intraspecific competition following the riker model
  bl<-aux_alphas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10)
  
  bl
})

# Loop to create the prediction for the number of offspring for TU considering growth rate + intraspecific + interspecific competition
pred_coex1Gen$predTu_ALL<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  # subset data to get the information for the selection regimes
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  #Predicting number of offspring with growth rate and the intra and interspecific competition following the riker model
  bl<-aux_alphas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10- aux_alphas$Tu_inter[1]*10)
  
  bl
})

# Same loop as above but now for Te
# Loop to create the prediction for the number of offspring for TE considering growth rate 
pred_coex1Gen$predTe_onlyLambda<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10
  
  bl
})

# Loop to create the prediction for the number of offspring for TE considering growth rate + intraspecific  competition
pred_coex1Gen$predTe_Lambda_INTRA<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10)
  
  bl
})

# Loop to create the prediction for the number of offspring for TE considering growth rate + intraspecific + interspecific competition
pred_coex1Gen$predTe_ALL<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10- aux_alphas$Te_inter[1]*10)
  
  bl
})

#' 
#' ### Predictions of offspring number considering lower estimates
#' Here we have the same approach but using the lower bounds.
#' 
## ---------------------------
# Loop to create the prediction for the number of offspring for TU considering only growth rate
pred_coex1Gen$predTu_onlyLambda_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Tu_lambda[1]*10
  
  bl
})
# Loop to create the prediction for the number of offspring for TU considering growth rate + intraspecific competition
pred_coex1Gen$predTu_Lambda_INTRA_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_lambdas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_lambdas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10)
  
  bl
})
# Loop to create the prediction for the number of offspring for TU considering growth rate + intraspecific + interspecific competition
pred_coex1Gen$predTu_ALL_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_lambdas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_lambdas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10- aux_alphas$Tu_inter[1]*10)
  
  bl
})

# Loop to create the prediction for the number of offspring for TE considering growth rate
pred_coex1Gen$predTe_onlyLambda_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10
  
  bl
})
# Loop to create the prediction for the number of offspring for TE considering growth rate + intraspecific competition
pred_coex1Gen$predTe_Lambda_INTRA_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_lambdas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_lambdas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10)
  
  bl
})
# Loop to create the prediction for the number of offspring for TE considering growth rate + intraspecific + interspecific competition
pred_coex1Gen$predTe_ALL_L<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_lambdas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_lambdas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10- aux_alphas$Te_inter[1]*10)
  
  bl
})


#' 
#' ###Predictions of offspring number considering upper estimates
#' Here we have the same approach but using the lower bounds.
#' 
## ---------------------------
# Loop to create the prediction for the number of offspring for TU considering only growth rate
pred_coex1Gen$predTu_onlyLambda_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Tu_lambda[1]*10
  
  bl
})
# Loop to create the prediction for the number of offspring for TU considering growth rate + intraspecific competition
pred_coex1Gen$predTu_lambda_INTRA_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_Lambdas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_Lambdas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10)
  
  bl
})
# Loop to create the prediction for the number of offspring for TU considering growth rate + intraspecific + interspecific competition
pred_coex1Gen$predTu_ALL_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_Lambdas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_Lambdas$Tu_lambda[1]*10*exp(-aux_alphas$Tu_intra[1]*10- aux_alphas$Tu_inter[1]*10)
  
  bl
})

# Loop to create the prediction for the number of offspring for TE considering growth rate 
pred_coex1Gen$predTe_onlyLambda_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_alphas$Te_lambda[1]*10
  
  bl
})

# Loop to create the prediction for the number of offspring for TE considering growth rate + intraspecific competition
pred_coex1Gen$predTe_lambda_INTRA_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_Lambdas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_Lambdas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10)
  
  bl
})

# Loop to create the prediction for the number of offspring for TE considering growth rate + intraspecific + interspecific competition
pred_coex1Gen$predTe_ALL_U<-sapply(c(1:length(pred_coex1Gen$Tu)), function(x){
  aux_alphas<-subset(param_all_REP_lower, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  aux_Lambdas<-subset(param_all_REP_upper, Tu_Regime==as.character(pred_coex1Gen$Tu[x]) & Te_Regime==as.character(pred_coex1Gen$Te[x]) & Environment==as.character(pred_coex1Gen$Environment[x]))
  
  bl<-aux_Lambdas$Te_lambda[1]*10*exp(-aux_alphas$Te_intra[1]*10- aux_alphas$Te_inter[1]*10)
  
  bl
})

#' 
#' #### Reshaping the data frame
#' This is to make it easier to plot the predictions
#' 
## ---------------------------
str(pred_coex1Gen)
# reshaping data of the mean estimates
pred_coex1Gen_long<-gather(pred_coex1Gen[,c("Te","Tu","Environment",  "predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL" )], parameter, value, c( "predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL" ))

# reshaping data of the lower bound estimates
pred_coex1Gen_long_L<-gather(pred_coex1Gen[,c("Te","Tu","Environment", "predTe_onlyLambda_L","predTe_Lambda_INTRA_L" ,"predTe_ALL_L" ,"predTu_onlyLambda_L","predTu_Lambda_INTRA_L","predTu_ALL_L" )], parameter, value_L, c( "predTe_onlyLambda_L","predTe_Lambda_INTRA_L" ,"predTe_ALL_L" ,"predTu_onlyLambda_L","predTu_Lambda_INTRA_L","predTu_ALL_L"))

#changing the names of the parameters so that they match between upper and lower bounds and the mean estimates
pred_coex1Gen_long_L$parameter2<-mapvalues(pred_coex1Gen_long_L$parameter,c("predTe_onlyLambda_L","predTe_Lambda_INTRA_L" ,"predTe_ALL_L" ,"predTu_onlyLambda_L","predTu_Lambda_INTRA_L","predTu_ALL_L"), c("predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL") )

# reshaping data of the upper estimates
pred_coex1Gen_long_U<-gather(pred_coex1Gen[,c("Te","Tu","Environment","predTe_onlyLambda_U","predTe_lambda_INTRA_U" ,"predTe_ALL_U" ,"predTu_onlyLambda_U","predTu_lambda_INTRA_U","predTu_ALL_U" )], parameter, value_U, c( "predTe_onlyLambda_U","predTe_lambda_INTRA_U" ,"predTe_ALL_U" ,"predTu_onlyLambda_U","predTu_lambda_INTRA_U","predTu_ALL_U") )

#changing the names of the parameters so that they match between upper and lower bounds and the mean estimates
pred_coex1Gen_long_U$parameter2<-mapvalues(pred_coex1Gen_long_U$parameter,c("predTe_onlyLambda_U","predTe_lambda_INTRA_U" ,"predTe_ALL_U" ,"predTu_onlyLambda_U","predTu_lambda_INTRA_U","predTu_ALL_U"), c("predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL") )

# creating the column to make it easier to plot
pred_coex1Gen_long$parameter2<-pred_coex1Gen_long$parameter

# Joining the lower and upper data frames to the mean data frame using the parameter column to join
pred_coex1Gen_long<-left_join(pred_coex1Gen_long, pred_coex1Gen_long_L, by=c("Te","Tu", "parameter2","Environment"))

pred_coex1Gen_long<-left_join(pred_coex1Gen_long, pred_coex1Gen_long_U, by=c("Te","Tu", "parameter2","Environment"))
# changing the names of the columns
colnames(pred_coex1Gen_long)<-c("Te", "Tu", "Environment", "parameter", "value", "parameter2","parameter_L", "value_L", "parameter_U", "value_U" )

# Transform to factor
pred_coex1Gen_long$parameter3<-factor(pred_coex1Gen_long$parameter, c("predTe_onlyLambda","predTe_Lambda_INTRA" ,"predTe_ALL" ,"predTu_onlyLambda","predTu_Lambda_INTRA","predTu_ALL"))



#' 
#' ## Plotting
#' 
#' Plot with predictions of the number of offspring produced by Te and Tu when considering different levels competition (no competition - only lambda, only intraspecific competition - lambda_Intra, intra and interspecific competition - ALL)
## ---------------------------
# This is for T. evansi
ggplot(subset(pred_coex1Gen_long, parameter=="predTe_onlyLambda" |  parameter=="predTe_Lambda_INTRA" |  parameter=="predTe_ALL"), aes(y=parameter3, x=value))+
  facet_grid(.~Environment, labeller=labeller(Environment=Env))+
  geom_errorbarh(aes(xmin=value_L, xmax=value_U, group=interaction(Te, Tu)), colour="black", height=0.5, position=position_dodge2(0.5))+
  geom_point(aes(fill=interaction(Te, Tu)),size=2.5, position=position_dodge2(0.5), stat="identity", shape=21)+
  geom_vline(xintercept = 1, colour="lightgray", linetype="dashed")+
  theme_bw()+
  theme_plots+
  scale_fill_brewer(palette = "Spectral", labels=c("Te no cadmium:Tu no cadmium", "Te cadmium:Tu no cadmium", "Te no cadmium:Tu cadmium", "Te cadmium:Tu cadmium"), name="")+
  guides(fill=guide_legend(nrow=2))+
  xlab(expression(paste("Predicted offspring production for ", italic("T. evansi"))))+
  scale_y_discrete(labels=c(expression(lambda+ alpha[ii] + alpha [ij]),expression(lambda+ alpha[ii]), expression(lambda)), limits=rev(levels(droplevels(subset(pred_coex1Gen_long, parameter=="predTe_onlyLambda" |  parameter=="predTe_Lambda_INTRA" |  parameter=="predTe_ALL"))$parameter3)))+
  theme(legend.position = "bottom", axis.text = element_text(size=12), axis.title = element_text(face="plain", size=12))+
  ylab("")
 save_plot("./Plots/Fig1A.pdf", width=17.5, height=10)

#' 
#' 
## ---------------------------
# This is for T. urticae
ggplot(subset(pred_coex1Gen_long, parameter=="predTu_onlyLambda" |  parameter=="predTu_Lambda_INTRA" |  parameter=="predTu_ALL"), aes(y=parameter3, x=value))+
  facet_grid(.~Environment, labeller=labeller(Environment=Env))+
  geom_errorbarh(aes(xmin=value_L, xmax=value_U, group=interaction(Te, Tu)), colour="black", height=0.5, position=position_dodge2(0.5))+
  geom_point(aes(fill=interaction(Te, Tu)),size=2.5, position=position_dodge2(0.5), stat="identity", shape=21)+
  geom_vline(xintercept = 1, colour="lightgray", linetype="dashed")+
  theme_bw()+
  theme_plots+
  scale_fill_brewer(palette = "Spectral", labels=c("Te no cadmium:Tu no cadmium", "Te cadmium:Tu no cadmium", "Te no cadmium:Tu cadmium", "Te cadmium:Tu cadmium"), name="")+
  xlab(expression(paste("Predicted offspring production for ", italic("T. urticae"))))+
  guides(fill=guide_legend(nrow=2))+
  scale_y_discrete(labels=c(expression(lambda+ alpha[ii] + alpha [ij]),expression(lambda+ alpha[ii]), expression(lambda)), limits=rev(levels(droplevels(subset(pred_coex1Gen_long, parameter=="predTu_onlyLambda" |  parameter=="predTu_Lambda_INTRA" |  parameter=="predTu_ALL"))$parameter3)))+
  theme(legend.position = "bottom", axis.text = element_text(size=12), axis.title = element_text(face="plain", size=12))+
  ylab("")
save_plot("./Plots/Fig1B.pdf", width=17.5, height=10)


#' # Figure 2
#' 
#' ### Importing parameters
#' 
## ---------------------------
# Importing parameters
test_struct<-read.csv("./Analyses/structural_REP_new.csv")
# removing the first column
test_struct<-test_struct[,-1]


#' 
## ---------------------------
# Since we are only focusing on positive growth rates, we transform all the negative slopes that would go to y <0 into 0 and all negative slopes that would go to x<0 into 90
test_struct[which(test_struct$a21_a11<0),"a21_a11"]<-0
test_struct[which(test_struct$a21_a11_lower<0),"a21_a11_lower"]<-0
test_struct[which(test_struct$a21_a11_upper<0),"a21_a11_upper"]<-0

test_struct[which(test_struct$a22_a12<0),"a22_a12"]<-90
test_struct[which(test_struct$a22_a12_lower<0),"a22_a12_lower"]<-90
test_struct[which(test_struct$a22_a12_upper<0),"a22_a12_upper"]<-90

#Since we will draw polygons with ggplot we need to have x and y coordinates for the 4 vertex

# This is just to guide the plotting. The values 0 and 3 correspond to the same point where the feasibility cone starts and ends. Then we have to calculate the coordinates for the middle
values_to_plot<-c(0,1, 2,3)
max_value<-8 # this value should be larger than the growth rate

  #All of the x start the same, but after calculating the y values, we will need to adjust the x so that the negative slopes in Tu are set to 0
# Creating a data frame to store the values
  CI_rep2<-as.data.frame(expand.grid(values=values_to_plot, Te_Regime=c("SR4", "SR5"), Tu_Regime=c("SR1", "SR2"), Environment=c("N", "Cd")))

  # Adding the x values corresponding to the coordinates for the mean, lower and upper intervals
CI_rep2$x<-mapvalues(CI_rep2$values, c(0,1,2,3), c(0,max_value,max_value,0))
CI_rep2$x_lower<-mapvalues(CI_rep2$values, c(0,1,2,3), c(0,max_value,max_value,0))
CI_rep2$x_upper<-mapvalues(CI_rep2$values, c(0,1,2,3), c(0,max_value,max_value,0))

# Loop to create the data frame where we calculate the vertex by multiplying the value of x by the slope (i.e. the edge of the feasibility domain)
aux_CI2<-as.data.frame(t(sapply(c(1:nrow(CI_rep2)), function(x){
  a<-subset(test_struct, Environment==CI_rep2$Environment[x] &
              Te_Regime==CI_rep2$Te_Regime[x] & Tu_Regime==CI_rep2$Tu_Regime[x])
  
  # If values are 0 or 3 then the points are always 0
  if(CI_rep2$values[x]==0 | CI_rep2$values[x]==3){
    b<-c(0,0,0)
    # If its 2 then we need to use the urticae edge
  }else if(CI_rep2$values[x]==2){
    vect<-c(a$a22_a12_lower[1],a$a22_a12[1],a$a22_a12_upper[1])
      # ordering to be sure that it matches (because facilitation switches values)
    vect<-vect[order(vect, decreasing=TRUE)]
    b<-c(CI_rep2$x_lower[x],CI_rep2$x[x], CI_rep2$x_upper[x])*vect
  }else if(CI_rep2$values[x]==1){ # If its 1 then we need to use the evansi edge
    vect<-c(a$a21_a11_lower[1],a$a21_a11[1],a$a21_a11_upper[1])
    # ordering to be sure that it matches (because facilitation switches values)
    vect<-vect[order(vect, decreasing=FALSE)]
    b<-c(CI_rep2$x_lower[x],CI_rep2$x[x], CI_rep2$x_upper[x])*vect
  }
  
  b
})))
# addin column names
colnames(aux_CI2)<-c("y_lower", "y", "y_upper")

#Joininig the data frame
CI_rep2<-as.data.frame(cbind(CI_rep2, aux_CI2))


#' 
#' #### no cadmium
#' 
## ---------------------------
normal_feas<-ggplot(subset(test_struct, Environment=="N"), aes(x=Tu_lambda, y=Te_lambda, colour=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Tu_Regime~ Te_Regime, labeller=labeller(Tu_Regime=regimeTu, Te_Regime=regimeTe) )+
  geom_polygon(data=subset(CI_rep2,Environment=="N"),aes(x=x_upper, y=y_upper), fill="#FADADD",  linewidth=0.85, colour=NA)+
  geom_polygon(data=subset(CI_rep2,Environment=="N"),aes(x=x, y=y), fill="#969696",  linewidth=0.85, colour="black")+
  geom_polygon(data=subset(CI_rep2,Environment=="N"),aes(x=x_lower, y=y_lower), fill="lightskyblue1",  linewidth=0.85, colour=NA)+
  geom_polygon(data=subset(CI_rep2,Environment=="N" & Te_Regime=="SR5" & Tu_Regime=="SR2"),aes(x=x, y=y), fill="#969696", alpha=0.75, linewidth=0.85, colour="black")+
  geom_point(colour="black")+
  geom_segment(data=subset(test_struct, Environment=="N" ), aes(xend=Tu_lambda, yend=Te_lambda,x=0, y=0),  arrow=arrow(length = unit(0.3, "cm")), colour="black", linewidth=1)+
  geom_segment(data=subset(test_struct, Environment=="N" ), aes(xend=Tu_lambda_lower, yend=Te_lambda_lower,x=0, y=0),  colour="black", linetype="dashed", linewidth=0.75)+
  geom_segment(data=subset(test_struct, Environment=="N" ), aes(xend=Tu_lambda_upper, yend=Te_lambda_upper,x=0, y=0), colour="black", linetype="dashed", linewidth=0.75)+
  theme_bw()+
  theme_plots+
  ylab(c("Intrinsic growth rate T. evansi"))+
  xlab(c("Intrinsic growth rate T. urticae"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  ylim(c(-200,200))+
  xlim(c(-200,200))+
  coord_cartesian(xlim =c(0.01,4), ylim=c(0.01,8), expand = TRUE)+
  ggtitle("No cadmium environment")
normal_feas
save_plot("./Plots/Fig2B_B.pdf", width=20, height=15)

#' 
#' #### cadmium
## ---------------------------
cadmium_feas<-ggplot(subset(test_struct, Environment=="Cd"), aes(x=Tu_lambda, y=Te_lambda, colour=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Tu_Regime~ Te_Regime, labeller=labeller(Tu_Regime=regimeTu, Te_Regime=regimeTe) )+
  geom_polygon(data=subset(CI_rep2,Environment=="Cd"),aes(x=x_lower, y=y_lower), fill="#FADADD",linewidth=0.85, colour=NA)+
  geom_polygon(data=subset(CI_rep2,Environment=="Cd"),aes(x=x, y=y), fill="#969696", linewidth=0.85,colour="black")+
  geom_polygon(data=subset(CI_rep2,Environment=="Cd"),aes(x=x_upper, y=y_upper), fill="lightskyblue1",  linewidth=0.85, colour=NA)+
  geom_point(colour="black")+
  geom_segment(data=subset(test_struct, Environment=="Cd" ), aes(xend=Tu_lambda, yend=Te_lambda,x=0, y=0),  arrow=arrow(length = unit(0.3, "cm")), colour="black", linewidth=1)+
  geom_segment(data=subset(test_struct, Environment=="Cd" ), aes(xend=Tu_lambda_lower, yend=Te_lambda_lower,x=0, y=0),  colour="black", linetype="dashed", linewidth=0.75)+
  geom_segment(data=subset(test_struct, Environment=="Cd" ), aes(xend=Tu_lambda_upper, yend=Te_lambda_upper,x=0, y=0), colour="black", linetype="dashed", linewidth=0.75)+
  theme_bw()+
  theme_plots+
  ylab(c("Intrinsic growth rate T. evansi"))+
  xlab(c("Intrinsic growth rate T. urticae"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  #ylim(c(0,8))+
  #xlim(c(0,5))+
  coord_cartesian(xlim =c(0.1,2), ylim=c(0.1,2.4), expand = TRUE)+
  ggtitle("Cadmium environment")

cadmium_feas
save_plot("./Plots/Fig2A_B.pdf", width=20, height=15)

#' 
#' #### both together
## ---------------------------
# joining the two plots
plot_grid(cadmium_feas, normal_feas, ncol=2, labels=c("A", "B") )
# save plots
save_plot("./Plots/Fig2.pdf", width=30, height=15)
save_plot("./Plots/Fig2.tiff", width=30, height=15)

#' 
#' # Figure 3
#' 
#' ### Importing parameters
#' 
## ---------------------------
struct_mat_REP_final2<-read.csv("./Analyses/min_distance_pooled.csv", header=TRUE)
struct_mat_w0_final<-read.csv("./Analyses/min_distance_per_replicate.csv", header=TRUE)


#' 
## ---------------------------
ggplot(struct_mat_w0_final, aes(x=interaction(Te_Regime, Tu_Regime), y=minDistance, fill=interaction(Te_Regime, Tu_Regime)))+
  facet_grid(.~Environment, labeller=labeller(Environment=Env))+
  geom_hline(yintercept = 0, colour="lightgrey", linetype="dashed")+
  geom_errorbar(data=struct_mat_w0_final, aes(ymin=minDistance_L, ymax=minDistance_U), alpha=0.35, position=position_dodge2(0.5), width=0.5)+
  geom_point(size=2, shape=24,position=position_dodge2(0.5), alpha=0.65)+
  geom_errorbar(data=struct_mat_REP_final2, aes(ymin=minDistance_L, ymax=minDistance_U), width=0.25)+
  geom_point(data=struct_mat_REP_final2,aes(x=interaction(Te_Regime, Tu_Regime), y=minDistance, fill=interaction(Te_Regime, Tu_Regime)), colour="black", size=2.5, shape=21 )+
  theme_bw()+
  theme_plots+
  scale_x_discrete(labels=c("Te no cadmium\nTu no cadmium","Te cadmium\nTu no cadmium","Te no cadmium\nTu cadmium","Te cadmium\nTu cadmium"))+
  #scale_color_manual(values=c("darkred", "darkblue"), labels=c("Cadmium", "Water"))+
  ylab("Distance to \n the closest edge")+
  scale_fill_brewer(palette = "Spectral", labels=c("Te no cadmium:Tu no cadmium", "Te cadmium:Tu no cadmium", "Te no cadmium:Tu cadmium", "Te cadmium:Tu cadmium"), name="")+
  xlab(expression(paste("Predicted offspring production for ", italic("T. urticae"))))+
  guides(fill=guide_legend(nrow=2))+
  xlab("Selection Regimes")+
  theme(legend.position = "bottom", strip.background =element_rect(colour="white"), panel.border = element_rect(colour="black"), axis.text.x = element_text(size=10, angle = 45, vjust = 0.63), axis.title = element_text(size=10, face="bold"), strip.text = element_text(size=10), axis.text.y = element_text(size=10))
save_plot("./Plots/Fig3.pdf", width=15, height=10)
save_plot("./Plots/Fig3.png", width=15, height=10)

#' 
#' # Figure 4
#' 
#' ### Importing data
#' 
## ---------------------------
sum_observed_coex_rep<-read.csv(file="./Analyses/popDyn_pooledData.csv", header=TRUE)

sum_observed_coex_ALL2<-read.csv(file="./Analyses/popDyn_replicates.csv", header=TRUE)

#' 
#' ### slopes and figure
## ---------------------------
# This is from the summary of the model so we can put the correct slope and CI in the plot, and not ggplot's model
slope_all<-log(1.91) #slope m7 obtained from emmeans
slope_ci_L<-log(1.91 - 0.156) 
slope_ci_U<-log(1.91 + 0.156) 

# Predicting the data for the slope
slopes_data<-data.frame(x=seq(0,1, 0.05))
slopes_data$ymin<-slope_ci_L*slopes_data$x
slopes_data$ymax<-slope_ci_U*slopes_data$x

# Figure
ggplot(sum_observed_coex_ALL2)+
  geom_abline(intercept=0, slope=1, color="black", linewidth=0.75, linetype="dashed")+
  geom_abline(intercept = 0, slope=slope_all)+
  geom_abline(intercept = 0, slope=slope_ci_L, colour="grey",linewidth=0.75)+
  geom_abline(intercept = 0, slope=slope_ci_U,  colour="grey",linewidth=0.75)+
  geom_ribbon(data=slopes_data, aes(ymin=ymin,ymax=ymax, x=x), fill="grey", alpha=0.5) +
  geom_errorbar(data=sum_observed_coex_ALL2,aes(x=meanRatio,y=mean_pred,ymin=mean_pred-sdPred, ymax=mean_pred+sdPred), width=0.01, colour="black")+
  geom_errorbarh(data=sum_observed_coex_ALL2,aes(y=mean_pred, xmin=meanRatio-sdRatio, xmax=meanRatio+sdRatio), height=0.01, colour="black")+
  geom_point(size=3, aes(x=meanRatio,y=mean_pred,fill=interaction(SRTu, SRTe), shape=Env))+
  scale_fill_manual(values=c("#D7191C", "#FDAE61" ,"#ABDDA4", "#2B83BA"), labels=c("Te no cadmium:Tu no cadmium", "Te cadmium: Tu no cadmium", "Te no cadmium: Tu cadmium", "Te cadmium: Tu cadmium"))+
  scale_shape_manual(values=c(22,23))+
  theme_plots+
  xlab("Observed ratio")+
  ylab("Predicted ratio")+
  ylim(c(-2,2))+
  xlim(c(-2,2))+
  theme(legend.position = "none")+
  coord_cartesian(xlim=c(0.3,0.96), ylim=c(0.3, 0.96))

save_plot("./Plots/Fig4.pdf", width=10, height=10)

save_plot("./Plots/Fig4.png", width=15, height=12)


#' 
#' # Supplementary figures
#' 
#' ## Figures S3 - S8
#' 
#' ### Importing data
#' Importing the data to get the parameter estimates for pooled data
## ---------------------------
param_all_REP<-read.csv("./Analyses/cxr_normal_REP_allEqual/parameters_cxr_normal_REP_best.csv")
param_all_REP_upper<-read.csv("./Analyses/cxr_normal_REP_allEqual/parameters_cxr_normal_REP_upper_best.csv")
param_all_REP_lower<-read.csv( "./Analyses/cxr_normal_REP_allEqual/parameters_cxr_normal_REP_lower_best.csv")
param_all_REP<-param_all_REP[,-1]
param_all_REP_upper<-param_all_REP_upper[,-1]
param_all_REP_lower<-param_all_REP_lower[,-1]

#' 
#' Importing the data to get the parameter estimates for each replicate
## ---------------------------
param_all_w0<-read.csv("./Analyses/cxr_normal_allEqual/parameters_cxr_normal_new.csv")
param_all_w0_upper<-read.csv("./Analyses/cxr_normal_allEqual/parameters_cxr_normal_upper_new.csv")
param_all_w0_lower<-read.csv( "./Analyses/cxr_normal_allEqual/parameters_cxr_normal_lower_new.csv")

param_all_w0<-param_all_w0[,-1]
param_all_w0_upper<-param_all_w0_upper[,-1]
param_all_w0_lower<-param_all_w0_lower[,-1]

#' 
#' ### Data wrangling pooled data
## ---------------------------
param_all_REP_long<-gather(param_all_REP, parameter, value,Tu_lambda:Te_inter )

param_all_REP_long$category<-mapvalues(param_all_REP_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

param_all_REP_lower_long<-gather(param_all_REP_lower, parameter, value,Tu_lambda:Te_inter )

param_all_REP_lower_long$category<-mapvalues(param_all_REP_lower_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

param_all_REP_upper_long<-gather(param_all_REP_upper, parameter, value,Tu_lambda:Te_inter )

param_all_REP_upper_long$category<-mapvalues(param_all_REP_upper_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

colnames(param_all_REP_lower_long)[5]<-"lower"
colnames(param_all_REP_upper_long)[5]<-"upper"

str(param_all_REP_long)

param_all_REP_long<-cbind(param_all_REP_long[,1:6],param_all_REP_lower_long$lower, param_all_REP_upper_long$upper)

colnames(param_all_REP_long)[7:8]<-c("lower","upper")

#' 
#' ### Data wrangling replicate data
## ---------------------------
# data wrangling

param_all_w0_long<-gather(param_all_w0, parameter, value,Tu_lambda:Te_inter )

param_all_w0_long$category<-mapvalues(param_all_w0_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

param_all_w0_lower_long<-gather(param_all_w0_lower, parameter, value,Tu_lambda:Te_inter )

param_all_w0_lower_long$category<-mapvalues(param_all_w0_lower_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

param_all_w0_upper_long<-gather(param_all_w0_upper, parameter, value,Tu_lambda:Te_inter )

param_all_w0_upper_long$category<-mapvalues(param_all_w0_upper_long$parameter, c("Tu_lambda", "Te_lambda", "Tu_intra", "Te_intra","Tu_inter", "Te_inter"), c("lambda", "lambda", "intra", "intra", "inter", "inter"))

colnames(param_all_w0_lower_long)[6]<-"lower"
colnames(param_all_w0_upper_long)[6]<-"upper"

str(param_all_w0_long)

param_all_w0_long<-cbind(param_all_w0_long[,1:6],param_all_w0_lower_long$lower, param_all_w0_upper_long$upper)

colnames(param_all_w0_long)[7:8]<-c("lower","upper")
str(param_all_w0_long)

#' 
#' ## Cadmium env
#' 
#' #### Growth rate
#' 
## ---------------------------
### Tu
Tu_gr_ev<-ggplot(data=subset(param_all_REP_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="Cd")), aes(fill=Tu_Regime, y=value, x=Tu_Regime))+
  geom_hline(yintercept = 1, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="Cd")), aes(y=value, x=Tu_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="Cd")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=Tu_Regime, y=value, x=Tu_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="Cd")), aes(ymin=lower, ymax=upper), colour="black", width=0.1)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="Cd")), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("lightgrey", "black"), labels=c("No cadmium", "Cadmium"), name="")+
  ylab("Intrinsic growth rate\n(T. urticae)")+
  xlab("T. urticae selection regime")+theme(legend.position = "None")

## Te
Te_gr_ev<-ggplot(data=subset(param_all_REP_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="Cd")), aes(fill=Te_Regime, y=value, x=Te_Regime))+
  geom_hline(yintercept = 1, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="Cd")), aes(y=value, x=Te_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="Cd")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=Te_Regime, y=value, x=Te_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="Cd")), aes(ymin=lower, ymax=upper), colour="black", width=0.1)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="Cd")), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("lightgrey", "black"), labels=c("No cadmium", "Cadmium"), name="")+
  ylab("Intrinsic growth rate\n(T. evansi)")+
  xlab("T. evansi selection regime")+ theme(legend.position = "None")

#' 
#' #### Intraspecific
## ---------------------------
### Tu
Tu_intra_ev<-ggplot(data=subset(param_all_REP_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="Cd")), aes(fill=Tu_Regime, y=value, x=Tu_Regime))+
  geom_hline(yintercept = 0, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="Cd")), aes(y=value, x=Tu_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="Cd")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=Tu_Regime, y=value, x=Tu_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="Cd")), aes(ymin=lower, ymax=upper), colour="black", width=0.1)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="Cd")), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("lightgrey", "black"), labels=c("No cadmium", "Cadmium"), name="")+
  ylab("Strength of intraspecific\ncompetition (T. urticae)")+
  xlab("T. urticae selection regime")+theme(legend.position = "None")+
  ylim(c(-0.05,0.15))

### Te
Te_intra_ev<-ggplot(data=subset(param_all_REP_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="Cd")), aes(fill=Te_Regime, y=value, x=Te_Regime))+
  geom_hline(yintercept = 0, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="Cd")), aes(y=value, x=Te_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="Cd")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=Te_Regime, y=value, x=Te_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="Cd")), aes(ymin=lower, ymax=upper), colour="black", width=0.1)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="Cd")), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("lightgrey", "black"), labels=c("No cadmium", "Cadmium"), name="")+
  ylab("Strength of intraspecific\ncompetition (T. evansi)")+
  xlab("T. evansi selection regime")+theme(legend.position = "None")+
  ylim(c(-0.05,0.15))


#' 
#' #### Interspecific
## ---------------------------
Tu_inter_ev<-ggplot(data=subset(param_all_REP_long, (parameter=="Tu_inter" & Environment=="Cd")), aes(fill=interaction(Te_Regime,Tu_Regime), y=value, x=Tu_Regime))+
  geom_hline(yintercept = 0, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Tu_inter" & Environment=="Cd")), aes(y=value, x=Tu_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Tu_inter"  & Environment=="Cd")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=interaction(Te_Regime,Tu_Regime), y=value, x=Tu_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Tu_inter" & Environment=="Cd")), aes(ymin=lower, ymax=upper),position=position_dodge2(0.3), colour="black", width=0.3)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Tu_inter" & Environment=="Cd")),position=position_dodge2(0.3), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("#D7191C", "#FDAE61" ,"#ABDDA4", "#2B83BA"), labels=c("Te no cadmium:Tu no cadmium", "Te cadmium: Tu no cadmium", "Te no cadmium: Tu cadmium", "Te cadmium: Tu cadmium"), name="")+
  ylab("Strength of interspecific\ncompetition (T. urticae)")+
  xlab("T. urticae selection regime")+
  guides(fill=guide_legend(nrow = 4))+
  theme(legend.position = "None", legend.text = element_text(size=10), legend.background = element_rect(fill=NA), legend.key.size = unit(0.2, 'cm'))+
  ylim(c(-0.15,0.2))

### Te
Te_inter_ev<-ggplot(data=subset(param_all_REP_long, (parameter=="Te_inter" & Environment=="Cd")), aes(fill=interaction(Te_Regime,Tu_Regime), y=value, x=Te_Regime))+
  geom_hline(yintercept = 0, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Te_inter" & Environment=="Cd")), aes(y=value, x=Te_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Te_inter"  & Environment=="Cd")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=interaction(Te_Regime,Tu_Regime), y=value, x=Te_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Te_inter" & Environment=="Cd")), aes(ymin=lower, ymax=upper),position=position_dodge2(0.3), colour="black", width=0.3)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Te_inter" & Environment=="Cd")),position=position_dodge2(0.3), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("#D7191C", "#FDAE61" ,"#ABDDA4", "#2B83BA"), labels=c("Te no cadmium:Tu no cadmium", "Te cadmium: Tu no cadmium", "Te no cadmium: Tu cadmium", "Te cadmium: Tu cadmium"), name="")+
  guides(fill=guide_legend(nrow = 4))+
  theme(legend.position = "None", legend.text = element_text(size=10), legend.background = element_rect(fill=NA), legend.key.size = unit(0.2, 'cm'))+
  ylab("Strength of interspecific\ncompetition (T. evansi)")+
  xlab("T. evansi selection regime")+
  ylim(c(-0.15,0.2))

#' 
#' ### Plotting
## ---------------------------
plot_grid(Te_gr_ev, Tu_gr_ev, labels=c("A", "B"))

save_plot("./Plots/FigS7.pdf", width=30, height=15)
save_plot("./Plots/FigS7.png", width=30, height=15)

plot_grid(Te_intra_ev, Tu_intra_ev, labels=c("A", "B"))

save_plot("./Plots/FigS3.pdf", width=30, height=15)
save_plot("./Plots/FigS3.png", width=30, height=15)

plot_grid(Te_inter_ev, Tu_inter_ev, labels=c("A", "B"))

save_plot("./Plots/FigS4.pdf", width=30, height=15)
save_plot("./Plots/FigS4.png", width=30, height=15)

#' 
#' ## Ancestral env
#' 
#' #### Growth rate
## ---------------------------
#Tu
Tu_gr_ev_N<-ggplot(data=subset(param_all_REP_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="N")), aes(fill=Tu_Regime, y=value, x=Tu_Regime))+
  geom_hline(yintercept = 1, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="N")), aes(y=value, x=Tu_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="N")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=Tu_Regime, y=value, x=Tu_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="N")), aes(ymin=lower, ymax=upper), colour="black", width=0.1)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Tu_lambda" & Te_Regime=="SR4" & Environment=="N")), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("lightgrey", "black"), labels=c("No cadmium", "Cadmium"), name="")+
  ylab("Intrinsic growth rate\n(T. urticae)")+
  xlab("T. urticae selection regime")+theme(legend.position = "None")+
  ylim(c(0,17))

# Te

Te_gr_ev_N<-ggplot(data=subset(param_all_REP_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="N")), aes(fill=Te_Regime, y=value, x=Te_Regime))+
  geom_hline(yintercept = 1, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="N")), aes(y=value, x=Te_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="N")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=Te_Regime, y=value, x=Te_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="N")), aes(ymin=lower, ymax=upper), colour="black", width=0.1)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Te_lambda" & Tu_Regime=="SR1" & Environment=="N")), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("lightgrey", "black"), labels=c("No cadmium", "Cadmium"), name="")+
  ylab("Intrinsic growth rate\n(T. evansi)")+
  xlab("T. evansi selection regime")+ theme(legend.position = "None")+
  ylim(c(0,17))


#' 
#' #### Intraspecific competition
#' 
## ---------------------------
#Tu
Tu_intra_ev_N<-ggplot(data=subset(param_all_REP_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="N")), aes(fill=Tu_Regime, y=value, x=Tu_Regime))+
  geom_hline(yintercept = 0, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="N")), aes(y=value, x=Tu_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="N")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=Tu_Regime, y=value, x=Tu_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="N")), aes(ymin=lower, ymax=upper), colour="black", width=0.1)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Tu_intra" & Te_Regime=="SR4" & Environment=="N")), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("lightgrey", "black"), labels=c("No cadmium", "Cadmium"), name="")+
  ylab("Strength of intraspecific\ncompetition (T. urticae)")+
  xlab("T. urticae selection regime")+theme(legend.position = "None")+
  ylim(c(-0.01,0.16))

# Te
Te_intra_ev_N<-ggplot(data=subset(param_all_REP_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="N")), aes(fill=Te_Regime, y=value, x=Te_Regime))+
  geom_hline(yintercept = 0, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="N")), aes(y=value, x=Te_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="N")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=Te_Regime, y=value, x=Te_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="N")), aes(ymin=lower, ymax=upper), colour="black", width=0.1)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Te_intra" & Tu_Regime=="SR1" & Environment=="N")), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("lightgrey", "black"), labels=c("No cadmium", "Cadmium"), name="")+
  ylab("Strength of intraspecific\ncompetition (T. evansi)")+
  xlab("T. evansi selection regime")+theme(legend.position = "None")+
  ylim(c(-0.01,0.16))

#' 
#' #### Interspecific competition
## ---------------------------
Tu_inter_ev_N<-ggplot(data=subset(param_all_REP_long, (parameter=="Tu_inter" & Environment=="N")), aes(fill=interaction(Te_Regime,Tu_Regime), y=value, x=Tu_Regime))+
  geom_hline(yintercept = 0, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Tu_inter" & Environment=="N")), aes(y=value, x=Tu_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Tu_inter"  & Environment=="N")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=interaction(Te_Regime,Tu_Regime), y=value, x=Tu_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Tu_inter" & Environment=="N")), aes(ymin=lower, ymax=upper),position=position_dodge2(0.3), colour="black", width=0.3)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Tu_inter" & Environment=="N")),position=position_dodge2(0.3), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("#D7191C", "#FDAE61" ,"#ABDDA4", "#2B83BA"), labels=c("Te no cadmium:Tu no cadmium", "Te cadmium: Tu no cadmium", "Te no cadmium: Tu cadmium", "Te cadmium: Tu cadmium"), name="")+
  ylab("Strength of interspecific\ncompetition (T. urticae)")+
  xlab("T. urticae selection regime")+
  guides(fill=guide_legend(nrow = 4))+
  theme(legend.position = "None", legend.text = element_text(size=10), legend.background = element_rect(fill=NA), legend.key.size = unit(0.2, 'cm'))+
  ylim(c(-0.05,0.25))

# Te
Te_inter_ev_N<-ggplot(data=subset(param_all_REP_long, (parameter=="Te_inter" & Environment=="N")), aes(fill=interaction(Te_Regime,Tu_Regime), y=value, x=Te_Regime))+
  geom_hline(yintercept = 0, colour="lightgray", linetype="dashed")+
  geom_errorbar(data=subset(param_all_w0_long, (parameter=="Te_inter" & Environment=="N")), aes(y=value, x=Te_Regime, ymin=lower, ymax=upper), position=position_dodge2(width=0.3), colour="grey", width=0.3)+
  geom_point(data=subset(param_all_w0_long, (parameter=="Te_inter"  & Environment=="N")), position=position_dodge2(width=0.3), size=2, shape=24,aes(fill=interaction(Te_Regime,Tu_Regime), y=value, x=Te_Regime), alpha=0.85, colour="darkgrey")+
  geom_errorbar(data=subset(param_all_REP_long, (parameter=="Te_inter" & Environment=="N")), aes(ymin=lower, ymax=upper),position=position_dodge2(0.3), colour="black", width=0.3)+
  geom_point(data=subset(param_all_REP_long, (parameter=="Te_inter" & Environment=="N")),position=position_dodge2(0.3), size=4, shape=21)+
  theme_plots+
  scale_x_discrete(labels=c("No cadmium","Cadmium"))+
  scale_fill_manual(values=c("#D7191C", "#FDAE61" ,"#ABDDA4", "#2B83BA"), labels=c("Te no cadmium:Tu no cadmium", "Te cadmium: Tu no cadmium", "Te no cadmium: Tu cadmium", "Te cadmium: Tu cadmium"), name="")+
  guides(fill=guide_legend(nrow = 4))+
  theme(legend.position = "None", legend.text = element_text(size=10), legend.background = element_rect(fill=NA), legend.key.size = unit(0.2, 'cm'))+
  ylab("Strength of interspecific\ncompetition (T. evansi)")+
  xlab("T. evansi selection regime")

#' 
#' ### Plotting
## ---------------------------
plot_grid(Te_gr_ev_N, Tu_gr_ev_N, labels=c("A", "B"))

save_plot("./Plots/FigS8.pdf", width=30, height=15)
save_plot("./Plots/FigS8.png", width=30, height=15)

plot_grid(Te_intra_ev_N, Tu_intra_ev_N, labels=c("A", "B"))

save_plot("./Plots/FigS5.pdf", width=30, height=15)
save_plot("./Plots/FigS5.png", width=30, height=15)

plot_grid(Te_inter_ev_N, Tu_inter_ev_N, labels=c("A", "B"))

save_plot("./Plots/FigS6.pdf", width=30, height=15)
save_plot("./Plots/FigS6.png", width=30, height=15)

#' 
#' ### Figure S9
#' 
#' ### Importing parameters
#' 
## ---------------------------
struct_mat_REP_final2<-read.csv("./Analyses/min_distance_pooled.csv")
struct_mat_w0_final<-read.csv("./Analyses/min_distance_per_replicate.csv")


#' 
## ---------------------------
# first part of the plot
ggplot_distTe<-ggplot(struct_mat_REP_final2, aes(x=interaction(Tu_Regime, Te_Regime), fill=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Environment~., labeller=labeller(Environment=Env))+
  geom_hline(yintercept = 0, colour="lightgrey", linetype="dashed")+
  geom_errorbar(data=struct_mat_w0_final, aes(ymin=Te_distance_lower, ymax=Te_distance_upper), width=0.5, position=position_dodge2(0.5), alpha=0.25)+
  geom_point(data=struct_mat_w0_final, size=1.5, shape=24, position=position_dodge2(0.5), colour="black", alpha=0.55, aes(y=Te_distance))

# second part of the plot
ggplot_distTe<-ggplot_distTe +  geom_errorbar(data=struct_mat_REP_final2, aes(ymin=Te_distance_lower, ymax=Te_distance_upper, y=Te_distance), width=0.1, position=position_dodge(1))+
  geom_point(data=struct_mat_REP_final2, size=2, position=position_dodge(1), shape=21, colour="black", aes(y=Te_distance))+
  theme_bw()+
  theme_plots+
  scale_x_discrete(labels=c("Te no cadmium\nTu no cadmium","Te cadmium\nTu no cadmium","Te no cadmium\nTu cadmium","Te cadmium\nTu cadmium"))+
  ylab("Distance to \n the T. evansi edge")+
  scale_fill_brewer(palette = "Spectral", labels=c("Te no cadmium:Tu no cadmium", "Te cadmium:Tu no cadmium", "Te no cadmium:Tu cadmium", "Te cadmium:Tu cadmium"), name="")+
  guides(fill=guide_legend(nrow=2))+
  xlab("")+
  theme(legend.position = "bottom", strip.background =element_rect(colour="white"), panel.border = element_rect(colour="black"), axis.text.x = element_text(size=10, angle = 45, vjust = 0.63), strip.text = element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10, face="bold"))

ggplot_distTe

#' 
#' 
## ---------------------------
# first part of the plot
ggplot_distTu<-ggplot(struct_mat_REP_final2, aes(x=interaction(Tu_Regime, Te_Regime), fill=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Environment~., labeller=labeller(Environment=Env))+
  geom_hline(yintercept = 0, colour="lightgrey", linetype="dashed")+
  geom_errorbar(data=struct_mat_w0_final, aes(ymin=Tu_distance_lower, ymax=Tu_distance_upper), width=0.5, position=position_dodge2(width=0.5), alpha=0.35)+
  geom_point(data=struct_mat_w0_final, size=1.5, shape=24, position=position_dodge2(width=0.5), colour="black", alpha=0.55, aes(y=Tu_distance))

# second part of the plot
ggplot_distTu<-ggplot_distTu+ geom_errorbar(data=struct_mat_REP_final2, aes(ymin=Tu_distance_lower, ymax=Tu_distance_upper), width=0.15)+
  geom_point(data=struct_mat_REP_final2, size=2, shape=21, colour="black", aes(y=Tu_distance))+
  theme_bw()+
  theme_plots+
  scale_x_discrete(labels=c("Te no cadmium\nTu no cadmium","Te cadmium\nTu no cadmium","Te no cadmium\nTu cadmium","Te cadmium\nTu cadmium"))+
  ylab("Distance to \n the T. urticae edge")+
  scale_fill_brewer(palette = "Spectral", labels=c("Te no cadmium:Tu no cadmium", "Te cadmium:Tu no cadmium", "Te no cadmium:Tu cadmium", "Te cadmium:Tu cadmium"), name="")+
  guides(fill=guide_legend(nrow=2))+
  xlab("Te Selection Regime : Tu Selection Regime")+
  theme(legend.position = "bottom", strip.background =element_rect(colour="white"), panel.border = element_rect(colour="black"), axis.text.x = element_text(size=10, angle = 45, vjust = 0.63), strip.text = element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10, face="bold"))

ggplot_distTu

#' 
## ---------------------------
# joining the two plots
plot_grid(ggplot_distTe + theme(legend.position="none"),ggplot_distTu + theme(legend.position="none"), ncol=2, labels=c("A", "B") )

save_plot("./Plots/FigS9.pdf", width=25, height=15)
save_plot("./Plots/FigS9.png", width=25, height=15)

#' 
