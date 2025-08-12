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

summary(m4)

summary(m3)
summary(m4)

emtrends(m3, var="pred_T1", type="response")
emtrends(m4, var="pred_T1", type="response")

#This forces the line to pass by the 0,0
m5<-glm(cbind(meanTe, meanTu)~0+pred_T1, data=sum_observed_coex_rep, family="binomial")
summary(m5)
emtrends(m5, var="pred_T1", type="response")
sum_observed_coex_rep<-sum_observed_coex_rep %>% mutate(yhat=predict(m5))

#### Testing per replicate

sum_observed_coex_ALL2<-sum2_observed_coex_ALL %>%
  group_by(SRTe, SRTu, Env) %>%
  summarize(meanRatio=mean(meanTeRatio, na.rm=TRUE), mean_pred=mean(predTeRatio, na.rm=TRUE), sdRatio=sd(meanTeRatio, na.rm=TRUE)/sqrt(5), sdPred=sd(predTeRatio, na.rm=TRUE)/sqrt(5), meanTe=mean(sum_Te, na.rm=TRUE), meanTu=mean(sum_Tu, na.rm=TRUE))

m4<- glmmTMB(cbind(meanTe, meanTu)~mean_pred*Env, data=sum_observed_coex_ALL2, family=binomial(link="logit"))

summary(m4)


# This is from the summary of the model so we can put the correct slope and CI in the plot, and not ggplot's model
slope_all<-invlogit(1.4211)
slope_ci_L<-invlogit(1.4211)-0.1861
slope_ci_U<-invlogit(1.4211)+0.1861

slopes_data<-data.frame(x=seq(0,1, 0.05))
slopes_data$ymin<-slope_ci_L*slopes_data$x
slopes_data$ymax<-slope_ci_U*slopes_data$x


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
  theme_plot+
  xlab("Observed ratio")+
  ylab("Predicted ratio")+
  ylim(c(0,1.2))+
  xlim(c(0,1.2))+
  theme(legend.position = "none")+
  coord_cartesian(xlim=c(0.3,0.96), ylim=c(0.3, 0.96))
save_plot("./Plots/Fig4.pdf", width=10, height=10)
save_plot("./Plots/Fig4.png", width=15, height=12)


#write.csv(pred_coex_RK_w0, "./PredictedPerReplicate.csv")
#write.csv(pred_coex_RK_REP, "./PredictedPooledReplicate.csv")