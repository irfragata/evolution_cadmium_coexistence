rm(list=ls())
library(plyr)
library(tidyverse)
library(fitdistrplus)
library(dplyr)
library(tidyr)
library(ggtext)
library(ggbreak)
library(effects)
library(MASS)
library(mvtnorm)
library(cowplot)
library(ggeffects)
library(ggtext)
library(gridExtra)
library(RColorBrewer)
library(LSAfun)
library(arm)
library(car)
library(fitdistrplus)
library(marginaleffects)

theme_plots<-theme(axis.text = element_text(size=14), axis.title = element_text(size=14, face="bold"), legend.text = element_text(size=12), strip.text = element_text(size=14), plot.title = element_text(size=14, face="bold"), panel.grid=element_line(colour="white"), panel.background = element_rect(fill="white") , axis.line = element_line(linewidth = 0.5, linetype = "solid",colour = "black"), strip.background = element_rect(fill="white"))

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


# 1 - Figure parameters

###### Importing data files

#To use the data sets already available in the repository, we can simply read the csv.

## Importing
param_all_REP<-read.csv("./Analyses/cxr_normal_REP/parameters_cxr_normal_REP.csv")
param_all_REP_upper<-read.csv("./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_upper.csv")
param_all_REP_lower<-read.csv( "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_lower.csv")
param_all_REP<-param_all_REP[,-1]
param_all_REP_upper<-param_all_REP_upper[,-1]
param_all_REP_lower<-param_all_REP_lower[,-1]


#Data wrangling to make it easier to plot the figures

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


###### Importing data files for replicate estimation

## Importing
param_all_w0<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal.csv")
param_all_w0_upper<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal_upper.csv")
param_all_w0_lower<-read.csv( "./Analyses/cxr_normal/parameters_cxr_normal_lower.csv")

param_all_w0<-param_all_w0[,-1]
param_all_w0_upper<-param_all_w0_upper[,-1]
param_all_w0_lower<-param_all_w0_lower[,-1]


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



colors_comb<-brewer.pal(name = "Spectral", 4)

#"#D7191C" "#FDAE61" "#ABDDA4" "#2B83BA"

### Cadmium env

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
  xlab("T. urticae selection regime")+theme(legend.position = "None")


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
  theme(legend.position = "None", legend.text = element_text(size=10), legend.background = element_rect(fill=NA), legend.key.size = unit(0.2, 'cm'))

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
  xlab("T. evansi selection regime")+theme(legend.position = "None")


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
  xlab("T. evansi selection regime")

plot_grid(Te_gr_ev, Tu_gr_ev, labels=c("A", "B"))


save_plot("./Plots2/FigS7.pdf", width=30, height=15)
save_plot("./Plots2/FigS7.png", width=30, height=15)

plot_grid(Te_intra_ev, Tu_intra_ev, labels=c("A", "B"))

save_plot("./Plots2/FigS3.pdf", width=30, height=15)
save_plot("./Plots2/FigS3.png", width=30, height=15)

plot_grid(Te_inter_ev, Tu_inter_ev, labels=c("A", "B"))

save_plot("./Plots2/FigS4.pdf", width=30, height=15)
save_plot("./Plots2/FigS4.png", width=30, height=15)


#### Ancestral env

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
  xlab("T. urticae selection regime")+theme(legend.position = "None")

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
  xlab("T. urticae selection regime")+theme(legend.position = "None")


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
  theme(legend.position = "None", legend.text = element_text(size=10), legend.background = element_rect(fill=NA), legend.key.size = unit(0.2, 'cm'))

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
  xlab("T. evansi selection regime")+ theme(legend.position = "None")

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
  xlab("T. evansi selection regime")+theme(legend.position = "None")


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

plot_grid(Te_gr_ev_N, Tu_gr_ev_N, labels=c("A", "B"))

save_plot("./Plots2/FigS8.pdf", width=30, height=15)
save_plot("./Plots2/FigS8.png", width=30, height=15)

plot_grid(Te_intra_ev_N, Tu_intra_ev_N, labels=c("A", "B"))

save_plot("./Plots2/FigS5.pdf", width=30, height=15)
save_plot("./Plots2/FigS5.png", width=30, height=15)

plot_grid(Te_inter_ev_N, Tu_inter_ev_N, labels=c("A", "B"))

save_plot("./Plots2/FigS6.pdf", width=30, height=15)
save_plot("./Plots2/FigS6.png", width=30, height=15)


# 2 - Figure dstances

## Estimating distance to the edge

struct_mat_REP_final<-read.csv("./Analyses/structural_REP.csv")
struct_mat_REP_final<-struct_mat_REP_final[,-1]

struct_mat_w0<-read.csv("./Analyses/structural_REP_w0.csv")
struct_mat_w0<-struct_mat_w0[,-1]


# CI is estimated using as reference the vector of the growth rate and then estimating the difference to the lower or upper edges. Negative distances indicate unfeasible systems

struct_mat_REP_final$distanceTu<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  if(struct_mat_REP_final$Feasibility[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11[x]))))
  }
  a
})

struct_mat_REP_final$distanceTe<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){ 
  if(struct_mat_REP_final$Feasibility[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12[x],struct_mat_REP_final$Te_lambda[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12[x],struct_mat_REP_final$Te_lambda[x]))))
  }
  a
})

#If the growth rate vector is within the feasibility domain, it will be closer to the smaller feasibility domain, if the growth rate vector is outside of the feasibility domain, the smallest distance will be to the wider feasibility 

struct_mat_REP_final$distanceTu_lower<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  if(struct_mat_REP_final$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11_lower[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11_lower[x]))))
  }
  a
})

struct_mat_REP_final$distanceTu_upper<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  if(struct_mat_REP_final$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11_upper[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11_upper[x]))))
  }
  a
})

struct_mat_REP_final$distanceTe_lower<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  if(struct_mat_REP_final$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12_lower[x],struct_mat_REP_final$Te_lambda[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12_lower[x],struct_mat_REP_final$Te_lambda[x]))))
  }
  a
})

struct_mat_REP_final$distanceTe_upper<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  if(struct_mat_REP_final$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12_upper[x],struct_mat_REP_final$Te_lambda[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12_upper[x],struct_mat_REP_final$Te_lambda[x]))))
  }
  a
})

# 
struct_mat_REP_final$distanceTu_lower
struct_mat_REP_final$distanceTu
struct_mat_REP_final$distanceTu_upper

struct_mat_REP_final$distanceTe_lower
struct_mat_REP_final$distanceTe
struct_mat_REP_final$distanceTe_upper

##### per replicate

# calculating y for the x corresponding to the lambda Tu, in the vector slope
struct_mat_w0$Tu_lambda[1]*struct_mat_w0$a21_a11[1]

# just to be sure, it is equivalent to use the Te or Tu lambda
acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[1], struct_mat_w0$Te_lambda[1]))
         ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[1],struct_mat_w0$Tu_lambda[1]*struct_mat_w0$a21_a11[1]))))

acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[1], struct_mat_w0$Te_lambda[1]))
         ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[1]/struct_mat_w0$a21_a11[1], struct_mat_w0$Te_lambda[1]))))

# per replicate
struct_mat_w0$distanceTu<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11[x]))))
  }
  a
})

struct_mat_w0$distanceTe<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12[x],struct_mat_w0$Te_lambda[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12[x],struct_mat_w0$Te_lambda[x]))))
  }
  a
})

struct_mat_w0$distanceTu_lower<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_lower[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_lower[x]))))
  }
  a
})


struct_mat_w0$distanceTu_upper<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_upper[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_upper[x]))))
  }
  a
} )

struct_mat_w0$distanceTe_lower<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12_lower[x],struct_mat_w0$Te_lambda[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12_lower[x],struct_mat_w0$Te_lambda[x]))))
  }
  
  a
} )

struct_mat_w0$distanceTe_upper<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12_upper[x],struct_mat_w0$Te_lambda[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12_upper[x],struct_mat_w0$Te_lambda[x]))))
  }
})


# Because of facilitation (specially when there is double facilitation estimated) the distance of the upper and lower estimates is switched. 
aux_Tu_orders<-as.data.frame(t(sapply(c(1:nrow(struct_mat_REP_final)), function(x){
  vect<-c(struct_mat_REP_final$distanceTu[x],struct_mat_REP_final$distanceTu_lower[x], struct_mat_REP_final$distanceTu_upper[x])
  ord<-order(vect,decreasing = FALSE)
  
  vect[ord]
})))

colnames(aux_Tu_orders)<-c("Tu_distance_lower","Tu_distance", "Tu_distance_upper" )

aux_Te_orders<-as.data.frame(t(sapply(c(1:nrow(struct_mat_REP_final)), function(x){
  vect<-c(struct_mat_REP_final$distanceTe[x],struct_mat_REP_final$distanceTe_lower[x], struct_mat_REP_final$distanceTe_upper[x])
  ord<-order(vect,decreasing = FALSE)
  
  vect[ord]
})))

colnames(aux_Te_orders)<-c("Te_distance_lower","Te_distance", "Te_distance_upper" )

aux_Tu_orders_w0<-as.data.frame(t(sapply(c(1:nrow(struct_mat_w0)), function(x){
  vect<-c(struct_mat_w0$distanceTu[x],struct_mat_w0$distanceTu_lower[x], struct_mat_w0$distanceTu_upper[x])
  ord<-order(vect,decreasing = FALSE)
  
  vect[ord]
})))

colnames(aux_Tu_orders_w0)<-c("Tu_distance_lower","Tu_distance", "Tu_distance_upper" )

aux_Te_orders_w0<-as.data.frame(t(sapply(c(1:nrow(struct_mat_w0)), function(x){
  vect<-c(struct_mat_w0$distanceTe[x],struct_mat_w0$distanceTe_lower[x], struct_mat_w0$distanceTe_upper[x])
  ord<-order(vect,decreasing = FALSE)
  
  vect[ord]
})))

colnames(aux_Te_orders_w0)<-c("Te_distance_lower","Te_distance", "Te_distance_upper" )

struct_mat_REP_final2<-as.data.frame(cbind(struct_mat_REP_final[,c(1,2,3)], aux_Tu_orders, aux_Te_orders))

struct_mat_w0_final<-as.data.frame(cbind(struct_mat_w0[,c(1,2,3,4)], aux_Tu_orders_w0, aux_Te_orders_w0))


### Estimating minimum distance between Te or Tu

aux_min_distance<-as.data.frame(t(sapply(c(1:nrow(struct_mat_REP_final2)), function(x){
  a<-min(c(abs(struct_mat_REP_final2$Tu_distance[x]),abs(struct_mat_REP_final2$Te_distance[x])))
  
  a_name<-colnames(struct_mat_REP_final2)[which(struct_mat_REP_final2[x,]==a | struct_mat_REP_final2[x,]==-a)]
  
  if(a_name=="Tu_distance"){
    vecF<-c(struct_mat_REP_final2$Tu_distance_lower[x], struct_mat_REP_final2$Tu_distance[x], struct_mat_REP_final2$Tu_distance_upper[x])
  }else if(a_name=="Te_distance"){
    vecF<- c(struct_mat_REP_final2$Te_distance_lower[x], struct_mat_REP_final2$Te_distance[x], struct_mat_REP_final2$Te_distance_upper[x])
  }
  
  vecF
})))

# This is to use the replicate values of the species that showed the lower distance (and not the minimal between species, otherwise it would create some weird patterns) 
colnames(aux_min_distance)<-c("minDistance_L", "minDistance", "minDistance_U")


aux_min_distance_w0<-as.data.frame(t(sapply(c(1:nrow(struct_mat_w0_final)), function(x){
  a<-min(c(abs(struct_mat_w0_final$Tu_distance[x]),abs(struct_mat_w0_final$Te_distance[x])))
  
  a_name<-colnames(struct_mat_w0_final)[which(struct_mat_w0_final[x,]==a | struct_mat_w0_final[x,]==-a)][1]
  
  if(a_name=="Tu_distance" |a_name=="Tu_distance_lower" | a_name=="Tu_distance_upper" ){
    vecF<-c(struct_mat_w0_final$Tu_distance_lower[x],struct_mat_w0_final$Tu_distance[x] , struct_mat_w0_final$Tu_distance_upper[x])
  }else if(a_name=="Te_distance" |a_name=="Te_distance_lower" | a_name=="Te_distance_upper"){
    vecF<- c(struct_mat_w0_final$Te_distance_lower[x], struct_mat_w0_final$Te_distance[x], struct_mat_w0_final$Te_distance_upper[x])
  }
  
  vecF
})))

str(struct_mat_w0_final)

colnames(aux_min_distance_w0)<-c("minDistance_L", "minDistance", "minDistance_U")

struct_mat_REP_final2<-as.data.frame(cbind(struct_mat_REP_final2, aux_min_distance))

struct_mat_w0_final<-as.data.frame(cbind(struct_mat_w0_final, aux_min_distance_w0))

##### Final figure 3 and S9

ggplot_distTe<-ggplot(struct_mat_REP_final2, aes(x=interaction(Tu_Regime, Te_Regime), fill=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Environment~., labeller=labeller(Environment=Env))+
  geom_hline(yintercept = 0, colour="lightgrey", linetype="dashed")+
  geom_errorbar(data=struct_mat_w0_final, aes(ymin=Te_distance_lower, ymax=Te_distance_upper), width=0.5, position=position_dodge2(0.5), alpha=0.25)+
  geom_point(data=struct_mat_w0_final, size=1.5, shape=24, position=position_dodge2(0.5), colour="black", alpha=0.55, aes(y=Te_distance))


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

ggplot_distTu<-ggplot(struct_mat_REP_final2, aes(x=interaction(Tu_Regime, Te_Regime), fill=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Environment~., labeller=labeller(Environment=Env))+
  geom_hline(yintercept = 0, colour="lightgrey", linetype="dashed")+
  geom_errorbar(data=struct_mat_w0_final, aes(ymin=Tu_distance_lower, ymax=Tu_distance_upper), width=0.5, position=position_dodge2(width=0.5), alpha=0.35)+
  geom_point(data=struct_mat_w0_final, size=1.5, shape=24, position=position_dodge2(width=0.5), colour="black", alpha=0.55, aes(y=Tu_distance))

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

plot_grid(ggplot_distTe + theme(legend.position="none"),ggplot_distTu + theme(legend.position="none"), ncol=2, labels=c("A", "B") )

save_plot("./Plots2/FigS9.pdf", width=25, height=15)
save_plot("./Plots2/FigS9.png", width=25, height=15)


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
save_plot("./Plots2/Fig3.pdf", width=15, height=10)
save_plot("./Plots2/Fig3.png", width=15, height=10)


# Figure feasibility domains
eval=FALSE # if you want to rerun things again you need to put to TRUE

test_struct<-read.csv("./Analyses/structural_REP.csv")# Same data frame as above, added in case the code for other figures above was not run
test_struct<-test_struct[,-1]

CI_rep<-read.csv("./Analyses/CI_rep.csv")
CI_rep<-CI_rep[,-1]

#test_struct<-struct_mat_REP_final
str(test_struct)
str(CI_rep)

ggplot(subset(test_struct, Environment=="N"), aes(x=Tu_lambda, y=Te_lambda, colour=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Tu_Regime~ Te_Regime, labeller=labeller(Tu_Regime=regimeTu, Te_Regime=regimeTe) )+
  #geom_polygon(data=subset(CI_rep,Environment=="N"),aes(x=x_upper, y=y_upper), fill="#fb6a4a", alpha=0.35, linewidth=0.85, colour=NA)+
  #geom_polygon(data=subset(CI_rep,Environment=="N"),aes(x=x, y=y), fill="#969696", alpha=0.75, linewidth=0.85, colour="black")+
  #geom_polygon(data=subset(CI_rep,Environment=="N"),aes(x=x_lower, y=y_lower), fill="lightskyblue1", alpha=0.9,  linewidth=0.85, colour=NA)+
  geom_point(colour="black")+
  geom_segment(data=subset(test_struct, Environment=="N" ), aes(xend=Tu_lambda, yend=Te_lambda,x=0, y=0),  arrow=arrow(length = unit(0.3, "cm")), colour="black", linewidth=1)+
  geom_segment(data=subset(test_struct, Environment=="N" ), aes(xend=Tu_lambda_lower, yend=Te_lambda_lower,x=0, y=0),  colour="black", linetype="dashed", linewidth=0.75)+
  geom_segment(data=subset(test_struct, Environment=="N" ), aes(xend=Tu_lambda_upper, yend=Te_lambda_upper,x=0, y=0), colour="black", linetype="dashed", linewidth=0.75)+
  geom_abline(aes(slope=a22_a12, intercept=0), colour="black", linewidth=2)+
  geom_abline(aes(slope=a21_a11, intercept=0), colour="darkgrey", linewidth=2)+
  geom_abline(aes(slope=a22_a12_lower, intercept=0), colour="darkorange", linewidth=2)+
  geom_abline(aes(slope=a21_a11_lower, intercept=0), colour="darkgreen", linewidth=2)+
  geom_abline(aes(slope=a22_a12_upper, intercept=0), colour="darkorange", linewidth=2, linetype="dotted")+
  geom_abline(aes(slope=a21_a11_upper, intercept=0), colour="darkgreen", linewidth=2, linetype="dotted")+
  #geom_abline(intercept=0, slope=360, linewidth=2)+
  theme_bw()+
  theme_plots+
  ylab(c("Intrinsic growth rate T. evansi"))+
  xlab(c("Intrinsic growth rate T. urticae"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  ylim(c(-100,100))+
  xlim(c(-100,100))+
  coord_cartesian(xlim =c(0.1,2.8), ylim=c(0.1,6), expand = TRUE)+
  ggtitle("No cadmium environment")


ggplot(subset(test_struct, Environment=="Cd"), aes(x=Tu_lambda, y=Te_lambda, colour=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Tu_Regime~ Te_Regime, labeller=labeller(Tu_Regime=regimeTu, Te_Regime=regimeTe) )+
  #geom_polygon(data=subset(CI_rep,Environment=="Cd"),aes(x=x_upper, y=y_upper), fill="#fb6a4a", alpha=0.35, linewidth=0.85, colour=NA)+
  #geom_polygon(data=subset(CI_rep,Environment=="Cd"),aes(x=x, y=y), fill="#969696", alpha=0.75, linewidth=0.85, colour="black")+
  #geom_polygon(data=subset(CI_rep,Environment=="Cd"),aes(x=x_lower, y=y_lower), fill="lightskyblue1", alpha=0.9,  linewidth=0.85, colour=NA)+
  geom_point(colour="black")+
  geom_segment(data=subset(test_struct, Environment=="Cd" ), aes(xend=Tu_lambda, yend=Te_lambda,x=0, y=0),  arrow=arrow(length = unit(0.3, "cm")), colour="black", linewidth=1)+
  geom_segment(data=subset(test_struct, Environment=="Cd" ), aes(xend=Tu_lambda_lower, yend=Te_lambda_lower,x=0, y=0),  colour="black", linetype="dashed", linewidth=0.75)+
  geom_segment(data=subset(test_struct, Environment=="Cd" ), aes(xend=Tu_lambda_upper, yend=Te_lambda_upper,x=0, y=0), colour="black", linetype="dashed", linewidth=0.75)+
  geom_abline(aes(slope=a22_a12, intercept=0), colour="black", linewidth=2)+
  geom_abline(aes(slope=a21_a11, intercept=0), colour="darkgrey", linewidth=2)+
  geom_abline(aes(slope=a22_a12_lower, intercept=0), colour="darkorange", linewidth=2)+
  geom_abline(aes(slope=a21_a11_lower, intercept=0), colour="darkgreen", linewidth=2)+
  geom_abline(aes(slope=a22_a12_upper, intercept=0), colour="darkorange", linewidth=2, linetype="dotted")+
  geom_abline(aes(slope=a21_a11_upper, intercept=0), colour="darkgreen", linewidth=2, linetype="dotted")+
  #geom_abline(intercept=0, slope=360, linewidth=2)+
  theme_bw()+
  theme_plots+
  ylab(c("Intrinsic growth rate T. evansi"))+
  xlab(c("Intrinsic growth rate T. urticae"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  ylim(c(-100,100))+
  xlim(c(-100,100))+
  coord_cartesian(xlim =c(0.1,2.8), ylim=c(0.1,6), expand = TRUE)+
  ggtitle("Cadmium environment")


##### Final figure

normal_feas<-ggplot(subset(test_struct, Environment=="N"), aes(x=Tu_lambda, y=Te_lambda, colour=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Tu_Regime~ Te_Regime, labeller=labeller(Tu_Regime=regimeTu, Te_Regime=regimeTe) )+
  geom_polygon(data=subset(CI_rep,Environment=="N"),aes(x=x_upper, y=y_upper), fill="#fb6a4a", alpha=0.35, linewidth=0.85, colour=NA)+
  geom_polygon(data=subset(CI_rep,Environment=="N"),aes(x=x, y=y), fill="#969696", alpha=0.75, linewidth=0.85, colour="black")+
  geom_polygon(data=subset(CI_rep,Environment=="N"),aes(x=x_lower, y=y_lower), fill="lightskyblue1", alpha=0.9,  linewidth=0.85, colour=NA)+
  geom_point(colour="black")+
  geom_segment(data=subset(test_struct, Environment=="N" ), aes(xend=Tu_lambda, yend=Te_lambda,x=0, y=0),  arrow=arrow(length = unit(0.3, "cm")), colour="black", linewidth=1)+
  geom_segment(data=subset(test_struct, Environment=="N" ), aes(xend=Tu_lambda_lower, yend=Te_lambda_lower,x=0, y=0),  colour="black", linetype="dashed", linewidth=0.75)+
  geom_segment(data=subset(test_struct, Environment=="N" ), aes(xend=Tu_lambda_upper, yend=Te_lambda_upper,x=0, y=0), colour="black", linetype="dashed", linewidth=0.75)+
  theme_bw()+
  theme_plots+
  ylab(c("Intrinsic growth rate T. evansi"))+
  xlab(c("Intrinsic growth rate T. urticae"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  ylim(c(-100,100))+
  xlim(c(-100,100))+
  coord_cartesian(xlim =c(0.1,2.8), ylim=c(0.1,6), expand = TRUE)+
  ggtitle("No cadmium environment")
normal_feas
#save_plot("./Plots/Fig2B_B.pdf", width=20, height=15)

cd_feas<-ggplot(subset(test_struct, Environment=="Cd"), aes(x=Tu_lambda, y=Te_lambda, colour=interaction(Tu_Regime, Te_Regime)))+
  facet_grid(Tu_Regime~ Te_Regime, labeller=labeller(Tu_Regime=regimeTu, Te_Regime=regimeTe) )+
  geom_polygon(data=subset(CI_rep,Environment=="Cd"),aes(x=x_upper, y=y_upper), fill="#fb6a4a", alpha=0.35,linewidth=0.85, colour=NA)+
  geom_polygon(data=subset(CI_rep,Environment=="Cd"),aes(x=x, y=y), fill="#969696", alpha=0.75, linewidth=0.85,colour="black")+
  geom_polygon(data=subset(CI_rep,Environment=="Cd"),aes(x=x_lower, y=y_lower), fill="lightskyblue1", alpha=0.9, linewidth=0.85, colour=NA)+
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
cd_feas
#save_plot("./Plots/Fig2A_B.pdf", width=20, height=15)


plot_grid(cd_feas, normal_feas, ncol=2, labels=c("A", "B") )
save_plot("./Plots/Fig2_all.pdf", width=40, height=20)
save_plot("./Plots/Fig2_all.tiff", width=40, height=20)


# Figure 4
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

#This forces the line to pass by 0,0
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
  theme_plots+
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


# Figure 1 - predictions

# Importing data if it was not imported before
if(!eval){
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
  theme_plots+
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
  theme_plots+
  scale_fill_brewer(palette = "Spectral", labels=c("Te no cadmium:Tu no cadmium", "Te cadmium:Tu no cadmium", "Te no cadmium:Tu cadmium", "Te cadmium:Tu cadmium"), name="")+
  xlab(expression(paste("Predicted offspring production for ", italic("T. urticae"))))+
  guides(fill=guide_legend(nrow=2))+
  scale_y_discrete(labels=c(expression(lambda+ alpha[ii] + alpha [ij]),expression(lambda+ alpha[ii]), expression(lambda)), limits=rev(levels(droplevels(subset(pred_coex1Gen_long, parameter=="predTu_onlyLambda" |  parameter=="predTu_Lambda_INTRA" |  parameter=="predTu_ALL"))$parameter3)))+
  theme(legend.position = "bottom", axis.text = element_text(size=12), axis.title = element_text(face="plain", size=12))+
  ylab("")
save_plot("./Plots/Fig1B.pdf", width=17.5, height=10)
