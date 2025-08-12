rm(list=ls())
library(plyr)
library(tidyverse)
library(fitdistrplus)
library(tidyr)
library(ggtext)
library(ggbreak)
library(effects)
library(cowplot)
library(ggeffects)
library(ggtext)
library(gridExtra)
library(RColorBrewer)

theme_plots<-theme(axis.text = element_text(size=14), axis.title = element_text(size=14, face="bold"), legend.text = element_text(size=12), strip.text = element_text(size=14), plot.title = element_text(size=14, face="bold"), panel.grid=element_line(colour="white"), panel.background = element_rect(fill="white") , axis.line = element_line(linewidth = 0.5, linetype = "solid",colour = "black"), strip.background = element_rect(fill="white"))

save_plot<-function(dir, width=15, height=10, ...){
  ggsave(dir, width = width, height = height, units = c("cm"))
}

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


save_plot("./Plots/FigS7.pdf", width=30, height=15)
save_plot("./Plots/FigS7.png", width=30, height=15)

plot_grid(Te_intra_ev, Tu_intra_ev, labels=c("A", "B"))

save_plot("./Plots/FigS3.pdf", width=30, height=15)
save_plot("./Plots/FigS3.png", width=30, height=15)

plot_grid(Te_inter_ev, Tu_inter_ev, labels=c("A", "B"))

save_plot("./Plots/FigS4.pdf", width=30, height=15)
save_plot("./Plots/FigS4.png", width=30, height=15)


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

save_plot("./Plots/FigS8.pdf", width=30, height=15)
save_plot("./Plots/FigS8.png", width=30, height=15)

plot_grid(Te_intra_ev_N, Tu_intra_ev_N, labels=c("A", "B"))

save_plot("./Plots/FigS5.pdf", width=30, height=15)
save_plot("./Plots/FigS5.png", width=30, height=15)

plot_grid(Te_inter_ev_N, Tu_inter_ev_N, labels=c("A", "B"))

save_plot("./Plots/FigS6.pdf", width=30, height=15)
save_plot("./Plots/FigS6.png", width=30, height=15)



