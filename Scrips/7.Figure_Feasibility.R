rm(list=ls())

library(plyr)
library(tidyverse)
library(tidyr)
library(ggtext)
library(ggbreak)
library(MASS)
library(cowplot)
library(ggeffects)
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

# Figure 2

test_struct<-read.csv("./Analyses/test_structural.csv")
test_struct<-test_struct[,-1]



str(test_struct)


#Code to create the confidence interval, not needed if you imported the CI_rep
# values_to_plot<-c(0,1, 2,3)
# 
# CI_rep<-as.data.frame(expand.grid(values=values_to_plot, Te_Regime=c("SR4", "SR5"), Tu_Regime=c("SR1", "SR2"), Environment=c("N", "Cd")))
# 
# CI_rep$x<-mapvalues(CI_rep$values, c(0,1,2,3), c(0,8,8,0))
# CI_rep$x_lower<-mapvalues(CI_rep$values, c(0,1,2,3), c(0,8,8,0))
# CI_rep$x_upper<-mapvalues(CI_rep$values, c(0,1,2,3), c(0,8,8,0))
# 
# CI_rep$y<-sapply(c(1:dim(CI_rep)[1]), function(x){
#   a<-subset(test_struct, Environment==CI_rep$Environment[x] &
#               Te_Regime==CI_rep$Te_Regime[x] & Tu_Regime==CI_rep$Tu_Regime[x])
#   
#   if(CI_rep$values[x]==0 | CI_rep$values[x]==3){
#     b<-0
#   }else{
#     if(CI_rep$values[x]==1){ 
#       b<-a$a22_a12[1]*CI_rep$x[x]
#     }else if(CI_rep$values[x]==2){
#       b<-a$a21_a11[1]*CI_rep$x[x]
#     }
#     
#   }
#   b
# })
# 
# 
# CI_rep$y_lower<-sapply(c(1:dim(CI_rep)[1]), function(x){
#   a<-subset(test_struct, Environment==CI_rep$Environment[x] &
#               Te_Regime==CI_rep$Te_Regime[x] & Tu_Regime==CI_rep$Tu_Regime[x])
#   
#   if(CI_rep$values[x]==0 | CI_rep$values[x]==3){
#     b<-0
#   }else{
#     if(a$Environment[1]=="N"){
#       if(a$Tu_Regime[1]=="SR2" & a$Te_Regime[1]=="SR5"){
#         if(CI_rep$values[x]==1){ 
#           b<-a$min_a22_a12[1]*CI_rep$x_lower[x]
#         }else if(CI_rep$values[x]==2){
#           b<-a$min_a21_a11[1]*CI_rep$x_lower[x]
#         }
#         
#       }else{
#         if(CI_rep$values[x]==1){ 
#           b<-a$max_a22_a12[1]*CI_rep$x_lower[x]
#         }else if(CI_rep$values[x]==2){
#           b<-a$min_a21_a11[1]*CI_rep$x_lower[x]
#         }}}
#     else if(a$Environment[1]=="Cd"){
#       if( a$Te_Regime[1]=="SR4"){
#         if(CI_rep$values[x]==1){ 
#           b<-a$min_a22_a12[1]*CI_rep$x_lower[x]
#         }else if(CI_rep$values[x]==2){
#           b<-a$max_a21_a11[1]*CI_rep$x_lower[x]
#         }
#         
#       }else{
#         if(CI_rep$values[x]==1){ 
#           b<-a$max_a22_a12[1]*CI_rep$x_lower[x]
#         }else if(CI_rep$values[x]==2){
#           b<-a$max_a21_a11[1]*CI_rep$x_lower[x]
#         }}
#     }}
#   b
# })
# 
# CI_rep$y_upper<-sapply(c(1:dim(CI_rep)[1]), function(x){
#   a<-subset(test_struct, Environment==CI_rep$Environment[x] &
#               Te_Regime==CI_rep$Te_Regime[x] & Tu_Regime==CI_rep$Tu_Regime[x])
#   
#   if(CI_rep$values[x]==0 | CI_rep$values[x]==3){
#     b<-0
#   }else{
#     if(a$Environment[1]=="N"){
#       if(a$Tu_Regime[1]=="SR2" & a$Te_Regime[1]=="SR5"){
#         if(CI_rep$values[x]==1){ 
#           b<-a$max_a22_a12[1]*CI_rep$x_upper[x]
#         }else if(CI_rep$values[x]==2){
#           b<-a$max_a21_a11[1]*CI_rep$x_upper[x]
#         }
#         
#       }else{
#         if(CI_rep$values[x]==1){ 
#           b<-a$min_a22_a12[1]*CI_rep$x_upper[x]
#         }else if(CI_rep$values[x]==2){
#           b<-a$max_a21_a11[1]*CI_rep$x_upper[x]
#         }}}
#     else if(a$Environment[1]=="Cd"){
#       if( a$Te_Regime[1]=="SR4"){
#         if(CI_rep$values[x]==1){ 
#           b<-a$max_a22_a12[1]*CI_rep$x_upper[x]
#         }else if(CI_rep$values[x]==2){
#           b<-a$min_a21_a11[1]*CI_rep$x_upper[x]
#         }
#         
#       }else{
#         if(CI_rep$values[x]==1){ 
#           b<-a$min_a22_a12[1]*CI_rep$x_upper[x]
#         }else if(CI_rep$values[x]==2){
#           b<-a$min_a21_a11[1]*CI_rep$x_upper[x]
#         }}
#     }}
#   b
# })


CI_rep<-read.csv("./Analyses/CI_rep.csv")
CI_rep<-CI_rep[,-1]

str(CI_rep)

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
  theme_plot+
  ylab(c("Intrinsic growth rate T. evansi"))+
  xlab(c("Intrinsic growth rate T. urticae"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  #ylim(c(0,8))+
  #xlim(c(0,5))+
  coord_cartesian(xlim =c(0.1,2.8), ylim=c(0.1,6), expand = TRUE)+
  ggtitle("No cadmium environment")
normal_feas
save_plot("./Plots/Fig2B_B.pdf", width=20, height=15)

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
  theme_plot+
  ylab(c("Intrinsic growth rate T. evansi"))+
  xlab(c("Intrinsic growth rate T. urticae"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  #ylim(c(0,8))+
  #xlim(c(0,5))+
  coord_cartesian(xlim =c(0.1,2), ylim=c(0.1,2.4), expand = TRUE)+
  ggtitle("Cadmium environment")
cd_feas
save_plot("./Plots/Fig2A_B.pdf", width=20, height=15)


plot_grid(cd_feas, normal_feas, ncol=2, labels=c("A", "B") )
save_plot("./Plots/Fig2_all.pdf", width=40, height=20)
save_plot("./Plots/Fig2_all.tiff", width=40, height=20)

