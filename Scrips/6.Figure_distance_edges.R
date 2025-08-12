rm(list=ls())

# Packages and function
library(tidyverse)
library(plyr)
library(dplyr)
library(tidyr)
library(ggtext)
library(ggbreak)
library(MASS)
library(mvtnorm)
library(cowplot)
library(ggeffects)
library(ggtext)
library(gridExtra)
library(RColorBrewer)
library(LSAfun)
library(arm)

theme_plot<-theme(axis.text = element_text(size=14), axis.title = element_text(size=14, face="bold"), legend.text = element_text(size=12), strip.text = element_text(size=14), plot.title = element_text(size=14, face="bold"), panel.grid=element_line(colour="white"), panel.background = element_rect(fill="white") , axis.line = element_line(linewidth = 0.5, linetype = "solid",colour = "black"), strip.background = element_rect(fill="white"))

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


## Estimating distance to the edge

### Importing to prevent problems with reading information

struct_mat_REP_final<-read.csv("./Analyses/structural_REP.csv")
struct_mat_REP_final<-struct_mat_REP_final[,-1]

struct_mat_w0<-read.csv("./Analyses/structural_REP_w0.csv")
struct_mat_w0<-struct_mat_w0[,-1]


#Since we only take the first quadrant, we will put all negative slopes for  to a very small number

struct_mat_REP_final$a21_a11[which(struct_mat_REP_final$a21_a11<0)]<-0
struct_mat_REP_final$a22_a12[which(struct_mat_REP_final$a22_a12<0)]<-90

struct_mat_REP_final$a21_a11_upper[which(struct_mat_REP_final$a21_a11_upper<0)]<-0
struct_mat_REP_final$a22_a12_upper[which(struct_mat_REP_final$a22_a12_upper<0)]<-90

struct_mat_w0$a21_a11[which(struct_mat_w0$a21_a11<0)]<-0
struct_mat_w0$a22_a12[which(struct_mat_w0$a22_a12<0)]<-90

struct_mat_w0$a21_a11_upper[which(struct_mat_w0$a21_a11_upper<0)]<-0
struct_mat_w0$a22_a12_upper[which(struct_mat_w0$a22_a12_upper<0)]<-90

struct_mat_w0$a21_a11_lower[which(struct_mat_w0$a21_a11_lower<0)]<-0
struct_mat_w0$a22_a12_lower[which(struct_mat_w0$a22_a12_lower<0)]<-90


### Testing difference to one of the edges of the cone

##### all replicates

# calculating y for the x corresponding to the lambda Tu, in the vector slope
struct_mat_REP_final$Tu_lambda[1]*struct_mat_REP_final$a21_a11[1]

# just to be sure, it is equivalent to use the Te or Tu lambda
acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[1], struct_mat_REP_final$Te_lambda[1]))
         ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[1],struct_mat_REP_final$Tu_lambda[1]*struct_mat_REP_final$a21_a11[1]))))

acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[1], struct_mat_REP_final$Te_lambda[1]))
         ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[1]/struct_mat_REP_final$a21_a11[1], struct_mat_REP_final$Te_lambda[1]))))

#And if I use a random value of 10
acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[1], struct_mat_REP_final$Te_lambda[1])),
         LSAfun::normalize(c(10,10*struct_mat_REP_final$a21_a11[1]))))
# also ok, so I can just use the lambda's to do this
# CI is done of the vector of the growth rate to the lower or upper edges

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
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11_upper[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11_lower[x]))))
  }
  a
})

struct_mat_REP_final$distanceTu_upper<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  if(struct_mat_REP_final$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11_lower[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x],struct_mat_REP_final$Tu_lambda[x]*struct_mat_REP_final$a21_a11_upper[x]))))
  }
  a
})

struct_mat_REP_final$distanceTe_lower<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  if(struct_mat_REP_final$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12_upper[x],struct_mat_REP_final$Te_lambda[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12_lower[x],struct_mat_REP_final$Te_lambda[x]))))
  }
  a
})

struct_mat_REP_final$distanceTe_upper<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  if(struct_mat_REP_final$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12_lower[x],struct_mat_REP_final$Te_lambda[x]))))*-1
  } else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_REP_final$Tu_lambda[x], struct_mat_REP_final$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_REP_final$Te_lambda[x]/struct_mat_REP_final$a22_a12_upper[x],struct_mat_REP_final$Te_lambda[x]))))
  }
  a
})


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
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_upper[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_lower[x]))))
  }
  a
})


struct_mat_w0$distanceTu_upper<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_U[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_lower[x]))))*-1
  }else{
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x],struct_mat_w0$Tu_lambda[x]*struct_mat_w0$a21_a11_upper[x]))))
  }
  a
} )

struct_mat_w0$distanceTe_lower<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$Feasibility_L[x]==0){
    a<-acos(Dot(LSAfun::normalize(c(struct_mat_w0$Tu_lambda[x], struct_mat_w0$Te_lambda[x]))
                ,LSAfun::normalize(c(struct_mat_w0$Te_lambda[x]/struct_mat_w0$a22_a12_upper[x],struct_mat_w0$Te_lambda[x]))))*-1
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


#ordering

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

### minimum distance

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
  theme_plot+
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
  theme_plot+
  scale_x_discrete(labels=c("Te no cadmium\nTu no cadmium","Te cadmium\nTu no cadmium","Te no cadmium\nTu cadmium","Te cadmium\nTu cadmium"))+
  ylab("Distance to \n the T. urticae edge")+
  scale_fill_brewer(palette = "Spectral", labels=c("Te no cadmium:Tu no cadmium", "Te cadmium:Tu no cadmium", "Te no cadmium:Tu cadmium", "Te cadmium:Tu cadmium"), name="")+
  guides(fill=guide_legend(nrow=2))+
  xlab("Te Selection Regime : Tu Selection Regime")+
  theme(legend.position = "bottom", strip.background =element_rect(colour="white"), panel.border = element_rect(colour="black"), axis.text.x = element_text(size=10, angle = 45, vjust = 0.63), strip.text = element_text(size=10), axis.text.y=element_text(size=10), axis.title = element_text(size=10, face="bold"))

ggplot_distTu

plot_grid(ggplot_distTe + theme(legend.position="none"),ggplot_distTu + theme(legend.position="none"), ncol=2, labels=c("A", "B") )

save_plot("./Plots/FigS9.pdf", width=25, height=15)
save_plot("./Plots/FigS9.png", width=25, height=15)


ggplot(struct_mat_w0_final, aes(x=interaction(Te_Regime, Tu_Regime), y=minDistance, fill=interaction(Te_Regime, Tu_Regime)))+
  facet_grid(.~Environment, labeller=labeller(Environment=Env))+
  geom_hline(yintercept = 0, colour="lightgrey", linetype="dashed")+
  geom_errorbar(data=struct_mat_w0_final, aes(ymin=minDistance_L, ymax=minDistance_U), alpha=0.35, position=position_dodge2(0.5), width=0.5)+
  geom_point(size=2, shape=24,position=position_dodge2(0.5), alpha=0.65)+
  geom_errorbar(data=struct_mat_REP_final2, aes(ymin=minDistance_L, ymax=minDistance_U), width=0.25)+
  geom_point(data=struct_mat_REP_final2,aes(x=interaction(Te_Regime, Tu_Regime), y=minDistance, fill=interaction(Te_Regime, Tu_Regime)), colour="black", size=2.5, shape=21 )+
  theme_bw()+
  theme_plot+
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
