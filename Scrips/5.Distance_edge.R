rm(list=ls())

library(plyr)
library(tidyverse)
library(car)
library(fitdistrplus)
library(tidyr)
library(ggtext)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(MASS)
library(mvtnorm)
library(effects)
library(LSAfun)
library(arm)

## Estimating distance to the edge
struct_mat_REP_final<-read.csv("./Analyses/structural_REP.csv")
struct_mat_REP_final<-struct_mat_REP_final[,-1]

struct_mat_w0<-read.csv("./Analyses/structural_REP_w0.csv")
struct_mat_w0<-struct_mat_w0[,-1]


#Since we only consider the first quadrant in our analyses, we set a minimum and maximum slope

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

#CI is estimated using the vector of the growth rate to the lower or upper edges

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

##### per replicate

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

# Write for the figure
#write.csv(struct_mat_REP_final2,"./min_distance_pooled.csv")
#write.csv(struct_mat_w0_final,"./min_distance_per_replicate.csv")


#### checking which is the shortest distance

struct_mat_REP_final$minDistance<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  min(c(abs(struct_mat_REP_final$distanceTu[x]), abs(struct_mat_REP_final$distanceTe[x])))
  
})

struct_mat_w0$minDistance<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  min(c(abs(struct_mat_w0$distanceTu[x]), abs(struct_mat_w0$distanceTe[x])))
  
})
x<-2
struct_mat_REP_final$minDistance_L<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  
  if(struct_mat_REP_final$minDistance[x]==struct_mat_REP_final$distanceTu[x] | struct_mat_REP_final$minDistance[x]== -struct_mat_REP_final$distanceTu[x]){
    res<-struct_mat_REP_final$distanceTu_lower[x]
  }else if(struct_mat_REP_final$minDistance[x]==struct_mat_REP_final$distanceTe[x] | struct_mat_REP_final$minDistance[x]== -struct_mat_REP_final$distanceTe[x]){
    res<- struct_mat_REP_final$distanceTe_lower[x]
  }else
    res<-NA
  res
})


struct_mat_REP_final$minDistance_U<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){

  if(struct_mat_REP_final$minDistance[x]==struct_mat_REP_final$distanceTu[x] | struct_mat_REP_final$minDistance[x]== -struct_mat_REP_final$distanceTu[x]){
    res<-struct_mat_REP_final$distanceTu_upper[x]
  }else if(struct_mat_REP_final$minDistance[x]==struct_mat_REP_final$distanceTe[x] | struct_mat_REP_final$minDistance[x]== -struct_mat_REP_final$distanceTe[x]){
    res<- struct_mat_REP_final$distanceTe_upper[x]
  }else
    res<-NA
  res

})


struct_mat_w0$minDistance_L<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$minDistance[x]== struct_mat_w0$distanceTu[x] |  - struct_mat_w0$minDistance[x]== struct_mat_w0$distanceTu[x]){
    res<- struct_mat_w0$distanceTu_lower[x]
  }else if(struct_mat_w0$minDistance[x]== struct_mat_w0$distanceTe[x] |  - struct_mat_w0$minDistance[x]== struct_mat_w0$distanceTe[x]){
    res<- struct_mat_w0$distanceTe_lower[x]
  }else
    res<-NA
  
  res
})


struct_mat_w0$minDistance_U<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  if(struct_mat_w0$minDistance[x]== struct_mat_w0$distanceTu[x] |  - struct_mat_w0$minDistance[x]== struct_mat_w0$distanceTu[x]){
    res<- struct_mat_w0$distanceTu_upper[x]
  }else if(struct_mat_w0$minDistance[x]== struct_mat_w0$distanceTe[x] |  - struct_mat_w0$minDistance[x]== struct_mat_w0$distanceTe[x]){
    res<- struct_mat_w0$distanceTe_upper[x]
  }else
    res<-NA
  
  res
})



struct_mat_w0$minDistance2_lower<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  a<-struct_mat_w0$minDistance_L[x]
  
  if(struct_mat_w0$Feasibility[x]==0)
    a2<-a*-1
  else
    a2<-a
})

struct_mat_REP_final$minDistance2_lower<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  a<-struct_mat_REP_final$minDistance_L[x]
  if(struct_mat_REP_final$Feasibility[x]==0)
    a2<-a*-1
  else
    a2<-a
})

struct_mat_w0$minDistance2_upper<-sapply(c(1:length(struct_mat_w0$ND)), function(x){
  a<-struct_mat_w0$minDistance_U[x]
  
  if(struct_mat_w0$Feasibility[x]==0)
    a2<-a*-1
  else
    a2<-a
})

struct_mat_REP_final$minDistance2_upper<-sapply(c(1:length(struct_mat_REP_final$ND)), function(x){
  a<-struct_mat_REP_final$minDistance_U[x]
  if(struct_mat_REP_final$Feasibility[x]==0)
    a2<-a*-1
  else
    a2<-a
})

##### Stats

# descdist(subset(struct_mat_w0, Environment=="Cd")$distanceTu, discrete = FALSE, boot = 500)
# descdist(subset(struct_mat_w0, Environment=="Cd")$distanceTe , discrete = FALSE, boot = 500)
# 
# descdist(subset(struct_mat_w0, Environment=="N")$distanceTu , discrete = FALSE, boot = 500)
# descdist(subset(struct_mat_w0, Environment=="N")$distanceTe , discrete = FALSE, boot = 500)


dist_Cd_Tu<-glmmTMB(abs(distanceTu)~Te_Regime*Tu_Regime, data=subset(struct_mat_w0, Environment=="Cd"), family=Gamma(link="log"))
dist_Cd_Tu2<-glmmTMB(abs(distanceTu)~Te_Regime*Tu_Regime, data=subset(struct_mat_w0, Environment=="Cd"), family=Gamma(link="identity"))
dist_Cd_Tu3<-glmmTMB(abs(distanceTu)~Te_Regime*Tu_Regime, data=subset(struct_mat_w0, Environment=="Cd"), family=gaussian(link="log"))

dist_Cd_Te<-glmmTMB(abs(distanceTe)~Te_Regime*Tu_Regime, data=subset(struct_mat_w0, Environment=="Cd"), family=Gamma(link="log"))
dist_Cd_Te2<-glmmTMB(abs(distanceTe)~Te_Regime*Tu_Regime, data=subset(struct_mat_w0, Environment=="Cd"), family=Gamma(link="identity"))
dist_Cd_Te3<-glmmTMB(abs(distanceTe)~Te_Regime*Tu_Regime, data=subset(struct_mat_w0, Environment=="Cd"), family=gaussian(link="log"))


anova(dist_Cd_Tu, dist_Cd_Tu2, dist_Cd_Tu3)
anova(dist_Cd_Te, dist_Cd_Te2, dist_Cd_Te3)

summary(dist_Cd_Tu)
summary(dist_Cd_Te)

dist_N<-glmmTMB(minDistance~Te_Regime*Tu_Regime, data=subset(struct_mat_w0, Environment=="N"), family=Gamma(link="log"))
dist_Cd<-glmmTMB(minDistance~Te_Regime*Tu_Regime, data=subset(struct_mat_w0, Environment=="Cd"), family=Gamma(link="log"))
summary(dist_Cd)

summary(dist_Cd)
summary(dist_N)


dist_glm<-glmmTMB(minDistance~Te_Regime*Tu_Regime*Environment, data=struct_mat_w0, family=Gamma(link="log"))
summary(dist_glm)

dist_glm2<-glmmTMB(minDistance~Te_Regime+Tu_Regime+Environment, data=struct_mat_w0, family=Gamma(link="log"))
dist_glm3<-glmmTMB(minDistance~Te_Regime*Environment+Tu_Regime*Environment, data=struct_mat_w0, family=Gamma(link="log"))
summary(dist_glm3)
Anova(dist_glm2)

anova(dist_glm, dist_glm2, dist_glm3)

dist_Cd_Te<-glmmTMB(distanceTe~Environment, data=struct_mat_w0, family=Gamma(link="log"))
summary(dist_Cd_Te)

dist_Cd_Tu<-glmmTMB(distanceTu~Environment, data=struct_mat_w0, family=Gamma(link="log"))
summary(dist_Cd_Tu)


dist_Cd<-glmmTMB(minDistance~Environment, data=struct_mat_w0, family=Gamma(link="log"))
summary(dist_Cd)

