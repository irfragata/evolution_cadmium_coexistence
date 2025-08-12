rm(list=ls())

library(plyr)
library(tidyverse)
library(car)
library(tidyr)
library(MASS)
library(mvtnorm)
library(DescTools)
library(performance)
library(LSAfun)
library(arm)

## Importing pooled data
param_all_REP<-read.csv("./Analyses/cxr_normal_REP/parameters_cxr_normal_REP.csv")
param_all_REP_upper<-read.csv("./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_upper.csv")
param_all_REP_lower<-read.csv( "./Analyses/cxr_normal_REP/parameters_cxr_normal_REP_lower.csv")
param_all_REP<-param_all_REP[,-1]
param_all_REP_upper<-param_all_REP_upper[,-1]
param_all_REP_lower<-param_all_REP_lower[,-1]

## Importing per replicate
param_all_w0<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal.csv")
param_all_w0_upper<-read.csv("./Analyses/cxr_normal/parameters_cxr_normal_upper.csv")
param_all_w0_lower<-read.csv( "./Analyses/cxr_normal/parameters_cxr_normal_lower.csv")

param_all_w0<-param_all_w0[,-1]
param_all_w0_upper<-param_all_w0_upper[,-1]
param_all_w0_lower<-param_all_w0_lower[,-1]


# Predicting coexistence (structural stability)

### Structural stability

### Defining functions

#input parameters:
#alpha = competition strength matrix 
#r = vector of intrinsic growth rates

#structural niche difference (output on a log scale)
Omega <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha, tol = 1e-40)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}

#vector defining the centroid of the feasibility domain
r_centroid <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}


#structural fitness difference (in degree)
theta <- function(alpha,r){
  r_c <- r_centroid(alpha)
  out <- acos(sum(r_c*r, na.rm = TRUE)/(sqrt(sum(r^2, na.rm = TRUE))*sqrt(sum(r_c^2, na.rm = TRUE))))*90/pi
  return(out)
}


#test if a system (alpha and r) is feasible (output 1 = feasible, 0 = not feasible)
test_feasibility <- function(alpha,r){
  out <- prod(solve(alpha,r)>0)
  return(out)
}



### 3.1.1 - estimating structural stability for pooled replicates
#x<-2
struct_mat_REP<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_REP[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_REP$Te_intra[x],param_all_REP$Te_inter[x],param_all_REP$Tu_inter[x], param_all_REP$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_REP$Te_lambda[x],param_all_REP$Tu_lambda[x])
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))


colnames(struct_mat_REP)<-c("ND", "FD", "Feasibility")

#For the lower we use the higher alphas with lower lambda, and for upper the other way around
# Since we have facilitation we have to actually test what is the lowest value

struct_mat_REP_U<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_REP_lower[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_REP_upper$Te_intra[x],param_all_REP_upper$Te_inter[x], param_all_REP_upper$Tu_inter[x],param_all_REP_upper$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_REP$Te_lambda[x],param_all_REP$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))


struct_mat_REP_L<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_REP_upper[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_REP_lower$Te_intra[x],param_all_REP_lower$Te_inter[x],param_all_REP_lower$Tu_inter[x],param_all_REP_lower$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_REP$Te_lambda[x],param_all_REP$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

colnames(struct_mat_REP_U)<-c("ND_U", "FD_U", "Feasibility_U")
colnames(struct_mat_REP_L)<-c("ND_L", "FD_L", "Feasibility_L")

struct_mat_REP<-cbind(param_all_REP, struct_mat_REP,struct_mat_REP_L,struct_mat_REP_U)

# To create the boundaries
bound_struct_rk_w0<-data.frame(ND=seq(0,1, 0.01))
bound_struct_rk_w0$FD<-45*bound_struct_rk_w0$ND

struct_mat_REP3<-struct_mat_REP

### 3.1.2 - estimating structural stability per replicate
struct_mat_w0<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_w0[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_w0$Te_intra[x],param_all_w0$Te_inter[x],param_all_w0$Tu_inter[x], param_all_w0$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_w0$Te_lambda[x],param_all_w0$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

colnames(struct_mat_w0)<-c("ND", "FD", "Feasibility")

struct_mat_w0_L<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_w0_lower[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_w0_upper$Te_intra[x],param_all_w0_upper$Te_inter[x],param_all_w0_upper$Tu_inter[x], param_all_w0_upper$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_w0$Te_lambda[x],param_all_w0$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

colnames(struct_mat_w0_L)<-c("ND_L", "FD_L", "Feasibility_L")

struct_mat_w0_U<-as.data.frame(t(as.data.frame(sapply(c(1:length(param_all_w0_upper[,1])), function(x){
  #print(x)
  aux_alpha<-matrix(c(param_all_w0_lower$Te_intra[x],param_all_w0_lower$Te_inter[x],param_all_w0_lower$Tu_inter[x], param_all_w0_lower$Tu_intra[x]), ncol=2, byrow=TRUE)
  aux_lambda<-c(param_all_w0$Te_lambda[x],param_all_w0$Tu_lambda[x] )
  
  om<-Omega(aux_alpha)
  tta<-theta(aux_alpha, aux_lambda)
  feas<- test_feasibility(aux_alpha, aux_lambda)
  
  c(om, tta, feas)
}))))

colnames(struct_mat_w0_U)<-c("ND_U", "FD_U", "Feasibility_U")

struct_mat_w0<-cbind(param_all_w0, struct_mat_w0,struct_mat_w0_L,struct_mat_w0_U)

bound_struct_rk_w0<-data.frame(ND=seq(0,1, 0.01))
bound_struct_rk_w0$FD<-45*bound_struct_rk_w0$ND


### Distance to the edge and who wins

struct_mat_REP$a21_a11<-struct_mat_REP$Te_inter/struct_mat_REP$Tu_intra
struct_mat_REP$a22_a12<-struct_mat_REP$Te_intra/struct_mat_REP$Tu_inter

## Now to calculate the upper and lower bounds we have to see which is the lowest value of competition (and those create the upper bounds) or the highest values of competition (those create the lower boundaries)

struct_mat_REP$a21_a11_upper<-param_all_REP_lower$Te_inter/param_all_REP_lower$Tu_intra
struct_mat_REP$a22_a12_upper<-param_all_REP_lower$Te_intra/param_all_REP_lower$Tu_inter

struct_mat_REP$a21_a11_lower<-param_all_REP_upper$Te_inter/param_all_REP_upper$Tu_intra
struct_mat_REP$a22_a12_lower<-param_all_REP_upper$Te_intra/param_all_REP_upper$Tu_inter

struct_mat_REP$Tu_lambda_lower<-param_all_REP_lower$Tu_lambda
struct_mat_REP$Te_lambda_lower<-param_all_REP_lower$Te_lambda
struct_mat_REP$Tu_lambda_upper<-param_all_REP_upper$Tu_lambda
struct_mat_REP$Te_lambda_upper<-param_all_REP_upper$Te_lambda

#### order
struct_mat_REP$min_a21_a11<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  min(c(abs(struct_mat_REP$a21_a11_lower[x]),abs(struct_mat_REP$a21_a11_upper[x])))})

struct_mat_REP$min_a22_a12<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  min(c(abs(struct_mat_REP$a22_a12_lower[x]),abs(struct_mat_REP$a22_a12_upper[x])))})

struct_mat_REP$max_a21_a11<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  max(c(abs(struct_mat_REP$a21_a11_lower[x]),abs(struct_mat_REP$a21_a11_upper[x])))})

struct_mat_REP$max_a22_a12<-sapply(c(1:dim(struct_mat_REP)[1]), function(x){
  max(c(abs(struct_mat_REP$a22_a12_lower[x]),abs(struct_mat_REP$a22_a12_upper[x])))})

#write.csv(struct_mat_REP, "Analyses/structural_REP.csv")

### per replicate
struct_mat_w0$a21_a11<-struct_mat_w0$Te_inter/struct_mat_w0$Tu_intra
struct_mat_w0$a22_a12<-struct_mat_w0$Te_intra/struct_mat_w0$Tu_inter

struct_mat_w0$a21_a11_lower<-param_all_w0_lower$Te_inter/param_all_w0_lower$Tu_intra
struct_mat_w0$a22_a12_lower<-param_all_w0_lower$Te_intra/param_all_w0_lower$Tu_inter

struct_mat_w0$a21_a11_upper<-param_all_w0_upper$Te_inter/param_all_w0_upper$Tu_intra
struct_mat_w0$a22_a12_upper<-param_all_w0_upper$Te_intra/param_all_w0_upper$Tu_inter

#write.csv(struct_mat_w0, "Analyses/structural_REP_w0.csv")
