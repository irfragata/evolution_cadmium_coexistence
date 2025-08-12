### Currently this function estimates the alpha and lambda



# Function to modify the data frame structure to fit the cxr package design
mod_df<-function(df,  sizeDF=length(df[,1]), colToUse= "design"){ #this function creates a dataframe with the structure required for the alpha_lambda estimation, # colToUse is a parameter that allows to say whether we want to use the number of individuals and competitors as predicted by the initia design (= design) or the true number (=true) --> this last one is still not implemented
  
  mat<-as.data.frame(matrix(ncol=5, nrow=sizeDF))
  colnames(mat)<-c("Regime", "background", "Focal", "Comp", "Growthrate") # focal regime, competitor regime, # of focal individuals, # competitor individuals, growthrate
  
  mat[,1]<-df$FocalSR
  mat[,2]<-df$CompSR
  
  # Putting the competitor as itself in the intraspecific treatments
  mat[,2]<-sapply(c(1:length(mat[,2])),function(x) {
    if(is.na(df$CompSR[x]))
      a<-mat[x,1]
    else
      a<-mat[x,2]
    
    
    a
  })
  
  mat[,3]<- df$Dens
  
  mat[,3]<-sapply(c(1:length(mat[,2])),function(x) {
    if(colToUse=="design"){
      if(mat$Regime[x]==mat$background[x])
        a <-mat[x,3]
      else
        a<- 1
    }
    a
  })
  
  mat[,4]<- df$Dens
  
  mat[,4]<-sapply(c(1:length(mat[,2])),function(x) {
    if(colToUse=="design"){
      if(mat$Regime[x]==mat$background[x])
        a <- 0
      else
        a<- mat[x,4] -1
    }
    a
  }) 
  
  mat[,5]<-df$GrowthRateOA
  
  return(mat) 
}

######### Function for the optimization procedure
# Riker function

### Fitting with negative binomial
compmodel0<-function(par, data){
  a<-par[1]
  sigma<- par[2]
  t<-length(data)
  
  pred<-rep(a, times=t) 
  
  
  llik<-dnorm(data, log(pred), sd=sigma, log=TRUE)
  
  if(length(which(is.na(llik)))>0){
    llik<- -10000
  }
  
  return(sum(-1*llik))  
}

## Fitting a negative binomial with interaction regardless of the species identity
# Only estimating intra not estimating lambda

compmodel1<-function(par, data, dens2, l1){
  lambda<-l1
  alpha_intra<-par[1]
  sigma<-par[2]
  dens<-dens2-1 # removing one individual because it is the focal individual
  
  #red<-a/(1+alpha*(dens)) 
  pred<-lambda*exp(-alpha_intra*dens)
  
  llik<-dnorm(data, log(pred), sd=sigma, log=TRUE)
  
  if(length(which(is.na(llik)))>0){
    llik<- -10000
  }
  
  return(sum(-1*llik))  
}

#model 3 - all species have different competitive effects, but we only reestimate the interspecific competition
#The lambda and intraspecific competition comes from the estimation with the last model

compmodel2<-function(par, data1, dens1, dens2, l1, i1){
  lambda<-l1
  alpha_intra<-i1
  alpha_inter<-par[1]# Starting with the same value of intra and inter
  sigma<-par[2]
  
  #tot_density<-apply(dens1[1,], MARGIN=2, FUN=sum)
  
  #pred<- lambda/ (1+ alpha_intra*dens1 + alpha_inter* dens2)
  pred<-lambda *exp(-alpha_inter*dens2)
  
  llik<-dnorm(data1, mean=log(pred), sd=sigma, log=TRUE)
  
  if(length(which(is.na(llik)))>0){
    llik<- -10000
  }
  
  return(sum(-1*llik))  
}

magic_rk<-function(filepath2,  data2, reps2, env2, comparisons, ...){
  
# 
#   filepath2 = "./Analyses/magic_riker/"
#   data2=rep2
#   reps2=1
#   env2="N"
#   comparisons = comparison_mat

filepath<-filepath2 #where to keep the files
#lam<-lam2 # matrix with the initial estimates of lambda for the two regimes used
data<-data2 # matrix with the data, it has to have 5 columns: Regime, background, # focal ind, # comp ind, growth rate
reps <-reps2 #replicate
env<-env2 #environment
comparison_mat<-comparisons #matrix with the comparisons to use
# example
# comparison_mat<-matrix(nrow=3, ncol=3)
# comparison_mat[1,]<-c(1,4,5)
# comparison_mat[2,]<-c(4,1,NA)
# comparison_mat[3,]<-c(5,1,NA)


parameterList<-list(alpha_matrix<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # where we keep the alpha information
                    lower_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # matrix that keeps the alpha lower CI
                    upper_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))),# matrix that keeps the alpha upper CI
                    lambda_est<-NULL, #store the lambda estimates
                    lower_lambda <-NULL, # store the lower lambda CI
                    upper_lambda <-NULL, # store the upper lambda CI
                    sigma_est<-NULL, #store the estimates for the sigma (error)
                    convergence_code<-NULL # convergence code to be sure that everyhing is ok
)
names(parameterList)<-c("alpha", "lower_alpha", "upper_alpha", "lambda", "lower_lambda", "upper_lambda", "sigma", "convergence_code")
aux_vec_end<-NULL
# we have to check which data points are NA
toRemove<-which(is.na(data$Growthrate))
if(length(toRemove)!=0){
  data<-data[-toRemove,]
}

compmodel0<-function(par, data){
  a<-par[1]
  sigma<- par[2]
  t<-length(data)
  
  pred<-rep(a, times=t) 
  
  
  llik<-dnorm(data, log(pred), sd=sigma, log=TRUE)
  
  if(length(which(is.na(llik)))>0){
    llik<- -10000
  }
  
  return(sum(-1*llik))  
}

## Fitting a negative binomial with interaction regardless of the species identity
# Only estimating intra not estimating lambda

compmodel1<-function(par, data, dens2, l1){
  lambda<-l1
  alpha_intra<-par[1]
  sigma<-par[2]
  dens<-dens2-1 # removing one individual because it is the focal individual
  
  #red<-a/(1+alpha*(dens)) 
  pred<-lambda*exp(-alpha_intra*dens)
  
  llik<-dnorm(data, log(pred), sd=sigma, log=TRUE)
  
  if(length(which(is.na(llik)))>0){
    llik<- -10000
  }
  
  return(sum(-1*llik))  
}

#model 3 - all species have different competitive effects, but we only reestimate the interspecific competition
#The lambda and intraspecific competition comes from the estimation with the last model

compmodel2<-function(par, data1, dens1, dens2, l1, i1){
  lambda<-l1
  alpha_intra<-i1
  alpha_inter<-par[1]# Starting with the same value of intra and inter
  sigma<-par[2]
  
  #tot_density<-apply(dens1[1,], MARGIN=2, FUN=sum)
  
  #pred<- lambda/ (1+ alpha_intra*dens1 + alpha_inter* dens2)
  pred<-lambda *exp(-alpha_inter*dens2)
  
  llik<-dnorm(data1, mean=log(pred), sd=sigma, log=TRUE)
  
  if(length(which(is.na(llik)))>0){
    llik<- -10000
  }
  
  return(sum(-1*llik))  
}

#Since the regimes are named as numbers, then we need to put them as chars
data$Regime<-as.character(data$Regime)
data$background<-as.character(data$background)
#print(data)

pdf(file=paste(filepath,"plot_",  reps,env,".pdf", sep = ""), width=11, height=8)

for(i in 1:length(comparison_mat[,1])){
  # subset out the rows that are needed from the competition matrix
  comp_points<-subset(data, data$Regime==as.character(comparison_mat[i,1]), drop=TRUE)
  
  ## now need to build up a vector of nonzero mite production
  ## and corresponding density vectors (held in a matrix) for background species 
  ## to use in model fitting
  
  ## start with the lambda rate- will add on to this for each background:
  
  rate<-NULL
  
  ##build density matrix (each row will be density of a species)
  dens<-matrix(0, nrow=length(comparison_mat[i,]), ncol=(nrow(comp_points)))
  row.names(dens)<-comparison_mat[i,]
  
  #print("here2")
  #print(i)
  ##for each background species the target competes against:
  
  # if there isn't enough comparisons (because in replicate 2 there is one less comparison to do)
  if(is.na(comparison_mat[i,3])){
    background_list<-comparison_mat[i,1:2]
  }else{
    background_list<-comparison_mat[i,]
  }
  
  ## use this counter to keep track of which column is next to have data added to
  ## set it to begin after the lambda points:
  # print(comp_points)
  
  start <- 1
  #print( background_list)
  
  for(j in 1:length(background_list)){
    ## LOOP creates a rate vector and a corresponding dens ("density") matrix with 
    ## background species as rows, and columns corresponding to the seed production
    ## vector.
    ## beginning columns correspond to lambda
    
    #take just the rows pertaining to a specific background sp:
    bg_points<- subset(comp_points, comp_points$background==background_list[j])
    
    ## which row of the density matrix corresponds?
    #The inter
    rownum <- which(row.names(dens)==background_list[j]) #here is the row to put the competitor individuals (of other species)
    #The intra
    rownumIntra<-comparison_mat[i,1] # here is row where to keep the number of individuals of the focal species
    #column to end with:
    end<-start + nrow(bg_points)-1
    #print(start)
    #print(end)
    
    
    ## drop in density values into matrix:
    
    dens[rownum, start:end]<-bg_points$Comp
    dens[as.character(rownumIntra), start:end]<-bg_points$Focal #always do this one last, because in intraspecific competition it will add the number of individuals
    
    ## add seed numbers into the rate vector
    rate<-c(rate, bg_points$Growthrate)
    
    ##increase start counter so that next species is offset from this one
    start<-end+1
  }
  
  #print(dens)
  
  #print(paste("Before starting models, its rep ",reps, comparison_mat[i,], background_list[j], sep=" "))
  ## we'll be working with log rate (for giving a lognormal error structure):
  
  log_rate<-log(rate+1)
  
  # All start from the same point
  par1<-c(1, 0.1)
  
  aux_num<-length(dens[,1])
  if(aux_num==3){
    #subsetting to get only the intraspecific density
    aux_dens<-dens[1,which(dens[2,]==0 & dens[3,]==0 & dens[1,]==1)]
    aux_lograte<-log_rate[which(dens[2,]==0 & dens[3,]==0  & dens[1,]==1)]
  }else if(aux_num==2){
    aux_dens<-dens[1,which(dens[2,]==0 & dens[1,]==1)]
    aux_lograte<-log_rate[which(dens[2,]==0 & dens[1,]==1)]
  }
  
  for(k in 1:100){
    testcomp1<-optim(par1, compmodel0, data=aux_lograte,  control = list(maxit=10000000, gamma=5), hessian = TRUE)
    par1<-testcomp1$par
    if(testcomp1$convergence==0){
      print(paste(comparison_mat[i,j],  "model 0 converged on rep", k, sep=" "))
      break
    }
  }
  #print(testcomp1$par)
  #############################
  ## model 2, one alpha      ##
  #############################
  
  ## pars here are alpha and sigma- alpha starts at 0.1
  
  par2<-c(0.1, 0.1)
  
  aux_num<-length(dens[,1])
  if(aux_num==3){
    #subsetting to get only the intraspecific density
    aux_dens<-dens[1,which(dens[2,]==0 & dens[3,]==0)]
    aux_lograte<-log_rate[which(dens[2,]==0 & dens[3,]==0)]
  }else if(aux_num==2){
    aux_dens<-dens[1,which(dens[2,]==0)]
    aux_lograte<-log_rate[which(dens[2,]==0)]
  }
  
  ### estimating intraspecific interactions
  for(k in 1:100){
    testcomp2<-optim(par2, compmodel1, dens=aux_dens, data=aux_lograte, l1=testcomp1$par[1], control = list(maxit=10000000, gamma=5), hessian = TRUE)
    par2<-testcomp2$par
    if(testcomp2$convergence==0){
      print(paste(comparison_mat[i,j],  "model 2 converged on rep", k, sep=" "))
      break
    }
  }
  #testcomp2$par
  #############################
  ## model 3, unique alphas  ##
  #############################
  
  
  # From here we have the alpha and lambda fixed, now we need to estimate with each of the competitors in separate
  # So we run again model three but now subset the data for each competitor
  
  # when it reestimates everything all the time
  # for(c in 2:length(background_list)){
  #   
  #   par3<-c(testcomp2$par[1],testcomp2$par[2],testcomp2$par[2], testcomp2$par[3]) #starting from the same value as intraspecific competition
  #   data_c<-log_rate[which(comp_points[,2]==comparison_mat[i,1] | comp_points[,2]==comparison_mat[i,c] )] 
  #   dens_intra<-dens[1,]
  #   dens_inter<-dens[c,]
  #   
  #   
  #  ## estimating interspecific
  #   for(k in 1:100){
  #     testcomp3<-optim(par3,data1=data_c, dens1=dens_intra,dens2=dens_inter, fn=compmodel2, control = list(maxit=100000,gamma=5), lower=c(0,-1,-1,0.0001), hessian=T)
  #     par3<-testcomp3$par
  #     if(testcomp3$convergence==0){
  #       print(paste(comparison_mat[i,j], "model 3 converged on rep", k, sep=" "))
  #       break
  #     }
  #   }
  
  
  for(c in 2:length(background_list)){
    
    par3<-c(0.1, testcomp2$par[2]) #starting from 0.1
    dens_intra<-dens[1,which(dens[c,]!=0) ] #here we only care about the interspecific competition, so on this row we choose only the data from the focal species but for the interspecific ocompteition
    dens_inter<-dens[c,which(dens[c,]!=0) ]
    data_c<-log_rate[which(dens[c,]!=0)] # since the columns match we can subset by the columns
    
    
    ## estimating interspecific
    for(k in 1:100){
      testcomp3<-optim(par3,data1=data_c, dens1=dens_intra,dens2=dens_inter, l1=testcomp1$par[1],i1=testcomp2$par[1], fn=compmodel2, control = list(maxit=100000, gamma=5), hessian=T)
      par3<-testcomp3$par
      if(testcomp3$convergence==0){
        print(paste(comparison_mat[i,j], "model 3 converged on rep", k, sep=" "))
        break
      }
    }
    
    
    ##################################
    ### save estimates from model 2 and 3 ##
    ##################################
    
    
    parameterList$lambda_est<-c(parameterList$lambda_est, testcomp2$par[1])
    
    parameterList$sigma_est<-c(parameterList$sigma_est, testcomp1$par[2])#Taking the sigma from the intraspecific + lambda estimation
    
    parameterList$convergence_code<-c(parameterList$convergence_code, testcomp3$convergence)
    
    # inter
    tryCatch(fisher_info<-solve(testcomp3$hessian), error=function(cond){fisher_info<-matrix(NA, ncol=1, nrow=1)
    print("Problem with hessian model2")
      message(cond)} )
    #lambda and intra
    tryCatch(fisher_info_intra<-solve(testcomp2$hessian), error=function(cond){print("Problem with hessian model3")
      message(cond)
      fisher_info_intra<-matrix(NA, ncol=2, nrow=2) } )
    tryCatch(fisher_info_lambda<-solve(testcomp1$hessian), error=function(cond){print("Problem with hessian model1")
      message(cond)
      fisher_info_lambda<-matrix(NA, ncol=2, nrow=2) } )
    
    #print("three")
    prop_sigma_inter<-sqrt(diag(fisher_info))
    upper_inter<-testcomp3$par[1]+1.96*prop_sigma_inter[1]
    lower_inter<-testcomp3$par[1]-1.96*prop_sigma_inter[1]
    
    prop_sigma_intra<-sqrt(diag(fisher_info_intra))
    upper_intra<-testcomp2$par[1]+1.96*prop_sigma_intra[1]
    lower_intra<-testcomp2$par[1]-1.96*prop_sigma_intra[1]
    
    prop_sigma_lambda<-sqrt(diag(fisher_info_lambda))
    parameterList$upper_lambda<- c(parameterList$upper_lambda, testcomp1$par[1]+1.96*prop_sigma_lambda[1])
    parameterList$lower_lambda <- c(parameterList$lower_lambda, testcomp1$par[1]-1.96*prop_sigma_lambda[1])                 
    
    
    ## in keeping with Lotka Volterra notation, we'll use alpha1_2 to indicate effect of
    ## sp 2 on growth of 1.  Following convention, i refers to rows and j to cols in a matrix
    ## so each step of the loop here (for a target sp) corresponds to one row of this matrix:
    
    parameterList$alpha[i,1]<-testcomp2$par[1]
    parameterList$alpha[i,c]<-testcomp3$par[1]
    
    parameterList$lower_alpha[i,1] <- lower_intra[1]
    parameterList$upper_alpha[i,1] <- upper_intra[1]
    
    parameterList$lower_alpha[i,c] <- lower_inter[1]
    parameterList$upper_alpha[i,c] <- upper_inter[1]
    
    ##note that in cases where there is no data for a particular species the alpha estimate for that species ends up as the starting value- we need to be careful of these as they are basically gargbage numbers.  Keeping them in up to now to keep the structure of the data constant, but will set them to NA here:
    
    #identify which species have no data in this fit:
    which(apply(dens, MARGIN=1, FUN=mean)==0)->no_data
    
    ## set their alphas to NA in the matrix:
    
    alpha_matrix[i,no_data]<-NA
    
    ###############################
    ## some diagnostics +  plots ##
    ###############################
    
    ## print an error to the console if any one of the three models failed to converge:
    
    if( testcomp2$convergence + testcomp3$convergence !=0){
      print(paste("at least one model did not converge for", comparison_mat[i,j], sep=" "))
    }
    
    splist<-c(comparison_mat[i,c(1,c)])
    ####### plots##########
    
    ##################################
    ## plot observed vs predicted:
    ##################
    
    par<-c(testcomp1$par[1], testcomp2$par[1], testcomp3$par[1])
    
    #from model 3 code:
    
    lambda<-par[1] #same as model 2
    a_intra<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
    a_inter<-par[3]
    
    pred <- lambda *  exp(-a_intra* dens[1,] - a_inter * dens[c,])
    
    min_max<-(c(min(c(log(rate+1), log(pred+1))), max(c(log(rate+1), log(pred+1)))))
    
    plot(log(rate+1), log(pred+1), xlim=min_max, ylim=min_max, xlab="log(observed rate)", ylab="log(predicted rate)", main=comparison_mat[i,1] )
    
    abline(a=0, b=1, lwd=2)
    
    
    #####################
    #### plot each fit for intra and inter
    ##########
    alphas<-c(a_intra, a_inter)
    names(alphas)<-splist
    
    plotlist<-as.character(splist)
    for(l in 1:length(plotlist)){
      
    lam_points<-dens[c(1,c),]
      
    ## which columns in the density dataframe have nozero values for species l ?
    

    if(l==1){
      cols<-which(dens[2,]==0 & dens[3,]==0)
    }else if(l>1){
      cols<-which(dens[l,]>0)
    }
      x<-dens[l,cols]
      y<-rate[cols]
      
      
      ##add lambdas:
      x<-c(x, rep(0,ncol(dens)))
      y<-c(y, rate)
      
      x_det<-seq(min(x), max(x), by=((max(x)-min(x))/1000  ))
      
      alpha_temp<-alphas[which(names(alphas)==plotlist[l])]
      y_pred<-lambda*exp(-alpha_temp*x_det)
      
      if(l==1){
        name_list<-paste("intra", plotlist[l], sep=" ")
      }else if(l==2){
        name_list<-paste("inter", plotlist[l], sep=" ")
      }
      
      if(length(unique(x))>1){
        plot(y~x, xlab="density", ylab="rate", main=name_list)
      } else {
        
        plot(x=0, y=0, main=name_list, type='n')
        text(0, 0, "no data")
        
      }
      
      lines(y_pred~x_det, col="red", lwd=2)	
    }
    aux_vec_end<- c(aux_vec_end, comparison_mat[i,1])
  }  # end of background list loop (c)
  
  #######
  
  #vector to store regime for lambda
  
  
} # end of big for loop
  
  dev.off()
    
results_same_time<-data.frame(aux_vec_end, parameterList$lambda_est, parameterList$lower_lambda, parameterList$upper_lambda, parameterList$sigma_est, parameterList$convergence_code)

names(results_same_time)<-c("species", "lambda", "lower_error", "upper error", "sigma", "convergence_code")
alpha_matrix<-parameterList$alpha
lower_alpha<-parameterList$lower_alpha
upper_alpha<-parameterList$upper_alpha

write.csv(results_same_time, paste(filepath, "lambda_estimates_", reps,env, ".csv", sep = ""))
write.csv(alpha_matrix, paste(filepath,"alpha_estimates_row_is_target_",  reps,env,".csv", sep = ""))
write.csv(lower_alpha, paste(filepath,"alpha_lower_errors_", reps,env,".csv", sep = ""))
write.csv(upper_alpha, paste(filepath,"alpha_upper_errors_", reps,env,".csv", sep = ""))

}



####################################################################################################################################################
### Riker model using lambda directly estimated from the data set                                                                               ####
####################################################################################################################################################

magic_rk_lambda<-function(filepath2,  data2, reps2, env2, comparisons, lam2, ...){
  ## Fitting a negative binomial with interaction regardless of the species identity
  # Only estimating intra not estimating lambda
  
  compmodel1<-function(par, data, dens2, l1){
    lambda<-l1
    alpha_intra<-par[1]
    sigma<-par[2]
    dens<-dens2-1 # removing one individual because it is the focal individual
    
    #red<-a/(1+alpha*(dens)) 
    pred<-lambda*exp(-alpha_intra*dens)
    
    llik<-dnorm(data, log(pred), sd=sigma, log=TRUE)
    
    if(length(which(is.na(llik)))>0){
      llik<- -1000000
    }
    
    return(sum(-1*llik))  
  }
  
  #model 3 - all species have different competitive effects, but we only reestimate the interspecific competition
  #The lambda and intraspecific competition comes from the estimation with the last model
  
  compmodel2<-function(par, data1, dens2, l1){
    lambda<-l1
    alpha_inter<-par[1]# Starting with the same value of intra and inter
    sigma<-par[2]
    
    #tot_density<-apply(dens1[1,], MARGIN=2, FUN=sum)
    
    #pred<- lambda/ (1+ alpha_intra*dens1 + alpha_inter* dens2)
    pred<-lambda* exp(-alpha_inter*dens2)
    
    llik<-dnorm(data1, mean=log(pred), sd=sigma, log=TRUE)
    
    if(length(which(is.na(llik)))>0){
      llik<- -1000000
    }
    
    return(sum(-1*llik))  
  }
  
  # 
  #   filepath2 = "./Analyses/magic_riker/"
  #   data2=rep2
  #   reps2=1
  #   env2="N"
  #   comparisons = comparison_mat
  
  filepath<-filepath2 #where to keep the files
  lam<-lam2 # matrix with the initial estimates of lambda for the two regimes used
  data<-data2 # matrix with the data, it has to have 5 columns: Regime, background, # focal ind, # comp ind, growth rate
  reps <-reps2 #replicate
  env<-env2 #environment
  comparison_mat<-comparisons #matrix with the comparisons to use
  # example
  # comparison_mat<-matrix(nrow=3, ncol=3)
  # comparison_mat[1,]<-c(1,4,5)
  # comparison_mat[2,]<-c(4,1,NA)
  # comparison_mat[3,]<-c(5,1,NA)
  
  
  parameterList<-list(alpha_matrix<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # where we keep the alpha information
                      lower_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # matrix that keeps the alpha lower CI
                      upper_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))),# matrix that keeps the alpha upper CI
                      lambda_est<-NULL, #store the lambda estimates
                      lower_lambda <-NULL, # store the lower lambda CI
                      upper_lambda <-NULL, # store the upper lambda CI
                      sigma_est<-NULL, #store the estimates for the sigma (error)
                      convergence_code<-NULL # convergence code to be sure that everyhing is ok
  )
  names(parameterList)<-c("alpha", "lower_alpha", "upper_alpha", "lambda", "lower_lambda", "upper_lambda", "sigma", "convergence_code")
  aux_vec_end<-NULL
  # we have to check which data points are NA
  toRemove<-which(is.na(data$Growthrate))
  if(length(toRemove)!=0){
    data<-data[-toRemove,]
  }
  
  
  #Since the regimes are named as numbers, then we need to put them as chars
  data$Regime<-as.character(data$Regime)
  data$background<-as.character(data$background)
  #print(data)
  
  pdf(file=paste(filepath,"plot_",  reps,env,".pdf", sep = ""), width=11, height=8)
  
  for(i in 1:length(comparison_mat[,1])){
    # subset out the rows that are needed from the competition matrix
    comp_points<-subset(data, data$Regime==as.character(comparison_mat[i,1]), drop=TRUE)
    
    ## now need to build up a vector of nonzero mite production
    ## and corresponding density vectors (held in a matrix) for background species 
    ## to use in model fitting
    
    ## start with the lambda rate- will add on to this for each background:
    
    rate<-NULL
    
    ##build density matrix (each row will be density of a species)
    dens<-matrix(0, nrow=length(comparison_mat[i,]), ncol=(nrow(comp_points)))
    row.names(dens)<-comparison_mat[i,]
    
    #print("here2")
    #print(i)
    ##for each background species the target competes against:
    
    # if there isn't enough comparisons (because in replicate 2 there is one less comparison to do)
    if(is.na(comparison_mat[i,3])){
      background_list<-comparison_mat[i,1:2]
    }else{
      background_list<-comparison_mat[i,]
    }
    
    ## use this counter to keep track of which column is next to have data added to
    ## set it to begin after the lambda points:
    # print(comp_points)
    
    start <- 1
    #print( background_list)
    
    for(j in 1:length(background_list)){
      ## LOOP creates a rate vector and a corresponding dens ("density") matrix with 
      ## background species as rows, and columns corresponding to the seed production
      ## vector.
      ## beginning columns correspond to lambda
      
      #take just the rows pertaining to a specific background sp:
      bg_points<- subset(comp_points, comp_points$background==background_list[j])
      
      ## which row of the density matrix corresponds?
      #The inter
      rownum <- which(row.names(dens)==background_list[j]) #here is the row to put the competitor individuals (of other species)
      #The intra
      rownumIntra<-comparison_mat[i,1] # here is row where to keep the number of individuals of the focal species
      #column to end with:
      end<-start + nrow(bg_points)-1
      #print(start)
      #print(end)
      
      
      ## drop in density values into matrix:
      
      dens[rownum, start:end]<-bg_points$Comp
      dens[as.character(rownumIntra), start:end]<-bg_points$Focal #always do this one last, because in intraspecific competition it will add the number of individuals
      
      ## add seed numbers into the rate vector
      rate<-c(rate, bg_points$Growthrate)
      
      ##increase start counter so that next species is offset from this one
      start<-end+1
    }
    
    #print(dens)
    
    #print(paste("Before starting models, its rep ",reps, comparison_mat[i,], background_list[j], sep=" "))
    ## we'll be working with log rate (for giving a lognormal error structure):
    
    log_rate<-log(rate+1)
    
    
    ### 
    ### Skipping estimating only lambda, lambda will be estimated directly from the data set
    ###
   
    #############################
    ## model 2, one alpha      ##
    #############################
    
    ## pars here are alpha and sigma- alpha starts at 0.1
    
    par2<-c(0.1, 0.1)
    
    aux_num<-length(dens[,1])
    if(aux_num==3){
      #subsetting to get only the intraspecific density
      aux_dens<-dens[1,which(dens[2,]==0 & dens[3,]==0)]
      aux_lograte<-log_rate[which(dens[2,]==0 & dens[3,]==0)]
    }else if(aux_num==2){
      aux_dens<-dens[1,which(dens[2,]==0)]
      aux_lograte<-log_rate[which(dens[2,]==0)]
    }
    
    ### estimating intraspecific interactions
    for(k in 1:100){
      testcomp2<-optim(par2, compmodel1, dens=aux_dens, data=aux_lograte, l1=subset(lam, SR==comparison_mat[i,1])$lambda, control = list(maxit=10000000), hessian = TRUE)
      par2<-testcomp2$par
      if(testcomp2$convergence==0){
        print(paste(comparison_mat[i,j],  "model 2 converged on rep", k, sep=" "))
        break
      }
    }

    for(c in 2:length(background_list)){
      
      par3<-c(0.1, 0.1) #starting from 0.1
      #here we only care about the interspecific competition, so on this row we choose only the data from the focal species but for the interspecific ocompteition
      dens_inter<-dens[c,which(dens[c,]!=0) ]
      data_c<-log_rate[which(dens[c,]!=0)] # since the columns match we can subset by the columns
      
      
      ## estimating interspecific
      for(k in 1:100){
        testcomp3<-optim(par3,fn=compmodel2,data1=data_c,dens2=dens_inter, l1=subset(lam, SR==comparison_mat[i,1])$lambda,  control = list(maxit=100000, gamma=5), hessian=T)
        par3<-testcomp3$par
        if(testcomp3$convergence==0){
          print(paste(comparison_mat[i,j], "model 3 converged on rep", k, sep=" "))
          break
        }
      }
      
      
      ##################################
      ### save estimates from model 2 and 3 ##
      ##################################
      
      
      parameterList$lambda_est<-c(parameterList$lambda_est, subset(lam, SR==comparison_mat[i,1])$lambda)
      
      parameterList$sigma_est<-c(parameterList$sigma_est, subset(lam, SR==comparison_mat[i,1])$sd_lambda)#Taking the sigma from the intraspecific + lambda estimation
      
      parameterList$convergence_code<-c(parameterList$convergence_code, testcomp3$convergence)
      
      # inter
      tryCatch(fisher_info<-solve(testcomp3$hessian), error=function(cond){fisher_info<-matrix(NA, ncol=1, nrow=1)
      print("Problem with hessian model2")
      message(cond)} )
      #lambda and intra
      tryCatch(fisher_info_intra<-solve(testcomp2$hessian), error=function(cond){print("Problem with hessian model3")
        message(cond)
        fisher_info_intra<-matrix(NA, ncol=2, nrow=2) } )
      #tryCatch(fisher_info_lambda<-solve(testcomp1$hessian), error=function(cond){print("Problem with hessian model1")
        #message(cond)
        #fisher_info_lambda<-matrix(NA, ncol=2, nrow=2) } )
      
      #print("three")
      prop_sigma_inter<-sqrt(diag(fisher_info))
      upper_inter<-testcomp3$par[1]+1.96*prop_sigma_inter[1]
      lower_inter<-testcomp3$par[1]-1.96*prop_sigma_inter[1]
      
      prop_sigma_intra<-sqrt(diag(fisher_info_intra))
      upper_intra<-testcomp2$par[1]+1.96*prop_sigma_intra[1]
      lower_intra<-testcomp2$par[1]-1.96*prop_sigma_intra[1]
      
      #prop_sigma_lambda<-sqrt(diag(fisher_info_lambda))
      parameterList$upper_lambda<- c(parameterList$upper_lambda, subset(lam, SR==comparison_mat[i,1])$lambda+1.96*subset(lam, SR==comparison_mat[i,1])$sd_lambda)
      parameterList$lower_lambda <- c(parameterList$lower_lambda, subset(lam, SR==comparison_mat[i,1])$lambda-1.96*subset(lam, SR==comparison_mat[i,1])$sd_lambda)                 
      
      
      ## in keeping with Lotka Volterra notation, we'll use alpha1_2 to indicate effect of
      ## sp 2 on growth of 1.  Following convention, i refers to rows and j to cols in a matrix
      ## so each step of the loop here (for a target sp) corresponds to one row of this matrix:
      
      parameterList$alpha[i,1]<-testcomp2$par[1]
      parameterList$alpha[i,c]<-testcomp3$par[1]
      
      parameterList$lower_alpha[i,1] <- lower_intra[1]
      parameterList$upper_alpha[i,1] <- upper_intra[1]
      
      parameterList$lower_alpha[i,c] <- lower_inter[1]
      parameterList$upper_alpha[i,c] <- upper_inter[1]
      
      ##note that in cases where there is no data for a particular species the alpha estimate for that species ends up as the starting value- we need to be careful of these as they are basically gargbage numbers.  Keeping them in up to now to keep the structure of the data constant, but will set them to NA here:
      
      #identify which species have no data in this fit:
      which(apply(dens, MARGIN=1, FUN=mean)==0)->no_data
      
      ## set their alphas to NA in the matrix:
      
      alpha_matrix[i,no_data]<-NA
      
      ###############################
      ## some diagnostics +  plots ##
      ###############################
      
      ## print an error to the console if any one of the three models failed to converge:
      
      if( testcomp2$convergence + testcomp3$convergence !=0){
        print(paste("at least one model did not converge for", comparison_mat[i,j], sep=" "))
      }
      
      splist<-c(comparison_mat[i,c(1,c)])
      ####### plots##########
      
      ##################################
      ## plot observed vs predicted:
      ##################
      
      par<-c(subset(lam, SR==comparison_mat[i,1])$lambda, testcomp2$par[1], testcomp3$par[1])
      
      #from model 3 code:
      
      lambda<-par[1] #same as model 2
      a_intra<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
      a_inter<-par[3]
      
      pred <- lambda* exp(-a_intra* (dens[1,]-1) - a_inter * dens[c,])
      
      min_max<-c(min(c(log(rate+1), log(pred+1))), max(c(log(rate+1), log(pred+1))))
      
      plot(log(rate+1), log(pred+1), xlim=min_max, ylim=min_max, xlab="log(observed rate)", ylab="log(predicted rate)", main=comparison_mat[i,1] )
      
      abline(a=0, b=1, lwd=2)
      
      
      #####################
      #### plot each fit for intra and inter
      ##########
      alphas<-c(a_intra, a_inter)
      names(alphas)<-splist
      
      plotlist<-as.character(splist)
      for(l in 1:length(plotlist)){
        
        lam_points<-dens[c(1,c),]
        
        ## which columns in the density dataframe have nozero values for species l ?
        
        
        if(l==1){
          cols<-which(dens[2,]==0 & dens[3,]==0)
        }else if(l>1){
          cols<-which(dens[l,]>0)
        }
        x<-dens[l,cols]
        y<-rate[cols]
        
        
        ##add lambdas:
        x<-c(x, rep(0,ncol(dens)))
        y<-c(y, rate)
        
        x_det<-seq(min(x), max(x), by=((max(x)-min(x))/1000))
        
        alpha_temp<-alphas[which(names(alphas)==plotlist[l])]
        y_pred<-lambda* exp(-alpha_temp*x_det)
        
        if(l==1){
          name_list<-paste("intra", plotlist[l], sep=" ")
        }else if(l==2){
          name_list<-paste("inter", plotlist[l], sep=" ")
        }
        
        if(length(unique(x))>1){
          plot(y~x, xlab="density", ylab="rate", main=name_list)
        } else {
          
          plot(x=0, y=0, main=name_list, type='n')
          text(0, 0, "no data")
          
        }
        
        lines(y_pred~x_det, col="red", lwd=2)	
      }
      aux_vec_end<- c(aux_vec_end, comparison_mat[i,1])
    }  # end of background list loop (c)
    
    #######
    
    #vector to store regime for lambda
    
    
  } # end of big for loop
  
  dev.off()
  
  results_same_time<-data.frame(aux_vec_end, parameterList$lambda_est, parameterList$lower_lambda, parameterList$upper_lambda, parameterList$sigma_est, parameterList$convergence_code)
  
  names(results_same_time)<-c("species", "lambda", "lower_error", "upper error", "sigma", "convergence_code")
  alpha_matrix<-parameterList$alpha
  lower_alpha<-parameterList$lower_alpha
  upper_alpha<-parameterList$upper_alpha
  
  write.csv(results_same_time, paste(filepath, "lambda_estimates_", reps,env, ".csv", sep = ""))
  write.csv(alpha_matrix, paste(filepath,"alpha_estimates_row_is_target_",  reps,env,".csv", sep = ""))
  write.csv(lower_alpha, paste(filepath,"alpha_lower_errors_", reps,env,".csv", sep = ""))
  write.csv(upper_alpha, paste(filepath,"alpha_upper_errors_", reps,env,".csv", sep = ""))
  
}



magic_rk_lambda_bounded<-function(filepath2,  data2, reps2, env2, comparisons, lam2, ...){
  ## Fitting a negative binomial with interaction regardless of the species identity
  # Only estimating intra not estimating lambda
  
  compmodel1<-function(par, data, dens2, l1){
    lambda<-l1
    alpha_intra<-par[1]
    sigma<-par[2]
    dens<-dens2-1 # removing one individual because it is the focal individual
    
    #red<-a/(1+alpha*(dens)) 
    pred<-lambda*exp(-alpha_intra*dens)
    
    llik<-dnorm(data, log(pred), sd=sigma, log=TRUE)
    
    if(length(which(is.na(llik)))>0){
      llik<- -1000000
    }
    
    return(sum(-1*llik))  
  }
  
  #model 3 - all species have different competitive effects, but we only reestimate the interspecific competition
  #The lambda and intraspecific competition comes from the estimation with the last model
  
  compmodel2<-function(par, data1, dens2, l1){
    lambda<-l1
    alpha_inter<-par[1]# Starting with the same value of intra and inter
    sigma<-par[2]
    
    #tot_density<-apply(dens1[1,], MARGIN=2, FUN=sum)
    
    #pred<- lambda/ (1+ alpha_intra*dens1 + alpha_inter* dens2)
    pred<-lambda* exp(-alpha_inter*dens2)
    
    llik<-dnorm(data1, mean=log(pred), sd=sigma, log=TRUE)
    
    if(length(which(is.na(llik)))>0){
      llik<- -1000000
    }
    
    return(sum(-1*llik))  
  }
  
  # 
  #   filepath2 = "./Analyses/magic_riker/"
  #   data2=rep2
  #   reps2=1
  #   env2="N"
  #   comparisons = comparison_mat
  
  filepath<-filepath2 #where to keep the files
  lam<-lam2 # matrix with the initial estimates of lambda for the two regimes used
  data<-data2 # matrix with the data, it has to have 5 columns: Regime, background, # focal ind, # comp ind, growth rate
  reps <-reps2 #replicate
  env<-env2 #environment
  comparison_mat<-comparisons #matrix with the comparisons to use
  # example
  # comparison_mat<-matrix(nrow=3, ncol=3)
  # comparison_mat[1,]<-c(1,4,5)
  # comparison_mat[2,]<-c(4,1,NA)
  # comparison_mat[3,]<-c(5,1,NA)
  
  
  parameterList<-list(alpha_matrix<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # where we keep the alpha information
                      lower_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # matrix that keeps the alpha lower CI
                      upper_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))),# matrix that keeps the alpha upper CI
                      lambda_est<-NULL, #store the lambda estimates
                      lower_lambda <-NULL, # store the lower lambda CI
                      upper_lambda <-NULL, # store the upper lambda CI
                      sigma_est<-NULL, #store the estimates for the sigma (error)
                      convergence_code<-NULL # convergence code to be sure that everyhing is ok
  )
  names(parameterList)<-c("alpha", "lower_alpha", "upper_alpha", "lambda", "lower_lambda", "upper_lambda", "sigma", "convergence_code")
  aux_vec_end<-NULL
  # we have to check which data points are NA
  toRemove<-which(is.na(data$Growthrate))
  if(length(toRemove)!=0){
    data<-data[-toRemove,]
  }
  
  
  #Since the regimes are named as numbers, then we need to put them as chars
  data$Regime<-as.character(data$Regime)
  data$background<-as.character(data$background)
  #print(data)
  
  pdf(file=paste(filepath,"plot_",  reps,env,".pdf", sep = ""), width=11, height=8)
  
  for(i in 1:length(comparison_mat[,1])){
    # subset out the rows that are needed from the competition matrix
    comp_points<-subset(data, data$Regime==as.character(comparison_mat[i,1]), drop=TRUE)
    
    ## now need to build up a vector of nonzero mite production
    ## and corresponding density vectors (held in a matrix) for background species 
    ## to use in model fitting
    
    ## start with the lambda rate- will add on to this for each background:
    
    rate<-NULL
    
    ##build density matrix (each row will be density of a species)
    dens<-matrix(0, nrow=length(comparison_mat[i,]), ncol=(nrow(comp_points)))
    row.names(dens)<-comparison_mat[i,]
    
    #print("here2")
    #print(i)
    ##for each background species the target competes against:
    
    # if there isn't enough comparisons (because in replicate 2 there is one less comparison to do)
    if(is.na(comparison_mat[i,3])){
      background_list<-comparison_mat[i,1:2]
    }else{
      background_list<-comparison_mat[i,]
    }
    
    ## use this counter to keep track of which column is next to have data added to
    ## set it to begin after the lambda points:
    # print(comp_points)
    
    start <- 1
    #print( background_list)
    
    for(j in 1:length(background_list)){
      ## LOOP creates a rate vector and a corresponding dens ("density") matrix with 
      ## background species as rows, and columns corresponding to the seed production
      ## vector.
      ## beginning columns correspond to lambda
      
      #take just the rows pertaining to a specific background sp:
      bg_points<- subset(comp_points, comp_points$background==background_list[j])
      
      ## which row of the density matrix corresponds?
      #The inter
      rownum <- which(row.names(dens)==background_list[j]) #here is the row to put the competitor individuals (of other species)
      #The intra
      rownumIntra<-comparison_mat[i,1] # here is row where to keep the number of individuals of the focal species
      #column to end with:
      end<-start + nrow(bg_points)-1
      #print(start)
      #print(end)
      
      
      ## drop in density values into matrix:
      
      dens[rownum, start:end]<-bg_points$Comp
      dens[as.character(rownumIntra), start:end]<-bg_points$Focal #always do this one last, because in intraspecific competition it will add the number of individuals
      
      ## add seed numbers into the rate vector
      rate<-c(rate, bg_points$Growthrate)
      
      ##increase start counter so that next species is offset from this one
      start<-end+1
    }
    
    #print(dens)
    
    #print(paste("Before starting models, its rep ",reps, comparison_mat[i,], background_list[j], sep=" "))
    ## we'll be working with log rate (for giving a lognormal error structure):
    
    log_rate<-log(rate+1)
    
    
    ### 
    ### Skipping estimating only lambda, lambda will be estimated directly from the data set
    ###
    
    #############################
    ## model 2, one alpha      ##
    #############################
    
    ## pars here are alpha and sigma- alpha starts at 0.1
    
    par2<-c(0.1, 0.1)
    
    aux_num<-length(dens[,1])
    if(aux_num==3){
      #subsetting to get only the intraspecific density
      aux_dens<-dens[1,which(dens[2,]==0 & dens[3,]==0)]
      aux_lograte<-log_rate[which(dens[2,]==0 & dens[3,]==0)]
    }else if(aux_num==2){
      aux_dens<-dens[1,which(dens[2,]==0)]
      aux_lograte<-log_rate[which(dens[2,]==0)]
    }
    
    ### estimating intraspecific interactions
    for(k in 1:100){
      testcomp2<-optim(par2, compmodel1, dens=aux_dens, data=aux_lograte, l1=subset(lam, SR==comparison_mat[i,1])$lambda, control = list(maxit=10000000, gamma=5), hessian = TRUE, method="L-BFGS-B", lower=c(-10,0.0001), upper=c(10,10))
      par2<-testcomp2$par
      if(testcomp2$convergence==0){
        print(paste(comparison_mat[i,j],  "model 2 converged on rep", k, sep=" "))
        break
      }
    }
    
    for(c in 2:length(background_list)){
      
      par3<-c(0.1, 0.1) #starting from 0.1
      #here we only care about the interspecific competition, so on this row we choose only the data from the focal species but for the interspecific ocompteition
      dens_inter<-dens[c,which(dens[c,]!=0) ]
      data_c<-log_rate[which(dens[c,]!=0)] # since the columns match we can subset by the columns
      
      
      ## estimating interspecific
      for(k in 1:100){
        testcomp3<-optim(par3,fn=compmodel2,data1=data_c,dens2=dens_inter, l1=subset(lam, SR==comparison_mat[i,1])$lambda,  control = list(maxit=100000, gamma=5), hessian=T, method="L-BFGS-B", lower=c(-10,0.0001), upper=c(10,10))
        par3<-testcomp3$par
        if(testcomp3$convergence==0){
          print(paste(comparison_mat[i,j], "model 3 converged on rep", k, sep=" "))
          break
        }
      }
      
      
      ##################################
      ### save estimates from model 2 and 3 ##
      ##################################
      
      
      parameterList$lambda_est<-c(parameterList$lambda_est, subset(lam, SR==comparison_mat[i,1])$lambda)
      
      parameterList$sigma_est<-c(parameterList$sigma_est, subset(lam, SR==comparison_mat[i,1])$sd_lambda)#Taking the sigma from the intraspecific + lambda estimation
      
      parameterList$convergence_code<-c(parameterList$convergence_code, testcomp3$convergence)
      
      # inter
      tryCatch(fisher_info<-solve(testcomp3$hessian), error=function(cond){fisher_info<-matrix(NA, ncol=1, nrow=1)
      print("Problem with hessian model2")
      message(cond)} )
      #lambda and intra
      tryCatch(fisher_info_intra<-solve(testcomp2$hessian), error=function(cond){print("Problem with hessian model3")
        message(cond)
        fisher_info_intra<-matrix(NA, ncol=2, nrow=2) } )
      #tryCatch(fisher_info_lambda<-solve(testcomp1$hessian), error=function(cond){print("Problem with hessian model1")
      #message(cond)
      #fisher_info_lambda<-matrix(NA, ncol=2, nrow=2) } )
      
      #print("three")
      prop_sigma_inter<-sqrt(diag(fisher_info))
      upper_inter<-testcomp3$par[1]+1.96*prop_sigma_inter[1]
      lower_inter<-testcomp3$par[1]-1.96*prop_sigma_inter[1]
      
      prop_sigma_intra<-sqrt(diag(fisher_info_intra))
      upper_intra<-testcomp2$par[1]+1.96*prop_sigma_intra[1]
      lower_intra<-testcomp2$par[1]-1.96*prop_sigma_intra[1]
      
      #prop_sigma_lambda<-sqrt(diag(fisher_info_lambda))
      parameterList$upper_lambda<- c(parameterList$upper_lambda, subset(lam, SR==comparison_mat[i,1])$lambda+1.96*subset(lam, SR==comparison_mat[i,1])$sd_lambda)
      parameterList$lower_lambda <- c(parameterList$lower_lambda, subset(lam, SR==comparison_mat[i,1])$lambda-1.96*subset(lam, SR==comparison_mat[i,1])$sd_lambda)                 
      
      
      ## in keeping with Lotka Volterra notation, we'll use alpha1_2 to indicate effect of
      ## sp 2 on growth of 1.  Following convention, i refers to rows and j to cols in a matrix
      ## so each step of the loop here (for a target sp) corresponds to one row of this matrix:
      
      parameterList$alpha[i,1]<-testcomp2$par[1]
      parameterList$alpha[i,c]<-testcomp3$par[1]
      
      parameterList$lower_alpha[i,1] <- lower_intra[1]
      parameterList$upper_alpha[i,1] <- upper_intra[1]
      
      parameterList$lower_alpha[i,c] <- lower_inter[1]
      parameterList$upper_alpha[i,c] <- upper_inter[1]
      
      ##note that in cases where there is no data for a particular species the alpha estimate for that species ends up as the starting value- we need to be careful of these as they are basically gargbage numbers.  Keeping them in up to now to keep the structure of the data constant, but will set them to NA here:
      
      #identify which species have no data in this fit:
      which(apply(dens, MARGIN=1, FUN=mean)==0)->no_data
      
      ## set their alphas to NA in the matrix:
      
      alpha_matrix[i,no_data]<-NA
      
      ###############################
      ## some diagnostics +  plots ##
      ###############################
      
      ## print an error to the console if any one of the three models failed to converge:
      
      if( testcomp2$convergence + testcomp3$convergence !=0){
        print(paste("at least one model did not converge for", comparison_mat[i,j], sep=" "))
      }
      
      splist<-c(comparison_mat[i,c(1,c)])
      ####### plots##########
      
      ##################################
      ## plot observed vs predicted:
      ##################
      
      par<-c(subset(lam, SR==comparison_mat[i,1])$lambda, testcomp2$par[1], testcomp3$par[1])
      
      #from model 3 code:
      
      lambda<-par[1] #same as model 2
      a_intra<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
      a_inter<-par[3]
      
      pred <- lambda* exp(-a_intra* (dens[1,]-1) - a_inter * dens[c,])
      
      min_max<-c(min(c(log(rate+1), log(pred+1))), max(c(log(rate+1), log(pred+1))))
      
      plot(log(rate+1), log(pred+1), xlim=min_max, ylim=min_max, xlab="log(observed rate)", ylab="log(predicted rate)", main=comparison_mat[i,1] )
      
      abline(a=0, b=1, lwd=2)
      
      
      #####################
      #### plot each fit for intra and inter
      ##########
      alphas<-c(a_intra, a_inter)
      names(alphas)<-splist
      
      plotlist<-as.character(splist)
      for(l in 1:length(plotlist)){
        
        lam_points<-dens[c(1,c),]
        
        ## which columns in the density dataframe have nozero values for species l ?
        
        
        if(l==1){
          cols<-which(dens[2,]==0 & dens[3,]==0)
        }else if(l>1){
          cols<-which(dens[l,]>0)
        }
        x<-dens[l,cols]
        y<-rate[cols]
        
        
        ##add lambdas:
        x<-c(x, rep(0,ncol(dens)))
        y<-c(y, rate)
        
        x_det<-seq(min(x), max(x), by=((max(x)-min(x))/1000))
        
        alpha_temp<-alphas[which(names(alphas)==plotlist[l])]
        y_pred<-lambda* exp(-alpha_temp*x_det)
        
        if(l==1){
          name_list<-paste("intra", plotlist[l], sep=" ")
        }else if(l==2){
          name_list<-paste("inter", plotlist[l], sep=" ")
        }
        
        if(length(unique(x))>1){
          plot(y~x, xlab="density", ylab="rate", main=name_list)
        } else {
          
          plot(x=0, y=0, main=name_list, type='n')
          text(0, 0, "no data")
          
        }
        
        lines(y_pred~x_det, col="red", lwd=2)	
      }
      aux_vec_end<- c(aux_vec_end, comparison_mat[i,1])
    }  # end of background list loop (c)
    
    #######
    
    #vector to store regime for lambda
    
    
  } # end of big for loop
  
  dev.off()
  
  results_same_time<-data.frame(aux_vec_end, parameterList$lambda_est, parameterList$lower_lambda, parameterList$upper_lambda, parameterList$sigma_est, parameterList$convergence_code)
  
  names(results_same_time)<-c("species", "lambda", "lower_error", "upper error", "sigma", "convergence_code")
  alpha_matrix<-parameterList$alpha
  lower_alpha<-parameterList$lower_alpha
  upper_alpha<-parameterList$upper_alpha
  
  write.csv(results_same_time, paste(filepath, "lambda_estimates_", reps,env, ".csv", sep = ""))
  write.csv(alpha_matrix, paste(filepath,"alpha_estimates_row_is_target_",  reps,env,".csv", sep = ""))
  write.csv(lower_alpha, paste(filepath,"alpha_lower_errors_", reps,env,".csv", sep = ""))
  write.csv(upper_alpha, paste(filepath,"alpha_upper_errors_", reps,env,".csv", sep = ""))
  
}


####################################################################################################################################################
### Riker model using lambda directly estimated from the data set                                                                               ####
####################################################################################################################################################

magic_rk_lambda<-function(filepath2,  data2, reps2, env2, comparisons, lam2, ...){
  ## Fitting a negative binomial with interaction regardless of the species identity
  # Only estimating intra not estimating lambda
  
  compmodel1<-function(par, data, dens2, l1){
    lambda<-l1
    alpha_intra<-par[1]
    sigma<-par[2]
    dens<-dens2-1 # removing one individual because it is the focal individual
    
    #red<-a/(1+alpha*(dens)) 
    pred<-lambda*exp(-alpha_intra*dens)
    
    llik<-dnorm(data, log(pred), sd=sigma, log=TRUE)
    
    if(length(which(is.na(llik)))>0){
      llik<- -1000000
    }
    
    return(sum(-1*llik))  
  }
  
  #model 3 - all species have different competitive effects, but we only reestimate the interspecific competition
  #The lambda and intraspecific competition comes from the estimation with the last model
  
  compmodel2<-function(par, data1, dens2, l1){
    lambda<-l1
    alpha_inter<-par[1]# Starting with the same value of intra and inter
    sigma<-par[2]
    
    #tot_density<-apply(dens1[1,], MARGIN=2, FUN=sum)
    
    #pred<- lambda/ (1+ alpha_intra*dens1 + alpha_inter* dens2)
    pred<-lambda* exp(-alpha_inter*dens2)
    
    llik<-dnorm(data1, mean=log(pred), sd=sigma, log=TRUE)
    
    if(length(which(is.na(llik)))>0){
      llik<- -1000000
    }
    
    return(sum(-1*llik))  
  }
  
  # 
  #   filepath2 = "./Analyses/magic_riker/"
  #   data2=rep2
  #   reps2=1
  #   env2="N"
  #   comparisons = comparison_mat
  
  filepath<-filepath2 #where to keep the files
  lam<-lam2 # matrix with the initial estimates of lambda for the two regimes used
  data<-data2 # matrix with the data, it has to have 5 columns: Regime, background, # focal ind, # comp ind, growth rate
  reps <-reps2 #replicate
  env<-env2 #environment
  comparison_mat<-comparisons #matrix with the comparisons to use
  # example
  # comparison_mat<-matrix(nrow=3, ncol=3)
  # comparison_mat[1,]<-c(1,4,5)
  # comparison_mat[2,]<-c(4,1,NA)
  # comparison_mat[3,]<-c(5,1,NA)
  
  
  parameterList<-list(alpha_matrix<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # where we keep the alpha information
                      lower_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # matrix that keeps the alpha lower CI
                      upper_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))),# matrix that keeps the alpha upper CI
                      lambda_est<-NULL, #store the lambda estimates
                      lower_lambda <-NULL, # store the lower lambda CI
                      upper_lambda <-NULL, # store the upper lambda CI
                      sigma_est<-NULL, #store the estimates for the sigma (error)
                      convergence_code<-NULL # convergence code to be sure that everyhing is ok
  )
  names(parameterList)<-c("alpha", "lower_alpha", "upper_alpha", "lambda", "lower_lambda", "upper_lambda", "sigma", "convergence_code")
  aux_vec_end<-NULL
  # we have to check which data points are NA
  toRemove<-which(is.na(data$Growthrate))
  if(length(toRemove)!=0){
    data<-data[-toRemove,]
  }
  
  
  #Since the regimes are named as numbers, then we need to put them as chars
  data$Regime<-as.character(data$Regime)
  data$background<-as.character(data$background)
  #print(data)
  
  pdf(file=paste(filepath,"plot_",  reps,env,".pdf", sep = ""), width=11, height=8)
  
  for(i in 1:length(comparison_mat[,1])){
    # subset out the rows that are needed from the competition matrix
    comp_points<-subset(data, data$Regime==as.character(comparison_mat[i,1]), drop=TRUE)
    
    ## now need to build up a vector of nonzero mite production
    ## and corresponding density vectors (held in a matrix) for background species 
    ## to use in model fitting
    
    ## start with the lambda rate- will add on to this for each background:
    
    rate<-NULL
    
    ##build density matrix (each row will be density of a species)
    dens<-matrix(0, nrow=length(comparison_mat[i,]), ncol=(nrow(comp_points)))
    row.names(dens)<-comparison_mat[i,]
    
    #print("here2")
    #print(i)
    ##for each background species the target competes against:
    
    # if there isn't enough comparisons (because in replicate 2 there is one less comparison to do)
    if(is.na(comparison_mat[i,3])){
      background_list<-comparison_mat[i,1:2]
    }else{
      background_list<-comparison_mat[i,]
    }
    
    ## use this counter to keep track of which column is next to have data added to
    ## set it to begin after the lambda points:
    # print(comp_points)
    
    start <- 1
    #print( background_list)
    
    for(j in 1:length(background_list)){
      ## LOOP creates a rate vector and a corresponding dens ("density") matrix with 
      ## background species as rows, and columns corresponding to the seed production
      ## vector.
      ## beginning columns correspond to lambda
      
      #take just the rows pertaining to a specific background sp:
      bg_points<- subset(comp_points, comp_points$background==background_list[j])
      
      ## which row of the density matrix corresponds?
      #The inter
      rownum <- which(row.names(dens)==background_list[j]) #here is the row to put the competitor individuals (of other species)
      #The intra
      rownumIntra<-comparison_mat[i,1] # here is row where to keep the number of individuals of the focal species
      #column to end with:
      end<-start + nrow(bg_points)-1
      #print(start)
      #print(end)
      
      
      ## drop in density values into matrix:
      
      dens[rownum, start:end]<-bg_points$Comp
      dens[as.character(rownumIntra), start:end]<-bg_points$Focal #always do this one last, because in intraspecific competition it will add the number of individuals
      
      ## add seed numbers into the rate vector
      rate<-c(rate, bg_points$Growthrate)
      
      ##increase start counter so that next species is offset from this one
      start<-end+1
    }
    
    #print(dens)
    
    #print(paste("Before starting models, its rep ",reps, comparison_mat[i,], background_list[j], sep=" "))
    ## we'll be working with log rate (for giving a lognormal error structure):
    
    log_rate<-log(rate+1)
    
    
    ### 
    ### Skipping estimating only lambda, lambda will be estimated directly from the data set
    ###
    
    #############################
    ## model 2, one alpha      ##
    #############################
    
    ## pars here are alpha and sigma- alpha starts at 0.1
    
    par2<-c(0.1, 0.1)
    
    aux_num<-length(dens[,1])
    if(aux_num==3){
      #subsetting to get only the intraspecific density
      aux_dens<-dens[1,which(dens[2,]==0 & dens[3,]==0)]
      aux_lograte<-log_rate[which(dens[2,]==0 & dens[3,]==0)]
    }else if(aux_num==2){
      aux_dens<-dens[1,which(dens[2,]==0)]
      aux_lograte<-log_rate[which(dens[2,]==0)]
    }
    
    ### estimating intraspecific interactions
    for(k in 1:100){
      testcomp2<-optim(par2, compmodel1, dens=aux_dens, data=aux_lograte, l1=subset(lam, SR==comparison_mat[i,1])$lambda, control = list(maxit=10000000), hessian = TRUE)
      par2<-testcomp2$par
      if(testcomp2$convergence==0){
        print(paste(comparison_mat[i,j],  "model 2 converged on rep", k, sep=" "))
        break
      }
    }
    
    for(c in 2:length(background_list)){
      
      par3<-c(0.1, 0.1) #starting from 0.1
      #here we only care about the interspecific competition, so on this row we choose only the data from the focal species but for the interspecific ocompteition
      dens_inter<-dens[c,which(dens[c,]!=0) ]
      data_c<-log_rate[which(dens[c,]!=0)] # since the columns match we can subset by the columns
      
      
      ## estimating interspecific
      for(k in 1:100){
        testcomp3<-optim(par3,fn=compmodel2,data1=data_c,dens2=dens_inter, l1=subset(lam, SR==comparison_mat[i,1])$lambda,  control = list(maxit=100000, gamma=5), hessian=T)
        par3<-testcomp3$par
        if(testcomp3$convergence==0){
          print(paste(comparison_mat[i,j], "model 3 converged on rep", k, sep=" "))
          break
        }
      }
      
      
      ##################################
      ### save estimates from model 2 and 3 ##
      ##################################
      
      
      parameterList$lambda_est<-c(parameterList$lambda_est, subset(lam, SR==comparison_mat[i,1])$lambda)
      
      parameterList$sigma_est<-c(parameterList$sigma_est, subset(lam, SR==comparison_mat[i,1])$sd_lambda)#Taking the sigma from the intraspecific + lambda estimation
      
      parameterList$convergence_code<-c(parameterList$convergence_code, testcomp3$convergence)
      
      # inter
      tryCatch(fisher_info<-solve(testcomp3$hessian), error=function(cond){fisher_info<-matrix(NA, ncol=1, nrow=1)
      print("Problem with hessian model2")
      message(cond)} )
      #lambda and intra
      tryCatch(fisher_info_intra<-solve(testcomp2$hessian), error=function(cond){print("Problem with hessian model3")
        message(cond)
        fisher_info_intra<-matrix(NA, ncol=2, nrow=2) } )
      #tryCatch(fisher_info_lambda<-solve(testcomp1$hessian), error=function(cond){print("Problem with hessian model1")
      #message(cond)
      #fisher_info_lambda<-matrix(NA, ncol=2, nrow=2) } )
      
      #print("three")
      prop_sigma_inter<-sqrt(diag(fisher_info))
      upper_inter<-testcomp3$par[1]+1.96*prop_sigma_inter
      lower_inter<-testcomp3$par[1]-1.96*prop_sigma_inter
      
      prop_sigma_intra<-sqrt(diag(fisher_info_intra))
      upper_intra<-testcomp2$par[1]+1.96*prop_sigma_intra
      lower_intra<-testcomp2$par[1]-1.96*prop_sigma_intra
      
      #prop_sigma_lambda<-sqrt(diag(fisher_info_lambda))
      parameterList$upper_lambda<- c(parameterList$upper_lambda, subset(lam, SR==comparison_mat[i,1])$lambda+1.96*subset(lam, SR==comparison_mat[i,1])$sd_lambda)
      parameterList$lower_lambda <- c(parameterList$lower_lambda, subset(lam, SR==comparison_mat[i,1])$lambda-1.96*subset(lam, SR==comparison_mat[i,1])$sd_lambda)                 
      
      
      ## in keeping with Lotka Volterra notation, we'll use alpha1_2 to indicate effect of
      ## sp 2 on growth of 1.  Following convention, i refers to rows and j to cols in a matrix
      ## so each step of the loop here (for a target sp) corresponds to one row of this matrix:
      
      parameterList$alpha[i,1]<-testcomp2$par[1]
      parameterList$alpha[i,c]<-testcomp3$par[1]
      
      parameterList$lower_alpha[i,1] <- lower_intra[1]
      parameterList$upper_alpha[i,1] <- upper_intra[1]
      
      parameterList$lower_alpha[i,c] <- lower_inter[1]
      parameterList$upper_alpha[i,c] <- upper_inter[1]
      
      ##note that in cases where there is no data for a particular species the alpha estimate for that species ends up as the starting value- we need to be careful of these as they are basically gargbage numbers.  Keeping them in up to now to keep the structure of the data constant, but will set them to NA here:
      
      #identify which species have no data in this fit:
      which(apply(dens, MARGIN=1, FUN=mean)==0)->no_data
      
      ## set their alphas to NA in the matrix:
      
      alpha_matrix[i,no_data]<-NA
      
      ###############################
      ## some diagnostics +  plots ##
      ###############################
      
      ## print an error to the console if any one of the three models failed to converge:
      
      if( testcomp2$convergence + testcomp3$convergence !=0){
        print(paste("at least one model did not converge for", comparison_mat[i,j], sep=" "))
      }
      
      splist<-c(comparison_mat[i,c(1,c)])
      ####### plots##########
      
      ##################################
      ## plot observed vs predicted:
      ##################
      
      par<-c(subset(lam, SR==comparison_mat[i,1])$lambda, testcomp2$par[1], testcomp3$par[1])
      
      #from model 3 code:
      
      lambda<-par[1] #same as model 2
      a_intra<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
      a_inter<-par[3]
      
      pred <- lambda* exp(-a_intra* (dens[1,]-1) - a_inter * dens[c,])
      
      min_max<-c(min(c(log(rate+1), log(pred+1))), max(c(log(rate+1), log(pred+1))))
      
      plot(log(rate+1), log(pred+1), xlim=min_max, ylim=min_max, xlab="log(observed rate)", ylab="log(predicted rate)", main=comparison_mat[i,1] )
      
      abline(a=0, b=1, lwd=2)
      
      
      #####################
      #### plot each fit for intra and inter
      ##########
      alphas<-c(a_intra, a_inter)
      names(alphas)<-splist
      
      plotlist<-as.character(splist)
      for(l in 1:length(plotlist)){
        
        lam_points<-dens[c(1,c),]
        
        ## which columns in the density dataframe have nozero values for species l ?
        
        
        if(l==1){
          cols<-which(dens[2,]==0 & dens[3,]==0)
        }else if(l>1){
          cols<-which(dens[l,]>0)
        }
        x<-dens[l,cols]
        y<-rate[cols]
        
        
        ##add lambdas:
        x<-c(x, rep(0,ncol(dens)))
        y<-c(y, rate)
        
        x_det<-seq(min(x), max(x), by=((max(x)-min(x))/1000))
        
        alpha_temp<-alphas[which(names(alphas)==plotlist[l])]
        y_pred<-lambda* exp(-alpha_temp*x_det)
        
        if(l==1){
          name_list<-paste("intra", plotlist[l], sep=" ")
        }else if(l==2){
          name_list<-paste("inter", plotlist[l], sep=" ")
        }
        
        if(length(unique(x))>1){
          plot(y~x, xlab="density", ylab="rate", main=name_list)
        } else {
          
          plot(x=0, y=0, main=name_list, type='n')
          text(0, 0, "no data")
          
        }
        
        lines(y_pred~x_det, col="red", lwd=2)	
      }
      aux_vec_end<- c(aux_vec_end, comparison_mat[i,1])
    }  # end of background list loop (c)
    
    #######
    
    #vector to store regime for lambda
    
    
  } # end of big for loop
  
  dev.off()
  
  results_same_time<-data.frame(aux_vec_end, parameterList$lambda_est, parameterList$lower_lambda, parameterList$upper_lambda, parameterList$sigma_est, parameterList$convergence_code)
  
  names(results_same_time)<-c("species", "lambda", "lower_error", "upper error", "sigma", "convergence_code")
  alpha_matrix<-parameterList$alpha
  lower_alpha<-parameterList$lower_alpha
  upper_alpha<-parameterList$upper_alpha
  
  write.csv(results_same_time, paste(filepath, "lambda_estimates_", reps,env, ".csv", sep = ""))
  write.csv(alpha_matrix, paste(filepath,"alpha_estimates_row_is_target_",  reps,env,".csv", sep = ""))
  write.csv(lower_alpha, paste(filepath,"alpha_lower_errors_", reps,env,".csv", sep = ""))
  write.csv(upper_alpha, paste(filepath,"alpha_upper_errors_", reps,env,".csv", sep = ""))
  
}

###########################################################################
################ BH model ################################################
###########################################################################

magic_bh_lambda<-function(filepath2,  data2, reps2, env2, comparisons, lam2, ...){
  ## Fitting a negative binomial with interaction regardless of the species identity
  # Only estimating intra not estimating lambda
  
  compmodel3<-function(par, data, dens2, l1){
    lambda<-l1
    alpha_intra<-par[1]
    sigma<-par[2]
    dens<-dens2-1 # removing one individual because it is the focal individual
    
    pred<-lambda/(1+alpha_intra*(dens)) 
    
    llik<-dnorm(data, log(pred), sd=sigma, log=TRUE)
    
    if(length(which(is.na(llik)))>0){
      llik<- -1000000
    }
    
    return(sum(-1*llik))  
  }
  
  #model 3 - all species have different competitive effects, but we only reestimate the interspecific competition
  #The lambda and intraspecific competition comes from the estimation with the last model
  
  compmodel4<-function(par, data1, dens2, l1){
    lambda<-l1
    alpha_inter<-par[1]# Starting with the same value of intra and inter
    sigma<-par[2]
    
    #tot_density<-apply(dens1[1,], MARGIN=2, FUN=sum)
    
    pred<- lambda/ (1+ alpha_inter* dens2)
    #pred<-lambda* exp(-alpha_inter*dens2)
    
    llik<-dnorm(data1, mean=log(pred), sd=sigma, log=TRUE)
    
    if(length(which(is.na(llik)))>0){
      llik<- -1000000
    }
    
    return(sum(-1*llik))  
  }
  
  # 
  #   filepath2 = "./Analyses/magic_riker/"
  #   data2=rep2
  #   reps2=1
  #   env2="N"
  #   comparisons = comparison_mat
  
  filepath<-filepath2 #where to keep the files
  lam<-lam2 # matrix with the initial estimates of lambda for the two regimes used
  data<-data2 # matrix with the data, it has to have 5 columns: Regime, background, # focal ind, # comp ind, growth rate
  reps <-reps2 #replicate
  env<-env2 #environment
  comparison_mat<-comparisons #matrix with the comparisons to use
  # example
  # comparison_mat<-matrix(nrow=3, ncol=3)
  # comparison_mat[1,]<-c(1,4,5)
  # comparison_mat[2,]<-c(4,1,NA)
  # comparison_mat[3,]<-c(5,1,NA)
  
  
  parameterList<-list(alpha_matrix<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # where we keep the alpha information
                      lower_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))), # matrix that keeps the alpha lower CI
                      upper_alpha<- matrix(0, nrow=length(comparison_mat[,1]), ncol=length(comparison_mat[1,]), dimnames = list(comparison_mat[,1],c("intra","inter1","inter2"))),# matrix that keeps the alpha upper CI
                      lambda_est<-NULL, #store the lambda estimates
                      lower_lambda <-NULL, # store the lower lambda CI
                      upper_lambda <-NULL, # store the upper lambda CI
                      sigma_est<-NULL, #store the estimates for the sigma (error)
                      convergence_code<-NULL # convergence code to be sure that everyhing is ok
  )
  names(parameterList)<-c("alpha", "lower_alpha", "upper_alpha", "lambda", "lower_lambda", "upper_lambda", "sigma", "convergence_code")
  aux_vec_end<-NULL
  # we have to check which data points are NA
  toRemove<-which(is.na(data$Growthrate))
  if(length(toRemove)!=0){
    data<-data[-toRemove,]
  }
  
  
  #Since the regimes are named as numbers, then we need to put them as chars
  data$Regime<-as.character(data$Regime)
  data$background<-as.character(data$background)
  #print(data)
  
  pdf(file=paste(filepath,"plot_",  reps,env,".pdf", sep = ""), width=11, height=8)
  
  for(i in 1:length(comparison_mat[,1])){
    # subset out the rows that are needed from the competition matrix
    comp_points<-subset(data, data$Regime==as.character(comparison_mat[i,1]), drop=TRUE)
    
    ## now need to build up a vector of nonzero mite production
    ## and corresponding density vectors (held in a matrix) for background species 
    ## to use in model fitting
    
    ## start with the lambda rate- will add on to this for each background:
    
    rate<-NULL
    
    ##build density matrix (each row will be density of a species)
    dens<-matrix(0, nrow=length(comparison_mat[i,]), ncol=(nrow(comp_points)))
    row.names(dens)<-comparison_mat[i,]
    
    #print("here2")
    #print(i)
    ##for each background species the target competes against:
    
    # if there isn't enough comparisons (because in replicate 2 there is one less comparison to do)
    if(is.na(comparison_mat[i,3])){
      background_list<-comparison_mat[i,1:2]
    }else{
      background_list<-comparison_mat[i,]
    }
    
    ## use this counter to keep track of which column is next to have data added to
    ## set it to begin after the lambda points:
    # print(comp_points)
    
    start <- 1
    #print( background_list)
    
    for(j in 1:length(background_list)){
      ## LOOP creates a rate vector and a corresponding dens ("density") matrix with 
      ## background species as rows, and columns corresponding to the seed production
      ## vector.
      ## beginning columns correspond to lambda
      
      #take just the rows pertaining to a specific background sp:
      bg_points<- subset(comp_points, comp_points$background==background_list[j])
      
      ## which row of the density matrix corresponds?
      #The inter
      rownum <- which(row.names(dens)==background_list[j]) #here is the row to put the competitor individuals (of other species)
      #The intra
      rownumIntra<-comparison_mat[i,1] # here is row where to keep the number of individuals of the focal species
      #column to end with:
      end<-start + nrow(bg_points)-1
      #print(start)
      #print(end)
      
      
      ## drop in density values into matrix:
      
      dens[rownum, start:end]<-bg_points$Comp
      dens[as.character(rownumIntra), start:end]<-bg_points$Focal #always do this one last, because in intraspecific competition it will add the number of individuals
      
      ## add seed numbers into the rate vector
      rate<-c(rate, bg_points$Growthrate)
      
      ##increase start counter so that next species is offset from this one
      start<-end+1
    }
    
    #print(dens)
    
    #print(paste("Before starting models, its rep ",reps, comparison_mat[i,], background_list[j], sep=" "))
    ## we'll be working with log rate (for giving a lognormal error structure):
    
    log_rate<-log(rate+1)
    
    
    ### 
    ### Skipping estimating only lambda, lambda will be estimated directly from the data set
    ###
    
    #############################
    ## model 2, one alpha      ##
    #############################
    
    ## pars here are alpha and sigma- alpha starts at 0.1
    
    par2<-c(0.1, 0.1)
    
    aux_num<-length(dens[,1])
    if(aux_num==3){
      #subsetting to get only the intraspecific density
      aux_dens<-dens[1,which(dens[2,]==0 & dens[3,]==0)]
      aux_lograte<-log_rate[which(dens[2,]==0 & dens[3,]==0)]
    }else if(aux_num==2){
      aux_dens<-dens[1,which(dens[2,]==0)]
      aux_lograte<-log_rate[which(dens[2,]==0)]
    }
    
    ### estimating intraspecific interactions
    for(k in 1:100){
      testcomp2<-optim(par2, compmodel3, dens=aux_dens, data=aux_lograte, l1=subset(lam, SR==comparison_mat[i,1])$lambda, control = list(maxit=10000000, gamma=5), hessian = TRUE, method="L-BFGS-B", lower=c(0,0.0001), upper=c(10,10))
      #testcomp2<-optim(par2, compmodel3, dens=aux_dens, data=aux_lograte, l1=subset(lam, SR==comparison_mat[i,1])$lambda, control = list(maxit=10000000, gamma=5), hessian = TRUE)
      par2<-testcomp2$par
      if(testcomp2$convergence==0){
        print(paste(comparison_mat[i,j],  "model 2 converged on rep", k, sep=" "))
        break
      }
    }
    
    for(c in 2:length(background_list)){
      
      par3<-c(0.1, 0.1) #starting from 0.1
      #here we only care about the interspecific competition, so on this row we choose only the data from the focal species but for the interspecific ocompteition
      dens_inter<-dens[c,which(dens[c,]!=0) ]
      data_c<-log_rate[which(dens[c,]!=0)] # since the columns match we can subset by the columns
      
      
      ## estimating interspecific
      for(k in 1:100){
        testcomp3<-optim(par3,fn=compmodel4,data1=data_c,dens2=dens_inter, l1=subset(lam, SR==comparison_mat[i,1])$lambda,  control = list(maxit=100000, gamma=5), hessian=T, method="L-BFGS-B", lower=c(0,0.0001), upper=c(10,10))
        #testcomp3<-optim(par3,fn=compmodel4,data1=data_c,dens2=dens_inter, l1=subset(lam, SR==comparison_mat[i,1])$lambda,  control = list(maxit=100000, gamma=5), hessian=T)
        par3<-testcomp3$par
        if(testcomp3$convergence==0){
          print(paste(comparison_mat[i,j], "model 3 converged on rep", k, sep=" "))
          break
        }
      }
      
      
      ##################################
      ### save estimates from model 2 and 3 ##
      ##################################
      
      
      parameterList$lambda_est<-c(parameterList$lambda_est, subset(lam, SR==comparison_mat[i,1])$lambda)
      
      parameterList$sigma_est<-c(parameterList$sigma_est, subset(lam, SR==comparison_mat[i,1])$sd_lambda)#Taking the sigma from the intraspecific + lambda estimation
      
      parameterList$convergence_code<-c(parameterList$convergence_code, testcomp3$convergence)
      
      # inter
      tryCatch(fisher_info<-solve(testcomp3$hessian), error=function(cond){fisher_info<-matrix(NA, ncol=1, nrow=1)
      print("Problem with hessian model2")
      message(cond)} )
      #lambda and intra
      tryCatch(fisher_info_intra<-solve(testcomp2$hessian), error=function(cond){print("Problem with hessian model3")
        message(cond)
        fisher_info_intra<-matrix(NA, ncol=2, nrow=2) } )
      #tryCatch(fisher_info_lambda<-solve(testcomp1$hessian), error=function(cond){print("Problem with hessian model1")
      #message(cond)
      #fisher_info_lambda<-matrix(NA, ncol=2, nrow=2) } )
      
      #print("three")
      prop_sigma_inter<-sqrt(diag(fisher_info))
      upper_inter<-testcomp3$par[1]+1.96*prop_sigma_inter[1]
      lower_inter<-testcomp3$par[1]-1.96*prop_sigma_inter[1]
      
      prop_sigma_intra<-sqrt(diag(fisher_info_intra))
      upper_intra<-testcomp2$par[1]+1.96*prop_sigma_intra[1]
      lower_intra<-testcomp2$par[1]-1.96*prop_sigma_intra[1]
      
      #prop_sigma_lambda<-sqrt(diag(fisher_info_lambda))
      parameterList$upper_lambda<- c(parameterList$upper_lambda, subset(lam, SR==comparison_mat[i,1])$lambda+1.96*subset(lam, SR==comparison_mat[i,1])$sd_lambda)
      parameterList$lower_lambda <- c(parameterList$lower_lambda, subset(lam, SR==comparison_mat[i,1])$lambda-1.96*subset(lam, SR==comparison_mat[i,1])$sd_lambda)                 
      
      
      ## in keeping with Lotka Volterra notation, we'll use alpha1_2 to indicate effect of
      ## sp 2 on growth of 1.  Following convention, i refers to rows and j to cols in a matrix
      ## so each step of the loop here (for a target sp) corresponds to one row of this matrix:
      
      parameterList$alpha[i,1]<-testcomp2$par[1]
      parameterList$alpha[i,c]<-testcomp3$par[1]
      
      parameterList$lower_alpha[i,1] <- lower_intra[1]
      parameterList$upper_alpha[i,1] <- upper_intra[1]
      
      parameterList$lower_alpha[i,c] <- lower_inter[1]
      parameterList$upper_alpha[i,c] <- upper_inter[1]
      
      ##note that in cases where there is no data for a particular species the alpha estimate for that species ends up as the starting value- we need to be careful of these as they are basically gargbage numbers.  Keeping them in up to now to keep the structure of the data constant, but will set them to NA here:
      
      #identify which species have no data in this fit:
      which(apply(dens, MARGIN=1, FUN=mean)==0)->no_data
      
      ## set their alphas to NA in the matrix:
      
      alpha_matrix[i,no_data]<-NA
      
      ###############################
      ## some diagnostics +  plots ##
      ###############################
      
      ## print an error to the console if any one of the three models failed to converge:
      
      if( testcomp2$convergence + testcomp3$convergence !=0){
        print(paste("at least one model did not converge for", comparison_mat[i,j], sep=" "))
      }
      
      splist<-c(comparison_mat[i,c(1,c)])
      ####### plots##########
      
      ##################################
      ## plot observed vs predicted:
      ##################
      
      par<-c(subset(lam, SR==comparison_mat[i,1])$lambda, testcomp2$par[1], testcomp3$par[1])
      
      #from model 3 code:
      
      lambda<-par[1] #same as model 2
      a_intra<-par[2]	## new parameters- use alpha estimate from model 2 as start value for fitting
      a_inter<-par[3]
      
      pred <- lambda/(1+ a_intra* (dens[1,]-1) + a_inter * dens[c,])
      
      min_max<-c(min(c(log(rate+1), log(pred+1)), na.rm=TRUE), max(c(log(rate+1), log(pred+1)), na.rm=TRUE))
      
      plot(log(rate), log(pred), xlim=min_max, ylim=min_max, xlab="log(observed rate)", ylab="log(predicted rate)", main=comparison_mat[i,1] )
      
      abline(a=0, b=1, lwd=2)
      
      
      #####################
      #### plot each fit for intra and inter
      ##########
      alphas<-c(a_intra, a_inter)
      names(alphas)<-splist
      
      plotlist<-as.character(splist)
      for(l in 1:length(plotlist)){
        
        lam_points<-dens[c(1,c),]
        
        ## which columns in the density dataframe have nozero values for species l ?
        
        
        if(l==1){
          cols<-which(dens[2,]==0 & dens[3,]==0)
        }else if(l>1){
          cols<-which(dens[l,]>0)
        }
        x<-dens[l,cols]
        y<-rate[cols]
        
        
        ##add lambdas:
        x<-c(x, rep(0,ncol(dens)))
        y<-c(y, rate)
        
        x_det<-seq(min(x), max(x), by=((max(x)-min(x))/1000))
        
        alpha_temp<-alphas[which(names(alphas)==plotlist[l])]
        y_pred<-lambda/(1+alpha_temp*x_det)
        
        if(l==1){
          name_list<-paste("intra", plotlist[l], sep=" ")
        }else if(l==2){
          name_list<-paste("inter", plotlist[l], sep=" ")
        }
        
        if(length(unique(x))>1){
          plot(y~x, xlab="density", ylab="rate", main=name_list)
        } else {
          
          plot(x=0, y=0, main=name_list, type='n')
          text(0, 0, "no data")
          
        }
        
        lines(y_pred~x_det, col="red", lwd=2)	
      }
      aux_vec_end<- c(aux_vec_end, comparison_mat[i,1])
    }  # end of background list loop (c)
    
    #######
    
    #vector to store regime for lambda
    
    
  } # end of big for loop
  
  dev.off()
  
  results_same_time<-data.frame(aux_vec_end, parameterList$lambda_est, parameterList$lower_lambda, parameterList$upper_lambda, parameterList$sigma_est, parameterList$convergence_code)
  
  names(results_same_time)<-c("species", "lambda", "lower_error", "upper error", "sigma", "convergence_code")
  alpha_matrix<-parameterList$alpha
  lower_alpha<-parameterList$lower_alpha
  upper_alpha<-parameterList$upper_alpha
  
  write.csv(results_same_time, paste(filepath, "lambda_estimates_", reps,env,"BH_", ".csv", sep = ""))
  write.csv(alpha_matrix, paste(filepath,"alpha_estimates_row_is_target_",  reps,env,"BH_",".csv", sep = ""))
  write.csv(lower_alpha, paste(filepath,"alpha_lower_errors_", reps,env,"BH_",".csv", sep = ""))
  write.csv(upper_alpha, paste(filepath,"alpha_upper_errors_", reps,env,"BH_",".csv", sep = ""))
  
}