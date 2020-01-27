library(coda)
HMC_helper <- function(x, y, U, ddall, epsilon, L, current_q, m=1, which) 
{
  q <- current_q
  p <- rnorm(length(q),0,1)
  current_p <- p
  
  # Make a half step for momentum at the beginning
  params_q          <- assign_paramater_values(params,which,v=q)
  p                 <- p - epsilon * ddall(x,y,params=params_q,which)/2
  
  # Alterate full steps for position and momentum
  
 allq<-c()
 allp<-c()
  
  for(i in 1:L)
  {
    # Make full step for position
    q <- q+epsilon*p/m
    allq<-rbind(allq,c(i,q))
    
    if( i!= L){ 
      params_q          <- assign_paramater_values(params_q,which,v=q)
      p                 <- p - epsilon * ddall(x,y,params=params_q,which)/2
      allp <-rbind(allp,c(i,p))
    }
  }
  
  
  #Negate momentum at end of trajectory to make proposal symmetric
  
  p = -p
  
  #Evaluate potential and kinetic energies at start and end of trajectory
  
  params_current          <- assign_paramater_values(params,which,v=current_q)
  current_U               <- U(x,y,params_current)
  
  
  params_proposed          <- assign_paramater_values(params,which,v=q)
  proposed_U               <- U(x,y,params_proposed)
  
  
  proposed_K <- sum(p^2)/2
  current_K  <- sum(current_p^2)/2
  #Accept or reject
  
  r <- current_U-proposed_U+current_K-proposed_K

  records <- list()
  
  records$q         <- q
  records$current_q <- current_q
  records$accept    <- ifelse(log(runif(1))< r,1,0)
  records$leapfrog  <- allq
  
  if(records$accept==1){
    records$final     <- records$q
    }else{
      records$final     <- records$current_q
    }
  
  # if(log(runif(1))< r)
  # {
  #   return(q)  #accept
  # }else{
  #   return(current_q) #reject
  # }
  
return(records)
  
}

HMC <- function(x, y, N, burnin, U, ddall, epsilon, L, current_q, which=c("gamma","beta1","sigma1","beta2","sigma2"), do.plot=TRUE, actual = NA){
  
  HMC_result        <- list()
  
  HMC_result$output <- matrix(nrow=N,ncol=length(current_q))
  HMC_result$all_q  <- matrix(nrow=N,ncol=length(current_q))
  
  HMC_result$output[1,]<- current_q
  HMC_result$all_q[1,] <- current_q
  
  HMC_result$accept <- c()
  HMC_result$ESS    <- c()
  HMC_result$leapfrog<- c()
  
  pb <- txtProgressBar(min=0, max=N, style=3)
  
  set.seed(12)
  for(j in 2:N){
    setTxtProgressBar(pb, j)
    
    a <- HMC_helper(x, y, U,ddall,epsilon,L,current_q =  HMC_result$output[j-1,],which=which)
    
    HMC_result$output[j,]<- a$final
    HMC_result$all_q[j,] <- a$q
    HMC_result$accept[j] <- a$accept
    HMC_result$leapfrog  <- rbind(HMC_result$leapfrog,a$leapfrog)
  }
  
  for(i in 1:ncol( HMC_result$output)){
    HMC_result$ESS <- c(HMC_result$ESS,effectiveSize(as.mcmc(HMC_result$output[burnin:length(HMC_result$output[,i]),i])))
  }
  
  ##plotting### 
  
  my_labels <- translate_which(which,text=TRUE)

  if(do.plot==TRUE){
    par(mfrow=c(2,2))
    for(i in 1:ncol(HMC_result$output)){
      hist(HMC_result$output[burnin:nrow(HMC_result$output),i],xlab=my_labels[i],main=my_labels[i])
      if (!is.na(actual[i])) abline(v=actual[i],col=2)
      plot(HMC_result$output[burnin:nrow(HMC_result$output),i],type="l",ylab=my_labels[i],main=my_labels[i])
      if (!is.na(actual[i])) abline(h=actual[i],col=2)
    }
  }
  
  return(HMC_result)
}

