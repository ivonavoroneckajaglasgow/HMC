###HMC Functions###

HMC_helper <- function(what, x, y, U, grad_U, epsilon, L, current_q) 
{
  q <- current_q
  p <- rnorm(length(q),0,1)
  current_p <- p
  
  # Make a half step for momentum at the beginning
  if(what=="mu1"){
    p <- p - epsilon * ddmu(which=1,y,prob,mu1=q,sigma1,mu2,sigma2)/2
  }
  if(what=="mu2"){
    p <- p - epsilon * ddmu(which=2,y,prob,mu1,sigma1,mu2=q,sigma2)/2   
  }
  if(what=="sigma1"){
    p <- p - epsilon * ddsigma(which=1,y,prob,mu1,sigma1=q,mu2,sigma2)/2
  }
  if(what=="sigma2"){
    p <- p - epsilon * ddsigma(which=2,y,prob,mu1,sigma1,mu2,sigma2=q)/2
  }
  if(what=="gamma"){
    prob_q <- p_calculator(gamma=q,x)
    p      <- p - epsilon * ddgamma(x,y,prob=prob_q,mu1,sigma1,mu2,sigma2)/2
  }
  
  # Alterate full steps for position and momentum
  
  for(i in 1:L)
  {
    #print(q)
    
    # Make full step for position
    q <- q+epsilon*p/m
    
    # Make a full step for the momentum, except at the end of trajectory
    if(what=="mu1"){
      if( i!= L) p <- p - epsilon * ddmu(which=1,y,prob,mu1=q,sigma1,mu2,sigma2) /2
    }
    if(what=="mu2"){
      if( i!= L) p <- p - epsilon * ddmu(which=2,y,prob,mu1,sigma1,mu2=q,sigma2) /2
    }
    if(what=="sigma1"){
      if( i!= L) p <- p - epsilon * ddsigma(which=1,y,prob,mu1,sigma1=q,mu2,sigma2)/2
    }
    if(what=="sigma2"){
      if( i!= L) p <- p - epsilon * ddsigma(which=2,y,prob,mu1,sigma1,mu2,sigma2=q)/2
    }
    if(what=="gamma"){
      if( i!= L){
      prob_q <- p_calculator(gamma=q,x)
      p      <- p - epsilon * ddgamma(x,y,prob=prob_q,mu1,sigma1,mu2,sigma2)/2
      }
    }
  }
  
  
  #Negate momentum at end of trajectory to make proposal symmetric
  
  p = -p
  
  #Evaluate potential and kinetic energies at start and end of trajectory
  if(what=="mu1"){
    current_U  <- U(y,prob,mu1=current_q,sigma1,mu2,sigma2)
    proposed_U <- U(y,prob,mu1=q,sigma1,mu2,sigma2)
  }
  if(what=="mu2"){
    current_U  <- U(y,prob,mu1,sigma1,mu2=current_q,sigma2)
    proposed_U <- U(y,prob,mu1,sigma1,mu2=q,sigma2)
  }
  if(what=="sigma1"){
    current_U  <- U(y,prob,mu1,sigma1=current_q,mu2,sigma2)
    proposed_U <- U(y,prob,mu1,sigma1=q,mu2,sigma2)
  }
  if(what=="sigma2"){
    current_U  <- U(y,prob,mu1,sigma1,mu2,sigma2=current_q)
    proposed_U <- U(y,prob,mu1,sigma1,mu2,sigma2=q)
  }
  if(what=="gamma"){
    prob_current<- p_calculator(gamma=current_q,x)
    current_U   <- U(y,prob_current,mu1,sigma1,mu2,sigma2)
    prob_q      <- p_calculator(gamma=q,x)
    proposed_U  <- U(y,prob_q,mu1,sigma1,mu2,sigma2)
  }
  
  proposed_K <- sum(p^2)/2
  current_K  <- sum(current_p^2)/2
  #Accept or reject
  
  r <- current_U-proposed_U+current_K-proposed_K
  
  if(log(runif(1))< r)
  {
    return(q)  #accept
  }else{
    return(current_q) #reject
  }
  
}

HMC <- function(what, x, y, N, burnin, U, grad_U, epsilon, L, current_q, do.plot=TRUE, actual = NA){
  
  HMC_result        <- list()
  
  result      <- matrix(nrow=N,ncol=length(current_q))
  result[1,]  <- current_q
  accept      <- numeric(N)
    
    pb <- txtProgressBar(min=0, max=N, style=3)
    for(j in 2:N){
      setTxtProgressBar(pb, j)
      result[j,]<-HMC_helper(what, x, y, U,grad_U,epsilon,L,current_q = result[j-1,])
      if(all(result[j-1,]==result[j,])) accept[[j]] <- 1
    }
    
    if(do.plot==TRUE){
      par(mfrow=c(1,2))
      for(i in 1:ncol(result)){
      hist(result[burnin:nrow(result),i],xlab=what,main=what)
      if (!is.na(actual)) abline(v=actual[i],col=2)
      plot(result[burnin:nrow(result),i],type="l",ylab=what,main=what)
      if (!is.na(actual)) abline(h=actual[i],col=2)
      }
    }
    
    ESS<-c()
    
    for(i in 1:ncol(result)){
      ESS <- c(ESS,effectiveSize(as.mcmc(result[burnin:length(result[,i])])))
    }
    
  HMC_result$ESS    <- ESS
  HMC_result$output <- result
  HMC_result$accept <- accept
  
  return(HMC_result)
}

