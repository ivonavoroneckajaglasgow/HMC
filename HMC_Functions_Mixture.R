library(coda)
HMC_helper <- function(x, y, U, ddall, epsilon, L, current_q, m=1, which) 
{
  q <- current_q
  p <- rnorm(length(q),0,1)
  current_p <- p
  
  # Make a half step for momentum at the beginning
  params_q          <- params
  params_q[[which]] <- q
  p                 <- p - epsilon * ddall(x,y,params=params_q,which)/2
  
  # Alterate full steps for position and momentum
  
  for(i in 1:L)
  {
    # Make full step for position
    q <- q+epsilon*p/m
    
    if( i!= L){ 
      params_q[[which]] <- q
      p                 <- p - epsilon * ddall(x,y,params=params_q,which)/2
    }
  }
  
  
  #Negate momentum at end of trajectory to make proposal symmetric
  
  p = -p
  
  #Evaluate potential and kinetic energies at start and end of trajectory
  
  params_current          <- params
  params_current[[which]] <- current_q
  current_U               <- U(x,y,params_current)
  
  
  params_proposed          <- params
  params_proposed[[which]] <- q
  proposed_U               <- U(x,y,params_proposed)
  
  
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

HMC <- function(x, y, N, burnin, U, ddall, epsilon, L, current_q, which, do.plot=TRUE, actual = NA){
  
  HMC_result        <- list()
  
  result      <- matrix(nrow=N,ncol=length(current_q))
  result[1,]  <- current_q
  accept      <- numeric(N)
  ESS         <-c()
  
  pb <- txtProgressBar(min=0, max=N, style=3)
  for(j in 2:N){
    setTxtProgressBar(pb, j)
    
    result[j,]<-HMC_helper(x, y, U,ddall,epsilon,L,current_q = result[j-1,],which=which)
    
    if(all(result[j-1,]!=result[j,])) accept[[j]] <- 1
  }
  
  for(i in 1:ncol(result)){
    ESS <- c(ESS,effectiveSize(as.mcmc(result[burnin:length(result[,i]),i])))
  }
  
  HMC_result$ESS    <- ESS
  HMC_result$output <- result
  HMC_result$accept <- accept
  
  ##plotting### 
  main_labels <- list(c("gamma0","gamma1"),
                      c("beta01","beta11"),
                      "sigma1",
                      c("beta02","beta12"),
                      "sigma2")
  
  my_labels   <- main_labels[[which]] 
  
  if(do.plot==TRUE){
    par(mfrow=c(1,2))
    for(i in 1:ncol(result)){
      hist(result[burnin:nrow(result),i],xlab=my_labels[i],main=my_labels[i])
      if (!is.na(actual[i])) abline(v=actual[i],col=2)
      plot(result[burnin:nrow(result),i],type="l",ylab=my_labels[i],main=my_labels[i])
      if (!is.na(actual[i])) abline(h=actual[i],col=2)
    }
  }
  
  return(HMC_result)
}

