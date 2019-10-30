mu1    <- 2
mu2    <- 5
sigma1 <- 0.1
sigma2 <- 2
prob   <- 0.5
m      <- 1

dmix <- function(x,prob,mu1,sigma1,mu2,sigma2){
  prob*dnorm(x,mu1,sqrt(sigma1))+(1-prob)*dnorm(x,mu2,sqrt(sigma2))
}

U <- function(x,prob,mu1,sigma1,mu2,sigma2){
  -log(dmix(x,prob,mu1,sigma1,mu2,sigma2))
}

grad_U <- function(x,prob,mu1,sigma1,mu2,sigma2){
    dens1 <- prob*dnorm(x, mu1, sqrt(sigma1))
    dens2 <- (1-prob)*dnorm(x, mu2, sqrt(sigma2))
    dens1/(dens1+dens2) * (x-mu1)/sigma1 + dens2/(dens1+dens2) * (x-mu2)/sigma2
}




#U  <- function(x,prob,mu1,sigma1,mu2,sigma2){
#  prob*dnorm(x,mu1,sqrt(sigma1))+(1-prob)*dnorm(x,mu2,sqrt(sigma2))
#}

#ddnorm <-  function(x,mu,sigma){
#  -dnorm(x,mu,sqrt(sigma))*(x-mu)/sigma
#}

#grad_U <-  function(x,prob,mu1,sigma1,mu2,sigma2){
#  prob*ddnorm(x,mu1,sigma1)+(1-prob)*ddnorm(x,mu2,sigma2)
#}


HMC <- function(U, grad_U, epsilon, L, current_q) 
{
  q <- current_q
  p <- rnorm(length(q),0,1)
  current_p <- p
  
  # Make a half step for momentum at the beginning
  
  p <- p - epsilon * grad_U(x=q,prob,mu1,sigma1,mu2,sigma2)/2
  
  # Alterate full steps for position and momentum
  
  for(i in 1:L)
  {
    # Make full step for position
    q <- q+epsilon*p/m
    
    # Make a full step for the momentum, except at the end of trajectory
    if( i!= L) p <- p - epsilon * grad_U(x=q,prob,mu1,sigma1,mu2,sigma2) /2
  }
    #Negate momentum at end of trajectory to make proposal symmetric
    
    p = -p
    
    #Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U  <- U(x=current_q,prob,mu1,sigma1,mu2,sigma2)
    current_K  <- sum(current_p^2)/2
    proposed_U <- U(x=q,prob,mu1,sigma1,mu2,sigma2)
    proposed_K <- sum(p^2)/2
    
    #Accept or reject
    
    r <- current_U-proposed_U+current_K-proposed_K
    #r <- (current_U+current_K)/(proposed_U+proposed_K)
    
    if(log(runif(1))< r)
    {
      return(q)  #accept
    }else{
      return(current_q) #reject
    }
    
}


HMC(U,grad_U,epsilon = 1,L=10,current_q = 3)

N       <- 5000
my.x    <- c()
my.x[1] <- 3

for(j in 2:N){
my.x[j]<-HMC(U,grad_U,epsilon = 0.5,L=10,current_q = my.x[j-1])
}

hist(my.x)

