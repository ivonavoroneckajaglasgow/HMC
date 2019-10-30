source("helper_functions.R")
load("data1AR.RData")
x  <-loc$x
y  <-loc$y
plot(x,y)
mu1    <- 3
mu2    <- 7
sigma1 <- 2
sigma2 <- 2

gamma  <- c(1,5)
prob   <- p_calculator(gamma,x)


dmix <- function(y,prob,mu1,sigma1,mu2,sigma2){
  prob*dnorm(y,mu1,sqrt(sigma1))+(1-prob)*dnorm(y,mu2,sqrt(sigma2))
}

U <- function(y,prob,mu1,sigma1,mu2,sigma2){
  -log(dmix(y,prob,mu1,sigma1,mu2,sigma2))
}

l <- function(y,prob,mu1,sigma1,mu2,sigma2){
  sum(log(dmix(y,prob,mu1,sigma1,mu2,sigma2)))
}

ddgamma <- function(x,y,prob,mu1,sigma1,mu2,sigma2){
-1/dmix(x,prob,mu1,sigma1,mu2,sigma2)*prob*(1-prob)*x*(dnorm(y,mu1,sqrt(sigma1))-dnorm(y,mu2,sqrt(sigma2)))
}

ddmu <- function(which,y,prob,mu1,sigma1,mu2,sigma2){ 
 first <- -1/dmix(y,prob,mu1,sigma1,mu2,sigma2)
 if(which==1) return( first*prob*dnorm(y,mu1,sqrt(sigma1))*(y-mu1)/sigma1 )
 if(which==2) return( first*(1-prob)*dnorm(y,mu2,sqrt(sigma2))*(y-mu2)/sigma2)
}

ddsigma <- function(which,y,prob,mu1,sigma1,mu2,sigma2){
 first <- -1/dmix(y,prob,mu1,sigma1,mu2,sigma2)
 if(which==1){
   second <- -dnorm(y,mu1,sqrt(sigma1))*prob
   result <- first*second*1/(2*sigma1)*(1-(y-mu1)^2/sigma1)
   return(result)
 }
 if(which==2){
   second <- -dnorm(y,mu2,sqrt(sigma2))*(1-prob)
   result <- first*second*1/(2*sigma2)*(1-(y-mu2)^2/sigma2)
   return(result)
 }
}


U(x,prob,mu1,sigma1,mu2,sigma2)
l(x,prob,mu1,sigma1,mu2,sigma2)
ddgamma(x,y,prob,mu1,sigma1,mu2,sigma2)
ddmu(which=1,x,prob,mu1,sigma1,mu2,sigma2)
ddmu(which=2,x,prob,mu1,sigma1,mu2,sigma2)
ddsigma(which=1,y,prob,mu1,sigma1,mu2,sigma2)
ddsigma(which=2,y,prob,mu1,sigma1,mu2,sigma2)

################################################################

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
