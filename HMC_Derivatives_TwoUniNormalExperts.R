###Derivatives###

dmix <- function(y,prob,mu1,sigma1,mu2,sigma2){
  prob*dnorm(y,mu1,sqrt(sigma1))+(1-prob)*dnorm(y,mu2,sqrt(sigma2))
}

U <- function(x,y,params){
  mu1    <- params$mu1
  mu2    <- params$mu2
  sigma1 <- params$sigma1
  sigma2 <- params$sigma2
  gamma  <- params$gamma
  prob   <- p_calculator(gamma,x)
  sum(-log(dmix(y,prob,mu1,sigma1,mu2,sigma2)))
}

l <- function(y,prob,mu1,sigma1,mu2,sigma2){
  sum(log(dmix(y,prob,mu1,sigma1,mu2,sigma2)))
}

ddall <- function(x,y,params,which){
  
mu1    <- params$mu1
mu2    <- params$mu2
sigma1 <- params$sigma1
sigma2 <- params$sigma2
gamma  <- params$gamma
prob   <- p_calculator(gamma,x)
    
result<-c(ddmu(y,prob,mu1,sigma1,mu2,sigma2),
          ddsigma(y,prob,mu1,sigma1,mu2,sigma2),
          ddgamma(x,y,prob,mu1,sigma1,mu2,sigma2))

return(result[which])
}

ddgamma <- function(x,y,prob,mu1,sigma1,mu2,sigma2){
  colSums(-1/dmix(y,prob,mu1,sigma1,mu2,sigma2)*prob*(1-prob)*(dnorm(y,mu1,sqrt(sigma1))-dnorm(y,mu2,sqrt(sigma2)))*cbind(1,x))
}

ddmu <- function(y,prob,mu1,sigma1,mu2,sigma2){ 
  first <- -1/dmix(y,prob,mu1,sigma1,mu2,sigma2)
  mu1   <- sum(first*prob*dnorm(y,mu1,sqrt(sigma1))*(y-mu1)/sigma1)
  mu2   <- sum(first*(1-prob)*dnorm(y,mu2,sqrt(sigma2))*(y-mu2)/sigma2)
  return(c(mu1,mu2))
}

ddsigma <- function(y,prob,mu1,sigma1,mu2,sigma2){
  
  first <- -1/dmix(y,prob,mu1,sigma1,mu2,sigma2)
  
  second_sigma1 <- -dnorm(y,mu1,sqrt(sigma1))*prob
  result_sigma1 <- sum(first*second_sigma1*1/(2*sigma1)*(1-(y-mu1)^2/sigma1))
  

  second_sigma2 <- -dnorm(y,mu2,sqrt(sigma2))*(1-prob)
  result_sigma2 <- sum(first*second_sigma2*1/(2*sigma2)*(1-(y-mu2)^2/sigma2))
  
    
  return(c(result_sigma1, result_sigma2))
}

