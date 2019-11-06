###Derivatives###

dmix <- function(y,prob,mu1,sigma1,mu2,sigma2){
  prob*dnorm(y,mu1,sqrt(sigma1))+(1-prob)*dnorm(y,mu2,sqrt(sigma2))
}

U <- function(y,prob,mu1,sigma1,mu2,sigma2){
  sum(-log(dmix(y,prob,mu1,sigma1,mu2,sigma2)))
}

l <- function(y,prob,mu1,sigma1,mu2,sigma2){
  sum(log(dmix(y,prob,mu1,sigma1,mu2,sigma2)))
}

ddgamma <- function(x,y,prob,mu1,sigma1,mu2,sigma2){
  colSums(-1/dmix(y,prob,mu1,sigma1,mu2,sigma2)*prob*(1-prob)*(dnorm(y,mu1,sqrt(sigma1))-dnorm(y,mu2,sqrt(sigma2)))*cbind(1,x))
}

ddmu <- function(which,y,prob,mu1,sigma1,mu2,sigma2){ 
  first <- -1/dmix(y,prob,mu1,sigma1,mu2,sigma2)
  if(which==1) return( sum(first*prob*dnorm(y,mu1,sqrt(sigma1))*(y-mu1)/sigma1) )
  if(which==2) return( sum(first*(1-prob)*dnorm(y,mu2,sqrt(sigma2))*(y-mu2)/sigma2))
}

ddsigma <- function(which,y,prob,mu1,sigma1,mu2,sigma2){
  first <- -1/dmix(y,prob,mu1,sigma1,mu2,sigma2)
  if(which==1){
    second <- -dnorm(y,mu1,sqrt(sigma1))*prob
    result <- first*second*1/(2*sigma1)*(1-(y-mu1)^2/sigma1)
  }
  if(which==2){
    second <- -dnorm(y,mu2,sqrt(sigma2))*(1-prob)
    result <- first*second*1/(2*sigma2)*(1-(y-mu2)^2/sigma2)
  }
  return(sum(result))
}