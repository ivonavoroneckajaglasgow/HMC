load("data1AR.RData")
x<-loc$x
y<-loc$y

###parameters estimated by sampler 1
params<-list(gamma=-c(-11.95,2.42),
             beta1=c(-1.47,2.37),
             sigma1=log(3.158193),      
             beta2=c(23.59,-2.75),
             sigma2=log(4.0552)
)

p_calculator = function(gamma,x){
  
  a_gamma       <- exp(gamma[1]+gamma[2]*x)
  p             <- a_gamma/(1+a_gamma)
  return(p)
}

transform_sigma <- function(params){
  params_new <- params
  params_new$sigma1 <- exp(params$sigma1)
  params_new$sigma2 <- exp(params$sigma2)
  params_new
}

getmu         <- function(x,params){
  mu1 <- params$beta1[1]+params$beta1[2]*x 
  mu2 <- params$beta2[1]+params$beta2[2]*x
  return(cbind(mu1,mu2))
}

dmix          <- function(x,y,params,transform_sigma=FALSE){
  if(transform_sigma==TRUE) params <- transform_sigma(params)
  
  prob   <- p_calculator(params$gamma,x)
  mu     <- getmu(x, params)
  return(prob*dnorm(y,mu[,1],sqrt(params$sigma1))+(1-prob)*dnorm(y,mu[,2],sqrt(params$sigma2)))
}

dmix_log <- function(x,y,params,transform_sigma=FALSE){
  log(dmix(x,y,params))
}

ddbeta_helper_log <- function(x,y,params,mu){
  
  prob    <- p_calculator(params$gamma,x)
  
  ###first component on a log scale
  second1 <- log(prob)+dnorm(y,mean=mu[,1],sd=sqrt(params$sigma1),log=TRUE)
  ###second component on a log scale
  second2 <- log(1-prob)+dnorm(y,mean=mu[,2],sd=sqrt(params$sigma2),log=TRUE)
  
  ###for each point find out which component is the largest###
  largest <- apply(cbind(second1,second2),1,max)
  
  ###calculate the density of the mixture on a log scale and subtract the largest###
  first   <- -1/(dmix_log(x,y,params)-largest)
  
  ###subtract the largest from both component densities and multiply by the above step###
  adjusted <- cbind(first*(second1-largest),first*(second2-largest))
  
  ###return back to the normal scale###
  return(exp(adjusted))
}

ddbeta        <- function(x,y,params){
  params        <- transform_sigma(params)
  mu            <- getmu(x, params)
  first         <- ddbeta_helper_log(x,y,params,mu)
  comp1_beta0   <- first[,1]*(y-mu[,1])/params$sigma1
  comp1_beta1   <- comp1_beta0*x
  comp2_beta0   <- first[,2]*(y-mu[,2])/params$sigma2
  comp2_beta1   <- comp2_beta0*x
  return(c(sum(comp1_beta0),sum(comp1_beta1),
           sum(comp2_beta0),sum(comp2_beta1)))
}
