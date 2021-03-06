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

U             <- function(x,y,params,transform_sigma=FALSE){
  sum(-log(dmix(x,y,params,transform_sigma)))
}

U_HMC         <- function(x,y,params,transform_sigma=TRUE){
  sum(-log(dmix(x,y,params,transform_sigma)))
}

l             <- function(x,y,params,transform_sigma=FALSE){
  sum(log(dmix(x,y,params,transform_sigma)))
}

ddgamma       <- function(x,y,params){
  params     <- transform_sigma(params)
  prob       <- p_calculator(params$gamma,x)
  mu         <- getmu(x, params)
  colSums(-1/dmix(x,y,params)*
            prob*(1-prob)*
            (dnorm(y,mu[,1],sqrt(params$sigma1))-dnorm(y,mu[,2],sqrt(params$sigma2)))
          *cbind(1,x))
  
}


ddbeta_helper <- function(x,y,params,mu){
  
  prob    <- p_calculator(params$gamma,x)

  ###first component 
  second1 <- prob*dnorm(y,mean=mu[,1],sd=sqrt(params$sigma1))
  ###second component 
  second2 <- (1-prob)*dnorm(y,mean=mu[,2],sd=sqrt(params$sigma2))

  first   <- -1/dmix(x,y,params)
  return(cbind(first*second1,first*second2))
}

ddbeta_helper_log <- function(x,y,params,mu){
  
  prob    <- p_calculator(params$gamma,x)
  
  ###first component on a log scale
  second1 <- log(prob)+dnorm(y,mean=mu[,1],sd=sqrt(params$sigma1),log=TRUE)
  ###second component on a log scale
  second2 <- log(1-prob)+dnorm(y,mean=mu[,2],sd=sqrt(params$sigma2),log=TRUE)
  
  ###for each point find out which component is the largest###
  largest <- apply(cbind(second1,second2),1,max)
  
  adjusted1 <- exp(second1-largest)
  adjusted2 <- exp(second2-largest)
  
  adjusted <- cbind(-adjusted1/(adjusted1+adjusted2),
                    -adjusted2/(adjusted1+adjusted2))

  return(adjusted)
  
}

ddbeta        <- function(x,y,params){
  params        <- transform_sigma(params)
  mu            <- getmu(x, params)
  #first         <- ddbeta_helper(x,y,params,mu)
  first         <- ddbeta_helper_log(x,y,params,mu)
  comp1_beta0   <- first[,1]*(y-mu[,1])/params$sigma1
  comp1_beta1   <- comp1_beta0*x
  comp2_beta0   <- first[,2]*(y-mu[,2])/params$sigma2
  comp2_beta1   <- comp2_beta0*x
  return(c(sum(comp1_beta0),sum(comp1_beta1),
           sum(comp2_beta0),sum(comp2_beta1)))
}



ddsigma       <- function(x,y,params){
  
  params  <- transform_sigma(params) 
  
  prob    <- p_calculator(params$gamma,x)
  mu      <- getmu(x, params)
  first   <- -1/dmix(x,y,params) 
  
  second_sigma1 <- -dnorm(y,mu[,1],sqrt(params$sigma1))*prob
  result_sigma1 <- sum(first*second_sigma1*1/(2*params$sigma1)*(1-(y-mu[,1])^2/params$sigma1))
  
  second_sigma2 <- -dnorm(y,mu[,2],sqrt(params$sigma2))*(1-prob)
  result_sigma2 <- sum(first*second_sigma2*1/(2*params$sigma2)*(1-(y-mu[,2])^2/params$sigma2))
  

  return(c(params$sigma1,params$sigma2)*c(result_sigma1, result_sigma2))
}



ddall         <- function(x,y,params,which=c("gamma","beta1","sigma1","beta2","sigma2")){
  
  which_helper <- translate_which(which)
  
  result<-c(   ddgamma=ddgamma(x,y,params),
               ddbeta1=ddbeta(x,y,params)[1:2],
               ddsigma1=ddsigma(x,y,params)[1],
               ddbeta2=ddbeta(x,y,params)[3:4],
               ddsigma2=ddsigma(x,y,params)[2]
  )
  
  result<- result[which_helper]
  names(result) <-translate_which(which, text=TRUE)
  
  return(result)
}


