### Functions

marginal.dist.y <- function(x, y, A, B, V, mu_beta){
  
  if(length(x)<=1) return(NA)
  
  n <- length(x)
  
  X <- cbind(1,x)
  
  V_inv      <- solve(V)
  
  V_star     <- solve(V_inv+crossprod(X))
  
  mu_star    <- V_star%*%(V_inv%*%mu_beta+t(X)%*%y)
  
  const      <- A*log(B)+lgamma(A+n/2)+0.5*determinant(V_star, logarithm=TRUE)$modulus[1]-(n/2)*log(2*pi)-lgamma(A)-0.5*
                determinant(V, logarithm=TRUE)$modulus[1]
  
  helper     <- B+0.5*(t(mu_beta)%*%V_inv%*%mu_beta+crossprod(y)-t(mu_star)%*%solve(V_star)%*%mu_star)
  
  log_helper <- -(A+n/2)*log(helper)
  
  answer     <- const+log_helper
  
  return(answer)
}


p_calculator = function(gamma,x){
  
  a_gamma       <- exp(gamma[1]+gamma[2]*x)
  p             <- a_gamma/(1+a_gamma)
  return(p)
}

assign_paramater_values <- function(params, which, v) {
  positions <- cumsum(c(0, sapply(params[which], length)))
  positions_start <- positions[-length(positions)]+1
  positions_end <- positions[-1]
  names(positions_start) <- which
  names(positions_end) <- which
  for (i in seq_along(which))
    params[[which[i]]] <- v[positions_start[i]:positions_end[i]] 
  params
}

translate_which <- function(which,text=FALSE){
  which_numeric <-c()
  
  for(i in 1:length(which)){
    which_numeric <-c(which_numeric,which(c("gamma","beta1","sigma1","beta2","sigma2")==which[i]))
  }
  
  if(text==TRUE){
    my_list<-list(c("gamma0","gamma1"),c("beta01","beta11"),c("sigma1"),c("beta02","beta12"),"sigma2")
  }else{
  my_list<-list(c(1,2),c(3,4),5,c(6,7),8)
  }
  return(unlist(my_list[which_numeric]))
}

