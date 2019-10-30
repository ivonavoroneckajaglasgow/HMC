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