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

assign_paramater_values <- function(params, which=c("gamma","beta1","sigma1","beta2","sigma2"), v) {
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

plot_HMC <- function(result, acceptance=NA, burnin, xlab, ylab, plot_every=1){
  
  result <- result[burnin:N,]
  acceptance <-acceptance[burnin:N]
  
  logic <- all(is.na(acceptance))==FALSE
  
  N<- nrow(result)
  
  if(logic){
    main <- paste("Acceptance rate=",round(mean(acceptance,na.rm = TRUE),2))
  }else{
    main <- c("Output of the HMC")
  }
  
  par(mfrow=c(1,1))
  plot(NULL,ylim=range(result[,2]),xlim=range(result[,1]),
       ylab=ylab,xlab=xlab)
  title(main)
  if(logic) legend("topright",c("accepted","rejected"),lty=1,col=c(1,2)) 
  
  
  
  my_seq <- seq(from=2,to=N,by=plot_every)
  
  for(i in my_seq){
    if(logic){
    colour <- c(2,1)[acceptance[i]+1]
    }else{
      colour<-1
    }
    
    points(c(result[i-1,1],result[i,1]),c(result[i-1,2],result[i,2]))
    lines(c(result[i-1,1],result[i,1]),c(result[i-1,2],result[i,2]),col=colour)
  
  }

}

plot_HMC_sigma <- function(result, acceptance=NA, burnin, plot_every=1){
  
  result <- result[burnin:N,]
  acceptance <-acceptance[burnin:N]
  
  logic <- all(is.na(acceptance))==FALSE
  
  if(logic){
    main <- paste("Acceptance rate=",round(mean(acceptance,na.rm = TRUE),2))
  }else{
    main <- c("Output of the HMC")
  }
  
  par(mfrow=c(1,2))
  
  for(i in 1:2){
  plot(NULL,ylim=range(result[,i]),xlim=c(1,(N-burnin)),
       ylab=paste("sigma",i),xlab="Iteration")
  title(main)
    
  my_seq <- seq(from=2,to=nrow(result),by=plot_every)
    
  for(j in my_seq){
    if(logic){
        colour <- c(2,1)[acceptance[j]+1]
      }else{
        colour<-1
      }
    
    points(j-1,result[j-1,i])
    points(j,result[j,i])
    lines(c(j-1,j),c(result[j-1,i],result[j,i]),col=colour)
  }
  }
}
