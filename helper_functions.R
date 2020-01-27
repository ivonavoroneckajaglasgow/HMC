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

plot_HMC <- function(result, acceptance=NA, burnin, xlab, ylab, plot_every=1,arrows=FALSE,lines=TRUE){
  
  ###trim burnin off###
  N          <- nrow(result)
  result     <- result[(burnin+1):N,]
  if(all(is.na(acceptance))) acceptance <- rep(1,times=N)
  acceptance <- acceptance[(burnin+1):N]
  
  ###has the first point been accepted?###
  accept     <- acceptance[1]
  col_temp   <- c(2,1)[accept+1]
  
  ###set up plotting area###
  par(mfrow=c(1,1))
  plot(NULL,ylim=range(result[,2]),xlim=range(result[,1]),ylab=ylab,xlab=xlab)
  if(all(acceptance==1)) {title("Output of HMC")}
  else{title(paste("Acceptance rate=",round(mean(acceptance,na.rm = TRUE),2)))}
  legend("topright",c("accepted","rejected"),lty=1,col=c(1,2)) 
  
  ###print the first point###
  current_point <- c(result[1,1],result[1,2])
  points(current_point[1],current_point[2], col=col_temp)
  
  ###set up plotting sequence###
  my_seq <- seq(from=2,to=nrow(result),by=plot_every)
  
  
  ###plot all the points###
  
  for(i in my_seq){
  
  ###if previous point has been accepted then set it as the current point###
    if(accept==1){
      current_point  <- c(result[i-1,1],result[i-1,2])
    }
    
  ###record the proposed point and if it has been accepted or not###
  ###if accepted-plot it in black, if rejected - in red  
    accept         <- acceptance[i]  
    proposed_point <- c(result[i,1],result[i,2])
    colour <- c(2,1)[acceptance[i]+1]
  
  ###plot the proposed point###
    points(proposed_point[1],proposed_point[2],col=colour)
    
  ###join the proposed point with the current point in relevant colour###
  ###if want the direction arrows to be shown, then set arrows=TRUE###
  if(lines==TRUE){
    if(arrows==FALSE){
      segments(current_point[1],current_point[2], proposed_point[1],proposed_point[2], col= colour)
    }else{
     arrows(current_point[1],current_point[2], proposed_point[1],proposed_point[2], col= colour)
    }
  }
  }
}

plot_HMC_sigma <- function(result, acceptance=NA, burnin, plot_every=1,arrows=FALSE,lines=TRUE){
  
  ###trim the burnin off###
  N          <- nrow(result)
  result     <- result[(burnin+1):N,]
  if(all(is.na(acceptance))) acceptance <- rep(1,times=N)
  acceptance <- acceptance[(burnin+1):N]
  
  ###has the first point been accepted?###
  accept     <- acceptance[1]
  col_temp   <- c(2,1)[accept+1]
  
  ###plot this for both sigmas###
  par(mfrow=c(1,2))
  
  for(i in 1:2){
  ###set up plotting area###
  plot(NULL,ylim=range(result[,i]),xlim=c(1,(N-burnin)),ylab=paste("sigma",i),xlab="Iteration")
  if(all(acceptance==1)) {title("Output of HMC")}
  else{title(paste("Acceptance rate=",round(mean(acceptance,na.rm = TRUE),2)))}
  
  ###set up plotting sequence###  
  my_seq <- seq(from=2,to=nrow(result),by=plot_every)
  
  
  ###set the first point as the current point###
  current_point <-c(1,result[1,i])
  
  ###plot the first point###
  points(current_point[1],current_point[2], col=col_temp)
  
    for(j in my_seq){
    
    ###if previous point has been accepted, update the current point###
    if(accept==1){
      current_point  <- c(j-1,result[j-1,i])
    }
      
    ###update the proposed point and its acceptance###  
    accept         <- acceptance[j]  
    proposed_point <- c(j,result[j,i])
    
    ###if point is accepted - plot in black, if rejected - in red###
    colour <- c(2,1)[accept+1]
    
    ###plot the proposed point in the relevant colour###
    points(proposed_point[1],proposed_point[2],col=colour)
   
    ###join the proposed point with the current point in relevant colour###
    ###if want the direction arrows to be shown, then set arrows=TRUE###
    if(lines==TRUE){
    if(arrows==FALSE){
      segments(current_point[1],current_point[2], proposed_point[1],proposed_point[2], col= colour)
    }else{
      arrows(current_point[1],current_point[2], proposed_point[1],proposed_point[2], col= colour)
    }
   }
  }
 }
}

plot_HMC_leapfrog <- function(result, leapfrog, acceptance=NA, burnin, xlab, ylab, plot_every=1, points=TRUE, arrows=FALSE){
  
  ###trim burnin off###
  N          <- nrow(result)
  L          <- max(unique(leapfrog[,1]))
  result     <- result[(burnin+1):N,]
  if(all(is.na(acceptance))) acceptance <- rep(1,times=N)
  acceptance <- acceptance[(burnin+1):N]
  leapfrog   <- leapfrog[-seq(from=1,to=burnin*L,by=1),]
  
  ###has the first point been accepted?###
  accept     <- acceptance[1]
  col_temp   <- c(2,1)[accept+1]
  
  ###set up plotting area###
  par(mfrow=c(1,1))
  plot(NULL,ylim=range(leapfrog[,3],result[,2]),xlim=range(leapfrog[,2],result[,1]),ylab=ylab,xlab=xlab)
  if(all(acceptance==1)) {title("Output of HMC")
  }else{title(paste("Acceptance rate=",round(mean(acceptance,na.rm = TRUE),2)))}
  legend("topright",c("accepted","rejected"),lty=1,col=c(1,2)) 
  
  ###print the first point###
  current_point <- c(result[1,1],result[1,2])
  points(current_point[1],current_point[2],col=col_temp,pch=20)
  
  ###set up plotting sequence###
  my_seq <- seq(from=2,to=nrow(result),by=plot_every)
  
  ###plot all the points###
  
  for(i in my_seq){
  ###if previous point has been accepted then set it as the current point###
  if(accept==1){
    current_point  <- c(result[i-1,1],result[i-1,2])
  }
  
  ###record the proposed point and if it has been accepted or not###
  ###if accepted-plot it in black, if rejected - in red  
  accept         <- acceptance[i]  
  proposed_point <- c(result[i,1],result[i,2])
  colour    <- c(2,1)[acceptance[i]+1]
  arrow_col <- c(2,3)[acceptance[i]+1]
  
  ###plot the proposed point###
  points(proposed_point[1],proposed_point[2],col=colour,pch=20)
  
  ###subset the corresponding leapfrog steps###
  my_leapfrog <- leapfrog[seq(from=(i-2)*L+1,to=(i-1)*L,by=1),]
  
  ###plot intermediate leapfrog steps###
  if(arrows==TRUE){
  arrows(current_point[1],current_point[2],my_leapfrog[1,2],my_leapfrog[1,3], col= arrow_col)
  }else{
    segments(current_point[1],current_point[2],my_leapfrog[1,2],my_leapfrog[1,3], col= arrow_col)
  }
  
  for(j in 2:L){
    if(points==TRUE)
    points(my_leapfrog[j,2],my_leapfrog[j,3],col=arrow_col)
    if(arrows==TRUE){
    arrows(my_leapfrog[j-1,2],my_leapfrog[j-1,3], my_leapfrog[j,2],my_leapfrog[j,3], col=arrow_col)
    }else{
      segments(my_leapfrog[j-1,2],my_leapfrog[j-1,3], my_leapfrog[j,2],my_leapfrog[j,3], col=arrow_col)
    }
  }

} 
  
  ###join the proposed point with the current point in relevant colour###
  ###if want the direction arrows to be shown, then set arrows=TRUE###
  #if(lines==TRUE){
   # if(arrows==FALSE){
    #  segments(current_point[1],current_point[2], proposed_point[1],proposed_point[2], col= colour)
    #}else{
    #  arrows(current_point[1],current_point[2], proposed_point[1],proposed_point[2], col= colour)
    #}
  
}
  
  
