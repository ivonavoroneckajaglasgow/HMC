###ADAM####

source("helper_functions.R")
#source("HMC_Functions_Mixture.R")
source("HMC_Derivatives_Mixture.R")
load("data1AR.RData")
x<-loc$x
y<-loc$y


alpha   <- 0.001
beta1   <- 0.9
beta2   <- 0.999
epsilon <- 1e-8

params <- list(gamma=-c(-11.95,2.42),
               beta1=c(-1.47,2.37),
               sigma1=1.15,          ###on log scale already
               beta2=c(23.59,1-2.75),
               sigma2=1.4
)

theta  <- params

m <-1e-8
v <-1e-8
t <-1e-8

old.theta <- list(gamma=-c(0,0),
                  beta1=c(0,0),
                  sigma1=0,          
                  beta2=c(0,0),
                  sigma2=0
)

while(all(abs(theta-old.theta)<1e-03)){ ###have to make this cond work for lists
 
#for(i in 1:50){
  
  old.theta <- theta
  
  g <- ddall(x,y,theta) 
  
  m <- beta1 * m + (1-beta1)*g
  v <- beta2 * v + (1-beta2)*g^2
  
  m_hat <- m/(1-beta1^t)
  v_hat <- v/(1-beta2^t)
  
  theta_helper        <- unlist(theta) - alpha * m_hat/(sqrt(v_hat)+epsilon) 
  names(theta_helper) <- NULL
  theta               <- assign_paramater_values(params = theta,v = theta_helper)
  
}

  theta

