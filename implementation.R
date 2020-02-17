source("helper_functions.R")
source("HMC_Functions_Mixture.R")
source("HMC_Derivatives_Mixture.R")
load("data1AR.RData")
x<-loc$x
y<-loc$y

x11()
par(mfrow=c(1,1))
plot(x,y)

###parameters estimated by sampler 1
params<-list(gamma=-c(-11.95,2.42),
             beta1=c(-1.47,2.37),
             sigma1=log(3.158193),      
             beta2=c(23.59,-2.75),
             sigma2=log(4.0552)
)




###estimate and visualise kernel density 

#library(MASS)
#library(plotly)
#den3d <- kde2d(x, y)

#plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()

###evaluate dmix on a grid and check that it looks ok               

#x_<-seq(from=1,to=8,length.out = 100)
#y_<-seq(from=1,to=10,length.out = 100)       
#f<-matrix(nrow=length(x_),ncol=length(y_))

#for(i in 1:length(x_)){
#  for(j in 1:length(y_)){
#    f[i,j]<-dmix(x_[i],y_[j],params)
#  }
#}

#plot_ly(x=x_, y=y_, z=f) %>% add_surface()

###test individual derivatives in HMC_NumericalTesting_Mixture.R

###test ddall returns what it is meant to

#test<-c(ddgamma(x,y,params),ddbeta(x,y,params)[1:2],ddsigma(x,y,params)[1],ddbeta(x,y,params)[3:4],ddsigma(x,y,params)[2])
#compare<-ddall(x,y,params)
#test==compare

###test estimating stuff###
###let's do gamma###

epsilon <- 0.27
L <- 20
N <- 1000
burnin <- 200

set.seed(7)
gamma_result<-HMC(x,y,N,burnin, U=U_HMC, ddall, epsilon, L, current_q=c(0,0), which="gamma")
plot_HMC(gamma_result$all_q,acceptance = gamma_result$accept,burnin, xlab="gamma0",ylab="gamma1")
plot_HMC(gamma_result$all_q,acceptance = gamma_result$accept,burnin, xlab="gamma0",ylab="gamma1",arrows = TRUE)

plot_HMC(gamma_result$output,acceptance = NA,burnin, xlab="gamma0",ylab="gamma1")
plot_HMC(gamma_result$output,acceptance = NA,burnin, xlab="gamma0",ylab="gamma1",lines=FALSE)

plot_HMC_leapfrog(result=gamma_result$all_q,
                  leapfrog=gamma_result$leapfrog,
                  acceptance=gamma_result$accept, 
                  burnin,xlab="gamma0",ylab="gamma1",points=FALSE)


est_gamma<-colMeans(gamma_result$output[burnin:N,])
est_mix<-p_calculator(est_gamma,x)

rbPal <- colorRampPalette(c('red','blue'))
Col <- rbPal(20)[as.numeric(cut(est_mix,breaks = 20))]
plot(x,y,pch = 20,col = Col)

mean(gamma_result$accept[burnin:N])

###beta###
epsilon <- 0.04
L <- 20
N <- 1000
burnin <- 200

set.seed(7)
beta_result<-HMC(x,y,N,burnin, U_HMC, ddall, epsilon, L, current_q=c(mean(y),0,mean(y),0), which=c("beta1","beta2"))
plot_HMC(result=beta_result$all_q[,1:2], acceptance= beta_result$accept,burnin=burnin, xlab="beta01", ylab="beta11")
#plot_HMC(result=beta_result$all_q[,1:2], acceptance= beta_result$accept,burnin=burnin, xlab="beta01", ylab="beta11",arrows = TRUE)
plot_HMC(result=beta_result$output[,1:2], acceptance= NA, burnin=burnin, xlab="beta01",ylab="beta11")
plot_HMC(result=beta_result$output[,1:2], acceptance= NA, burnin=burnin, xlab="beta01",ylab="beta11",lines = FALSE)

plot_HMC_leapfrog(result=beta_result$all_q[,1:2],
                  leapfrog=beta_result$leapfrog[,c(1,2,3)],
                  acceptance=beta_result$accept, 
                  burnin, xlab="beta01",ylab="beta11",points=FALSE)


plot_HMC(result=beta_result$all_q[,3:4], acceptance= beta_result$accept,burnin=burnin, xlab="beta01", ylab="beta11")
#plot_HMC(result=beta_result$all_q[,3:4], acceptance= beta_result$accept,burnin=burnin, xlab="beta01", ylab="beta11",arrows = TRUE)
plot_HMC(result=beta_result$output[,3:4], acceptance= NA, burnin=burnin, xlab="beta01",ylab="beta11")
plot_HMC(result=beta_result$output[,3:4], acceptance= NA, burnin=burnin, xlab="beta01",ylab="beta11",lines = FALSE)

plot_HMC_leapfrog(result=beta_result$all_q[,3:4],
                  leapfrog=beta_result$leapfrog[,c(1,4,5)],
                  acceptance=beta_result$accept, 
                  burnin, xlab="beta02",ylab="beta12",points=FALSE)


beta_est <- colMeans(beta_result$output[burnin:N,])
plot(x,y)
lines(x,beta_est[1]+beta_est[2]*x)
lines(x,beta_est[3]+beta_est[4]*x,col=2)

mean(beta_result$accept[burnin:N])

###sigma###
epsilon <- 0.05
L <- 20
N <- 1000
burnin <- 200

set.seed(7)
sigma_result<-HMC(x,y,N,burnin, U_HMC, ddall, epsilon, L, current_q=c(1,1), which=c("sigma1","sigma2"))
plot_HMC(sigma_result$all_q,acceptance =sigma_result$accept,burnin, xlab="sigma1",ylab="sigma2")
#plot_HMC(sigma_result$all_q,acceptance =sigma_result$accept,burnin, xlab="sigma1",ylab="sigma2",arrows = TRUE)

plot_HMC_sigma(sigma_result$all_q,acceptance = sigma_result$accept,burnin = burnin)
plot_HMC_sigma(sigma_result$all_q,acceptance = sigma_result$accept,burnin = burnin,arrows = TRUE)
plot_HMC_sigma(sigma_result$all_q,acceptance = NA,burnin = burnin)
plot_HMC_sigma(sigma_result$all_q,acceptance = NA,burnin = burnin,lines=FALSE)

plot_HMC_leapfrog(result=sigma_result$all_q,
                  acceptance = sigma_result$accept,
                  leapfrog = sigma_result$leapfrog,
                  burnin, xlab="sigma1",ylab="sigma2",points=FALSE)


sigma_est <- exp(colMeans(sigma_result$output[burnin:N,]))
par(mfrow=c(1,1))
plot(x,y)
lines(x,beta_est[1]+beta_est[2]*x)
lines(x,beta_est[1]+beta_est[2]*x-2*sqrt(sigma_est[1]),lty=2)
lines(x,beta_est[1]+beta_est[2]*x+2*sqrt(sigma_est[1]),lty=2)
lines(x,beta_est[3]+beta_est[4]*x,col=2)
lines(x,beta_est[3]+beta_est[4]*x-2*sqrt(sigma_est[2]),lty=2,col=2)
lines(x,beta_est[3]+beta_est[4]*x+2*sqrt(sigma_est[2]),lty=2,col=2)

mean(sigma_result$accept[burnin:N])

###all together? DOES NOT WORK

epsilon <- 0.05
L <- 20
N <- 1000
burnin <- 200

set.seed(7)
all_result<-HMC(x,y,N,burnin, U_HMC, ddall, epsilon, L, 
                  current_q=c(1,1,mean(y),0,1,mean(y),0,1))
