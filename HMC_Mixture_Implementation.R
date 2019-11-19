source("helper_functions.R")
source("HMC_Functions_Mixture.R")
source("HMC_Derivatives_Mixture.R")
load("data1AR.RData")
x<-loc$x
y<-loc$y

plot(x,y)

###actual parameters estimated using sampler I below###

params<-list(gamma=-c(-11.95,2.42),
             beta1=c(-1.47,2.37),
             sigma1=1.15,          ###on log scale already
             beta2=c(23.59,1-2.75),
             sigma2=1.4
)


###Derivatives###

ddall(x,y,params)


###Set parameters###
N      <- 3000
burnin <- 1000
epsilon<- 0.01
L      <- 10


result_gamma<-HMC(x,y,N,burnin,U,ddall,epsilon,L,current_q = c(1,1), which="gamma")
par(mfrow=c(1,1))
plot(result_gamma$output[seq(burnin,N,by=5),],type="l")
plot(x,p_calculator(result_gamma$output[N,],x))

result_betas<-HMC(x,y,N,burnin,U,ddall,epsilon,L,current_q = rep(c(mean(y),0),2), which=c("beta1","beta2"))
par(mfrow=c(1,1))
plot(result_betas$output[seq(burnin,N,by=5),c(1,3)],type="l")
plot(result_betas$output[seq(burnin,N,by=5),c(2,4)],type="l")
plot(x,y,main="Last Iteration Estimates")
abline(a=result_betas$output[N,1],b=result_betas$output[N,2],col=2)
abline(a=result_betas$output[N,3],b=result_betas$output[N,4],col=3)

# result_beta1<-HMC(x,y,N,burnin,U,ddall,epsilon,L,current_q = c(mean(y),0), which="beta1")
# par(mfrow=c(1,1))
# plot(result_beta1$output[seq(burnin,N,by=5),],type="l")

result_sigmas<-HMC(x,y,N,burnin,U,ddall,epsilon,L,current_q = c(0.5,0.5), which=c("sigma1","sigma2"))
par(mfrow=c(1,1))
plot(result_sigmas$output[seq(burnin,N,by=5),],type="l")



