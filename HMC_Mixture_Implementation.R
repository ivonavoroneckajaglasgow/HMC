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
             sigma1=0.14,
             beta2=c(23.59,1-2.75),
             sigma2=0.34
)


N      <- 3000
burnin <- 1000
epsilon<- 0.01
L      <- 10


result_gamma<-HMC(x,y,N,burnin,U,ddall,epsilon,L,current_q = c(1,1), which=1)
par(mfrow=c(1,1))
plot(result_gamma$output[seq(burnin,N,by=5),],type="l")
plot(x,p_calculator(result_gamma$output[N,],x))

result_beta1<-HMC(x,y,N,burnin,U,ddall,epsilon,L,current_q = c(1,1), which=2)
par(mfrow=c(1,1))
plot(result_beta1$output[seq(burnin,N,by=5),],type="l")

