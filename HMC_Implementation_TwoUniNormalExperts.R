source("helper_functions.R")
source("HMC_Functions_TwoUniNormalExperts.R")
source("HMC_Derivatives_TwoUniNormalExperts.R")
library(coda)

###Generate a mixture of two univariate normal experts and test for sampling mu and sigma (gamma constant for now)

mu1    <- 2
mu2    <- 5
sigma1 <- 0.1
sigma2 <- 0.1
m      <- 1

N     <- 1000

set.seed(5)

p_1   <- 0.5
alloc <- sample(c(1,2),N,replace=TRUE,prob = c(p_1,1-p_1))
y     <- c()

for(i in 1:length(alloc)){
if(alloc[i]==1)
   y[i]<- rnorm(1,mean=mu1,sd=sqrt(sigma1))
if(alloc[i]==2)
   y[i]<- rnorm(1,mean=mu2,sd=sqrt(sigma2))
}


gamma           <- -c(-8,2.5)
prob            <- p_calculator(gamma,y)

###Data###
hist(y)
plot(y,prob)

###Test sampling###
###mu1###

result_mu1<-HMC(what="mu1",y, y, N=2000,burnin=500,U=U,grad_U=ddmu,epsilon = 0.009,L=10,current_q = mean(y),actual=mu1)
result_mu1$ESS
mean(result_mu1$accept)

###TESTING A RANGE OF EPSILONS CODE BELOW###
# 
# try_epsilon   <- seq(from=0.001,to=0.01,by=0.001)
# ESS_record    <- c()
# accept_record <-c()
# for(i in 1:length(try_epsilon)){
# set.seed(3)
# result_mu1<-HMC(what="mu1",y, y, N=2000,burnin=500,U=U,grad_U=ddmu,epsilon = try_epsilon[i],L=10,current_q = mean(y),actual=mu1)
# ESS_record<-c(ESS_record,result_mu1$ESS)
# accept_record<-c(accept_record,mean(result_mu1$accept))
# }
# plot(try_epsilon,ESS_record,pch=c(1:length(try_epsilon)),main="mu test")
# abline(h=1500,col=1)
# plot(try_epsilon,accept_record,pch=c(1:length(try_epsilon)),main="mu test")

###mu2###
result_mu2<-HMC(what="mu2",y, y, N=2000,burnin=500,U=U,grad_U=ddmu,epsilon = 0.009,L=10,current_q = mean(y),actual=mu2)
result_mu2$ESS
mean(result_mu2$accept)
###sigma1###
result_sigma1<-HMC(what="sigma1", y, y, N=2000,burnin=500,U=U,grad_U=ddsigma,epsilon = 0.001,L=10,current_q =1, actual=sigma1)
result_sigma1$ESS
mean(result_sigma1$accept)
###sigma2###
result_sigma2<-HMC(what="sigma2", y, y, N=2000,burnin=500,U=U,grad_U=ddsigma,epsilon = 0.001,L=10,current_q =1, actual=sigma2)
result_sigma2$ESS
mean(result_sigma2$accept)


############Import data with covariate x to test gamma too################
load("data1AR.RData")
x<-loc$x
y<-loc$y

plot(x,y)

result_gamma <- HMC(what="gamma", x, y, N=2000,burnin=500,U=U,grad_U=ddgamma,epsilon = 0.08,L=10,current_q =c(0,0))
result_gamma$ESS
mean(result_gamma$accept)

# try_epsilon   <- seq(from=0,to=0.15,by=0.01)[-1]
# ESS_record    <- c()
# accept_record <-c()
# for(i in 1:length(try_epsilon)){
#  set.seed(3)
#  result_gamma <- HMC(what="gamma", x, y, N=2000,burnin=500,U=U,grad_U=ddgamma,epsilon = try_epsilon[i],L=10,current_q =c(0,0))
#  ESS_record<-c(ESS_record,result_gamma$ESS[1])
#  accept_record<-c(accept_record,mean(result_gamma$accept))
# }
#  plot(try_epsilon,ESS_record,pch=c(1:length(try_epsilon)),main="Gamma Test")
#  abline(h=1500,col=1)
#  plot(try_epsilon,accept_record,pch=c(1:length(try_epsilon)),main="Gamma Test")

 