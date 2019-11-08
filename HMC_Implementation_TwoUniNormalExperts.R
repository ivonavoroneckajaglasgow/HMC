source("helper_functions.R")
source("HMC_Functions_TwoUniNormalExperts.R")
source("HMC_Derivatives_TwoUniNormalExperts.R")
library(coda)

###Generate a mixture of two univariate normal experts and test for sampling mu and sigma (gamma constant for now)

params        <- list()
params$mu1    <- 3
params$mu2    <- 5
params$sigma1 <- 0.1
params$sigma2 <- 0.5
params$gamma  <- -c(-8,2.5)

N     <- 1000
m     <- 1

set.seed(5)

p_1   <- 0.5
alloc <- sample(c(1,2),N,replace=TRUE,prob = c(p_1,1-p_1))
y     <- c()

for(i in 1:length(alloc)){
if(alloc[i]==1)
   y[i]<- rnorm(1,mean=params$mu1,sd=sqrt(params$sigma1))
if(alloc[i]==2)
   y[i]<- rnorm(1,mean=params$mu2,sd=sqrt(params$sigma2))
}


gamma           <- -c(-8,2.5)
prob            <- p_calculator(gamma,y)

###Data###
hist(y)
plot(y,prob)

###Get all derivatives###

ddall(y,y,params,which=c(1,2,3,4,5,6)) ### needs to be (x,y) for gamma to make sense

###Get derivatives for mus only###

ddall(y,y,params,which=c(1,2))


###Test sampling mus only###
###mu1 & mu2###

result<-HMC(y, y, N=3000,burnin=1000,U=U,ddall,epsilon = 0.01,L=10,
            current_q = c(mean(y),mean(y)), which = c(1,2), 
            actual=c(params$mu1,params$mu2))

mean(result$accept)
result$ESS

