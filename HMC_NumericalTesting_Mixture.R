source("HMC_Derivatives_Mixture.R")
source("helper_functions.R")
load("data1AR.RData")
x<-loc$x
y<-loc$y

params<-list(gamma=c(8,-2.5),
             beta1=c(5,2),
             sigma1=log(3.158193), ### ls=log(sigma)
             beta2=c(3,1),
             sigma2=log(4.0552)
)

h <- 1e-6


###test betas####
params_test <- params
params_test$beta1[1] <- params_test$beta1[1]+h

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/h,
  ddbeta = ddbeta(x,y,params)[1]
)

params_test <- params
params_test$beta1[2] <- params_test$beta1[2]+h

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/h,
  ddbeta = ddbeta(x,y,params)[2]
)

params_test <- params
params_test$beta2[1] <- params_test$beta2[1]+h

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/h,
  ddbeta = ddbeta(x,y,params)[3]
)


params_test <- params
params_test$beta2[2] <- params_test$beta2[2]+h

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/h,
  ddbeta = ddbeta(x,y,params)[4]
)

###test gamma###
params_test <- params
params_test$gamma[1] <- params_test$gamma[1]+h

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/h,
  ddgamma = ddgamma(x,y,params)[1]
)

params_test <- params
params_test$gamma[2] <- params_test$gamma[2]+h

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/h,
  ddgamma = ddgamma(x,y,params)[2]
)

###test sigma###

params_test        <- params
params_test$sigma1 <- params_test$sigma1+h

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/h,
  ddsigma1 = ddsigma(x,y,params)[1]
)

params_test        <- params
params_test$sigma2 <- params_test$sigma2+h

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/h,
  ddsigma2 = ddsigma(x,y,params)[2]
)

###test log sigma???###
params_new        <- params
params_new$sigma1 <- exp(params_new$sigma1)*exp(h) ### exp(params_new$sigma1+h)
params_new$sigma2 <- exp(params_new$sigma2)

params_compare        <- params
params_compare$sigma1 <- exp(params_compare$sigma1)
params_compare$sigma2 <- exp(params_compare$sigma2)

data.frame(
  LogLik= (l(x,y,params_new)-l(x,y,params_compare))/h, 
  U=  (sum(U(x,y,params_new))-sum(U(x,y,params_compare)))/h,
  ddlogsigma1 = ddlogsigma(x,y,params)[1]
)

params_new        <- params
params_new$sigma1 <- exp(params_new$sigma1)
params_new$sigma2 <- exp(params_new$sigma2)*exp(h) ### exp(params_new$sigma2+h)

params_compare        <- params
params_compare$sigma1 <- exp(params_compare$sigma1)
params_compare$sigma2 <- exp(params_compare$sigma2)

data.frame(
  LogLik= (l(x,y,params_new)-l(x,y,params_compare))/h, 
  U=  (sum(U(x,y,params_new))-sum(U(x,y,params_compare)))/h,
  ddlogsigma1 = ddlogsigma(x,y,params)[2]
)





