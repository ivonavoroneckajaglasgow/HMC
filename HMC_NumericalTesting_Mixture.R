load("data1AR.RData")
x<-loc$x
y<-loc$y

params<-list(gamma=c(8,-2.5),
             beta1=c(5,2),
             sigma1=1,
             beta2=c(3,1),
             sigma2=1
)

params_test <- params
params_test$beta1[1] <- params_test$beta1[1]+1e-5

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/1e-5, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/1e-5,
  ddbeta = ddbeta(x,y,params)[1]
)

params_test <- params
params_test$beta1[2] <- params_test$beta1[2]+1e-5

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/1e-5, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/1e-5,
  ddbeta = ddbeta(x,y,params)[2]
)

params_test <- params
params_test$beta2[1] <- params_test$beta2[1]+1e-5

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/1e-5, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/1e-5,
  ddbeta = ddbeta(x,y,params)[3]
)


params_test <- params
params_test$beta2[2] <- params_test$beta2[2]+1e-5

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/1e-5, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/1e-5,
  ddbeta = ddbeta(x,y,params)[4]
)

params_test <- params
params_test$gamma[1] <- params_test$gamma[1]+1e-5

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/1e-5, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/1e-5,
  ddgamma = ddgamma(x,y,params)[1]
)

params_test <- params
params_test$gamma[2] <- params_test$gamma[2]+1e-5

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/1e-5, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/1e-5,
  ddgamma = ddgamma(x,y,params)[2]
)

params_test <- params
params_test$sigma1 <- params_test$sigma1+1e-5

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/1e-5, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/1e-5,
  ddsigma1 = ddsigma(x,y,params)[1]
)

params_test        <- params
params_test$sigma2 <- params_test$sigma2+1e-5

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params))/1e-5, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params)))/1e-5,
  ddsigma2 = ddsigma(x,y,params)[2]
)
