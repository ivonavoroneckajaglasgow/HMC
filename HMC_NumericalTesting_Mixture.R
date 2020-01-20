h<-1e-6

params_test <- params
params_test$gamma[1] <- params_test$gamma[1]+h

data.frame(
  LogLik= (l(x,y,params_test,transform_sigma = TRUE)-l(x,y,params,transform_sigma = TRUE))/h, 
  U=  (sum(U(x,y,params_test,transform_sigma = TRUE))-sum(U(x,y,params,transform_sigma = TRUE)))/h,
  ddgamma = ddgamma(x,y,params)[1]
)

params_test <- params
params_test$gamma[2] <- params_test$gamma[2]+h

data.frame(
  LogLik= (l(x,y,params_test,transform_sigma = TRUE)-l(x,y,params,transform_sigma = TRUE))/h, 
  U=  (sum(U(x,y,params_test,transform_sigma = TRUE))-sum(U(x,y,params,transform_sigma = TRUE)))/h,
  ddgamma = ddgamma(x,y,params)[2]
)


params_test <- params
params_test$beta1[1] <- params_test$beta1[1]+h

data.frame(
  LogLik= (l(x,y,params_test,transform_sigma = TRUE)-l(x,y,params,transform_sigma = TRUE))/h, 
  U=  (sum(U(x,y,params_test,transform_sigma = TRUE))-sum(U(x,y,params,transform_sigma = TRUE)))/h,
  ddbeta = ddbeta(x,y,params)[1]
)

params_test <- params
params_test$beta1[2] <- params_test$beta1[2]+h

data.frame(
  LogLik= (l(x,y,params_test,transform_sigma = TRUE)-l(x,y,params,transform_sigma = TRUE))/h, 
  U=  (sum(U(x,y,params_test,transform_sigma = TRUE))-sum(U(x,y,params,transform_sigma = TRUE)))/h,
  ddbeta = ddbeta(x,y,params)[2]
)

params_test <- params
params_test$beta2[1] <- params_test$beta2[1]+h

data.frame(
  LogLik= (l(x,y,params_test,transform_sigma = TRUE)-l(x,y,params,transform_sigma = TRUE))/h, 
  U=  (sum(U(x,y,params_test,transform_sigma = TRUE))-sum(U(x,y,params,transform_sigma = TRUE)))/h,
  ddbeta = ddbeta(x,y,params)[3]
)


params_test <- params
params_test$beta2[2] <- params_test$beta2[2]+h

data.frame(
  LogLik= (l(x,y,params_test,transform_sigma = TRUE)-l(x,y,params,transform_sigma = TRUE))/h, 
  U=  (sum(U(x,y,params_test,transform_sigma = TRUE))-sum(U(x,y,params,transform_sigma = TRUE)))/h,
  ddbeta = ddbeta(x,y,params)[4]
)

params_test           <- transform_sigma(params)
params_test$sigma1    <- params_test$sigma1*exp(h)

params_compare        <- transform_sigma(params)

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params_compare))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params_compare)))/h,
  ddsigma1 = ddsigma(x,y,params)[1]
)

params_test           <- transform_sigma(params)
params_test$sigma2    <- params_test$sigma2*exp(h)

params_compare        <- transform_sigma(params)

data.frame(
  LogLik= (l(x,y,params_test)-l(x,y,params_compare))/h, 
  U=  (sum(U(x,y,params_test))-sum(U(x,y,params_compare)))/h,
  ddsigma2 = ddsigma(x,y,params)[2]
)
