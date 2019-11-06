###Numerical Testing of the derivative functions###

###check for gamma derivative####
###need to load data1AR.RData for this as x is also required###
data.frame(
LogLik=
  c((l(y,p_calculator(gamma+c(1e-5,0), x), mu1,sigma1,mu2,sigma2)-l(y,p_calculator(gamma, x), mu1,sigma1,mu2,sigma2))/1e-5,
    (l(y,p_calculator(gamma+c(0,1e-5), x), mu1,sigma1,mu2,sigma2)-l(y,p_calculator(gamma, x), mu1,sigma1,mu2,sigma2))/1e-5),
U=
  c((sum(U(y,p_calculator(gamma+c(1e-5,0), x), mu1,sigma1,mu2,sigma2))-sum(U(y,p_calculator(gamma, x), mu1,sigma1,mu2,sigma2)))/1e-5,
    (sum(U(y,p_calculator(gamma+c(0,1e-5), x), mu1,sigma1,mu2,sigma2))-sum(U(y,p_calculator(gamma, x), mu1,sigma1,mu2,sigma2)))/1e-5),
ddgamma=
  ddgamma(x,y,p_calculator(gamma, x),mu1,sigma1,mu2,sigma2)
)

###check for mu1 derivarive###
data.frame(
  LogLik= (l(y,prob, mu1+1e-5,sigma1,mu2,sigma2)-l(y,prob, mu1,sigma1,mu2,sigma2))/1e-5, 
  U=  (sum(U(y,prob, mu1+1e-5,sigma1,mu2,sigma2))-sum(U(y,prob, mu1,sigma1,mu2,sigma2)))/1e-5,
  ddmu1= ddmu(which=1,y,prob,mu1,sigma1,mu2,sigma2)
)

###check for mu2 derivarive###
data.frame(
  LogLik= (l(y,prob, mu1,sigma1,mu2+1e-5,sigma2)-l(y,prob, mu1,sigma1,mu2,sigma2))/1e-5, 
  U=  (sum(U(y,prob, mu1,sigma1,mu2+1e-5,sigma2))-sum(U(y,prob, mu1,sigma1,mu2,sigma2)))/1e-5,
  ddmu2= ddmu(which=2,y,prob,mu1,sigma1,mu2,sigma2)
)


###check for sigma1 derivarive###
data.frame(
  LogLik= (l(y,prob, mu1,sigma1+1e-5,mu2,sigma2)-l(y,prob, mu1,sigma1,mu2,sigma2))/1e-5, 
  U=  (sum(U(y,prob, mu1,sigma1+1e-5,mu2,sigma2))-sum(U(y,prob, mu1,sigma1,mu2,sigma2)))/1e-5,
  ddsigma1= ddsigma(which=1,y,prob,mu1,sigma1,mu2,sigma2)
)

###check for sigma2 derivarive###
data.frame(
  LogLik= (l(y,prob, mu1,sigma1,mu2,sigma2+1e-5)-l(y,prob, mu1,sigma1,mu2,sigma2))/1e-5, 
  U=  (sum(U(y,prob, mu1,sigma1,mu2,sigma2+1e-5))-sum(U(y,prob, mu1,sigma1,mu2,sigma2)))/1e-5,
  ddsigma2= ddsigma(which=2,y,prob,mu1,sigma1,mu2,sigma2)
)

