###Numerical Testing of the derivative functions###

###check for gamma derivative####
data.frame(
LogLik=
  c((l(x,p_calculator(gamma+c(1e-4,0), x), mu1,sigma1,mu2,sigma2)-l(x,p_calculator(gamma, x), mu1,sigma1,mu2,sigma2))/1e-4,
    (l(x,p_calculator(gamma+c(0,1e-4), x), mu1,sigma1,mu2,sigma2)-l(x,p_calculator(gamma, x), mu1,sigma1,mu2,sigma2))/1e-4),
U=
  c((sum(U(x,p_calculator(gamma+c(1e-4,0), x), mu1,sigma1,mu2,sigma2))-sum(U(x,p_calculator(gamma, x), mu1,sigma1,mu2,sigma2)))/1e-4,
    (sum(U(x,p_calculator(gamma+c(0,1e-4), x), mu1,sigma1,mu2,sigma2))-sum(U(x,p_calculator(gamma, x), mu1,sigma1,mu2,sigma2)))/1e-4),
ddgamma=
  ddgamma(x,y,prob,mu1,sigma1,mu2,sigma2)
)

###check for mu1 derivarive###
data.frame(
  LogLik= (l(x,prob, mu1+1e-4,sigma1,mu2,sigma2)-l(x,prob, mu1,sigma1,mu2,sigma2))/1e-4, 
  U=  (sum(U(x,prob, mu1+1e-4,sigma1,mu2,sigma2))-sum(U(x,prob, mu1,sigma1,mu2,sigma2)))/1e-4,
  ddmu1= ddmu(which=1,x,prob,mu1,sigma1,mu2,sigma2)
)

###check for mu2 derivarive###
data.frame(
  LogLik= (l(x,prob, mu1,sigma1,mu2+1e-4,sigma2)-l(x,prob, mu1,sigma1,mu2,sigma2))/1e-4, 
  U=  (sum(U(x,prob, mu1,sigma1,mu2+1e-4,sigma2))-sum(U(x,prob, mu1,sigma1,mu2,sigma2)))/1e-4,
  ddmu1= ddmu(which=2,x,prob,mu1,sigma1,mu2,sigma2)
)


###check for sigma1 derivarive###
data.frame(
  LogLik= (l(x,prob, mu1,sigma1+1e-4,mu2,sigma2)-l(x,prob, mu1,sigma1,mu2,sigma2))/1e-4, 
  U=  (sum(U(x,prob, mu1,sigma1+1e-4,mu2,sigma2))-sum(U(x,prob, mu1,sigma1,mu2,sigma2)))/1e-4,
  ddsigma1= ddsigma(which=1,x,prob,mu1,sigma1,mu2,sigma2)
)

###check for sigma2 derivarive###
data.frame(
  LogLik= (l(x,prob, mu1,sigma1,mu2,sigma2+1e-4)-l(x,prob, mu1,sigma1,mu2,sigma2))/1e-4, 
  U=  (sum(U(x,prob, mu1,sigma1,mu2,sigma2+1e-4))-sum(U(x,prob, mu1,sigma1,mu2,sigma2)))/1e-4,
  ddsigma2= ddsigma(which=2,x,prob,mu1,sigma1,mu2,sigma2)
)
