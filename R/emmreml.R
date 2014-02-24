emmreml <-
function(y,X,Z,K){
  q=dim(X)[2]
  
  n=length(y)
  
  
  spI<- diag(n)
  S<-spI-tcrossprod(X%*%solve(crossprod(X)),X)
  ZK<-Z%*%K
  offset<-log(n)
  ZKZt<-tcrossprod(ZK,Z)
  ZKZtandoffset<-ZKZt+offset*spI
  SZKZtSandoffset<-{S%*%ZKZtandoffset}%*%S
  
  svdSZKZtSandspI<-eigen(SZKZtSandoffset, symmetric=TRUE)
  Ur<-svdSZKZtSandspI$vectors[,1:(n-q)]
  
  lambda<-svdSZKZtSandspI$values[1:(n-q)]-offset
  eta<-crossprod(Ur,y)
  minimfunc<-function(delta){(n-q)*log(sum(eta^2/{lambda+delta}))+sum(log(lambda+delta))}
  optimout<-optimize(minimfunc, lower=9^(-9), upper=9^9, tol=.000001)
  deltahat<-optimout$minimum
 
  Hinvhat<-solve(ZKZt+deltahat*spI)
  XtHinvhat<-crossprod(X,Hinvhat)
  betahat<-solve(XtHinvhat%*%X,XtHinvhat%*%y)
  ehat<-(y-{X%*%betahat})
  Hinvhatehat<-Hinvhat%*%ehat
  
  sigmausqhat<-sum(eta^2/{lambda+deltahat})/(n-q)
  
  sigmaesqhat<-deltahat*sigmausqhat
  
  uhat<-crossprod(ZK,Hinvhatehat)
  
  return(list(Vu=sigmausqhat,Ve=sigmaesqhat,betahat=betahat,uhat=uhat))
}
