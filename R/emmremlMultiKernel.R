emmremlMultiKernel <-
function(y, X, Zlist, Klist){
	q=dim(X)[2]
  	n=length(y)
   	lz<-length(Zlist)
  
  spI<- diag(n)
  S<-spI-tcrossprod(X%*%solve(crossprod(X)),X)
  
 
 Z<-c()
 for (i in 1:lz){
  Z<-cbind(Z, Zlist[[i]])
  }
  

 
  minimfunctionouter<-function(weights=rep(1/lz, lz)){
    weights=weights/sum(weights)
    Klistweighted<-Klist
  for (i in 1:lz){Klistweighted[[i]]<-weights[i]*Klist[[i]]}

 K<-.bdiag(Klistweighted)
 
  K<-as.matrix(K)
  

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
  optimout<-optimize(minimfunc, lower=0, upper=10000, tol=1e-9)
  return(optimout$objective)
  }
  weights<-optim(par=rep(1/lz, lz), fn=minimfunctionouter, 
      method = "L-BFGS-B", lower = rep(0, lz), upper = rep(1,lz))$par
  
  weights<-weights/sum(weights)
   Klistweighted<-Klist
  for (i in 1:lz){Klistweighted[[i]]<-weights[i]*Klist[[i]]}
 K<-.bdiag(Klistweighted)
 K<-as.matrix(K)
  
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
  optimout<-optimize(minimfunc, lower=0, upper=10000, tol=1e-9)
  deltahat<-optimout$minimum

  #optimout$minimum
  #print(deltahat)
  Hinvhat<-solve(ZKZt+deltahat*spI)
  XtHinvhat<-crossprod(X,Hinvhat)
  betahat<-solve(XtHinvhat%*%X,XtHinvhat%*%y)
  ehat<-(y-{X%*%betahat})
  Hinvhatehat<-Hinvhat%*%ehat
  
  sigmausqhat<-sum(eta^2/{lambda+deltahat})/(n-q)
  sigmaesqhat<-deltahat*sigmausqhat
  uhat<-crossprod(ZK,Hinvhatehat)
  
  return(list(Vu=sigmausqhat,Ve=sigmaesqhat,betahat=betahat,uhat=uhat, weights=weights))

}
