
emmremlMultivariate<-function(Y,X, Z, K,tolpar=1e-04, tolparinv=1e-04){
   Z<-t(Z)
 ECM1<-function(ytl, xtl, Vgt,Vet,Bt, deltal){
 Vlt=deltal*Vgt+Vet
 	invVlt<-solve(Vlt+tolparinv*diag(d))
 return(list(Vlt=Vlt, gtl=deltal*Vgt%*%invVlt%*%(ytl-Bt%*%xtl), Sigmalt=deltal*Vgt-deltal*Vgt%*%invVlt%*%(deltal*Vgt)))

}

wrapperECM1<-function(l){
ytl<-Yt[,l]
xtl<-Xt[,l]
 deltal<-eigZKZt$values[l]

return( ECM1(ytl=ytl, xtl=xtl, Vgt=Vgt,Vet=Vet,Bt=Bt, deltal=deltal))
 }



Vgfunc<-function(l){
	Vgl<-tcrossprod(outfromECM1[[l]]$gtl)
		return((1/n)*(1/eigZKZt$values[l])*(Vgl + outfromECM1[[l]]$Sigmalt))
		}
		
Vefunc<-function(l){
		etl <- Yt[,l] - Bt%*%Xt[,l] - outfromECM1[[l]]$gtl
		return((1/n)*((tcrossprod(etl)+ outfromECM1[[l]]$Sigmalt)))
		}



 if (sum(is.na(Y))==0){
  N<-nrow(K)
  KZt<-tcrossprod(K,Z)
  ZKZt<-Z%*%KZt
   
 eigZKZt = eigen(ZKZt)
  eigZKZt$values<-eigZKZt$values
 n<-nrow(ZKZt)
  d<-nrow(Y)
 Yt = Y%*%eigZKZt$vectors
 Xt = X%*%eigZKZt$vectors

Vgt =cov(t(Y))/2
 Vet =cov(t(Y))/2
  XttinvXtXtt<-t(Xt)%*%solve(tcrossprod(Xt))
Bt<-Yt%*%XttinvXtXtt
Vetm1<-Vet
repeat{
outfromECM1<-lapply(1:n, wrapperECM1)
Vetm1<-Vet

	Gt=sapply(outfromECM1, function(x) {cbind(x$gtl)})
	Bt = (Yt - Gt) %*% XttinvXtXtt
	
	listVgts <- lapply(1:n,Vgfunc)
	Vgt<-Reduce('+', listVgts)
	listVets <- lapply(1:n,Vefunc)
	Vet<-Reduce('+', listVets)
	
	if(abs(sum(diag(Vet - Vetm1)))/abs(sum(diag(Vetm1)))<tolpar){break}
}

HobsInve<-solve(kronecker(ZKZt,Vgt)+kronecker(diag(n),Vet)+tolparinv*diag(d*n), matrix(Y - Bt%*%X,ncol=1, byrow=F))
gpred<-kronecker(K,Vgt)%*%(kronecker(t(Z),diag(d)))%*%HobsInve
Gpred<-matrix(gpred, nrow=nrow(Y), byrow=F)
colnames(Gpred)<-rownames(K)
return(list(Bhat=Bt,Vg=Vgt,Ve=Vet, Gpred=Gpred))
}
}
