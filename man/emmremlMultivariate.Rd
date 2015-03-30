\name{emmremlMultivariate}
\alias{emmremlMultivariate}
\title{
Function to fit multivariate Gaussian mixed model with with known covariance structure.
}
\description{
This function estimates the parameters of the model \deqn{Y=BX+GZ+ E} where \eqn{Y} is the \eqn{d x n} matrix of response variable, \eqn{X} is a \eqn{q x n} known design matrix of fixed effects,  \eqn{Z} is a \eqn{l x n} known design matrix of random effects, \eqn{B} is \eqn{d x q} matrix of fixed effects coefficients and  \eqn{G} and \eqn{E} are independent matrix variate variables with \eqn{N_{d x l}(0, V_G, K)} and \eqn{N_{d x n}(0, V_E, I_n)} correspondingly. It also produces the BLUPs for the random effects G and some other statistics.
}
\usage{
emmremlMultivariate(Y, X, Z, K,varBhat=TRUE,varGhat=TRUE,
PEVGhat=TRUE, test=T,tolpar=1e-06, tolparinv=1e-06)
}
\arguments{
  \item{Y}{\eqn{d x n} matrix of response variable
}
  \item{X}{\eqn{q x n} known design matrix of fixed effects
}
  \item{Z}{\eqn{l x n} known design matrix of random effects
}
  \item{K}{\eqn{l x l} matrix of known relationships
}
\item{varBhat}{TRUE or FALSE
}
\item{varGhat}{ TRUE or FALSE
 	}
\item{PEVGhat}{  TRUE or FALSE
}
\item{test}{TRUE or FALSE
}
  \item{tolpar}{ tolerance parameter for convergence
}
  \item{tolparinv}{ tolerance parameter for matrix inverse
}
}

\value{
\item{Vg }{Estimate of \eqn{V_G}}
\item{Ve}{Estimate of \eqn{V_E}}
\item{Bhat}{BLUEs for \eqn{B}}
\item{Gpred}{BLUPs for \eqn{G}}
\item{XsqtestB}{\eqn{\chi^2} test statistics for testing whether the fixed effect coefficients are equal to zero.}
\item{pvalB}{pvalues obtained from large sample theory for the fixed effects. We report the pvalues adjusted by the "padjust" function for all fixed effect coefficients.}
\item{XsqtestG}{\eqn{\chi^2} test statistic values for testing whether the BLUPs are equal to zero.}
\item{pvalG}{pvalues obtained from large sample theory for the BLUPs. We report the pvalues adjusted by the "padjust" function.}
\item{varGhat}{Large sample variance for BLUPs.}
\item{varBhat}{Large sample variance for the elements of B.}
\item{PEVGhat}{Prediction error variance estimates for the BLUPs.}
\item{loglik}{loglikelihood for the model.}
}



\examples{
l=20
n<-15
m<-40

M<-matrix(rbinom(m*l,2,.2),nrow=l)
rownames(M)<-paste("l",1:nrow(M))
beta1<-rnorm(m)*exp(rbinom(m,5,.2))
beta2<-rnorm(m)*exp(rbinom(m,5,.1))
beta3<- rnorm(m)*exp(rbinom(m,5,.1))+beta2

g1<-M\%*\%beta1
g2<-M\%*\%beta2
g3<-M\%*\%beta3
e1<-sd(g1)*rnorm(l)
e2<-(-e1*2*sd(g2)/sd(g1)+.25*sd(g2)/sd(g1)*rnorm(l))
e3<-1*(e1*.25*sd(g2)/sd(g1)+.25*sd(g2)/sd(g1)*rnorm(l))

y1<-10+g1+e1
y2<--50+g2+e2
y3<--5+g3+e3

Y<-rbind(t(y1),t(y2), t(y3))

colnames(Y)<-rownames(M)
cov(t(Y))
Y[1:3,1:5]

K<-cov(t(M))
K<-K/mean(diag(K))
rownames(K)<-colnames(K)<-rownames(M)
X<-matrix(1,nrow=1,ncol=l)
colnames(X)<-rownames(M)
Z<-diag(l)
rownames(Z)<-colnames(Z)<-rownames(M)
SampleTrain<-sample(rownames(Z),n)
Ztrain<-Z[rownames(Z)\%in\%SampleTrain,]
Ztest<-Z[!(rownames(Z)\%in\%SampleTrain),]


##For a quick answer, tolpar is set to 1e-4.  Correct this in practice.
outfunc<-emmremlMultivariate(Y=Y\%*\%t(Ztrain),
X=X\%*\%t(Ztrain), Z=t(Ztrain),
K=K,tolpar=1e-4,varBhat=FALSE,
varGhat=FALSE, PEVGhat=FALSE, test=FALSE)

Yhattest<-outfunc$Gpred\%*\%t(Ztest)

cor(cbind(Ztest\%*\%Y[1,],Ztest\%*\%outfunc$Gpred[1,],
Ztest\%*\%Y[2,],Ztest\%*\%outfunc$Gpred[2,],Ztest\%*\%Y[3,],Ztest\%*\%outfunc$Gpred[3,]))


outfuncRidgeReg<-emmremlMultivariate(Y=Y\%*\%t(Ztrain),X=X\%*\%t(Ztrain), Z=t(Ztrain\%*\%M),
K=diag(m),tolpar=1e-5,varBhat=FALSE,varGhat=FALSE,
PEVGhat=FALSE, test=FALSE)

Gpred2<-outfuncRidgeReg$Gpred\%*\%t(M)
cor(Ztest\%*\%Y[1,],Ztest\%*\%Gpred2[1,])
cor(Ztest\%*\%Y[2,],Ztest\%*\%Gpred2[2,])
cor(Ztest\%*\%Y[3,],Ztest\%*\%Gpred2[3,])
}
