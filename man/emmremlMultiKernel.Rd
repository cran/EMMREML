\name{emmremlMultiKernel}
\alias{emmremlMultiKernel}
\title{
Function to fit Gaussian mixed model with multiple mixed effects with known covariances.
}
\description{
This function is a wrapper for the emmreml to fit Gaussian mixed model with multiple mixed effects with known covariances. The model fitted is \deqn{y=X\beta+Z_1 u_1 +Z_2 u_2+...Z_k u_k+ e} where \eqn{y} is the $n$ vector of response variable, \eqn{X} is a \eqn{n x q} known design matrix of fixed effects,  \eqn{Z_j} is a \eqn{n x l_j} known design matrix of random effects for \eqn{j=1,2,...,k}, \eqn{\beta} is \eqn{n x 1} vector of fixed effects coefficients and  \eqn{U=(u^t_1, u^t_2,..., u^t_k)^t} and \eqn{e} are independent variables with \eqn{N_L(0, blockdiag(\sigma^2_{u_1} K_1,\sigma^2_{u_2} K_2,...,\sigma^2_{u_k} K_k))} and \eqn{N_n(0, \sigma^2_e I_n)} correspondingly. It also produces the BLUPs for the L dimensional random effects \eqn{U}.

The variance parameters for random effects are estimated as  \eqn{(\hat{w}_1, \hat{w}_2,...,\hat{w}_k)*\hat{\sigma^2_u}} where \eqn{w=(w_1,w_2,..., w_k)} are the kernel weights.

}
\usage{
emmremlMultiKernel(y, X, Zlist, Klist)
}
\arguments{
  \item{y}{
\eqn{n x 1} numeric vector
}
  \item{X}{
\eqn{n x q} matrix
}
  \item{Zlist}{
list of random effects design matrices of  dimensions \eqn{n x l_1},...,\eqn{n x l_k}
}
  \item{Klist}{
list of known relationship matrices of  dimensions \eqn{l_1 x l_1},...,\eqn{l_k x l_k}
}
}

\value{
\item{Vu }{Estimate of \eqn{\sigma^2_u}}
\item{Ve}{Estimate of \eqn{\sigma^2_e}}
\item{betahat}{BLUEs for \eqn{\beta}}
\item{uhat}{BLUPs for \eqn{u}}
\item{weights}{Estimates of kernel weights}
}
\author{
Deniz Akdemir
}
\examples{

####example
#Data from Gaussian process with three 
#(total four, including residuals) independent 
#sources of variation 

n=100
M1<-matrix(rnorm(n*10), nrow=n)

M2<-matrix(rnorm(n*20), nrow=n)

M3<-matrix(rnorm(n*5), nrow=n)


#Relationship matrices
K1<-cov(t(M1))
K2<-cov(t(M2))
K3<-cov(t(M3))
K1=K1/mean(diag(K1))
K2=K2/mean(diag(K2))
K3=K3/mean(diag(K3))


#Generate data
covY<-2*(.2*K1+.7*K2+.1*K3)+diag(n)

Y<-10+crossprod(chol(covY),rnorm(n))

#training set
Trainsamp<-sample(1:100, 80)

funout<-emmremlMultiKernel(y=Y[Trainsamp], X=matrix(rep(1, n)[Trainsamp], ncol=1), 
Zlist=list(diag(n)[Trainsamp,], diag(n)[Trainsamp,], diag(n)[Trainsamp,]),
Klist=list(K1,K2, K3))
#weights
funout$weights
#Correlation of predictions with true values in test set
uhatmat<-matrix(funout$uhat, ncol=3)
uhatvec<-rowSums(uhatmat)

cor(Y[-Trainsamp], uhatvec[-Trainsamp])
}