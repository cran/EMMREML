\name{emmreml}
\alias{emmreml}
\title{
Solver for Gaussian mixed model with known covariance structure.
}
\description{
This function estimates the parameters of the model \deqn{y=X\beta+Z u+ e} where \eqn{y} is the \eqn{n} vector of response variable, \eqn{X} is a \eqn{n x q} known design matrix of fixed effects,  \eqn{Z} is a \eqn{n x l} known design matrix of random effects, \eqn{\beta} is \eqn{q x 1} vector of fixed effects coefficients and  \eqn{u} and \eqn{e} are independent variables with \eqn{N_l(0, \sigma^2_u K)} and \eqn{N_n(0, \sigma^2_e I_n)} correspondingly. It also produces the BLUPs for the random effects.
}
\usage{
emmreml(y, X, Z, K)
}
\arguments{
  \item{y}{
\eqn{n x 1} numeric vector
}
  \item{X}{
\eqn{n x q} matrix
}
  \item{Z}{
\eqn{n x l} matrix
}
  \item{K}{
\eqn{l x l} matrix of known relationships
}
}

\value{
\item{Vu }{Estimate of \eqn{\sigma^2_u}}
\item{Ve}{Estimate of \eqn{\sigma^2_e}}
\item{betahat}{BLUEs for \eqn{\beta}}
\item{uhat}{BLUPs for \eqn{u}}

}

\examples{
n=300
M1<-matrix(rnorm(n*100), nrow=n)
K1<-cov(t(M1))
K1=K1/mean(diag(K1))

covY<-2*K1+1*diag(n)

Y<-10+crossprod(chol(covY),rnorm(n))


#training set
Trainset<-sample(1:300, 200)

funout<-emmreml(y=Y[Trainset], X=matrix(rep(1, n)[Trainset], ncol=1),
 Z=diag(n)[Trainset,], K=K1)

cor(Y[-Trainset], funout$uhat[-Trainset])

}

