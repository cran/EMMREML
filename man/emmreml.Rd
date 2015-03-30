\name{emmreml}
\alias{emmreml}
\title{
Solver for Gaussian mixed model with known covariance structure.
}
\description{
This function estimates the parameters of the model \deqn{y=X\beta+Z u+ e} where \eqn{y} is the \eqn{n} vector of response variable, \eqn{X} is a \eqn{n x q} known design matrix of fixed effects,  \eqn{Z} is a \eqn{n x l} known design matrix of random effects, \eqn{\beta} is \eqn{q x 1} vector of fixed effects coefficients and  \eqn{u} and \eqn{e} are independent variables with \eqn{N_l(0, \sigma^2_u K)} and \eqn{N_n(0, \sigma^2_e I_n)} correspondingly. It also produces the BLUPs and some other useful statistics like large sample estimates of variances and PEV.}
\usage{
emmreml(y, X, Z, K,varbetahat=TRUE,varuhat=TRUE, PEVuhat=TRUE, test=T)
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
 \item{varbetahat}{
 TRUE or FALSE
}
 \item{varuhat}{
 TRUE or FALSE
}
 \item{PEVuhat}{
 TRUE or FALSE
}
 \item{test}{
 TRUE or FALSE
}

}

\value{
\item{Vu }{Estimate of \eqn{\sigma^2_u}}
\item{Ve}{Estimate of \eqn{\sigma^2_e}}
\item{betahat}{BLUEs for \eqn{\beta}}
\item{uhat}{BLUPs for \eqn{u}}
\item{Xsqtestbeta}{\eqn{\chi^2} test statistics for testing whether the fixed effect coefficients are equal to zero.}
\item{pvalbeta}{pvalues obtained from large sample theory for the fixed effects. We report the pvalues adjusted by the "padjust" function for all fixed effect coefficients.}

\item{Xsqtestu}{\eqn{\chi^2} test statistic values for testing whether the BLUPs are equal to zero.}
\item{pvalu}{pvalues obtained from large sample theory for the BLUPs. We report the pvalues adjusted by the "padjust" function.}
\item{varuhat}{Large sample variance for the BLUPs.}
\item{varbetahat}{Large sample variance for the \eqn{\beta}'s.}
\item{PEVuhat}{Prediction error variance estimates for the BLUPs.}
\item{loglik}{loglikelihood for the model.}
}

\examples{
n=200
M1<-matrix(rnorm(n*300), nrow=n)
K1<-cov(t(M1))
K1=K1/mean(diag(K1))

covY<-2*K1+1*diag(n)

Y<-10+crossprod(chol(covY),rnorm(n))


#training set
Trainset<-sample(1:n, 150)

funout<-emmreml(y=Y[Trainset], X=matrix(rep(1, n)[Trainset], ncol=1),
 Z=diag(n)[Trainset,], K=K1)

cor(Y[-Trainset], funout$uhat[-Trainset])

}

