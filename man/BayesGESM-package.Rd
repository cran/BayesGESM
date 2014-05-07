\name{BayesGESM-package}
\alias{BayesGESM-package}
\docType{package}
\title{Bayesian Analysis of Generalized Elliptical Semiparametric Models}
\description{This package allows obtain the statistical inference based on the Bayesian approach for regression models under the assumption that
independent additive errors follow normal, Student-t, slash, contaminated normal, Laplace or symmetric hyperbolic distribution,
where both location and dispersion parameters of the response variable distribution include nonparametric additive components
described by B-splines.}
\details{
\tabular{ll}{
Package: \tab BayesGESM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-05\cr
License: \tab GPL-2 | GPL-3\cr
}
}
\author{
Luz Marina Rondon <lumarp@gmail.com> and Heleno Bolfarine

Maintainer: Luz Marina Rondon
}
\references{Rondon, L.M. and Bolfarine, H. (2014) Bayesian Analysis of Generalized Elliptical Semiparametric Models. Computational Statistics and Data Analysis, (submitted).}

\examples{
library(ssym)
data(Erabbits)
attach(Erabbits)

fit <- BayesGESM(wlens, np.mu=age, family= "Hyperbolic", eta=0.8,
				 burn.in=10, post.sam.s=10, thin=1)			 
summary(fit)
}
