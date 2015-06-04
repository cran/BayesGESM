\name{BayesGESM-package}
\alias{BayesGESM-package}
\docType{package}
\title{Bayesian Analysis of Generalized Elliptical Semi-Parametric Models and Flexible Measurement Error Models}
\description{
This package allows to perform the statistical inference based on the Bayesian approach for regression models under the assumption
that independent additive errors follow normal, Student-t, slash, contaminated normal, Laplace or symmetric hyperbolic distributions, i.e.,
additive errors follow a scale mixtures of normal distributions.
The regression models considered in this package are: \emph{(i)} Generalized elliptical semi-parametric models, where both location
and dispersion parameters of the response variable distribution include non-parametric additive components described by using B-splines; and
\emph{(ii)} Flexible measurement error models under the presence of homoscedastic and heteroscedastic random errors, 
which admit explanatory variables with and without measurement additive errors as well as
the presence of non-parametric components approximated by using B-splines.}
\details{
\tabular{ll}{
Package: \tab BayesGESM\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2015-06-04\cr
License: \tab GPL-2 | GPL-3\cr
}
}
\author{
Luz Marina Rondon <lumarp@gmail.com> and Heleno Bolfarine

Maintainer: Luz Marina Rondon
}
\references{Rondon, L.M. and Bolfarine, H. (2015) Bayesian Analysis of Generalized Elliptical Semi-parametric Models. (submitted).

			Rondon, L.M. and Bolfarine, H. (2015). Bayesian analysis of flexible measurement error models.(submitted)
}

\examples{
######### Example for Generalized Elliptical Semi-parametric Models #####
#library(ssym)
#data(Erabbits)
#Erabbits2 <- Erabbits[order(Erabbits$age,Erabbits$wlens),]
#attach(Erabbits2)

#fit <- gesm(wlens ~ bsp(age) | bsp(age), family= "ContNormal", eta=c(0.8,0.9),
#				 burn.in=1000, post.sam.s=5000, thin=2)			 
#summary(fit)

######### Example for Flexible Measurement Error Models ################
#### Ragweed Pollen ####
#library(SemiPar)
#data(ragweed)
#attach(ragweed)
#ragweedn <- as.data.frame(ragweed[year==1993,])
#
#model <- fmem(sqrt(ragweed) ~ wind.speed | rain + temperature + bsp(day.in.seas),
#			   data=ragweedn,family="Normal", burn.in=500, post.sam.s=2000,
#			   thin=10, omeg=1)
#summary(model)
# bsp.graph.fmem(model, 1, xlab="day.in.seas", ylab="f(day.in.seas)")
#						
#
#### Boston Data set #########
#library(MASS)
#data(Boston)
#attach(Boston)
#model <- fmem(log(medv) ~ nox | crim + rm + bsp(lstat) + bsp(dis), data=Boston,
#              family="ContNormal", burn.in=10000, post.sam.s=5000, omeg=4, thin=10)
#summary(model)			   
#bsp.graph.fmem(model,1) ### for variable lstat
#bsp.graph.fmem(model,2) ### for variable dis
#								 
}