\name{BayesGESM-package}
\alias{BayesGESM-package}
\docType{package}
\title{Bayesian Analysis of Generalized Elliptical Semi-parametric Models and Flexible Measurement Error Models}
\description{
This package allows to perform the statistical inference based on the Bayesian approach for regression models under the assumption
that independent additive errors follow normal, Student-t, slash, contaminated normal, Laplace or symmetric hyperbolic distributions, i.e.,
additive errors follow a scale mixtures of normal distributions.
The regression models considered in this package are: \emph{(i)} Generalized elliptical semi-parametric models, where both location
and dispersion parameters of the response variable distribution include nonparametric additive components described by using B-splines; and
\emph{(ii)} Flexible measurement error models, which admit explanatory variables with and without measurement additive errors as well as
the presence of a nonparametric component approximated by using B-splines.}
\details{
\tabular{ll}{
Package: \tab BayesGESM\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2015-01-09\cr
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
###### Example for Generalized Elliptical Semi-parametric Models #####
#library(ssym)
#data(Erabbits)
#Erabbits2 <- Erabbits[order(Erabbits$age,Erabbits$wlens),]
#attach(Erabbits2)

#fit <- gesm(wlens ~ bsp(age) | bsp(age), family= "Hyperbolic", eta=0.8,
#				 burn.in=100, post.sam.s=100, thin=2)			 
#summary(fit)

###### Example for Flexible Measurement Error Models ################
#library(SemiPar)
#library(Formula)
#### Ragweed Pollen ####
#data(ragweed)
#attach(ragweed)
#ragweed2 <- ragweed[year==1993]
#day.in.seas  <- day.in.seas[year==1993]
#temperature <- temperature[year==1993]
#rain <- rain[year==1993]
#wind.speed <- wind.speed[year==1993]
#ragweedn <- data.frame(ragweed2,day.in.seas,temperature,rain,wind.speed)
#
#model <- fmem(sqrt(ragweed2) ~ wind.speed | rain + temperature + bsp(day.in.seas),
#			   data=ragweedn,family="Normal", burn.in=500, post.sam.s=2000,
#			   thin=10, omeg=1)
#summary(model)
#
}

