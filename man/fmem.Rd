\name{fmem}
\alias{fmem}
\title{Flexible Measurement Error Models}
\description{
\bold{fmem} is used to obtain the statistical inference based on the Bayesian approach for the structural version
of the flexible measurement error models under the presence of homoscedastic and heteroscedastic random errors.
These models admits vectors of explanatory variables with and without measurement error as well as the presence of nonlinear effects, which is approximated
by using B-splines. The error-prone variables and the random error follow scale mixtures of normal distributions.}
\usage{
fmem(formula, data, omeg, family, eta, burn.in, post.sam.s, thin, heter)
}
\arguments{
 \item{formula}{a symbolic description of the systematic component of the model to be fitted. See details for further information.}
 \item{data}{an optional data frame, list or environment containing the variables in the model.}
 \item{omeg}{(optional) the ratio \eqn{\omega=\sigma_{\epsilon}^2/\sigma_{\xi}^2}, this value must be specified when the model of interest is the homocedastic
             flexible measurement model. If this value is not specified is assumed to be 1, that is, \eqn{\sigma_y^2=\sigma_{\xi}^2}}
 \item{family}{a description of the error-prone variables and the random error distributions to be used in the model.
 			   Supported distributions include \emph{Normal}, \emph{Student-t}, \emph{Slash},\emph{Hyperbolic},
			   \emph{Laplace} and \emph{ContNormal}, which correspond to normal, Student-t, slash, symmetric hyperbolic, Laplace and
			   contaminated normal distributions, respectively.}
  \item{eta}{(optional) a numeric value or numeric vector that represents the extra parameter of the specified error distribution. This parameter can be assumed known or unknown.}
  \item{burn.in}{the number of burn-in iterations for the MCMC algorithm.}
  \item{post.sam.s}{the required size for the posterior sample of interest parameters.}
  \item{thin}{(optional) the thinning interval used in the simulation to obtain the required size for the posterior sample.}
  \item{heter}{(optional) An object type list that contains the values \eqn{\sigma_{\epsilon_i}^2} and \eqn{\Sigma_{\xi_i}} for all \eqn{i} (\eqn{i=1,ldots,n}).
  The objects have to be specified as \code{sigma2y} and \code{sigma2xi}, i.e. \code{heter <- list(sigma2y, sigma2xi)}.
  If this argument is not specified the adjusted model is the version homocedastic.}
}
\details{
The argument \emph{formula} comprises of three parts, namely: \emph{(i)} observed response variable;
\emph{(ii)} covariates with measurement error; and \emph{(iii)} covariates without measurement error including the non-parametric
components, which can be specified by using the function \code{bsp()}.
The first two parts are separated by the symbol "~" and the second and third parts are separated by the symbol "|".

This function allows to fit the measurement error model under the presence of homocedastic and heterocedastic random errors. These models
admits vectors of explanatory variables with and without measurement error as well as the presence of nonlinear effects
approximated by using B-splines. The model investigated is the structural
version, as the error-prone variables follow scale mixtures of normal distributions.
}
\value{
\item{chains}{A matrix that contains the posterior sample of interest parameters. Each column represents the marginal posterior sample of each parameter.}
\item{res}{a vector of quantile residuals, proposed by Dunn and Smyth (1996) in the context of classical inference, but suited here to the Bayesian case.}
\item{K-L}{a vector of case-deletion influence measures based on the Kullback-Leibler divergence.}
\item{X_2}{a vector of case-deletion influence measures based on the X2-Distance divergence.}
\item{DIC}{DIC criterion for model selection.}
\item{LMPL}{Log-marginal pseudo-likelihood for model selection.}
}

\references{Rondon, L.M. and Bolfarine, H. (2015). Bayesian analysis of flexible measurement error models. (submitted)}
\author{Luz Marina Rondon <lumarp@gmail.com> and Heleno Bolfarine}

\seealso{
\code{\link{bsp}}
\code{\link{bsp.graph.fmem}}
}


\examples{
#library(SemiPar)
#### Ragweed Pollen ####
#data(ragweed)
#attach(ragweed)
############ Example ragweed data
#ragweed2 <- ragweed[year==1993]
#day.in.seas  <- day.in.seas[year==1993]
#temperature <- temperature[year==1993]
#rain <- rain[year==1993]
#wind.speed <- wind.speed[year==1993]
#ragweedn <- data.frame(ragweed2,day.in.seas,temperature,rain,wind.speed)

#model <- fmem(sqrt(ragweed2) ~ wind.speed | rain + temperature + bsp(day.in.seas),
#			   data=ragweedn,family="Normal", burn.in=500, post.sam.s=2000,
#			   thin=10, omeg=1)
#summary(model)
#
###   Plot non-parametric component
#bsp.graph.fmem(model, which=1, xlab="Day", ylab="f(Day)")


############ Example Boston data 
#library(MASS)
#data(Boston)
#attach(Boston)

#model <- fmem(log(medv) ~ nox | crim + rm + bsp(lstat) + bsp(dis), data=Boston,
#              family="ContNormal", burn.in=10000, post.sam.s=5000, omeg=4, thin=10)
#summary(model)
#
###   Plot non-parametric components
#bsp.graph.fmem(model, which=1, xlab="lstat", ylab="f(lstat)") ### for variable lstat
#bsp.graph.fmem(model, which=2, xlab="dis", ylab="f(dis)") ### for variable dis
#

}

