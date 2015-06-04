\name{gesm}
\alias{gesm}
\title{Generalized Elliptical Semi-parametric Models}
\description{
\bold{gesm} is used to	obtain the statistical inference based on the Bayesian approach for regression models under the assumption that
independent additive errors follow a scale mixtures of normal distribution (i.e., normal, Student-t, slash, contaminated normal, Laplace and
symmetric hyperbolic distribution), where both location and dispersion parameters of the response variable distribution include nonparametric
additive components described by B-splines.}
\usage{
gesm(formula, data, family, eta, burn.in, post.sam.s, thin)
}
\arguments{
  \item{formula}{a symbolic description of the systematic component of the model to be fitted. This description allows parametric and nonparametric functions in
  				 the locacion and dispersion parameters. See details for further information.}
  \item{data}{an optional data frame, list or environment containing the variables in the model.}
  \item{family}{a description of the error distribution to be used in the model. Supported distributions include \emph{Normal}, \emph{Student-t}, \emph{Slash},
               \emph{Hyperbolic}, \emph{Laplace} and \emph{ContNormal}, which correspond to normal, Student-t, slash, symmetric hyperbolic, Laplace and
			   contaminated normal distributions, respectively.}
  \item{eta}{(optional) a numeric value or numeric vector that represents the extra parameter of the specified error distribution. This parameter can be assumed known or unknown.}
  \item{burn.in}{the number of burn-in iterations for the MCMC algorithm.}
  \item{post.sam.s}{the required size for the posterior sample of interest parameters.}
  \item{thin}{(optional) the thinning interval used in the simulation to obtain the required size for the posterior sample.}
}
\details{
The argument \emph{formula} comprises three parts, namely: \emph{(i)} observed response variable; \emph{(ii)} covariates for the
location parameter including the nonparametric components; and \emph{(iii)} covariates for the dispersion parameter including the
nonparametric components.
The first two parts are separated by the symbol "~" and the second and third parts are separated by the symbol "|".
Furthemore, the nonparametric components can be specified by using the function \code{bsp()} in the second and third parts of the argument
\emph{formula}.

We implemented an efficient MCMC algorithm by combining Gibbs sampler and Metropolis-Hastings algorithm, which is mainly based on the ability
of the B-splines to be expressed linearly and on the fact that the distribution of the model error can be obtained as scale mixture of
normal distributions. We assume that a priori, the four parameters vectors (parametric and nonparametric components on location and dispersion
submodels) are independent and normally distributed. The considered values for hyperparameters enable a direct comparison of the results with
those obtained under the classical approach.
}
\value{
\item{chains}{A matrix that contains the posterior sample of interest parameters. Each column represents the marginal posterior sample of each parameter.}
\item{res}{a vector of quantile residuals, proposed by Dunn and Smyth (1996) in the context of classical inference, but suited here to the Bayesian case.}
\item{K-L}{a vector of case-deletion influence measures based on the Kullback-Leibler divergence.}
\item{X_2}{a vector of case-deletion influence measures based on the X2-Distance divergence.}
\item{DIC}{DIC criterion for model selection.}
\item{LMPL}{Log-marginal pseudo-likelihood for model selection.}
}

\references{Rondon, L.M. and Bolfarine, H. (2015) Bayesian Analysis of Generalized Elliptical Semi-parametric Models. (submitted).

			Dunn, P.K. e Smyth, G.K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics. 5, 236-244. }

\author{Luz Marina Rondon <lumarp@gmail.com> and Heleno Bolfarine}

\seealso{
\code{\link{bsp}}
}

\examples{
##########   European Rabbit   ##########
#library(ssym)
#data(Erabbits)
#Erabbits2 <- Erabbits[order(Erabbits$age,Erabbits$wlens),]
#attach(Erabbits2)

#fit <- gesm(wlens ~ bsp(age) | bsp(age), family= "ContNormal", eta=c(0.8,0.9),
#				 burn.in=1000, post.sam.s=5000, thin=10)			 
#summary(fit)

##########   Plot nonparametric components for the location and dispersion parameters 
#par(mfrow=c(1,2))
#bsp.graph.gesm(fit, which=1, age, xlab="Rabbit age", ylab="f(age)", main="Location")
#bsp.graph.gesm(fit, which=2, age, xlab="Rabbit age", ylab="g(age)", main="Dispersion")
									
##########   Residual plot
#par(mfrow=c(1,2))
#plot(fit$res, ylim=c(-2,2), xlab="Index", ylab="", main="Residuals", cex=0.3,
#	 type="p", lwd=3)
#abline(h=0,lty=3)
#qqnorm(fit$res, xlim=c(-2,2), ylim=c(-2,2), xlab="Quantile", ylab="Residuals",
#	   cex=0.3, type="p", lwd=3)
#abline(0,1,lty=3)

##########   Influence measures plot
#par(mfrow=c(1,2))
#plot(fit$KL, xlab="Index", ylab="", main="K-L divergence", cex=0.3, type="p", lwd=3)
#abline(h=3*mean(fit$KL))
#plot(fit$X_2, xlab="Index", ylab="", main="X2 divergence", cex=0.3, type="p", lwd=3)
#abline(h=3*mean(fit$X_2))
}

