\name{BayesGESM}
\alias{BayesGESM}
\title{Bayesian Analysis of Generalized Elliptical Semiparametric Models}
\description{
\bold{BayesGESM} is used to	obtain the statistical inference based on the Bayesian approach for regression models under the assumption that
independent additive errors follow a scale mixtures of normal distribution (i.e., normal, Student-t, slash, contaminated normal, Laplace and
symmetric hyperbolic distribution), where both location and dispersion parameters of the response variable distribution include nonparametric
additive components described by B-splines.}
\usage{
BayesGESM(response, formula.mu, np.mu, formula.dis, np.dis, family, eta, burn.in,
		  post.sam.s, thin)
}
\arguments{
  \item{response}{the response variable.}
  \item{formula.mu}{(optional) a symbolic description of the parametric function to be used for describe the location parameter.}
  \item{np.mu}{(optional) an explanatory (continuous) variable associated to the nonparametric function to be used for describe the location parameter.}
  \item{formula.dis}{(optional) a symbolic description of the parametric function to be used for describe the dispersion parameter.}
  \item{np.dis}{(optional) an explanatory (continuous) variable associated to the nonparametric function to be used for describe the dispersion parameter.}
  \item{family}{a description of the error distribution to be used in the model. Supported distributions include \emph{Normal}, \emph{Student-t}, \emph{Slash},
               \emph{Hyperbolic}, \emph{Laplace} and \emph{ContNormal}, which correspond to normal, Student-t, slash, symmetric hyperbolic, Laplace and
			   contaminated normal distributions, respectively.}
  \item{eta}{a numeric value or numeric vector that represents the extra parameter of the specified error distribution.}
  \item{burn.in}{the number of burn-in iterations for the MCMC algorithm.}
  \item{post.sam.s}{the required size for the posterior sample of interest parameters.}
  \item{thin}{(optional) the thinning interval used in the simulation to obtain the required size for the posterior sample.}
}
\details{
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
\item{EAIC}{EAIC criterion for model selection.}
\item{EBIC}{EBIC criterion for model selection.}
\item{LMPL}{Log-marginal pseudo-likelihood for model selection.}
}

\references{Rondon, L.M. and Bolfarine, H. (2014) Bayesian Analysis of Generalized Elliptical Semiparametric Models. Computational Statistics and Data Analysis, (submitted).
			Dunn, P.K. e Smyth, G.K. (1996). Randomized quantile residuals. Journal of Computational and Graphical Statistics. 5, 236-244. }

\author{Luz Marina Rondon <lumarp@gmail.com> and Heleno Bolfarine}


\examples{
##########   European Rabbit   ##########
library(ssym)
data(Erabbits)
Erabbits2 <- Erabbits[order(Erabbits$age,Erabbits$wlens),]
attach(Erabbits2)


fit <- BayesGESM(wlens, np.mu=age, np.dis=age, family= "Hyperbolic", eta=0.8,
				 burn.in=10, post.sam.s=10, thin=1)			 
summary(fit)

##########   Plot of posterior mean of location and dispersion parameters 
par(mfrow=c(1,2))
np.mu <- fit$B\%*\%(apply(fit$chains, 2, mean)[(fit$p+1):(fit$p+fit$k1)])
plot(age, wlens, xlim=range(age), ylim=range(wlens), type="p", cex=0.3, lwd=3, ylab="",
	 xlab="", main="Location")
par(new=TRUE)
plot(age, np.mu, xlim=range(age), ylim=range(wlens), type="l", col="red",
	 ylab="dry weight of eye lens", xlab="Age")

xi <- besselK(fit$eta, 2)/(fit$eta*besselK(fit$eta,1))
np.dis <- fit$D\%*\%(apply(fit$chains, 2, mean)[(fit$p+fit$k1+fit$q+1):(fit$p+fit$k1+fit$q+fit$k2)])
aa <- (wlens - np.mu)^2/xi
plot(age, aa, xlim=range(age), ylim=range(aa), cex=0.3, type="p", xlab="", ylab="", lwd=3,
	 main="Dispersion")
par(new=TRUE)
plot(age, exp(np.dis), xlim=range(age), ylim=range(aa), type="l", col="red", ylab="",
	 xlab="Age")

##########   Residual plot
par(mfrow=c(1,2))
plot(fit$res, ylim=c(-2,2), xlab="Index", ylab="", main="Residuals", cex=0.3,
	 type="p", lwd=3)
abline(h=0,lty=3)
qqnorm(fit$res, xlim=c(-2,2), ylim=c(-2,2), xlab="Quantile", ylab="Residuals",
	   cex=0.3, type="p", lwd=3)
abline(0,1,lty=3)

##########   Influence measures plot
par(mfrow=c(1,2))
plot(fit$KL, xlab="Index", ylab="", main="K-L divergence", cex=0.3, type="p", lwd=3)
abline(h=3*mean(fit$KL))
plot(fit$X_2, xlab="Index", ylab="", main="X2 divergence", cex=0.3, type="p", lwd=3)
abline(h=3*mean(fit$X_2))

}

