\name{mcmc.fmem}
\alias{mcmc.fmem}
\title{MCMC algorithm for Flexible Measurement Error Models}
\description{
This function implements the MCMC algorithm derived in order to draw samples of
the posterior distribution of the interest parameters in flexible measurement error models.
}
\usage{
mcmc.fmem(params)
}
\arguments{
\item{params}{An object type list, which provides the setup (i.e., values of hyperparameters, initial values, model matrices, basis functions for
the B-splines, burn-in and posterior sample size) of the model requested  by the user.}
}
\author{Luz Marina Rondon <lumarp@gmail.com> and Heleno Bolfarine}
