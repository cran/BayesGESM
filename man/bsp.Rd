\name{bsp}
\alias{bsp}
\title{Tool to approximate smooth functions by B-splines.}
\description{
\bold{bsp} is used to approximate smooth functions by B-splines.}
\usage{
bsp(x,kn)
}
\arguments{
  \item{x}{values of the explanatory variable.}
  \item{kn}{(optional) the number of internal knots . The default value is \eqn{[n^{1/5}]}, where \emph{n} is the sample size}.
}
\value{
\item{x}{the B-spline basis matrix, which is cubic and has an intercept.}
}
\details{
This function uses the routine \code{bs()} of the R package \emph{splines}.
}

\references{De Boor, C. (1978). A practical Guide to Splines. Applied Mathematical Sciences. Springer-Verlag, New York.}
\author{Luz Marina Rondon <lumarp@gmail.com> and Heleno Bolfarine}



