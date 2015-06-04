\name{bsp.graph.gesm}
\alias{bsp.graph.gesm}
\title{Tool for plotting the nonlinear effects that are approximated by using
B-splines for generalized elliptical semi-parametric models}
\description{
\bold{bsp.graph} displays the graphs of the nonparametric effects from an object of
the class \code{gesm()}.
}
\usage{
bsp.graph.gesm(object, which, var, xlab, ylab, main)
}
\arguments{
  \item{object}{An object of the class \code{gesm()} or \code{fmem()}.}
  \item{which}{An integer value, where 1 indicates location submodel, and 2 indicates dispersion submodel.}
  \item{var}{The name of the variable that is approximate using B-splines.}
  \item{xlab}{(optional) A title for the \emph{x} axis}
  \item{ylab}{(optional) A title for the \emph{y} axis}
  \item{main}{(optional) An overall title for the graph.}
}
\author{Luz Marina Rondon <lumarp@gmail.com> and Heleno Bolfarine}


