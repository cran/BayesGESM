\name{TexasData}
\alias{TexasData}
\docType{data}
\title{Relationship between income and demographic composition in Texas}
\description{
This dataset consists of some demographic variables for each of the 254 municipalities of the state of Texas, USA, which were obtained from the US Census of 2010.
From the \emph{American Community Survey} (ACS), collected during the period 2009-2013, some variables of labor force and labor market of the households in each municipality
of the state of Texas, were also included. Because the variables obtained from the ACS are subjected to sampling error, their analysis is based on the assumption of
measurement error, where the sampling variances are considered as the variances of the measurement errors.
}
\usage{data(TexasData)}
\format{
  A data frame with 254 observations on the following 13 variables.
  \describe{
    \item{\code{County}}{name of municipality}
    \item{\code{MeanIng}}{Average family income in the last year, in hundreds of dollars}
    \item{\code{VarMeanIng}}{Variance of average family income}
    \item{\code{PopLabor16}}{Percentage of people in the Labor force (over 16 years).}
    \item{\code{VarPopLabor16}}{Variance of percentage of people in the Labor force }
    \item{\code{P_Desemp}}{Unemployment rate}
    \item{\code{VarP_Desemp}}{Variance of unemployment rate}
    \item{\code{PbMasc}}{Male population}
    \item{\code{PbFem}}{Female population}
    \item{\code{P_PbHisp}}{Percentage of Hispanic population }
    \item{\code{P_PbNegra}}{Percentage of black population}
    \item{\code{DensPobla}}{Population density}
    \item{\code{P_PbFem}}{Percentage of female population}
  }
}
\source{ \url{http://www.census.gov/}, \url{http://www.census.gov/acs/www/}
}
\examples{
#data(TexasData)
#nu <- 3
#zeta <- nu/(nu-1)
#heter <- list(sigma2y=VarMeanIng/zeta, sigma2xi=cbind(VarPopLabor16,
#          VarP_Desemp)/zeta)
#model <- fmem(MeanIng/100 ~  PopLabor16 + P_Desemp | log(DensPobla) +
#         P_PbFem +  P_PNegra + bsp(P_PbHisp), data=TexasData, family="Slash",
#         eta=nu, burn.in=10000, post.sam.s=10000, heter=heter,thin=10)
#summary(model)
}
\keyword{datasets}
