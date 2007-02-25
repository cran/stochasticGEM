\name{SIR.MLE}
\alias{SIR.MLE}
\title{Fit a fully observed standard SIR general epidemic model}
\description{ 
  \code{SIR.MLE} is used to fit the standard SIR general epidemic model. This is an internal
  function which can also be called directly by the user. It requires that the 
  epidemic is fully observed, implying that the infection times are also available.
}

\usage{
SIR.MLE(parameters, N, infectionTimes, removalTimes,
    infectious.density = "exponential")
}

\arguments{
  \item{parameters}{starting values for parameters. a vector of length 2. needed by \code{optim}. }
  \item{N}{initial susceptible individuals}
  \item{infectionTimes}{infection times}
  \item{removalTimes}{removal times}
  \item{infectious.density}{type of density. can be either \code{exponential}, \code{gamma} and \code{weibull} }
}

\seealso{
\code{\link{SIR.MH}}
}

\value{
  a \code{list} of components returned by \code{optim}.
}

\references{
O'Neill, P.D. & Roberts, G.O. (1999). 'Bayesian inference for partially observed stochastic epidemics' \emph{J.R. Statist. Soc. A.} \bold{162}, 121-129.

O'Neill, P.D. & Becker, N.G. (2001). 'Inference for an epidemic when susceptibility varies' \emph{Biostatistics} \bold{2}, 99-108.

Streftaris, G. & Gibson, G.J. (2004). 'Bayesian inference for stochastic epidemics in closed populations' \emph{Statistical Modelling} \bold{4}, 63-75.

}



\author{ Eugene Zwane \email{e.zwane@gmail.com}}}

\examples{
    data(smallpox)
    SIR.MLE(c(1/1000,1/10), 119, smallpox-16, smallpox,
    	infectious.density = "exponential")
}
\keyword{models}
