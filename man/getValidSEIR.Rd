\name{getValidSEIR}
\alias{getValidSEIR}
\title{Fit a fully observed standard SIR general epidemic model}
\description{ 
  \code{getValidSEIR} is used to fit the standard SIR general epidemic model. This is an internal
  function which can also be called directly by the user. 
}

\usage{
getValidSEIR(N,removalTimes,fLatent=NULL)
}

\arguments{
  \item{N}{initial susceptible individuals}
  \item{removalTimes}{removal times}
  \item{fLatent}{some assumed value of the latency duration }
}

\details{
    This is a simple and likely to be in-efficient way of getting valid SEIR starting infection and 'end of latency' times. It uses the maximum interremoval time and the assumed latency duration.
}


\seealso{
\code{\link{SEIR.MH}}, \code{\link{fLatentSEIR.MH}}
}


\author{ Eugene Zwane \email{e.zwane@gmail.com}}}

\examples{
    data(smallpox)
    getValidSEIR(119, smallpox,	10)
}
\keyword{models}
