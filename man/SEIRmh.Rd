\name{SEIR.MH}
\alias{SEIR.MH}
\title{Fit a partially observed SEIR general epidemic model.}
\description{ 
  \code{SEIR.MH} is used to fit the SEIR general epidemic model via Markov chain Monte Carlo. In this compartmental model only the removal times are observed. The computation is based on the algorithm that is developed by O`Neill & Becker (2001).
}

\usage{
SEIR.MH(N, infectionTimes = NULL, InfPeriod = c(1,20), 
    latencyTimes = NULL, LatPeriod = c(1,20), removalTimes,  START = NULL, 
    priorValues = NULL, bayesReps = 10000, burnIn = 0, bayesThin = 1, 
    verbose = FALSE, missingInfectionTimes = TRUE)
}

\arguments{
  \item{N}{initial susceptible individuals}
  \item{infectionTimes}{removal times}
  \item{InfPeriod}{A range of giving the possible infectious periods}
  \item{latencyTimes}{end of latency times}
  \item{LatPeriod}{A range of giving the possible latency durations}
  \item{removalTimes}{removal times}
  \item{START}{A \code{vector} with 2 elements, the infection and removal rate respectively. Defaults to NULL.}
  \item{priorValues}{A \code{list} with elements\code{infectionRate}, \code{removalRate} and \code{theta}, the first 2 being vectors of length 2 containing the gamma prior coefficients for the corresponding parameters, and \code{theta} being a scalar for the exponential prior of the infection time of the initial infective. Defaults to NULL.}
  \item{bayesReps}{A positive integer denoting the number of MCMC draws. The default is 10000}
  \item{burnIn}{A positive integer denoting the burn-in interval for the Markov chain, i.e., the number of initial draws that should not be stored. The default is 0.}
  \item{bayesThin}{A positive integer denoting the thinning interval for the Markov chain, i.e., the interval between successive values of the Markov chain. The default is 1.}
  \item{verbose}{Used to check activity of MCMC sampler. A dot is printed at every \code{bayesReps}/100 iteration.}
  \item{missingInfectionTimes}{Are missing values updated or fixed. By default the infection times are updated.}
}

\details{
    If certain elements of the starting values are missing an attempt is made to get suitable starting values.
}


\value{
  a \code{list} of components containing the folowing elements:
  \item{logLikelihood}{A Markov chain Monte Carlo object of the 'psuedo' log likelihood for each completed data set. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{infectionTimes}{Posterior mean of the infection times.}
  \item{removalTimes}{Removal times.}
  \item{infRateSEIR}{A Markov chain Monte Carlo object of the Gibbs draws for the infection rate. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{latRateSEIR}{A Markov chain Monte Carlo object of the Gibbs draws for the latency rate. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{remRateSEIR}{A Markov chain Monte Carlo object of the Gibbs draws for the removal rate. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{acceptRate}{Number of accepted draws for the infection times.}
  \item{bayesReps}{The number of MCMC draws}
  \item{burnIn}{The burn-in interval for the Markov chain.}
  \item{bayesThin}{The thinning interval for the Markov chain.}
  \item{bayesOut}{Number of saved iterations.}
  \item{infectiousPeriod}{A Markov chain Monte Carlo object of the infectious period. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{reproductionNumber}{A Markov chain Monte Carlo object of the reproduction number. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{initialSusceptible}{Initial susceptible individuals.}
  \item{initialInfective}{Initial infective individuals.}
}

\author{ Eugene Zwane \email{e.zwane@gmail.com}}}


\references{
Gibson, G.J. & Renshaw, E. (1998). 'Estimating the parameters in stochastic epidemic models using Markov chain models' \emph{IMA Journal of Mathematics Applied in Medicine & Biology} \bold{16}, 19-40.

Hoehle, M., Jorgensen, E. & O'Neill, P.D. (2005). 'Inference in disease transmission experiments by using stochastic epidemic models' \emph{Appl. Statist.} \bold{54}, 349-366.

O'Neill, P.D. & Roberts, G.O. (1999). 'Bayesian inference for partially observed stochastic epidemics' \emph{J.R. Statist. Soc. A.} \bold{162}, 121-129.

O'Neill, P.D. & Becker, N.G. (2001). 'Inference for an epidemic when susceptibility varies' \emph{Biostatistics} \bold{2}, 99-108.

O'Neill, P.D. (2002). 'A tutorial introduction to Bayesian inference for stochastic epidemic models using Markov chain Monte Carlo methods' \emph{Mathematical Biosciences} \bold{180}, 103-114.
}

\seealso{
\code{\link{fLatentSEIR.MH}}, \code{\link{getValidSEIR}}, \code{\link{SIR.MH}}, \code{\link{SIR.MLE}}
}

\examples{ 
	data(smallpox) 
	priors <- list(infectionRate = c(1.00,0.01),
		latencyRate = c(1.00,0.01),
		removalRate = c(1.00,0.01), 
		theta = c(0.01,0.01))
	validStartTimes <- getValidSEIR(119,smallpox,13)
	temp <- SEIR.MH(N=119,
		infectionTimes=validStartTimes$infectionTimes,
		latencyTimes=validStartTimes$latencyTimes,
		removalTimes=smallpox, priorValues=priors,
		bayesReps=1000,burnIn=500,bayesThin=1)
	summary(temp$infRateSEIR)
	summary(temp$remRateSEIR)
	summary(temp$infectiousPeriod)
	summary(temp$reproductionNumber)
}
\keyword{models}
\keyword{optimize}
