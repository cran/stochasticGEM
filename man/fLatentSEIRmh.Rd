\name{fLatentSEIR.MH}
\alias{fLatentSEIR.MH}
\title{Fit a partially observed SEIR general epidemic model with fixed latency duration.}
\description{ 
  \code{fLatentSEIR.MH} is used to fit the SEIR general epidemic model with fixed latency duration via Markov chain Monte Carlo. In this compartmental model only the removal times are observed. The computation is based on the algorithm that is developed by O`Neill & Becker (2001).
}

\usage{
fLatentSEIR.MH(N, infectionTimes = NULL, latencyTimes = NULL, removalTimes,  
    START = NULL, priorValues = NULL, fixedLatencyDuration = NULL, minInfPeriod = 1,
    bayesReps = 10000, burnIn = 0, bayesThin = 1, updateIniInfGibbs=TRUE,
    verbose = FALSE, missingInfectionTimes = TRUE,
    infectious.density = "exponential")
}

\arguments{
  \item{N}{initial susceptible individuals}
  \item{infectionTimes}{removal times}
  \item{latencyTimes}{end of latency times}
  \item{removalTimes}{removal times}
  \item{START}{A \code{vector} with 2 or 3 elements, the infection rate and the removal rate for exponential distributed infectious periods, or the infection rate and the scale and shape parameter of the infectious period respectively for gamma/Weibull distributed infectious periods. Defaults to NULL.}
  \item{priorValues}{A \code{list} with elements \code{infectionRate} and \code{removalRate} for an exponential density or \code{infectionRate}, \code{infectiousScale} and \code{infectiousShape} for other densities, being vectors of length 2 containing the gamma prior coefficients for the corresponding parameters. A prior, \code{theta}, for the initial infection time is needed for the \code{exponential} (a scalar) and \code{gamma} (a vector of \code{length} 2). Priors are necessary, for instance when using a gamma density the first element of \code{theta} should be at least 1.}
  \item{fixedLatencyDuration}{The 'assumed' fixed value of the latency duration. The default is NULL}
  \item{minInfPeriod}{The minimum possible infectious period. The default is 1}  
  \item{bayesReps}{A positive integer denoting the number of MCMC draws. The default is 10000}
  \item{burnIn}{A positive integer denoting the burn-in interval for the Markov chain, i.e., the number of initial draws that should not be stored. The default is 0.}
  \item{bayesThin}{A positive integer denoting the thinning interval for the Markov chain, i.e., the interval between successive values of the Markov chain. The default is 1.}
  \item{updateIniInfGibbs}{Update the initial infection time. The default is TRUE}  
  \item{verbose}{Used to check activity of MCMC sampler. A dot is printed at every \code{bayesReps}/100 iteration.}
  \item{missingInfectionTimes}{Are missing values updated or fixed. By default the infection times are updated.}
  \item{infectious.density}{type of density. can be either \code{exponential}, \code{gamma} and \code{weibull} }
}

\details{
    If certain elements of the starting values are missing an attempt is made to get suitable starting values. If 'infection.density' is \code{exponential} the \code{remRateSEIR} is returned instead of \code{infScaleSEIR} and \code{infShapeSEIR}.
}


\value{
  a \code{list} of components containing the folowing elements:
  \item{logLikelihood}{A Markov chain Monte Carlo object of the 'psuedo' log likelihood for each completed data set. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{fixedLatencyDuration}{The 'assumed' fixed value of the latency duration.}
  \item{infectionTimes}{Posterior mean of the infection times.}
  \item{removalTimes}{Removal times.}
  \item{infRateSEIR}{A Markov chain Monte Carlo object of the Gibbs draws for the infection rate. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{remRateSEIR}{A Markov chain Monte Carlo object of the Gibbs draws for the removal rate. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{infScaleSEIR}{A Markov chain Monte Carlo object of the Gibbs draws for the scale parameter of the gamma/Weibull distribution. The function \code{mcmc} in the \code{coda} library is used to create the object.}
  \item{infShapeSEIR}{A Markov chain Monte Carlo object of the Metropolis-Hastings draws for the shape parameter of the gamma/Weibull distribution. The function \code{mcmc} in the \code{coda} library is used to create the object.}
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
O'Neill, P.D. & Roberts, G.O. (1999). 'Bayesian inference for partially observed stochastic epidemics' \emph{J.R. Statist. Soc. A.} \bold{162}, 121-129.

O'Neill, P.D. & Becker, N.G. (2001). 'Inference for an epidemic when susceptibility varies' \emph{Biostatistics} \bold{2}, 99-108.

O'Neill, P.D. (2002). 'A tutorial introduction to Bayesian inference for stochastic epidemic models using Markov chain Monte Carlo methods' \emph{Mathematical Biosciences} \bold{180}, 103-114.
}

\seealso{
\code{\link{getValidSEIR}}, \code{\link{SEIR.MH}}, \code{\link{SIR.MH}}, \code{\link{SIR.MLE}}
}

\examples{ 
    data(smallpox) 
    priors <- list(infectionRate = c(1.0,0.01), removalRate = c(1.0,0.01), theta = c(0,0)) 
    temp <- fLatentSEIR.MH(N=119,
        removalTimes=smallpox, priorValues=priors,
        fixedLatencyDuration = 10, bayesReps=1000,
        burnIn=500,bayesThin=1)
    summary(temp$infRateSEIR)
    summary(temp$remRateSEIR)
    summary(temp$infectiousPeriod)
    summary(temp$reproductionNumber)
    #
    #
    priors <- list(infectionRate = c(1.00,0.01),
        infectiousScale = c(1.0,0.01), infectiousShape = c(1.0,0.01), 
        theta = c(1.0,0.01))
    temp <- fLatentSEIR.MH(N=119,
        removalTimes=smallpox, priorValues=priors,
        fixedLatencyDuration = 10, bayesReps=1000,
        burnIn=500,bayesThin=1,infectious.density = "gamma")
    summary(temp$infRateSEIR)
    summary(temp$infScaleSEIR)
    summary(temp$infShapeSEIR)
    summary(temp$infectiousPeriod)
    summary(temp$reproductionNumber)
    #
    #
    priors <- list(infectionRate = c(1.00,0.01),
        infectiousScale = c(1.0,0.01), infectiousShape = c(1.0,0.01))
    temp <- fLatentSEIR.MH(N=119,
        removalTimes=smallpox, priorValues=priors,
        fixedLatencyDuration = 10, bayesReps=1000,
        burnIn=500,bayesThin=1,infectious.density = "weibull")
    summary(temp$infRateSEIR)
    summary(temp$infScaleSEIR)
    summary(temp$infShapeSEIR)
    summary(temp$infectiousPeriod)
    summary(temp$reproductionNumber)
}
\keyword{models}
\keyword{optimize}
