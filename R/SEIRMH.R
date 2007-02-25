SEIR.MH <- function(N, infectionTimes = NULL, InfPeriod = c(1,20), 
    latencyTimes = NULL, LatPeriod = c(1,20), removalTimes,  START = NULL, 
    priorValues = NULL, bayesReps = 10000, burnIn = 0, bayesThin = 1, 
    verbose = FALSE, missingInfectionTimes = TRUE){
    if(is.null(bayesThin)){bayesThin <- 1}
    if(is.null(burnIn)){burnIn <- trunc(bayesReps/(bayesThin*2))}
    bayesStart <- burnIn+bayesThin
    bayesOut <- trunc((bayesReps - burnIn)/bayesThin)
    if(verbose){
        cat("\n")
        cat("Bayes replications        ==   ",bayesReps,"\n")
        cat("Start of Bayes analysis   ==   ",bayesStart,"\n")
        cat("Thin                      ==   ",bayesThin,"\n")
        cat("Reps used in analysis     ==   ",bayesOut,"\n")
    }
    removalTimes <- as.double(removalTimes)# removalTimes should be double
    minimum.infection.duration <- max(unlist(lapply(2:length(removalTimes), function(i){
        removalTimes[i]-removalTimes[(i-1)]})))
    minimum.infection.duration <- (minimum.infection.duration + 1)# force initial sample to be valid
    initialValues <- list()
    if(!is.null(latencyTimes)){initialValues$latencyTimes <- latencyTimes}
    else{initialValues$latencyTimes <- removalTimes-minimum.infection.duration}
    if(!is.null(infectionTimes)){initialValues$infectionTimes <- infectionTimes}
    else{initialValues$infectionTimes <- removalTimes-2*(minimum.infection.duration-1)}
    if(is.null(START)){
        initialValues$infectionRate <- length(removalTimes)/N
        initialValues$latencyRate <- 1/minimum.infection.duration
        initialValues$removalRate <- 1/minimum.infection.duration
    }
    else{
        initialValues$infectionRate <- START[1]
        initialValues$latencyRate <- START[2]
        initialValues$removalRate <- START[3]
    }
    infectedBeforeDay <- unlist(lapply(1:length(removalTimes), function(i){
        return(ifelse(removalTimes[i] == 0, 0,
        max(removalTimes[(removalTimes < removalTimes[i])]) ))
    }))
    initialInfected <- sum(removalTimes==0)
    minimumLikelyInfectionTime <- (-50)
    otherParameters <- list(infectedBeforeDay=infectedBeforeDay,
        minimum.infection.duration=minimum.infection.duration,
        InfPeriod=InfPeriod,
        LatPeriod=LatPeriod,
        initialInfected=initialInfected,
        minimumLikelyInfectionTime=minimumLikelyInfectionTime,
        missingInfectionTimes=missingInfectionTimes,
        verbose=verbose)
    result <- .Call("expMH_SEIR", N, removalTimes, otherParameters, 
        priorValues, initialValues, bayesReps, bayesStart, bayesThin, bayesOut)
    result$logLikelihood <- mcmc(data=result$logLikelihood, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$infRateSEIR <- mcmc(data=result$infRateSEIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$latRateSEIR <- mcmc(data=result$latRateSEIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$remRateSEIR <- mcmc(data=result$remRateSEIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$infectiousPeriod <- mcmc(data=1/result$remRateSEIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$latencyDuration <- mcmc(data=1/result$latRateSEIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$reproductionNumber <- mcmc(data=N*result$infRateSEIR/result$remRateSEIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$bayesReps <- bayesReps
    result$burnIn <- burnIn
    result$bayesThin <- bayesThin
    result$removalTimes <- removalTimes
    result$bayesOut <- bayesOut
    result$initialSusceptible <- N
    result$initialInfective <- sum(removalTimes==0)
    return(result)
}
