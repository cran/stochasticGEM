fLatentSEIR.MH <- function(N, infectionTimes = NULL, latencyTimes = NULL, removalTimes,  
    START = NULL, priorValues = NULL,fixedLatencyDuration = NULL, minInfPeriod = 1,
    bayesReps = 10000, burnIn = 0, bayesThin = 1, updateIniInfGibbs=TRUE,
    verbose = FALSE, missingInfectionTimes = TRUE,
    infectious.density = "exponential"){
    if(!(infectious.density %in% c("exponential","gamma","weibull"))){
        stop(" 'infectious.density' should be one of\n\t 'exponential','gamma' and 'weibull'")
    }
    if(infectious.density == "exponential"){
        if(START != 2 && !is.null(START)){
        stop("Starting values should be of length 2 or 'NULL' for exponential")
        }
    }
    else{
        if(START != 3 && !is.null(START)){
        stop("Starting values should be of length 3 or 'NULL' for gamma/Weibull")
        }
    }
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
    # infection times should be valid SEIR
    initialValues <- list()
    hypotheticalStartTimes <- getValidSEIR(N = N, removalTimes = removalTimes,
        fLatent = fixedLatencyDuration)
    if(!is.null(infectionTimes)){initialValues$infectionTimes <- infectionTimes}
    else{initialValues$infectionTimes <- hypotheticalStartTimes$infectionTimes}
    if(!is.null(latencyTimes)){initialValues$latencyTimes <- latencyTimes}
    else{initialValues$latencyTimes <- hypotheticalStartTimes$latencyTimes}
    minimum.infection.duration <- max(unlist(lapply(2:length(removalTimes), function(i){
        removalTimes[i]-removalTimes[(i-1)]})))
    minimum.infection.duration <- (minimum.infection.duration + 1)# force initial sample to be valid
    if(infectious.density == "exponential"){
        if(is.null(START)){
        initialValues$infectionRate <- length(removalTimes)/N
        initialValues$removalRate <- 1/minimum.infection.duration
        }   
        else{
        initialValues$infectionRate <- START[1]
        initialValues$removalRate <- START[2]
        }
    }
    else{
        if(is.null(START)){
        initialValues$infectionRate <- length(removalTimes)/N
        initialValues$infectiousScale <- 1/minimum.infection.duration
        initialValues$infectiousShape <- 1
        }   
        else{
        initialValues$infectionRate <- START[1]
        initialValues$infectiousScale <- START[2]
        initialValues$infectiousShape <- START[3]
        }
    }
    infectedBeforeDay <- unlist(lapply(1:length(removalTimes), function(i){
        return(ifelse(removalTimes[i] == 0, 0,
        max(removalTimes[(removalTimes < removalTimes[i])]) ))
    }))
    infectedBeforeDay <- removalTimes-2
    initialInfected <- sum(removalTimes==0)
    minimumLikelyInfectionTime <- (-50)
    otherParameters <- list(infectedBeforeDay=infectedBeforeDay,
        minimum.infection.duration=minimum.infection.duration,
        fixedLatent=fixedLatencyDuration,
        minInfPeriod=minInfPeriod,
        initialInfected=initialInfected,
        minimumLikelyInfectionTime=minimumLikelyInfectionTime,
        missingInfectionTimes=missingInfectionTimes,
        updateIniInfGibbs=updateIniInfGibbs,
        verbose=verbose)
    if(infectious.density == "exponential"){
        result <- .Call("expMH_fSEIR", N, removalTimes, otherParameters, 
            priorValues, initialValues, bayesReps, bayesStart, bayesThin, bayesOut)
        result$remRateSEIR <- mcmc(data=result$remRateSEIR, start = (burnIn+bayesThin), 
            end = bayesReps, thin = bayesThin)
        result$infectiousPeriod <- mcmc(data=1/result$remRateSEIR, start = (burnIn+bayesThin), 
            end = bayesReps, thin = bayesThin)
        result$reproductionNumber <- mcmc(data=N*result$infRateSEIR/result$remRateSEIR, start = (burnIn+bayesThin), 
            end = bayesReps, thin = bayesThin)
    }
    else if(infectious.density == "gamma"){
        result <- .Call("gammaMH_fSEIR", N, removalTimes, otherParameters, 
            priorValues, initialValues, bayesReps, bayesStart, bayesThin, bayesOut)
        result$infScaleSEIR <- mcmc(data=result$infScaleSEIR, start = (burnIn+bayesThin), 
            end = bayesReps, thin = bayesThin)
        result$infShapeSEIR <- mcmc(data=result$infShapeSEIR, start = (burnIn+bayesThin), 
            end = bayesReps, thin = bayesThin)
        result$infectiousPeriod <- mcmc(data=result$infShapeSEIR/result$infScaleSEIR,
            start = (burnIn+bayesThin), end = bayesReps, thin = bayesThin)
        result$reproductionNumber <- mcmc(data=N*result$infRateSEIR*result$infectiousPeriod, 
            start = (burnIn+bayesThin), end = bayesReps, thin = bayesThin)
    }
    else if(infectious.density == "weibull"){
        result <- .Call("weibullMH_fSEIR", N, removalTimes, otherParameters, 
            priorValues, initialValues, bayesReps, bayesStart, bayesThin, bayesOut)
        result$infScaleSEIR <- mcmc(data=result$infScaleSEIR, start = (burnIn+bayesThin), 
            end = bayesReps, thin = bayesThin)
        result$infShapeSEIR <- mcmc(data=result$infShapeSEIR, start = (burnIn+bayesThin), 
            end = bayesReps, thin = bayesThin)
        result$infectiousPeriod <- mcmc(data=result$infScaleSEIR^(-1/result$infShapeSEIR)*gamma(1+1/result$infShapeSEIR),
            start = (burnIn+bayesThin), end = bayesReps, thin = bayesThin)
        result$reproductionNumber <- mcmc(data=N*result$infRateSEIR*result$infectiousPeriod, 
            start = (burnIn+bayesThin), end = bayesReps, thin = bayesThin)
    }    
    result$logLikelihood <- mcmc(data=result$logLikelihood, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$infRateSEIR <- mcmc(data=result$infRateSEIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$bayesReps <- bayesReps
    result$burnIn <- burnIn
    result$bayesThin <- bayesThin
    result$removalTimes <- removalTimes
    result$fixedLatencyDuration <- fixedLatencyDuration
    result$bayesOut <- bayesOut
    result$initialSusceptible <- N
    result$initialInfective <- sum(removalTimes==0)
    return(result)
}
