SIR.MH <- function(N, infectionTimes = NULL, removalTimes, START = NULL, 
    priorValues = NULL, bayesReps = 10000, burnIn = 0, bayesThin = 1, 
    ARS = TRUE, verbose = FALSE, missingInfectionTimes = TRUE, 
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
    minimum.infection.duration <- max(unlist(lapply(2:length(removalTimes), function(i){
        removalTimes[i]-removalTimes[(i-1)]})))
    minimum.infection.duration <- (minimum.infection.duration + 1)# force initial sample to be valid
    # use homogeneous SIR model
    initialValues <- list()
    if(!is.null(infectionTimes)){initialValues$infectionTimes <- infectionTimes}
    else{initialValues$infectionTimes <- removalTimes-minimum.infection.duration}
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
    initialInfected <- sum(removalTimes==0)
    minimumLikelyInfectionTime <- (-50)
    otherParameters <- list(infectedBeforeDay=infectedBeforeDay,
        minimum.infection.duration=minimum.infection.duration,
        initialInfected=initialInfected,
        minimumLikelyInfectionTime=minimumLikelyInfectionTime,
        missingInfectionTimes=missingInfectionTimes,
        verbose=verbose)
    if(infectious.density == "exponential"){
        result <- .Call("expMH_SIR", N, removalTimes, otherParameters, 
        priorValues, initialValues, bayesReps, bayesStart, bayesThin,
        bayesOut)
        result$remRateSIR <- mcmc(data=result$remRateSIR, start = (burnIn+bayesThin), 
            end = bayesReps, thin = bayesThin)
        result$infectiousPeriod <- mcmc(data=1/result$remRateSIR, start = (burnIn+bayesThin), 
            end = bayesReps, thin = bayesThin)
        result$reproductionNumber <- mcmc(data=N*result$infRateSIR/result$remRateSIR, 
            start = (burnIn+bayesThin), end = bayesReps, thin = bayesThin)
    }
    else if(infectious.density == "gamma"){
        if(ARS){
        result <- .Call("gammaMH_SIRars", N, removalTimes, otherParameters, 
        priorValues, initialValues, bayesReps, bayesStart, bayesThin, 
        bayesOut)
        }
        else{
        result <- .Call("gammaMH_SIR", N, removalTimes, otherParameters, 
        priorValues, initialValues, bayesReps, bayesStart, bayesThin, 
        bayesOut)
        }
        result$infScaleSIR <- mcmc(data=result$infScaleSIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
        result$infShapeSIR <- mcmc(data=result$infShapeSIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
        result$infectiousPeriod <- mcmc(data=result$infShapeSIR/result$infScaleSIR,
        start = (burnIn+bayesThin), end = bayesReps, thin = bayesThin)
        result$reproductionNumber <- mcmc(data=N*result$infRateSIR*result$infectiousPeriod, 
        start = (burnIn+bayesThin), end = bayesReps, thin = bayesThin)
    }
    else if(infectious.density == "weibull"){
        if(ARS){
        result <- .Call("weibullMH_SIRars", N, removalTimes, otherParameters, priorValues, initialValues,
        bayesReps, bayesStart, bayesThin, bayesOut)
        }
        else{
        result <- .Call("weibullMH_SIR", N, removalTimes, otherParameters, priorValues, initialValues,
        bayesReps, bayesStart, bayesThin, bayesOut)
        }
        result$infScaleSIR <- mcmc(data=result$infScaleSIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
        result$infShapeSIR <- mcmc(data=result$infShapeSIR, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
        result$infectiousPeriod <- mcmc(data=result$infScaleSIR^(-1/result$infShapeSIR)*gamma(1+1/result$infShapeSIR),
        start = (burnIn+bayesThin), end = bayesReps, thin = bayesThin)
        result$reproductionNumber <- mcmc(data=N*result$infRateSIR*result$infectiousPeriod, 
        start = (burnIn+bayesThin), end = bayesReps, thin = bayesThin)
    }
    result$logLikelihood <- mcmc(data=result$logLikelihood, start = (burnIn+bayesThin), 
        end = bayesReps, thin = bayesThin)
    result$infRateSIR <- mcmc(data=result$infRateSIR, start = (burnIn+bayesThin), 
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
