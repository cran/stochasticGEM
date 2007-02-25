getValidSEIR <- function(N,removalTimes,fLatent = NULL){
    nInfected <- nRemoved <- length(removalTimes)
    infectionTimes <- latencyTimes <- rep(0, nRemoved)
    inter.removal.times <- unlist(lapply(2:length(removalTimes), function(i){
        removalTimes[i]-removalTimes[(i-1)]}))
    latencyTimes[2] <- inter.removal.times[1]
    for(i in 2:length(inter.removal.times)){
        latencyTimes[(i+1)] <-
            ifelse(inter.removal.times[i]==0,latencyTimes[i],inter.removal.times[i])
    }
    minimum.infection.duration <- (max(latencyTimes)+1)
    latencyTimes[-1] <- removalTimes[-1]-(latencyTimes[-1]+1)+
        seq(from=0.1,to=0.9,length.out=nRemoved-1)
    infectionTimes[-1] <- latencyTimes[-1]-fLatent
    # Change this in the presence of information on initial infective #
    latencyTimes[1] <- (infectionTimes[2])-1
    infectionTimes[1] <- latencyTimes[1]-fLatent
    temp <- .C("consistent_SEIR",as.double(infectionTimes),as.double(latencyTimes), 
        as.double(removalTimes), as.integer(N), as.integer(nInfected), as.integer(nRemoved), 
        as.double(rep(0, (nInfected+2*nRemoved))), as.integer(rep(0, (nInfected+2*nRemoved))),
        as.integer(rep(0, (nInfected+2*nRemoved+1))), as.integer(rep(0, (nInfected+2*nRemoved+1))), 
        as.integer(rep(0, (nInfected+2*nRemoved+1))) )
    return(list(infectionTimes=infectionTimes,latencyTimes=latencyTimes))
}
