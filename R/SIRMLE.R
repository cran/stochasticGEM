SIR.MLE <- function(parameters, N, infectionTimes, removalTimes, infectious.density = "exponential"){
	nInfected <- length(infectionTimes);
	nRemoved <- length(removalTimes);
	initialInfected <- sum(removalTimes==0)
	if(!(infectious.density %in% c("exponential","gamma","weibull"))){
		stop(" 'infectious.density' should be one of\n\t 'exponential','gamma' and 'weibull'")
	}
	if(nInfected != nRemoved){stop("Number infected not equal to number removed")}
	if(nRemoved > (N+initialInfected)){stop("Number of cases larger than actual population")}
	if(infectious.density == "exponential" && length(parameters) != 2){
		stop("exponential density should have 2 (starting) parameters")}
	if(infectious.density != "exponential" && length(parameters) != 3){
		stop("Weibull and gamma densities should have 3 (starting) parameters")}
	SIR.FIT <- function(parameters, N, infectionTimes, removalTimes, 
		nInfected, nRemoved, infectious.density = infectious.density){
		temp <- switch(infectious.density,
		exponential = .C("expLkhood_SIR",
			as.double(parameters),as.double(infectionTimes), as.double(removalTimes),
			as.integer(N), as.integer(nInfected), as.integer(nRemoved), as.double(0),
			as.double(rep(0, (nInfected+nRemoved))), as.integer(rep(0, (nInfected+nRemoved))),
			as.integer(rep(0, (nInfected+nRemoved+1))), as.integer(rep(0, (nInfected+nRemoved+1))) ),
		gamma = .C("gammaLkhood_SIR",
			as.double(parameters), as.double(infectionTimes),
			as.double(removalTimes), as.integer(N), as.integer(nInfected), as.integer(nRemoved),
			as.double(0),as.double(rep(0, (nInfected+nRemoved))), 
			as.integer(rep(0, (nInfected+nRemoved))),
			as.integer(rep(0, (nInfected+nRemoved+1))), as.integer(rep(0, (nInfected+nRemoved+1))) ),
		weibull = .C("weibullLkhood_SIR",
			as.double(parameters),as.double(infectionTimes),
			as.double(removalTimes), as.integer(N), as.integer(nInfected), as.integer(nRemoved),
			as.double(0),as.double(rep(0, (nInfected+nRemoved))), 
			as.integer(rep(0, (nInfected+nRemoved))),
			as.integer(rep(0, (nInfected+nRemoved+1))), as.integer(rep(0, (nInfected+nRemoved+1))) )
		)
		likelihood <- temp[[7]]
		return(-likelihood)
	}
	temp <- optim(fn = SIR.FIT, par = parameters, N = N, infectionTimes = infectionTimes, 
		nInfected = nInfected, removalTimes = removalTimes, nRemoved = nRemoved, 
		infectious.density = infectious.density)
	return(temp)
}