#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include "staticVariables.h"

/*
Susceptible-Latent-Infectious-Removed MCMC analysis:
	. Exponentially distributed latent periods
	. Exponentially distributed infectiousness periods
*/
static void expLikelihood_SEIR(double *parameters, double *infectionTimes,
	double *endLatency, double *removalTimes, int *N, int *nInfected, int *nRemoved, 
	double *sumSI, double *sumLatents, double *sumInfectious, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *LL, int *II)
{
	int i, k = 0, initialInfective=0, nEvents;
	double sumLogBeta=0, sumLogInfections=0, sumLatent=0, sumInfection=0, sumBetaSI=0;
	nEvents = *nInfected+ *nInfected+ *nRemoved;
	for(i = 0; i < *nInfected; ++i){
		allTimes[(i*3)] = infectionTimes[i];
		indicator[(i*3)] = 3;
		allTimes[(i*3)+1] = endLatency[i];
		indicator[(i*3)+1] = 2;
		allTimes[(i*3)+2] = removalTimes[i];
		indicator[(i*3)+2] = 1;
		if(removalTimes[i]==0){++initialInfective;}
	}
	rsort_with_index(allTimes,indicator,nEvents);
	SS[0] = *N+initialInfective;
	LL[0] = 0;
	II[0] = 0;
	for(i = 1; i < (nEvents+1); ++i){
		if(indicator[(i-1)] == 3){
			SS[i] = SS[(i-1)]-1;
			LL[i] = LL[(i-1)]+1;
			II[i] = II[(i-1)];
		}
		else if(indicator[(i-1)] == 2){
			SS[i] = SS[(i-1)];
			LL[i] = LL[(i-1)]-1;
			II[i] = II[(i-1)]+1;
		}
		else{
			SS[i] = SS[(i-1)];
			LL[i] = LL[(i-1)];
			II[i] = II[(i-1)]-1;
		}
	}
	*sumSI = 0;
	for(i = 1; i < nEvents; ++i){/* "0" is the start of observation	*/
		if(allTimes[i] != allTimes[i-1]){k = i;}
		sumBetaSI+=parameters[0]*II[i]*SS[i]*(allTimes[i]-allTimes[(i-1)]);
		if(indicator[i] == 1 && II[k] != 0){sumLogInfections+=log(II[k]);}
		if(indicator[i] == 2 && LL[k] != 0){sumLogInfections+=log(LL[k]);}
		if(indicator[i] == 3 && II[k] != 0){sumLogInfections+=log(parameters[0]*SS[k]*II[k]);}
		*sumSI+=SS[i]*II[i]*(allTimes[i]-allTimes[(i-1)]);
	}
	*sumLatents=0;
	*sumInfectious=0;	
	for(i = 0; i < *nInfected; ++i){
		sumLatent+=dexp((endLatency[i]-infectionTimes[i]),1/parameters[1],TRUE);
		sumInfection+=dexp((removalTimes[i]-endLatency[i]),1/parameters[2],TRUE);
		*sumLatents+=(endLatency[i]-infectionTimes[i]);
		*sumInfectious+=(removalTimes[i]-endLatency[i]);
	}
	*likelihood=sumLogBeta+sumLogInfections-sumBetaSI+sumLatent+sumInfection;
}

/*
Susceptible-Latent-Infectious-Removed MCMC analysis:
	. Exponentially distributed latent periods
	. Exponentially distributed infectiousness periods
*/
SEXP expMH_SEIR(SEXP N, SEXP removalTimes, SEXP otherParameters, SEXP priorValues,
	SEXP initialValues, SEXP bayesReps, SEXP bayesStart, SEXP bayesThin, SEXP bayesOut){
	/* Declarations  */
	int ii, jj, kk, ll, nInfected, nRemoved, nEvents, nProtected=0, initialInfected;
	SEXP infRateSEIR, latRateSEIR, remRateSEIR, logLikelihood;
	SEXP parameters;
	SEXP infectionTimes, infectionCandidate, infectedBeforeDay;
	SEXP latencyTimes, latencyCandidate;
	SEXP allTimes, indicator, SS, II, LL;
	double infRate, remRate, latRate, oldLkhood, newLkhood, minimumLikelyInfectionTime;
	double infRatePrior[2], remRatePrior[2], latRatePrior[2], thetaprior[2];
	double sumSI, sumLatents, sumInfectious, likelihood,logR;	
	double minInfPeriod,maxInfPeriod,minLatPeriod,maxLatPeriod;	
	int acceptRate=0, consistent=0, verbose;
	SEXP retParameters, parNames, acceptanceRate;
	SEXP infTimes, latTimes;
	/* Code  */
	GetRNGstate(); /* should be before a call to a random number generator */
	minInfPeriod = REAL(getListElement(otherParameters, "InfPeriod"))[0];
	maxInfPeriod = REAL(getListElement(otherParameters, "InfPeriod"))[1];
	minLatPeriod = REAL(getListElement(otherParameters, "LatPeriod"))[0];
	maxLatPeriod = REAL(getListElement(otherParameters, "LatPeriod"))[1];
	initialInfected = INTEGER(getListElement(otherParameters, "initialInfected"))[0];
	verbose = INTEGER(getListElement(otherParameters, "verbose"))[0];
	PROTECT(N = AS_INTEGER(N));
	++nProtected;
	PROTECT(removalTimes = AS_NUMERIC(removalTimes));
	++nProtected;
	/* priors and starting values */
	PROTECT(priorValues = AS_LIST(priorValues));
	++nProtected;
	PROTECT(initialValues = AS_LIST(initialValues));
	++nProtected;
	nRemoved = LENGTH(removalTimes); /* number of individuals removed */
	/* bayes replications, thin, etc */
	PROTECT(bayesReps = AS_INTEGER(bayesReps));
	++nProtected;
	PROTECT(bayesStart = AS_INTEGER(bayesStart));
	++nProtected;
	PROTECT(bayesThin = AS_INTEGER(bayesThin));
	++nProtected;
	PROTECT(bayesOut = AS_INTEGER(bayesOut));
	++nProtected;
	PROTECT(infRateSEIR = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(latRateSEIR = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(remRateSEIR = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(logLikelihood = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(parameters = allocVector(REALSXP,3));
	++nProtected;
	/* Infection times */
	PROTECT(infectionTimes = allocVector(REALSXP,nRemoved));
	++nProtected;
	PROTECT(infectionCandidate = allocVector(REALSXP,nRemoved));
	++nProtected;
	PROTECT(infectedBeforeDay = allocVector(REALSXP,nRemoved));
	++nProtected;
	PROTECT(infTimes = allocVector(REALSXP,nRemoved));
	++nProtected;
	/* Latency times */
	PROTECT(latencyTimes = allocVector(REALSXP,nRemoved));
	++nProtected;
	PROTECT(latencyCandidate = allocVector(REALSXP,nRemoved));
	++nProtected;
	PROTECT(latTimes = allocVector(REALSXP,nRemoved));
	++nProtected;
	for(jj = 0; jj < nRemoved; ++jj){
		REAL(infectionTimes)[jj] = REAL(getListElement(initialValues, "infectionTimes"))[jj];
		REAL(infectionCandidate)[jj] = REAL(infectionTimes)[jj];
		REAL(latencyTimes)[jj] = REAL(getListElement(initialValues, "latencyTimes"))[jj];
		REAL(latencyCandidate)[jj] = REAL(latencyTimes)[jj];
		REAL(infectedBeforeDay)[jj] = REAL(getListElement(otherParameters, "infectedBeforeDay"))[jj];
		REAL(infTimes)[jj] = 0;
		REAL(latTimes)[jj] = 0;
	}
	nInfected = LENGTH(infectionTimes);
	nEvents = 2*nInfected+nRemoved;
	PROTECT(allTimes = allocVector(REALSXP,nEvents));
	++nProtected;
	PROTECT(indicator = allocVector(INTSXP,nEvents));
	++nProtected;
	PROTECT(SS = allocVector(INTSXP,nEvents+1));
	++nProtected;
	PROTECT(II = allocVector(INTSXP,nEvents+1));
	++nProtected;
	PROTECT(LL = allocVector(INTSXP,nEvents+1));
	++nProtected;
	/* working variables */
	minimumLikelyInfectionTime = REAL(getListElement(otherParameters, "minimumLikelyInfectionTime"))[0];
	infRate = REAL(getListElement(initialValues, "infectionRate"))[0];
	latRate = REAL(getListElement(initialValues, "latencyRate"))[0];
	remRate = REAL(getListElement(initialValues, "removalRate"))[0];
	for(ii = 0; ii < 2; ++ii){
		infRatePrior[ii] = REAL(getListElement(priorValues, "infectionRate"))[ii];
		latRatePrior[ii] = REAL(getListElement(priorValues, "latencyRate"))[ii];
		remRatePrior[ii] = REAL(getListElement(priorValues, "removalRate"))[ii];
		thetaprior[ii] = REAL(getListElement(priorValues, "theta"))[ii];
	}
	REAL(parameters)[0] = infRate;
	REAL(parameters)[1] = latRate;
	REAL(parameters)[2] = remRate;
	expLikelihood_SEIR(REAL(parameters),REAL(infectionTimes),REAL(latencyTimes),
		REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
		&sumSI, &sumLatents, &sumInfectious, &likelihood,
		REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(LL),INTEGER(II));
	oldLkhood = likelihood;
	/*Rprintf("Likelihood before   =   %f\n",likelihood); */
	for(ii = 1; ii <= INTEGER(bayesReps)[0]; ++ii){
		infRate = rgamma(nInfected-1+infRatePrior[0],1/(sumSI+infRatePrior[1]));/* update infRate */
		latRate = rgamma(nRemoved+latRatePrior[0],1/(sumLatents+latRatePrior[1]));
		remRate = rgamma(nRemoved+remRatePrior[0],1/(sumInfectious+remRatePrior[1]));
		/*
		latRate = trgamma(1/maxLatPeriod,1/minLatPeriod,
			nRemoved+latRatePrior[0],(sumLatents+latRatePrior[1]));
		remRate = trgamma(1/maxInfPeriod,1/minInfPeriod,
			nRemoved+remRatePrior[0],(sumInfectious+remRatePrior[1]));
		*/
		REAL(parameters)[0] = infRate;
		REAL(parameters)[1] = latRate;
		REAL(parameters)[2] = remRate;
		expLikelihood_SEIR(REAL(parameters),REAL(infectionTimes),REAL(latencyTimes),
			REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
			&sumSI, &sumLatents, &sumInfectious, &likelihood,
			REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(LL),INTEGER(II));
		oldLkhood = likelihood;
		kk = ceil(unif_rand()*(nRemoved-1)); /*initial infection time excluded */
		consistent=0;
		while(consistent==0){
		/* Sample latency times */
		if(REAL(removalTimes)[kk]==0){REAL(latencyCandidate)[kk] =
			runif((minimumLikelyInfectionTime+minLatPeriod),REAL(infectionTimes)[kk+1]);}
		/* for first infection : infected after end of latency of initial case */
		else if(kk == 1 && REAL(latencyTimes)[kk+1] > (REAL(removalTimes)[kk]-minInfPeriod)){
			if((REAL(latencyTimes)[0]+minLatPeriod) > REAL(removalTimes)[kk]-maxInfPeriod){
				REAL(latencyCandidate)[kk] =
					runif((REAL(latencyTimes)[0]+minLatPeriod),
						(REAL(removalTimes)[kk]-minInfPeriod));}
			else{
				REAL(latencyCandidate)[kk] =
					runif((REAL(removalTimes)[kk]-maxInfPeriod),
						(REAL(removalTimes)[kk]-minInfPeriod));}
		}
		else if(kk == 1 && REAL(latencyTimes)[kk+1] <= (REAL(removalTimes)[kk]-minInfPeriod)){
			if((REAL(latencyTimes)[0]+minLatPeriod) > REAL(removalTimes)[kk]-maxInfPeriod){
				REAL(latencyCandidate)[kk] =
					runif((REAL(latencyTimes)[0]+minLatPeriod),
						REAL(latencyTimes)[kk+1]);}
			else{
				REAL(latencyCandidate)[kk] =
					runif((REAL(removalTimes)[kk]-maxInfPeriod),
						REAL(latencyTimes)[kk+1]);}
		}
		/* for last removal */
		else if(kk == nRemoved-1){
			REAL(latencyCandidate)[kk] =
			runif(REAL(latencyTimes)[kk-1],(REAL(removalTimes)[kk]-minInfPeriod));}
		/* for all other cases */
		else if(REAL(latencyTimes)[kk+1] > (REAL(removalTimes)[kk]-minInfPeriod)){
			if(REAL(latencyTimes)[kk-1] > REAL(removalTimes)[kk]-maxInfPeriod){
				REAL(latencyCandidate)[kk]=
					runif(REAL(latencyTimes)[kk-1],
						(REAL(removalTimes)[kk]-minInfPeriod));}
			else{
				REAL(latencyCandidate)[kk]=
					runif(REAL(removalTimes)[kk]-maxInfPeriod,
						(REAL(removalTimes)[kk]-minInfPeriod));}
		}
		else{
			if(REAL(latencyTimes)[kk-1] > REAL(removalTimes)[kk]-maxInfPeriod){
				REAL(latencyCandidate)[kk]=
					runif(REAL(latencyTimes)[kk-1],	REAL(latencyTimes)[kk+1]);}
			else{
				REAL(latencyCandidate)[kk]=
					runif(REAL(removalTimes)[kk]-maxInfPeriod,REAL(latencyTimes)[kk+1]);}
		}
		/* Sample infection times (make them ordered) */
		if(REAL(removalTimes)[kk]==0){REAL(infectionCandidate)[kk] =
			runif(minimumLikelyInfectionTime,REAL(latencyCandidate)[kk]-minLatPeriod);}
		/* for first infection : infected after end of latency of initial case */
		else if(kk == 1 && 
			REAL(infectionTimes)[kk+1] > (REAL(latencyCandidate)[kk]-minLatPeriod)){
				REAL(infectionCandidate)[kk] = runif(REAL(latencyTimes)[0],
					(REAL(latencyCandidate)[kk]-minLatPeriod));}
		else if(kk == 1 && 
			REAL(infectionTimes)[kk+1] <= (REAL(latencyCandidate)[kk]-minLatPeriod)){
				REAL(infectionCandidate)[kk] = runif(REAL(latencyTimes)[0],
					REAL(infectionTimes)[kk+1]);}
		/* for last removal */
		else if(kk == nRemoved-1){REAL(infectionCandidate)[kk] =
			runif(REAL(infectionTimes)[kk-1],(REAL(latencyCandidate)[kk]-minLatPeriod));}
		/* for all other cases */
		else if(REAL(infectionTimes)[kk+1] > (REAL(latencyCandidate)[kk]-minLatPeriod)){
			REAL(infectionCandidate)[kk]=
			runif(REAL(infectionTimes)[kk-1],(REAL(latencyCandidate)[kk]-minLatPeriod));}
		else{REAL(infectionCandidate)[kk] =
			runif(REAL(infectionTimes)[kk-1],REAL(infectionTimes)[kk+1]);}
		/* Check for consistency with SEIR models: if not repeat until consistent */
		consistent=consistent_SEIR(REAL(infectionCandidate),REAL(latencyCandidate),
			REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
			REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(LL),INTEGER(II));
		}
		expLikelihood_SEIR(REAL(parameters),REAL(infectionCandidate),REAL(latencyCandidate),
			REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
			&sumSI, &sumLatents, &sumInfectious, &likelihood,
			REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(LL),INTEGER(II));
		newLkhood = likelihood;
		logR = (newLkhood-oldLkhood);
		if(log(unif_rand()) < logR){
			REAL(infectionTimes)[kk] = REAL(infectionCandidate)[kk];
			REAL(latencyTimes)[kk] = REAL(latencyCandidate)[kk];
			++acceptRate;
		}
		REAL(latencyCandidate)[kk] = REAL(latencyTimes)[kk];
		REAL(infectionCandidate)[kk] = REAL(infectionTimes)[kk];
		/* update latency/infection times for index case */
		REAL(latencyTimes)[0] = REAL(infectionTimes)[1]-
			rexp(1/(latRate+remRate+thetaprior[0]));
		REAL(latencyCandidate)[0] = REAL(latencyTimes)[0];
		REAL(infectionTimes)[0] = REAL(latencyTimes)[0]-
			rexp(1/(infRate*INTEGER(N)[0]+latRate+thetaprior[1]));
		REAL(infectionCandidate)[0] = REAL(infectionTimes)[0];
		/*
		REAL(latencyTimes)[0] = REAL(infectionTimes)[1]-
			trexponential(minLatPeriod,maxLatPeriod,(latRate+remRate+thetaprior[0]));
		REAL(latencyCandidate)[0] = REAL(latencyTimes)[0];
		REAL(infectionTimes)[0] = REAL(latencyTimes)[0]-
			trexponential(minInfPeriod,maxInfPeriod,(infRate*INTEGER(N)[0]+latRate+thetaprior[1]));
		REAL(infectionCandidate)[0] = REAL(infectionTimes)[0];
		*/
		/* update latency/infection times for index case */
		expLikelihood_SEIR(REAL(parameters),REAL(infectionTimes),REAL(latencyTimes),
			REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
			&sumSI, &sumLatents, &sumInfectious, &likelihood,
			REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(LL),INTEGER(II));
		oldLkhood = likelihood;
		kk = ceil(INTEGER(bayesReps)[0]/100);
		ll = ceil(INTEGER(bayesReps)[0]/ 10);
		if(verbose){
			if((ii % kk) == 0){Rprintf(".");}
			if((ii % ll) == 0){Rprintf("   %d\n",ii);}
		}
		if((ii >= (INTEGER(bayesStart)[0])) &&
			((ii-INTEGER(bayesStart)[0]) % INTEGER(bayesThin)[0] == 0)){
			ll = (ii - (INTEGER(bayesStart)[0]))/INTEGER(bayesThin)[0];
			REAL(logLikelihood)[ll] = likelihood;
			REAL(infRateSEIR)[ll] = infRate;
			REAL(latRateSEIR)[ll] = latRate;
			REAL(remRateSEIR)[ll] = remRate;
			for(jj = 0; jj < nRemoved; ++jj){
				REAL(infTimes)[jj] += REAL(infectionTimes)[jj];
				REAL(latTimes)[jj] += REAL(latencyTimes)[jj];
			}
		}
	}
	PutRNGstate(); /* after using random number generators. */
	for(jj = 0; jj < nRemoved; ++jj){
		REAL(infTimes)[jj] = REAL(infTimes)[jj]/INTEGER(bayesOut)[0];
		REAL(latTimes)[jj] = REAL(latTimes)[jj]/INTEGER(bayesOut)[0];
	}
	/* Print infection times and removal times at last iteration */
	if(verbose){
		for(jj = 0; jj < nRemoved; ++jj){
			Rprintf("%2d  %8.4f   %8.4f   %2.0f\n",jj,
				REAL(infTimes)[jj],REAL(latTimes)[jj],REAL(removalTimes)[jj]);
		}
	}
	PROTECT(retParameters = NEW_LIST(5));
	++nProtected;
	PROTECT(acceptanceRate = allocVector(INTSXP,1));
	++nProtected;
	INTEGER(acceptanceRate)[0] = acceptRate;
	PROTECT(parNames = allocVector(STRSXP,5));
	++nProtected;
	SET_STRING_ELT(parNames, 0, mkChar("logLikelihood"));
	SET_STRING_ELT(parNames, 1, mkChar("infRateSEIR"));
	SET_STRING_ELT(parNames, 2, mkChar("latRateSEIR"));
	SET_STRING_ELT(parNames, 3, mkChar("remRateSEIR"));
	SET_STRING_ELT(parNames, 4, mkChar("acceptanceRate"));
	setAttrib(retParameters, R_NamesSymbol,parNames);
	
	SET_ELEMENT(retParameters, 0, logLikelihood);
	SET_ELEMENT(retParameters, 1, infRateSEIR);
	SET_ELEMENT(retParameters, 2, latRateSEIR);
	SET_ELEMENT(retParameters, 3, remRateSEIR);
	SET_ELEMENT(retParameters, 4, acceptanceRate);
	/*
	SET_ELEMENT(retParameters, 4, initialLat);
	SET_ELEMENT(retParameters, 5, initialInf);
	SET_ELEMENT(retParameters, 6, timeInfected);
	SET_ELEMENT(retParameters, 7, timeLatent);
	*/
	UNPROTECT(nProtected);
	return(retParameters);
}

/*
Susceptible-Latent-Infectious-Removed MCMC analysis:
	. Exponentially distributed latent periods
	. Exponentially distributed infectiousness periods
*/
void expLkhood_SEIR(double *parameters, double *infectionTimes,
	double *endLatency, double *removalTimes,
	int *N, int *nInfected, int *nRemoved, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *LL, int *II)
{
	int i, k = 0, initialInfective=0, nEvents;
	double sumLogBeta=0, sumLogInfections=0, sumLatent=0, sumInfectious=0, sumBetaSI=0;
	nEvents = *nInfected+ *nInfected+ *nRemoved;
	for(i = 0; i < *nInfected; ++i){
		allTimes[(i*3)] = infectionTimes[i];
		indicator[(i*3)] = 3;
		allTimes[(i*3)+1] = endLatency[i];
		indicator[(i*3)+1] = 2;
		allTimes[(i*3)+2] = removalTimes[i];
		indicator[(i*3)+2] = 1;
		if(removalTimes[i]==0){++initialInfective;}
	}
	rsort_with_index(allTimes,indicator,nEvents);
	SS[0] = *N+initialInfective;
	LL[0] = 0;
	II[0] = 0;
	for(i = 1; i < (nEvents+1); ++i){
		if(indicator[(i-1)] == 3){
			SS[i] = SS[(i-1)]-1;
			LL[i] = LL[(i-1)]+1;
			II[i] = II[(i-1)];
		}
		else if(indicator[(i-1)] == 2){
			SS[i] = SS[(i-1)];
			LL[i] = LL[(i-1)]-1;
			II[i] = II[(i-1)]+1;
		}
		else{
			SS[i] = SS[(i-1)];
			LL[i] = LL[(i-1)];
			II[i] = II[(i-1)]-1;
		}
	}
	/*sumLogBeta=(*nInfected-initialInfective)*log(parameters[0]);*/
	for(i = 1; i < nEvents; ++i){
		if(allTimes[i] != allTimes[i-1]){k = i;}
		sumBetaSI+=parameters[0]*II[i]*SS[i]*(allTimes[i]-allTimes[(i-1)]);
		if(indicator[i] == 1 && II[k] != 0){sumLogInfections+=log(II[k]);}
		if(indicator[i] == 2 && LL[k] != 0){sumLogInfections+=log(LL[k]);}
		if(indicator[i] == 3 && II[k] != 0){sumLogInfections+=log(parameters[0]*SS[k]*II[k]);}
	}
	for(i = 0; i < *nRemoved; ++i){
		sumLatent+=dexp((endLatency[i]-infectionTimes[i]),1/parameters[1],TRUE);
		sumInfectious+=dexp((removalTimes[i]-endLatency[i]),1/parameters[2],TRUE);
	}
	*likelihood=sumLogBeta+sumLogInfections-sumBetaSI+sumLatent+sumInfectious;
}

