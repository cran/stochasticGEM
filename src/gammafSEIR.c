#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include "staticVariables.h"

/*
Susceptible-Latent-Infectious-Removed MCMC analysis:
	. Fixed latent periods
	. Gamma distributed infectiousness periods
*/
static void gammaLikelihood_fSEIR(double *parameters, double *infectionTimes,
	double *endLatency, double *removalTimes, int *N, int *nInfected, int *nRemoved, 
	double *sumSI, double *sumInfectious, double *likelihood,
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
		if(indicator[i] == 3 && II[k] != 0){sumLogInfections+=log(parameters[0]*SS[k]*II[k]);}
		*sumSI+=SS[i]*II[i]*(allTimes[i]-allTimes[(i-1)]);
	}
	*sumInfectious=0;	
	for(i = 0; i < *nInfected; ++i){
		sumLatent+=dexp((endLatency[i]-infectionTimes[i]),1/parameters[1],TRUE);
		sumInfection+=dgamma((removalTimes[i]-endLatency[i]),parameters[3],1/parameters[2],TRUE);
		*sumInfectious+=(removalTimes[i]-endLatency[i]);
	}
	*likelihood=sumLogBeta+sumLogInfections-sumBetaSI+sumLatent+sumInfection;
}

/*
Susceptible-Latent-Infectious-Removed MCMC analysis:
	. Fixed latent periods
	. Gamma distributed infectiousness periods
*/
SEXP gammaMH_fSEIR(SEXP N, SEXP removalTimes, SEXP otherParameters, SEXP priorValues,
	SEXP initialValues, SEXP bayesReps, SEXP bayesStart, SEXP bayesThin, SEXP bayesOut){
	/* Declarations  */
	int ii, jj, kk, ll, nInfected, nRemoved, nEvents, nProtected=0, initialInfected;
	SEXP infRateSEIR, infScaleSEIR, infShapeSEIR, logLikelihood;
	SEXP parameters;
	SEXP infectionTimes, infectionCandidate, infectedBeforeDay;
	SEXP latencyTimes, latencyCandidate;
	SEXP allTimes, indicator, SS, II, LL;
	double infRate, latRate, infScale, infShape, oldLkhood, newLkhood, minimumLikelyInfectionTime,fixedLatent;
	double infRatePrior[2], infScalePrior[2], infShapePrior[2], thetaprior[2];/* priors values */
	double sumSI, sumInfectious, likelihood,logR,infPeriod,minInfPeriod,maxInfPeriod;	
	AuxStruct ext;
	int acceptRate=0, consistent=0, verbose, updateIniInfGibbs;
	SEXP retParameters, parNames, acceptanceRate;
	SEXP infTimes, latTimes;
	double slice,lowerLimit=0,upperLimit=0;
	/* Code  */
	GetRNGstate(); /* should be before a call to a random number generator */
	fixedLatent = REAL(getListElement(otherParameters, "fixedLatent"))[0];
	latRate = 1/fixedLatent;
	minInfPeriod = REAL(getListElement(otherParameters, "minInfPeriod"))[0];
	initialInfected = INTEGER(getListElement(otherParameters, "initialInfected"))[0];
	verbose = INTEGER(getListElement(otherParameters, "verbose"))[0];
	updateIniInfGibbs = INTEGER(getListElement(otherParameters, "updateIniInfGibbs"))[0];
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
	PROTECT(infScaleSEIR = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(infShapeSEIR = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(logLikelihood = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(parameters = allocVector(REALSXP,4));
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
		REAL(latencyTimes)[jj] = REAL(getListElement(initialValues, "latencyTimes"))[jj];
		REAL(infectionTimes)[jj] = REAL(latencyTimes)[jj]-fixedLatent;
		REAL(latencyCandidate)[jj] = REAL(latencyTimes)[jj];
		REAL(infectionCandidate)[jj] = REAL(infectionTimes)[jj];
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
	infRate = REAL(getListElement(initialValues, "infectionRate"))[0];
	infScale = REAL(getListElement(initialValues, "infectiousScale"))[0];
	infShape = REAL(getListElement(initialValues, "infectiousShape"))[0];
	minimumLikelyInfectionTime = REAL(getListElement(otherParameters, "minimumLikelyInfectionTime"))[0];
	for(ii = 0; ii < 2; ++ii){
		infRatePrior[ii] = REAL(getListElement(priorValues, "infectionRate"))[ii];
		infScalePrior[ii] = REAL(getListElement(priorValues, "infectiousScale"))[ii];
		infShapePrior[ii] = REAL(getListElement(priorValues, "infectiousShape"))[ii];
		thetaprior[ii] = REAL(getListElement(priorValues, "theta"))[ii];
	}
	REAL(parameters)[0] = infRate;
	REAL(parameters)[1] = latRate;
	REAL(parameters)[2] = infScale;
	REAL(parameters)[3] = infShape;
	gammaLikelihood_fSEIR(REAL(parameters),REAL(infectionTimes),REAL(latencyTimes),
		REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
		&sumSI, &sumInfectious, &likelihood,REAL(allTimes),INTEGER(indicator),
		INTEGER(SS),INTEGER(LL),INTEGER(II));
	oldLkhood = likelihood;
	/* Update shape parameters */
	ext = (AuxStruct) R_alloc(1, sizeof(aux_struct));
	for(ii = 1; ii <= INTEGER(bayesReps)[0]; ++ii){
		/*
		Rprintf("%5d  %7.4f  %7.4f  %7.4f\n",ii,infRate,infScale,infShape);
		*/
		infRate = rgamma(nInfected-1+infRatePrior[0],1/(sumSI+infRatePrior[1]));/* update infRate */
		infScale = rgamma(nRemoved*infShape+infScalePrior[0],1/(sumInfectious+infScalePrior[1]));
		/* Update shape parameter */
		ext->startTimes = REAL(latencyTimes);
		ext->endTimes = REAL(removalTimes);
		ext->nRemoved = nRemoved;
		ext->shapePrior = infShapePrior;
		ext->scale = infScale;
		slice = gammaShapelkhoodARS(infShape,ext) - exp_rand();
		lowerLimit=1;upperLimit=100;
		steppingOut(slice,infShape,ext,gammaShapelkhoodARS,&lowerLimit,&upperLimit);
		steppingOutShrinkage(slice,&infShape,ext,gammaShapelkhoodARS,lowerLimit,upperLimit);
		/* Update shape parameter */
		REAL(parameters)[0] = infRate;
		REAL(parameters)[1] = latRate;
		REAL(parameters)[2] = infScale;
		REAL(parameters)[3] = infShape;
		gammaLikelihood_fSEIR(REAL(parameters),REAL(infectionTimes),REAL(latencyTimes),
			REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
			&sumSI, &sumInfectious, &likelihood,REAL(allTimes),INTEGER(indicator),
			INTEGER(SS),INTEGER(LL),INTEGER(II));
		oldLkhood = likelihood;
		if(updateIniInfGibbs){kk = ceil(unif_rand()*(nRemoved-1));}
		else{kk = floor(unif_rand()*nRemoved);}
		consistent=0;
		if(kk == 0){maxInfPeriod=REAL(removalTimes)[kk]-(minimumLikelyInfectionTime+fixedLatent);}
		else{maxInfPeriod=REAL(removalTimes)[kk]-(REAL(latencyTimes)[0]+fixedLatent);}
		while(consistent==0){
			infPeriod = runif(minInfPeriod,maxInfPeriod);
			/*REAL(latencyCandidate)[kk] = REAL(removalTimes)[kk]-infPeriod;*/
			if(REAL(removalTimes)[kk]==0){REAL(latencyCandidate)[kk] =
				runif((minimumLikelyInfectionTime+fixedLatent),REAL(infectionTimes)[kk+1]);}
			else if(kk == 1 && REAL(latencyTimes)[kk+1] > (REAL(removalTimes)[kk]-minInfPeriod)){
				REAL(latencyCandidate)[kk] =
					runif((REAL(latencyTimes)[kk-1]+fixedLatent),
						(REAL(removalTimes)[kk]-minInfPeriod));}
			else if(kk == 1 && REAL(latencyTimes)[kk+1] <= (REAL(removalTimes)[kk]-minInfPeriod)){
				REAL(latencyCandidate)[kk] =
					runif((REAL(latencyTimes)[kk-1]+fixedLatent),REAL(latencyTimes)[kk+1]);}
			else if(kk == nRemoved-1){REAL(latencyCandidate)[kk] =
				runif(REAL(latencyTimes)[kk-1],(REAL(removalTimes)[kk]-minInfPeriod));}
			else if(REAL(latencyTimes)[kk+1] > (REAL(removalTimes)[kk]-minInfPeriod)){
				REAL(latencyCandidate)[kk]=
				runif(REAL(latencyTimes)[kk-1],(REAL(removalTimes)[kk]-minInfPeriod));}
			else{REAL(latencyCandidate)[kk] =
				runif(REAL(latencyTimes)[kk-1],REAL(latencyTimes)[kk+1]);}
			REAL(infectionCandidate)[kk] = REAL(latencyCandidate)[kk]-fixedLatent;
			consistent=consistent_SEIR(REAL(infectionCandidate),REAL(latencyCandidate),
				REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
				REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(LL),INTEGER(II));
		}
		gammaLikelihood_fSEIR(REAL(parameters),REAL(infectionCandidate),REAL(latencyCandidate),
			REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
			&sumSI, &sumInfectious, &likelihood,REAL(allTimes),INTEGER(indicator),
			INTEGER(SS),INTEGER(LL),INTEGER(II));
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
		if(updateIniInfGibbs){
			REAL(latencyTimes)[0]=0;
			while(REAL(latencyTimes)[0] > REAL(infectionTimes)[1]){
				REAL(latencyTimes)[0] = -rgamma(infShape-1+thetaprior[0],
					1/(latRate+infScale+thetaprior[1]));
			}
			REAL(latencyCandidate)[0] = REAL(latencyTimes)[0];
			REAL(infectionTimes)[0] = REAL(latencyTimes)[0]-fixedLatent;
			REAL(infectionCandidate)[0] = REAL(infectionTimes)[0];
		}
		/* update latency/infection times for index case */
		gammaLikelihood_fSEIR(REAL(parameters),REAL(infectionTimes),REAL(latencyTimes),
			REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
			&sumSI, &sumInfectious, &likelihood,REAL(allTimes),INTEGER(indicator),
			INTEGER(SS),INTEGER(LL),INTEGER(II));
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
			REAL(infScaleSEIR)[ll] = infScale;
			REAL(infShapeSEIR)[ll] = infShape;
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
	PROTECT(retParameters = NEW_LIST(7));
	++nProtected;
	PROTECT(acceptanceRate = allocVector(INTSXP,1));
	++nProtected;
	INTEGER(acceptanceRate)[0] = acceptRate;
	PROTECT(parNames = allocVector(STRSXP,7));
	++nProtected;
	SET_STRING_ELT(parNames, 0, mkChar("logLikelihood"));
	SET_STRING_ELT(parNames, 1, mkChar("infRateSEIR"));
	SET_STRING_ELT(parNames, 2, mkChar("infScaleSEIR"));
	SET_STRING_ELT(parNames, 3, mkChar("infShapeSEIR"));
	SET_STRING_ELT(parNames, 4, mkChar("infectionTimes"));
	SET_STRING_ELT(parNames, 5, mkChar("latencyTimes"));
	SET_STRING_ELT(parNames, 6, mkChar("acceptanceRate"));
	setAttrib(retParameters, R_NamesSymbol,parNames);
	SET_ELEMENT(retParameters, 0, logLikelihood);
	SET_ELEMENT(retParameters, 1, infRateSEIR);
	SET_ELEMENT(retParameters, 2, infScaleSEIR);
	SET_ELEMENT(retParameters, 3, infShapeSEIR);
	SET_ELEMENT(retParameters, 4, infTimes);
	SET_ELEMENT(retParameters, 5, latTimes);
	SET_ELEMENT(retParameters, 6, acceptanceRate);
	UNPROTECT(nProtected);
	return(retParameters);
}
