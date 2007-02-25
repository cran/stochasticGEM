#include "arms.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include "staticVariables.h"

/*
Susceptible-Infectious-Removed MCMC analysis:
	. Gamma distributed infectiousness periods
*/
SEXP gammaMH_SIRars(SEXP N, SEXP removalTimes, SEXP otherParameters, SEXP priorValues, SEXP initialValues,
		SEXP bayesReps, SEXP bayesStart, SEXP bayesThin, SEXP bayesOut){
	/* Declarations  */
	int ii, jj, kk, ll, nInfected, nRemoved, nProtected=0, initialInfected;
	SEXP infRateSIR, infScaleSIR, infShapeSIR, logLikelihood;/*, timeInfected, timeDim, initialInf ; */
	SEXP parameters, infectionTimes, candidateTimes, infectedBeforeDay;
	SEXP allTimes, indicator, SS, II;
	double infRate, infScale, infShape, oldLkhood, newLkhood, minimumLikelyInfectionTime;/*starting values */
	double infRatePrior[2], infScalePrior[2], infShapePrior[2], thetaprior[2];	 /* priors values */
	double sumSI, sumDurationInfectious, likelihood, logR;
	AuxStruct ext;
	int acceptRate[2], consistent=0, verbose, missingInfectionTimes;
	SEXP retParameters, parNames, acceptanceRate;
	SEXP infTimes;
	/* define ARS parameters */
	int err, ninit=10, dometrop=0;
	double xl=DOUBLE_EPS, xr=10, mlShape=1;
	/* Code  */
	initialInfected = INTEGER(getListElement(otherParameters, "initialInfected"))[0];
	verbose = INTEGER(getListElement(otherParameters, "verbose"))[0];
	missingInfectionTimes = INTEGER(getListElement(otherParameters, "missingInfectionTimes"))[0];
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
	PROTECT(infRateSIR = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(infScaleSIR = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(infShapeSIR = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(logLikelihood = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	/*
	PROTECT(timeInfected = allocVector(REALSXP, nRemoved * INTEGER(bayesOut)[0]));
	++nProtected;
	PROTECT(timeDim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(timeDim)[0] = nRemoved;
	INTEGER(timeDim)[1] = INTEGER(bayesOut)[0];
	setAttrib(timeInfected, R_DimSymbol, timeDim);
	PROTECT(initialInf = allocVector(REALSXP, INTEGER(bayesOut)[0]));
	++nProtected;
	*/
	PROTECT(parameters = allocVector(REALSXP,2));
	++nProtected;
	PROTECT(infectionTimes = allocVector(REALSXP,nRemoved));
	++nProtected;
	PROTECT(candidateTimes = allocVector(REALSXP,nRemoved));
	++nProtected;
	PROTECT(infectedBeforeDay = allocVector(REALSXP,nRemoved));
	++nProtected;
	PROTECT(infTimes = allocVector(REALSXP,nRemoved));
	++nProtected;
	for(jj = 0; jj < nRemoved; ++jj){
		REAL(infectionTimes)[jj] = REAL(getListElement(initialValues, "infectionTimes"))[jj];
		REAL(candidateTimes)[jj] = REAL(infectionTimes)[jj];
		REAL(infectedBeforeDay)[jj] = REAL(getListElement(otherParameters, "infectedBeforeDay"))[jj];
		REAL(infTimes)[jj] = 0;
	}
	nInfected = LENGTH(infectionTimes);
	PROTECT(allTimes = allocVector(REALSXP,nRemoved+nInfected));
	++nProtected;
	PROTECT(indicator = allocVector(INTSXP,nRemoved+nInfected));
	++nProtected;
	PROTECT(SS = allocVector(INTSXP,nRemoved+nInfected+1));
	++nProtected;
	PROTECT(II = allocVector(INTSXP,nRemoved+nInfected+1));
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
	REAL(parameters)[1] = infScale;
	REAL(parameters)[2] = infShape;
	gammaLikelihood_SIR(REAL(parameters),REAL(infectionTimes),
		REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
		&sumSI, &sumDurationInfectious, &likelihood,
		REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(II));
	oldLkhood = likelihood;
	ext = (AuxStruct) R_alloc(1, sizeof(aux_struct));
	for(ii = 0; ii < 2; ++ii){
		acceptRate[ii] = 0;}
	GetRNGstate(); /* should be before a call to a random number generator */
	for(ii = 1; ii <= INTEGER(bayesReps)[0]; ++ii){
		infRate = rgamma(nInfected-1+infRatePrior[0],1/(sumSI+infRatePrior[1]));/* update infRate */
		infScale = rgamma(nRemoved*infShape+infScalePrior[0],1/(sumDurationInfectious+infScalePrior[1]));
		/* Update shape parameter */
		ext->startTimes = REAL(infectionTimes);
		ext->endTimes = REAL(removalTimes);
		ext->nRemoved = nRemoved;
		ext->shapePrior = infShapePrior;
		ext->scale = infScale;
		/* Adaptive rejection sampling */
		err = arms_simple(ninit,&xl,&xr,gammaShapelkhoodARS,ext,dometrop,
			&infShape,&mlShape);
		infShape = mlShape;
		/* Adaptive rejection sampling */
		REAL(parameters)[0] = infRate;
		REAL(parameters)[1] = infScale;
		REAL(parameters)[2] = infShape;
		if(missingInfectionTimes){
			gammaLikelihood_SIR(REAL(parameters),REAL(infectionTimes),
				REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
				&sumSI, &sumDurationInfectious, &likelihood,
				REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(II));
			oldLkhood = likelihood;
			kk = ceil(unif_rand()*(nRemoved-1)); /* initial infection time excluded */
			consistent=0;
			if(kk == nRemoved-1){
				REAL(candidateTimes)[kk] =
					runif(REAL(infectionTimes)[kk-1], REAL(infectedBeforeDay)[kk]);}
			else if((REAL(infectionTimes)[kk+1] > REAL(infectedBeforeDay)[kk])){
				REAL(candidateTimes)[kk] =
					runif(REAL(infectionTimes)[kk-1], REAL(infectedBeforeDay)[kk]);}
			else{REAL(candidateTimes)[kk] =
					runif(REAL(infectionTimes)[kk-1], REAL(infectionTimes)[kk+1]);}
			gammaLikelihood_SIR(REAL(parameters),REAL(candidateTimes),
				REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
				&sumSI, &sumDurationInfectious, &likelihood,
				REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(II));
			newLkhood = likelihood;
			logR = (newLkhood-oldLkhood);
			if(log(unif_rand()) < logR){
				REAL(infectionTimes)[kk] = REAL(candidateTimes)[kk];
				++acceptRate[0];
			}
			REAL(infectionTimes)[0] = REAL(infectionTimes)[1]/* update initial infection time */
				-rgamma(infShape+thetaprior[0]-1,1/(infRate*INTEGER(N)[0]+infScale+thetaprior[1]));
			REAL(candidateTimes)[0] = REAL(infectionTimes)[0];  /* update candidate */
			REAL(candidateTimes)[kk] = REAL(infectionTimes)[kk];/* times */
		}
		gammaLikelihood_SIR(REAL(parameters),REAL(infectionTimes),
			REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
			&sumSI, &sumDurationInfectious, &likelihood,
			REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(II));
		oldLkhood = likelihood;
		kk = ceil(INTEGER(bayesReps)[0]/100);
		ll = ceil(INTEGER(bayesReps)[0]/ 10);
		if(verbose == 1){
			if((ii % kk) == 0){Rprintf(".");}
			if((ii % ll) == 0){Rprintf("   %d\n",ii);}
		}
		if((ii >= (INTEGER(bayesStart)[0])) &&
			((ii-INTEGER(bayesStart)[0]) % INTEGER(bayesThin)[0] == 0)){
			ll = (ii - (INTEGER(bayesStart)[0]))/INTEGER(bayesThin)[0];
			/*REAL(initialInf)[ll] = REAL(infectionTimes)[0]; */
			REAL(logLikelihood)[ll] = likelihood;
			REAL(infRateSIR)[ll] = infRate;
			REAL(infScaleSIR)[ll] = infScale;
			REAL(infShapeSIR)[ll] = infShape;
			for(jj = 0; jj < nRemoved; ++jj){
				REAL(infTimes)[jj] += REAL(infectionTimes)[jj];
			}
			/*
			for(jj = 0; jj < nRemoved; ++jj){
				REAL(timeInfected)[(nRemoved*ll+jj)] = REAL(infectionTimes)[jj];
			}
			*/
		}
	}
	PutRNGstate(); /* after using random number generators. */
	for(jj = 0; jj < nRemoved; ++jj){
		REAL(infTimes)[jj] = REAL(infTimes)[jj]/INTEGER(bayesOut)[0];
	}
	if(verbose){
		for(jj = 0; jj < nRemoved; ++jj){
			Rprintf("%2d  %8.4f   %2.0f\n",jj,
				REAL(infTimes)[jj],REAL(removalTimes)[jj]);
		}
	}
	PROTECT(retParameters = NEW_LIST(6));
	++nProtected;
	PROTECT(acceptanceRate = allocVector(INTSXP,2));
	++nProtected;
	for(ii = 0; ii < 2; ++ii){
		INTEGER(acceptanceRate)[ii] = acceptRate[ii];}
	PROTECT(parNames = allocVector(STRSXP,6));
	++nProtected;
	SET_STRING_ELT(parNames, 0, mkChar("logLikelihood"));
	SET_STRING_ELT(parNames, 1, mkChar("infRateSIR"));
	SET_STRING_ELT(parNames, 2, mkChar("infScaleSIR"));
	SET_STRING_ELT(parNames, 3, mkChar("infShapeSIR"));
	SET_STRING_ELT(parNames, 4, mkChar("infectionTimes"));
	SET_STRING_ELT(parNames, 5, mkChar("acceptanceRate"));
	setAttrib(retParameters, R_NamesSymbol,parNames);
	
	SET_ELEMENT(retParameters, 0, logLikelihood);
	SET_ELEMENT(retParameters, 1, infRateSIR);
	SET_ELEMENT(retParameters, 2, infScaleSIR);
	SET_ELEMENT(retParameters, 3, infShapeSIR);
	SET_ELEMENT(retParameters, 4, infTimes);
	SET_ELEMENT(retParameters, 5, acceptanceRate);
	UNPROTECT(nProtected);
	return(retParameters);
}
