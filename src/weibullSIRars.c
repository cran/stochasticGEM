#include "arms.h"
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include "staticVariables.h"

/*
Susceptible-Infectious-Removed MCMC analysis:
	. Infectiousness periods are distributed as a Weibull
*/
SEXP weibullMH_SIRars(SEXP N, SEXP removalTimes, SEXP otherParameters, SEXP priorValues, SEXP initialValues,
		SEXP bayesReps, SEXP bayesStart, SEXP bayesThin, SEXP bayesOut){
	/* Declarations  */
	int ii, jj, kk, ll, nInfected, nRemoved, nProtected=0, initialInfected;
	SEXP infRateSIR, infScaleSIR, infShapeSIR, logLikelihood;/*, timeInfected, timeDim, initialInf ; */
	SEXP parameters, infectionTimes, candidateTimes, infectedBeforeDay;
	SEXP allTimes, indicator, SS, II;
	double infRate, infScale, infShape, oldLkhood, newLkhood, minimumLikelyInfectionTime;/*starting values */
	double infRatePrior[2], infScalePrior[2], infShapePrior[2];	 /* priors values */
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
	}
	REAL(parameters)[0] = infRate;
	REAL(parameters)[1] = infScale;
	REAL(parameters)[2] = infShape;
	weibullLikelihood_SIR(REAL(parameters),REAL(infectionTimes),
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
		infScale = rgamma(nRemoved+infScalePrior[0],1/(sumDurationInfectious+infScalePrior[1]));
		/* Update nu */
		ext->startTimes = REAL(infectionTimes);
		ext->endTimes = REAL(removalTimes);
		ext->nRemoved = nRemoved;
		ext->shapePrior = infShapePrior;
		ext->scale = infScale;
		/* Adaptive rejection sampling */
		err = arms_simple(ninit,&xl,&xr,weiShapelkhoodARS,ext,dometrop,
			&infShape,&mlShape);
		infShape = mlShape;
		/* Adaptive rejection sampling */
		/* Update nu */
		REAL(parameters)[0] = infRate;
		REAL(parameters)[1] = infScale;
		REAL(parameters)[2] = infShape;
		if(missingInfectionTimes){
			weibullLikelihood_SIR(REAL(parameters),REAL(infectionTimes),
				REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
				&sumSI, &sumDurationInfectious, &likelihood,
				REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(II));
			oldLkhood = likelihood;
			kk = floor(unif_rand()*nRemoved);/* initial infection time included */
			consistent=0;
			if(REAL(removalTimes)[kk]==0){
				REAL(candidateTimes)[kk] =
					runif(minimumLikelyInfectionTime, REAL(infectionTimes)[kk+1]);}
			else if(kk == nRemoved-1){
				REAL(candidateTimes)[kk] =
					runif(REAL(infectionTimes)[kk-1], REAL(infectedBeforeDay)[kk]);}
			else if((REAL(infectionTimes)[kk+1] > REAL(infectedBeforeDay)[kk])){
				REAL(candidateTimes)[kk] =
					runif(REAL(infectionTimes)[kk-1], REAL(infectedBeforeDay)[kk]);}
			else{REAL(candidateTimes)[kk] =
					runif(REAL(infectionTimes)[kk-1], REAL(infectionTimes)[kk+1]);}
			weibullLikelihood_SIR(REAL(parameters),REAL(candidateTimes),
				REAL(removalTimes), &INTEGER(N)[0], &nInfected, &nRemoved,
				&sumSI, &sumDurationInfectious, &likelihood,
				REAL(allTimes),INTEGER(indicator),INTEGER(SS),INTEGER(II));
			newLkhood = likelihood;
			logR = newLkhood-oldLkhood;			
			if(log(runif(0,1)) < logR){
				REAL(infectionTimes)[kk] = REAL(candidateTimes)[kk];
				++acceptRate[0];
			}
			REAL(candidateTimes)[kk] = REAL(infectionTimes)[kk];
		}
		weibullLikelihood_SIR(REAL(parameters),REAL(infectionTimes),
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
	/*
	SET_ELEMENT(retParameters, 4, initialInf);
	SET_ELEMENT(retParameters, 5, timeInfected);
	*/
	UNPROTECT(nProtected);
	return(retParameters);
}
