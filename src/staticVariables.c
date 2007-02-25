#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>

#define SHAPE_SMALL     1.0e-2
#define SHAPE_LARGE     1.0e+2


SEXP getListElement(SEXP list, char *str)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;
	for (i = 0; i < length(list); i++)
	if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
		elmt = VECTOR_ELT(list, i);
		break;
	}
	return elmt;
}

double trexponential(double lower_cut_off, double upper_cut_off, double rate)
{
	double lower_p,upper_p,r;
	lower_p = pexp(lower_cut_off, 1/(rate),TRUE,FALSE);
	upper_p = pexp(upper_cut_off, 1/(rate),TRUE,FALSE);
	r = lower_cut_off+qexp(unif_rand()*(upper_p-lower_p), 1/(rate),TRUE,FALSE);
	return(r);
}

double trgamma(double lower_cut_off, double upper_cut_off, double shape, double rate)
{
	double lower_p,upper_p,r;
	lower_p = pgamma(lower_cut_off, shape, 1/(rate),TRUE,FALSE);
	upper_p = pgamma(upper_cut_off, shape, 1/(rate),TRUE,FALSE);
	r = lower_cut_off+qgamma(unif_rand()*(upper_p-lower_p), shape, 1/(rate),TRUE,FALSE);
	return(r);
}


/* Define structure used in shape likelihoods */
typedef struct aux_struct
{
	int nRemoved;
	double *startTimes;
	double *endTimes;
	double *shapePrior;
	double scale;
} aux_struct, *AuxStruct;

typedef struct aux_dm
{
	int nEvents;
	int *SS;
	int *II;
	int *ZZ;
	double *allTimes;
	double *Prior;
	double Rate;
} aux_dm, *AuxDM;

/* Check whether new infection time is consistent with an SIR model:       */
/* whenever there is an infection there should be infectious individuals   */
/* just before the infection time. Returns 1 for consistent and 0 otherwise*/
int consistent_SIR(double *infectionTimes, double *removalTimes,
	int *N, int *nInfected, int *nRemoved, double *allTimes,
	int *indicator, int *SS, int *II)
{
	int i,k=0,initialInfective=0, nEvents, consistent=1;
	nEvents = *nInfected+*nRemoved; 
	for(i = 0; i < *nInfected; ++i){
		allTimes[(i*2)] = infectionTimes[i];
		allTimes[(i*2)+1] = removalTimes[i];
		indicator[(i*2)] = 2;
		indicator[(i*2)+1] = 1;
		if(removalTimes[i]==0){++initialInfective;}
	}
	rsort_with_index(allTimes,indicator,nEvents);
	SS[0] = *N+initialInfective;
	II[0] = 0;
	for(i = 1; i < (nEvents+1); ++i){
		if(indicator[(i-1)] == 2){
			SS[i] = SS[(i-1)]-1;
			II[i] = II[(i-1)]+1;}
		else{
			SS[i] = SS[(i-1)];
			II[i] = II[(i-1)]-1;}
	}
	for(i = 1; i < nEvents; ++i){/* "0" is the start of observation */
		if(allTimes[i] != allTimes[i-1]){k = i;}
		if(indicator[i] == 2){consistent*=(II[k]>0);}
	}
	return(consistent);
}

/* Check whether new infection time is consistent with an SEIR model:      */
/* whenever there is an infection there should be infectious individuals   */
/* just before the infection time. Returns 1 for consistent and 0 otherwise*/
int consistent_SEIR(double *infectionTimes, double *endLatency,
	double *removalTimes, int *N, int *nInfected, int *nRemoved, 
	double *allTimes, int *indicator, int *SS, int *LL, int *II)
{
	int i, k = 0, initialInfective=0, nEvents, consistent=1;
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
	/*Rprintf("%d  %d  %d  %d  %f\n",SS[0],LL[0],II[0],indicator[0],allTimes[0]);*/
	for(i = 2; i < nEvents; ++i){/*SEIR starts when 1 individual is infective*/
		if(allTimes[i] != allTimes[i-1]){k = i;}
		if(indicator[i] == 3){consistent*=(II[k]>0);}
	}
	return(consistent);
}

/* Weibull distribution */
double weiShapelkhood(int n, double *par, void *ex)
{
	double *endTimes, *startTimes, *shapePrior, scale, parameter, lkhood;
	int ii, nRemoved;
	AuxStruct AS = (AuxStruct) ex;
	startTimes = AS->startTimes;
	endTimes = AS->endTimes;
	nRemoved = AS->nRemoved;
	shapePrior = AS->shapePrior;
	scale = AS->scale;
	parameter = *par;
	lkhood = (nRemoved+shapePrior[0]-1)*log(parameter)-shapePrior[1]*parameter;
	for(ii = 0; ii < nRemoved; ++ii){
		lkhood+=(parameter-1)*log(endTimes[ii]-startTimes[ii])-scale*
			R_pow((endTimes[ii]-startTimes[ii]),parameter);
	}
	return(-lkhood); /* minimizing the likelihood */
}

void weiFderivShape(int n, double *par, double *gr, void *ex)
{
	double *endTimes, *startTimes, *shapePrior, scale, firstDeriv,parameter;
	int ii, nRemoved;
	AuxStruct AS = (AuxStruct) ex;
	startTimes = AS->startTimes;
	endTimes = AS->endTimes;
	nRemoved = AS->nRemoved;
	shapePrior = AS->shapePrior;
	scale = AS->scale;
	parameter = *par;
	firstDeriv=(nRemoved+shapePrior[0]-1)/parameter-shapePrior[1];
	for(ii = 0; ii < nRemoved; ++ii){
		firstDeriv+=log(endTimes[ii]-startTimes[ii])-
			scale*log(endTimes[ii]-startTimes[ii])*
			R_pow((endTimes[ii]-startTimes[ii]),parameter);
	}
	*gr=-firstDeriv;
}

double weiSderivShape(double shapePost, void *ex)
{
	double *endTimes, *startTimes, *shapePrior, scale, secderiv;
	int ii, nRemoved;
	AuxStruct AS = (AuxStruct) ex;
	startTimes = AS->startTimes;
	endTimes = AS->endTimes;
	nRemoved = AS->nRemoved;
	shapePrior = AS->shapePrior;
	scale = AS->scale;
	secderiv=(-nRemoved-shapePrior[0]+1)/R_pow(shapePost,2);
	for(ii = 0; ii < nRemoved; ++ii){
		secderiv-=scale*R_pow(log(endTimes[ii]-startTimes[ii]),2)*
			R_pow((endTimes[ii]-startTimes[ii]),shapePost);
	}
	return(secderiv);
}

double weiShapelkhoodARS(double parameter, void *ex)
{
	double *endTimes, *startTimes, *shapePrior, scale, lkhood;
	int ii, nRemoved;
	AuxStruct AS = (AuxStruct) ex;
	startTimes = AS->startTimes;
	endTimes = AS->endTimes;
	nRemoved = AS->nRemoved;
	shapePrior = AS->shapePrior;
	scale = AS->scale;
	lkhood = (nRemoved+shapePrior[0]-1)*log(parameter)-shapePrior[1]*parameter;
	for(ii = 0; ii < nRemoved; ++ii){
		lkhood+=(parameter-1)*log(endTimes[ii]-startTimes[ii])-scale*
			R_pow((endTimes[ii]-startTimes[ii]),parameter);
	}
	return(lkhood); /* The likelihood */
}


/* Gamma distribution */
double gammaShapelkhood(int n, double *par, void *ex)
{
	double *endTimes, *startTimes, *shapePrior, scale, parameter, lkhood;
	int ii, nRemoved;
	AuxStruct AS = (AuxStruct) ex;
	startTimes = AS->startTimes;
	endTimes = AS->endTimes;
	nRemoved = AS->nRemoved;
	shapePrior = AS->shapePrior;
	scale = AS->scale;
	parameter = *par;
	lkhood = (nRemoved*parameter)*log(scale)+(shapePrior[1]-1)*log(parameter)-
		nRemoved*lgammafn(parameter)-shapePrior[0]*parameter;
	for(ii = 0; ii < nRemoved; ++ii){
		lkhood+=(parameter-1)*log(endTimes[ii]-startTimes[ii]);
	}
	return(-lkhood); /* minimizing the likelihood */
}

void gammaFderivShape(int n, double *par, double *gr, void *ex)
{
	double *endTimes, *startTimes, *shapePrior, scale, firstDeriv,parameter;
	int ii, nRemoved;
	AuxStruct AS = (AuxStruct) ex;
	startTimes = AS->startTimes;
	endTimes = AS->endTimes;
	nRemoved = AS->nRemoved;
	shapePrior = AS->shapePrior;
	scale = AS->scale;
	parameter = *par;
	firstDeriv=nRemoved*log(scale)+(shapePrior[1]-1)/parameter-
		nRemoved*digamma(parameter)-shapePrior[0];
	for(ii = 0; ii < nRemoved; ++ii){
		firstDeriv+=log(endTimes[ii]-startTimes[ii]);
	}
	*gr=-firstDeriv; /* minimizing the likelihood */
}

double gammaSderivShape(double shapePost, void *ex)
{
	double *endTimes, *startTimes, *shapePrior, scale, secderiv;
	int nRemoved;
	AuxStruct AS = (AuxStruct) ex;
	startTimes = AS->startTimes;
	endTimes = AS->endTimes;
	nRemoved = AS->nRemoved;
	shapePrior = AS->shapePrior;
	scale = AS->scale;
	secderiv=(-shapePrior[1]+1)/R_pow(shapePost,2)-nRemoved*trigamma(shapePost);
	return(secderiv);
}

/* Gamma distribution */
double gammaShapelkhoodARS(double parameter, void *ex)
{
	double *endTimes, *startTimes, *shapePrior, scale, lkhood;
	int ii, nRemoved;
	AuxStruct AS = (AuxStruct) ex;
	startTimes = AS->startTimes;
	endTimes = AS->endTimes;
	nRemoved = AS->nRemoved;
	shapePrior = AS->shapePrior;
	scale = AS->scale;
	lkhood = (nRemoved*parameter)*log(scale)+(shapePrior[1]-1)*log(parameter)-
		nRemoved*lgammafn(parameter)-shapePrior[0]*parameter;
	for(ii = 0; ii < nRemoved; ++ii){
		lkhood+=(parameter-1)*log(endTimes[ii]-startTimes[ii])-
			scale*(endTimes[ii]-startTimes[ii]);
	}
	return(lkhood); /* The likelihood */
}

/*
Incomplete data log-likelihood:
	Susceptible-Infectious-Removed MCMC analysis:
		. Exponentially distributed infectiousness periods
*/
void infRateMixLkhood_SIR(double *parameters, double *infectionTimes, double *removalTimes,
	int *N, int *nInfected, int *nRemoved, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *II)
{
	int i,j=0,k=0,initialInfective=0, nEvents, *infectedAfter;
	double sumLogInfections=0, sumDurationDensity=0, sumBetaSI=0;
	double *totalSI_at_infectionTime,*individualSI_at_infectionTime,*SI; 
	double *forceInfectiousBefore;
	infectedAfter=(int *) R_alloc(*nInfected, sizeof(int));
	totalSI_at_infectionTime=(double *) R_alloc(*nInfected, sizeof(double));
	individualSI_at_infectionTime=(double *) R_alloc(*nInfected, sizeof(double));
	SI=(double *) R_alloc(*nInfected, sizeof(double));
	forceInfectiousBefore=(double *) R_alloc(*nInfected, sizeof(double));
	nEvents = *nInfected+*nRemoved; 
	for(i = 0; i < *nInfected; ++i){
		allTimes[(i*2)] = infectionTimes[i];
		allTimes[(i*2)+1] = removalTimes[i];
		indicator[(i*2)] = 2;
		indicator[(i*2)+1] = 1;
		if(removalTimes[i]==0){++initialInfective;}
	}
	rsort_with_index(allTimes,indicator,nEvents);
	SS[0] = *N+initialInfective;
	II[0] = 0;
	for(i = 1; i < (nEvents+1); ++i){
		if(indicator[(i-1)] == 2){
			SS[i] = SS[(i-1)]-1;
			II[i] = II[(i-1)]+1;}
		else{
			SS[i] = SS[(i-1)];
			II[i] = II[(i-1)]-1;}
	}
	for(i = 0; i < *nRemoved; ++i){totalSI_at_infectionTime[i] = 0;}
	infectedAfter[0]= *nInfected;
	forceInfectiousBefore[0]=0;
	for(i = 1; i < nEvents; ++i){
		if(allTimes[i] != allTimes[i-1]){k = i;}
		if(indicator[i] == 1 && II[k] != 0){sumLogInfections+=log(II[k]);}
		if(indicator[i] == 2){
			++j;/* increment index for infected */
			/*infectedAfter[j]=infectedAfter[(j-1)]-1;*/
			infectedAfter[j]=(*nInfected)-(SS[0]-SS[k]);
			forceInfectiousBefore[j]=II[k]*SS[k];
		}
		totalSI_at_infectionTime[j] += SS[i]*II[i]*(allTimes[i]-allTimes[(i-1)]);
	}
	/* computing the force of infectious pressure on susceptibles  */
	SI[0]=0;
	individualSI_at_infectionTime[0]=0;
	for(i = 1; i < *nRemoved; ++i){
		/* individual contribution between infection times */
		individualSI_at_infectionTime[i] = totalSI_at_infectionTime[i]/infectedAfter[i];
		/* individual contribution at their infection time */
		SI[i] = individualSI_at_infectionTime[i]+SI[(i-1)];
		if(forceInfectiousBefore[i] != 0){
		sumBetaSI+=log(parameters[0]*parameters[1]*forceInfectiousBefore[i]*exp(-parameters[1]*SI[i])+
			(1-parameters[0])*parameters[2]*forceInfectiousBefore[i]*exp(-parameters[2]*SI[i]));
		}
	}
	/* computing the force of infectious on susceptibles  */
	for(i = 0; i < *nRemoved; ++i){
		sumDurationDensity+=dexp((removalTimes[i]-infectionTimes[i]),1/parameters[3],TRUE);
	}
	*likelihood=sumLogInfections+sumBetaSI+sumDurationDensity;
}


/*
Incomplete data log-likelihood:
	Susceptible-Infectious-Removed MCMC analysis:
		. Exponentially distributed infectiousness periods
*/
void remRateMixLkhood_SIR(double *parameters, double *infectionTimes, double *removalTimes,
	int *N, int *nInfected, int *nRemoved, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *II)
{
	int i,j=0,k=0,initialInfective=0, nEvents, *InfectiousBefore;
	double sumLogBeta=0, sumLogInfections=0, sumDurationDensity=0, sumBetaSI=0;
	InfectiousBefore=(int *) R_alloc(*nRemoved, sizeof(int));
	nEvents = *nInfected+*nRemoved; 
	for(i = 0; i < *nInfected; ++i){
		allTimes[(i*2)] = infectionTimes[i];
		allTimes[(i*2)+1] = removalTimes[i];
		indicator[(i*2)] = 2;
		indicator[(i*2)+1] = 1;
		if(removalTimes[i]==0){++initialInfective;}
	}
	rsort_with_index(allTimes,indicator,nEvents);
	SS[0] = *N+initialInfective;
	II[0] = 0;
	for(i = 1; i < (nEvents+1); ++i){
		if(indicator[(i-1)] == 2){
			SS[i] = SS[(i-1)]-1;
			II[i] = II[(i-1)]+1;}
		else{
			SS[i] = SS[(i-1)];
			II[i] = II[(i-1)]-1;}
	}
	/*sumLogBeta=(*nInfected-initialInfective)*log(parameters[1]);*/
	for(i = 1; i < nEvents; ++i){/* "0" is the start of observation */
		if(allTimes[i] != allTimes[i-1]){k = i;}
		sumBetaSI+=parameters[1]*II[i]*SS[i]*(allTimes[i]-allTimes[(i-1)]);
		if(indicator[i] == 1){
			InfectiousBefore[j]=II[k];
			++j;/* increment index for removed */
		}
		if(indicator[i] == 2 && II[i] != 0){sumLogInfections+=log(parameters[1]*SS[k]*II[k]);}
	}
	for(i = 0; i < *nRemoved; ++i){
		sumDurationDensity+=log(parameters[0]*parameters[2]*InfectiousBefore[i]*
				exp(-parameters[2]*(removalTimes[i]-infectionTimes[i]))+
			(1-parameters[0])*parameters[3]*InfectiousBefore[i]*
				exp(-parameters[3]*(removalTimes[i]-infectionTimes[i])));
	}
	*likelihood=sumLogBeta+sumLogInfections-sumBetaSI+sumDurationDensity;
}

void steppingOut(double slice, double shape, void *ext,  double (*likefun)(double,void *ext),
		double *lowerLimit, double *upperLimit){
	int j,k,m = 10;
	double w = 1.0,y;
	*lowerLimit = shape - w * unif_rand();
	*upperLimit = *lowerLimit + w;
	if(*lowerLimit<SHAPE_SMALL){*lowerLimit = SHAPE_SMALL;}
	if(*upperLimit>SHAPE_LARGE){*upperLimit = SHAPE_LARGE;}
	j = floor(m*unif_rand());
	k = (m - 1) - j;
	y = slice;
	while(j > 0 && y < likefun(*lowerLimit,ext)){
		*lowerLimit = *lowerLimit - w;
		if(*lowerLimit<SHAPE_SMALL){*lowerLimit = SHAPE_SMALL;}
		j = j - 1;
	}
	while(k > 0 && y < likefun(*upperLimit,ext)){
		*upperLimit = *upperLimit + w;
		if(*upperLimit>SHAPE_LARGE){*upperLimit = SHAPE_LARGE;}
		k = k - 1;
	}
}

void steppingOutShrinkage(double slice, double *shape, void *ext, double (*likefun)(double,void *ext),
		double lowerLimit, double upperLimit){
	double shapeCandidate, lowerLT, upperLT, y;
	lowerLT = lowerLimit;
	upperLT = upperLimit;
	shapeCandidate = lowerLT + unif_rand()*(upperLT - lowerLT);
	y = slice;
	while(y >= likefun(shapeCandidate,ext)){
		shapeCandidate = lowerLT + unif_rand()*(upperLT - lowerLT);
		if(shapeCandidate < *shape){lowerLT = shapeCandidate;}
		else{upperLT = shapeCandidate;}
	}
	*shape = shapeCandidate;
}

void doubling(double slice, double shape, void *ext,  double (*likefun)(double,void *ext),
		double *lowerLimit, double *upperLimit){
	int k,p = 3;
	double w = 0.1,y;
	*lowerLimit = shape - w * unif_rand();
	*upperLimit = *lowerLimit + w;
	k = p;
	y = slice;
	while(k > 0 && (y < likefun(*lowerLimit,ext) || y < likefun(*upperLimit,ext))){
		if(unif_rand() < 1/2){*lowerLimit = *lowerLimit - (*upperLimit - *lowerLimit);}
		else{*upperLimit = *upperLimit + (*upperLimit - *lowerLimit);}
		k = k - 1;
	}
}

void doublingShrinkage(double slice, double *shape, void *ext, double (*likefun)(double,void *ext),
		double lowerLimit, double upperLimit){
	double shapeCandidate, lowerLT, upperLT, w = 0.1, y, M;
	int D = 0;
	lowerLT = lowerLimit;
	upperLT = upperLimit;
	shapeCandidate = lowerLT + unif_rand()*(upperLT - lowerLT);
	y = slice;
	while(upperLT-lowerLT > 1.1*w){
		M = (lowerLT+upperLT)/2;
		if((*shape < M && shapeCandidate >= M) || (*shape >= M && shapeCandidate < M)){
			D = 1;}
		if(*shape < M){upperLT = M;}
		else{lowerLT = M;}
		if((D = 1) && (y >= likefun(lowerLT,ext) && y >= likefun(upperLT,ext))){
			shapeCandidate = lowerLT + unif_rand()*(upperLT - lowerLT);}
	}
	*shape = shapeCandidate;
}

double infectionDM(double parameter, void *ex)
{
	double *allTimes, *infPrior, remRate, infRate, lkhood;
	int i, k=0, *SS, *II, *ZZ, nEvents;
	AuxDM AS = (AuxDM) ex;
	SS = AS->SS;
	II = AS->II;
	ZZ = AS->ZZ;
	allTimes = AS->allTimes;
	nEvents = AS->nEvents;
	infPrior = AS->Prior;
	remRate = AS->Rate;
	infRate = parameter;
	lkhood = (infPrior[0]-1)*log(infRate)-infPrior[1]*infRate;
	for(i = 1; i < nEvents; ++i){
		if(allTimes[i] != allTimes[i-1]){k=i;}
		lkhood+=dexp((allTimes[k]-allTimes[(k-1)]),
			1/(SS[k]*II[k]*infRate+II[k]*remRate),TRUE);
	}
	return(lkhood); /* The likelihood */
}

double removalDM(double parameter, void *ex)
{
	double *allTimes, *remPrior, remRate, infRate, lkhood;
	int i, k = 0, *SS, *II, *ZZ, nEvents;
	AuxDM AS = (AuxDM) ex;
	SS = AS->SS;
	II = AS->II;
	ZZ = AS->ZZ;
	allTimes = AS->allTimes;
	nEvents = AS->nEvents;
	remPrior = AS->Prior;
	infRate = AS->Rate;
	remRate = parameter;
	lkhood = (remPrior[0]-1)*log(remRate)-remPrior[1]*remRate;
	for(i = 1; i < nEvents; ++i){
		if(allTimes[i]==allTimes[i-1]){k=i;}
		lkhood+=dexp((allTimes[k]-allTimes[(k-1)]),
			1/(SS[k]*II[k]*infRate+II[k]*remRate),TRUE);
	}
	return(lkhood); /* The likelihood */
}
