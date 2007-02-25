
extern double trexponential(double lower_cut_off, double upper_cut_off, double rate);

extern double trgamma(double lower_cut_off, double upper_cut_off, double shape, double rate);

extern double trweibull(double lower_cut_off, double upper_cut_off, double shape, double rate);

/* Static Variables */
extern SEXP getListElement(SEXP list, char *str);

/* consistency functions */
/* SIR model */
extern int consistent_SIR(double *infectionTimes, double *removalTimes,
	int *N, int *nInfected, int *nRemoved, double *allTimes,
	int *indicator, int *SS, int *II);
/* SEIR model */
extern int consistent_SEIR(double *infectionTimes,	double *endLatency,
	double *removalTimes, int *N, int *nInfected, int *nRemoved, 
	double *allTimes, int *indicator, int *SS, int *LL, int *II);

/* Shape likelihoods and derivatives for use with 'vmmin' */
extern double weiShapelkhood(int n, double *par, void *ex);

extern void weiFderivShape(int n, double *par, double *gr, void *ex);

extern double weiSderivShape(double shapePost, void *ex);

extern double weiShapelkhoodARS(double parameter, void *ex);

extern double gammaShapelkhood(int n, double *par, void *ex);

extern void gammaFderivShape(int n, double *par, double *gr, void *ex);

extern double gammaSderivShape(double shapePost, void *ex);

extern double gammaShapelkhoodARS(double parameter, void *ex);
/* Weibull functions called by ARS */
extern void weibullLikelihood_SIR(double *parameters, double *infectionTimes,
	double *removalTimes, int *N, int *nInfected, int *nRemoved,
	double *sumSI, double *sumDurationInfectious, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *II);

extern void weibullLikelihood_SEIR(double *parameters, double *infectionTimes,
	double *endLatency, double *removalTimes, int *N, int *nInfected, int *nRemoved, 
	double *sumSI, double *sumLatents, double *sumInfectious, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *LL, int *II);

/* Gamma functions called by ARS */
extern void gammaLikelihood_SIR(double *parameters, double *infectionTimes,
	double *removalTimes, int *N, int *nInfected, int *nRemoved,
	double *sumSI, double *sumDurationInfectious, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *II);

extern void gammaLikelihood_SEIR(double *parameters, double *infectionTimes,
	double *endLatency, double *removalTimes, int *N, int *nInfected, int *nRemoved, 
	double *sumSI, double *sumLatents, double *sumInfectious, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *LL, int *II);

/* Truncated gamma random variates */
extern int n_optimal(double scale);

extern double gamma_right(double shape, double scale, double t);

extern double gamma_left(double shape, double scale, double t);

/* Define structure used in shape likelihoods */
typedef struct aux_struct
{
	int nRemoved;
	double *startTimes;
	double *endTimes;
	double *shapePrior;
	double scale;
} aux_struct, *AuxStruct;

extern void infRateMixLkhood_SIR(double *parameters, double *infectionTimes, double *removalTimes,
	int *N, int *nInfected, int *nRemoved, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *II);

extern void remRateMixLkhood_SIR(double *parameters, double *infectionTimes, double *removalTimes,
	int *N, int *nInfected, int *nRemoved, double *likelihood,
	double *allTimes, int *indicator, int *SS, int *II);

/* Slice sampling */
extern void steppingOut(double slice, double shape, void *ext, double (*likefun)(double,void *ext),
	double *lowerLimit, double *upperLimit);
extern void steppingOutShrinkage(double slice, double *shape, void *ext, double (*likefun)(double,void *ext),
	double lowerLimit, double upperLimit);
extern void doubling(double slice, double shape, void *ext, double (*likefun)(double,void *ext),
	double *lowerLimit, double *upperLimit);
extern void doublingShrinkage(double slice, double *shape, void *ext, double (*likefun)(double,void *ext),
	double lowerLimit, double upperLimit);

typedef struct aux_dm
{
	int nEvents;
	int *SS;
	int *II;
	double *ZZ;
	double *allTimes;
	double *Prior;
	double Rate;
} aux_dm, *AuxDM;

extern double infectionDM(double parameter, void *ex);
extern double removalDM(double parameter, void *ex);
