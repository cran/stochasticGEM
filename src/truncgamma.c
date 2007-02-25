#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>

/* Adapted from "Simulation of right and left truncated Gamma distributions   */
/* by mixtures", Statistics and computing, 1997:173-181						  */
/******************************************************************************/
/*     algorithm to generate right truncated gamma distribution               */
/******************************************************************************/
/* the following function gives the optimal number of components              */
/* for p=0.95 fixed                                                           */
static int n_optimal(double scale)
{
    double nr,q;
    q=qnorm(0.95,0,1,TRUE,FALSE);/* quantile of standard normal */    
    nr=floor(R_pow(q+sqrt(q*q+4*scale),2)/4.0);
    return((int) nr);
}
/* the following function returns a random number from TG^-(shape,scale,1)    */ 
double gamma_right(double shape, double scale, double t)
{
	int n,i,j,nprotect=0;
	double x=0,M_x=1,rho_x=0;
	SEXP w_bar,w_tilde;
	scale= scale*t;
	n=n_optimal(scale);/* determine "n" */
	PROTECT(w_bar = allocVector(REALSXP, (n+1)));
	++nprotect;
	PROTECT(w_tilde = allocVector(REALSXP, (n+1)));
	++nprotect;
	REAL(w_bar)[1]=1; REAL(w_tilde)[1]=1;/* line 6 page 175 */
	for (i=1; i<n;i++)/* k \in 1,...,N-1 */
	{/* Equation 5 in article */
		REAL(w_bar)[i+1]=REAL(w_bar)[i]*scale/(shape+i);
		REAL(w_tilde)[i+1]=REAL(w_tilde)[i]+REAL(w_bar)[i+1];
	};
	REAL(w_tilde)[0]=0;/* line 7 page 175 */
	for(i=0; i<=n; i++) 
	{ 
		REAL(w_tilde)[i]=REAL(w_tilde)[i]/REAL(w_tilde)[n];
	};
	while (unif_rand()*M_x > rho_x)/* Eqn. A_1N: 3 */
	{
		j=0;
		while(unif_rand()>REAL(w_tilde)[j]){j=j+1;};
		x=rbeta(shape, (double) j+1);/* Eqn. A_2: 2 */
		M_x=0; rho_x=0;
		for (i=1; i<=n; i++)
		{
			M_x += R_pow(scale,i-1)/gammafn(i);/* Eqn. A_1N: 2, gammafn is gamma function */
			rho_x += exp(scale*x)*R_pow(scale*(1-x),i-1)/gammafn(i);/* Eqn. A_1N: 2 */
		};
		M_x = 1/M_x; rho_x = 1/rho_x;
	};
	UNPROTECT(nprotect);
	return(x*t);
}

/******************************************************************************/
/*algorithm to generate left  truncated gamma distribution                    */
/******************************************************************************/
/* the following function returns a random number from TG^+(shape,scale,1)    */
/* when "shape" a is an interger                                              */
/* See Devroye, L. (1985)  Non-Uniform Random Variate Generation              */
/* Springer-Verlag, New-York.                                                 */
static double gamma_int_shape(double shape, double scale)
{
	double x,nprotect=0;
	SEXP w_bar,w_tilde;
	int i;
	PROTECT(w_bar = allocVector(REALSXP, (int) shape));
	++nprotect;
	PROTECT(w_tilde = allocVector(REALSXP, (int) shape));
	++nprotect;
	REAL(w_bar)[0]=1.0;REAL(w_tilde)[0]=1.0;
	for(i=1; i< (int)shape ;i++)/* 2<=k<=A */ 
	{/* Proposition 3.1 */
		REAL(w_bar)[i]=REAL(w_bar)[(i-1)]*(shape-(i+1)+1)/scale;/* k=i+1 */
		REAL(w_tilde)[i]=REAL(w_tilde)[(i-1)]+REAL(w_bar)[i];
	};
	for(i=0; i< (int) shape; i++)
	{/* Proposition 3.1 */
		REAL(w_tilde)[i]=REAL(w_tilde)[i]/REAL(w_tilde)[(int)shape-1];
	};
	i=0;
	while(unif_rand()>REAL(w_tilde)[i]){i=i+1;};/* A3 : 2 */
	x=rgamma( (double) (i+1), scale); /* k=i+1 */
	UNPROTECT(nprotect);
	return(x+1);/* A3 : 3 */
}
/* the following function returns a random number from TG^+(shape,scale,1) */
double gamma_left(double shape, double scale, double t)
{
	double x=0,rho_x=0,M;
	scale= scale*t;
	if (shape<1.0){/* If shape is less than 1: A6 */
		while (rho_x < unif_rand()){
			x=1-log(1-unif_rand())/scale;
			/* x=gamma_int_shape(ceil(shape), scale); */
			rho_x=R_pow(x,shape-ceil(shape));
		};
	}
	else{
		if (scale<=shape){/* If shape is less than scale */
			M=exp(floor(shape)-shape);
			while (rho_x <= unif_rand()*M){/* Eqn. A5: 3 */
				x=gamma_int_shape(floor(shape), scale*floor(shape)/shape);/* Eqn. A5: 1 */
				rho_x = R_pow(x, shape-floor(shape))*exp(-x*scale*(1-floor(shape)/shape));
			};
		}
		else{/* If shape is greater than scale */
			M=R_pow(shape/scale,shape-floor(shape))*exp(floor(shape)-shape);
			while (rho_x <= unif_rand()*M){/* Eqn. A5: 3 */
				x=gamma_int_shape(floor(shape), scale-shape+floor(shape));/* Eqn. A5: 1 */
				rho_x=R_pow(x, shape-floor(shape))*exp(-x*(shape-floor(shape)));/* Eqn. A5: 2 */
			};
		};
	};
	return(x*t);
}


/* the following function returns a random number from TG^-(shape,scale,1)    */ 
SEXP gammaRight(SEXP n, SEXP a, SEXP b, SEXP t)
{
	int i,nprotected=0;
	SEXP truncGam;
	PROTECT(n = AS_INTEGER(n));
	++nprotected;
	PROTECT(a = AS_NUMERIC(a));
	++nprotected;
	PROTECT(b = AS_NUMERIC(b));
	++nprotected;
	PROTECT(t = AS_NUMERIC(t));
	++nprotected;
	PROTECT(truncGam = allocVector(REALSXP,INTEGER(n)[0]));
	++nprotected;
	GetRNGstate(); /* should be before a call to a random number generator */
	for(i = 1; i <= INTEGER(n)[0]; i++){
		REAL(truncGam)[(i-1)] = gamma_right(REAL(a)[0], REAL(b)[0], REAL(t)[0]);
	}
	PutRNGstate(); /* after using random number generators.	*/
	UNPROTECT(nprotected);
	return(truncGam);
}

/* the following function returns a random number from TG^+(shape,scale,1) */
SEXP gammaLeft(SEXP n, SEXP a, SEXP b, SEXP t)
{
	int i,nprotected=0;
	SEXP truncGam;
	PROTECT(n = AS_INTEGER(n));
	++nprotected;
	PROTECT(a = AS_NUMERIC(a));
	++nprotected;
	PROTECT(b = AS_NUMERIC(b));
	++nprotected;
	PROTECT(t = AS_NUMERIC(t));
	++nprotected;
	PROTECT(truncGam = allocVector(REALSXP,INTEGER(n)[0]));
	++nprotected;
	GetRNGstate(); /* should be before a call to a random number generator */
	for(i = 1; i <= INTEGER(n)[0]; i++){
		REAL(truncGam)[(i-1)] = gamma_left(REAL(a)[0], REAL(b)[0], REAL(t)[0]);
	}
	PutRNGstate(); /* after using random number generators.	*/
	UNPROTECT(nprotected);
	return(truncGam);
}
