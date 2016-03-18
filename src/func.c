/*
Description: 
	Auxilary functions

Revision History:
	2014

Functions:
	Auxilary functions
	* 
Author:
	Albert Oude Nijhuis <albertoudenijhuis@gmail.com>

Institute:
	Delft University of Technology
	
Zephyros version:
	0.4

Project:
	EU FP7 program, the UFO project

Dissemination:
	Confidential, only for members of the UFO project. Potentially public in the future.

Acknowledgement and citation:
	Whenever this code used for publication of scientific results,
	the code writer should be informed, acknowledged and referenced.

Note:
	If you have any suggestions for improvements or amendments, please inform the author of this code.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "func.h"



//n is even
//tau = [0, 1, ..., n/2-1, -n/2, ..., -1]
//n is odd
//tau = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1]
void f_tau(
	int 	*n,
	int 	*tau
	)
{
	int i;
	
	for ( i = 0; i < *n; i++ ) {
		if (i < (*n/2.)) {
			tau[i] = i;			
		} else {
			tau[i] = i - *n;			
		}
	}
}

//n is even
//f = [0, 1, ..., n/2-1, -n/2, ..., -1] / (d*n)
//n is odd
//f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)
void f_fftfreq(
	int 	*n,
	double 	*d,
	double	*freq
	)
{
	int i;
	int *tau	= malloc(*n * sizeof(int));

	f_tau(n, tau);

	for ( i = 0; i < *n; i++ ) {
		freq[i] = tau[i] / (*d * *n);
	}	
	
	free(tau);
}
	

int isnanorinf(
	double *x)
{
	return (	isnan(*x) |
				isinf(*x) |
				(fabs(*x) > (0.9 * DBL_MAX)) |		//do not trust huge values
				(fabs(*x - 999.9) < 1.e-5) |		//typical fill value
				(fabs(*x - 9999.9) < 1.e-5) |		//typical fill value
				(fabs(*x - -999.9) < 1.e-5) |		//typical fill value
				(fabs(*x - -9999.9) < 1.e-5) 		//typical fill value
			);
}

double pnan(
	double *x)
{
	if (isnanorinf(x)) {
		return -999.9;
	} else {
		return *x;
	}
}

void calc_avg(
	int		*n,
	double 	*dat,
	double	*avg,
	int		*nok
	)
{
	int 	i;
	double 	total;

	*nok = 0;
	total = 0.;
	
	for ( i = 0; i < *n; i++ ) {
		if (!isnanorinf(&dat[i])) {
			total = total + dat[i];
			*nok = *nok + 1;
		}
	}
		
	if (*nok == 0) {
		*avg 		= NAN;
	} else {
		*avg	= total / *nok;
	}
}

void calc_avg_ang(
	int		*n,
	double 	*dat,
	double	*avg,
	int		*nok
	)
{
	int i;
	
	double *dat_cos = malloc(*n * sizeof(double));
	double *dat_sin = malloc(*n * sizeof(double));
	double avg_cos;
	double avg_sin;
		
	for ( i = 0; i < *n; i++ ) {
		dat_cos[i] = cos(dat[i]);
		dat_sin[i] = sin(dat[i]);
	}

	calc_avg(n, dat_cos, &avg_cos, nok);
	calc_avg(n, dat_sin, &avg_sin, nok);
	
	*avg	= func_modulo(atan2(avg_sin, avg_cos), 2. * M_PI);
		
	free(dat_cos);
	free(dat_sin);
}

void calc_avg_variance(
	int		*n,
	double 	*dat,
	double	*avg,
	double	*variance,
	int		*nok
	)
{
	int 	i;

	calc_avg(n, dat, avg, nok);
	
	if (*nok == 0) {
		*variance 	= NAN;
	} else {
		*variance = 0.;
		for ( i = 0; i < *n; i++ ) {
			if (!isnanorinf(&dat[i])) {			
				*variance = *variance + pow(dat[i] - *avg, 2);
			}
		}
		*variance = *variance / *nok;
	}
}




void calc_avg_variance_ang(
	int		*n,
	double 	*dat,
	double	*avg,
	double	*variance,
	int		*nok
	)
{
	int i;
	
	double dummy;	
	double *dat2 = malloc(*n * sizeof(double));
		
	calc_avg_ang(n, dat, avg, nok);
	
	//calculate differences with respect to average angle
	for ( i = 0; i < *n; i++ ) {
		dat2[i] = func_modulo(dat[i] - *avg, 2. * M_PI) ;
		if (dat2[i] > M_PI) dat2[i] = (2. * M_PI) - dat2[i];
	}

	calc_avg_variance(n, dat2, &dummy, variance, nok);
	
	free(dat2);
}

void running_average(
	int		*np,			//number of points
	double 	*dat,
	int		*nav_in,		//number of points to do averaging
	double	*running_avg,
	double 	*running_variance,
	int		*calc_variance		//0 = do not calculate variance, 1 = do calculate variance
	)
{
	int 	i, nav, lag;
	
	double 	this_running_avg, this_running_avg_old;
	double 	this_running_var, this_running_var_old;
	int		nok, nok_old;
	
	nav 	= fmin(*nav_in, *np);
	lag 	= fmax(1, nav / 2.);
	if (*calc_variance) {
		calc_avg_variance(&nav, dat, &this_running_avg, &this_running_var, &nok);
	} else {
		calc_avg(&nav, dat, &this_running_avg, &nok);		
	}
	
	for ( i = 0; i < (nav - lag); i++ ) {
		running_avg[i]		= this_running_avg;
		if (*calc_variance) {		
			running_variance[i]	= this_running_var;
		} else {
			running_variance[i]	= NAN;			
		}
	}
	
	this_running_avg_old	= this_running_avg;
	this_running_var_old	= this_running_var;
	nok_old					= nok;
	
	for ( i = nav; i < *np; i++ ) {
		//update nok
		if (!isnanorinf(&dat[i-nav])) 	nok = nok - 1;
		if (!isnanorinf(&dat[i])) 		nok = nok + 1;
		
		//update running average with new point
		if (nok_old == 0) {
			this_running_avg = 0.;
		} else {
			this_running_avg = ((1. * nok_old) / nok) * this_running_avg_old;
		}
		if (!isnanorinf(&dat[i-nav])) 	this_running_avg = this_running_avg - (dat[i - nav] / nok);
		if (!isnanorinf(&dat[i])) 		this_running_avg = this_running_avg + (dat[i] / nok);
		
		if (*calc_variance) {
			//update running variance with new point
			if (isnanorinf(&this_running_var_old)) {
				this_running_var = 0.;
			} else {
				this_running_var = ((1. * nok_old) / nok) * (this_running_var_old + pow(this_running_avg_old - this_running_avg,2));
				if (!isnanorinf(&dat[i-nav])) 	this_running_var = this_running_var - (pow(dat[i-nav]-this_running_avg,2)	/ nok);
				if (!isnanorinf(&dat[i])) 		this_running_var = this_running_var + (pow(dat[i]-this_running_avg,2) 	/ nok);
			}
		}
		
		if (nok == 0) {
			running_avg[i-lag]		= NAN;
			running_variance[i-lag]	= NAN;
		} else {
			running_avg[i-lag]		= this_running_avg;
			if (*calc_variance) {		
				running_variance[i-lag]	= this_running_var;
			} else {
				running_variance[i-lag] = NAN;
			}
		}
						
		//set old running average and variance
		this_running_avg_old	= this_running_avg;
		this_running_var_old	= this_running_var;
		nok_old					= nok;
		}
		
		//the last points
		for ( i = *np - lag; i < *np; i++ ) {
			if (nok == 0) {
				running_avg[i]		= NAN;
				running_variance[i] = NAN;
			} else {
				running_avg[i]		= this_running_avg;
				if (*calc_variance) {					
					running_variance[i]	= this_running_var;
				} else {
					running_variance[i]	= NAN;
				}
			}
		}
		
		//if input was inf  or nan, then also the output is
		for ( i = 0; i < *np; i++ ) {
			if (isnanorinf(&dat[i])) {
				running_avg[i]		= NAN;
				running_variance[i]	= NAN;
			}
		}
}


void running_average_ang(
	int		*np,			//number of points
	double 	*dat,
	int		*nav_in,		//number of points to do averaging
	double	*running_avg,
	double 	*running_variance,
	int		*calc_variance		//0 = do not calculate variance, 1 = do calculate variance
	)
{
	int i;
	
	double *dat_cos = malloc(*np * sizeof(double));
	double *dat_sin = malloc(*np * sizeof(double));
	double *running_avg_cos = malloc(*np * sizeof(double));
	double *running_avg_sin = malloc(*np * sizeof(double));
	double *running_variance_cos = malloc(*np * sizeof(double));
	double *running_variance_sin = malloc(*np * sizeof(double));
	
	double *dat2 	= malloc(*np * sizeof(double));
	double *dummy	= malloc(*np * sizeof(double));
	
	int calc_var;
	
	for ( i = 0; i < *np; i++ ) {
		dat_cos[i] = cos(dat[i]);
		dat_sin[i] = sin(dat[i]);
	}

	calc_var = 0;
	running_average(np, dat_cos, nav_in, running_avg_cos, running_variance_cos, &calc_var);
	running_average(np, dat_sin, nav_in, running_avg_sin, running_variance_sin, &calc_var);
	
	for ( i = 0; i < *np; i++ ) {
		running_avg[i]		= func_modulo(atan2(running_avg_sin[i], running_avg_cos[i]), 2. * M_PI);
	}
	
	
	if (*calc_variance) {
		//calculate differences with average angle
		//This is not exact. I.e. if there is a large trent in the running average this will fail...
		for ( i = 0; i < *np; i++ ) {
			dat2[i] = func_modulo(dat[i], 2. * M_PI) - running_avg[i];
			if (dat2[i] > M_PI) dat2[i] = (2. * M_PI) - dat2[i];
		}

		calc_var = 1;
		running_average(np, dat2, nav_in, dummy, running_variance, &calc_var);
	} else {
		for ( i = 0; i < *np; i++ ) {		
			running_variance[i] = NAN;
		}
	}
	
	free(dat_cos);
	free(dat_sin);
	free(running_avg_cos);
	free(running_avg_sin);
	free(running_variance_cos);
	free(running_variance_sin);

	free(dat2);
	free(dummy);
}


void rpowspec(
	int 	*n,
	double 	*xreal,
	double 	*power,
	int		*periodic
	)
{
	double complex  *c 		= malloc(*n * sizeof(double complex));
	double complex *acfft 	= malloc(*n * sizeof(double complex));
	double *ac						= malloc(*n * sizeof(double));
	double complex *x2		= malloc(*n * sizeof(double complex));
	double mu;
	int	i, nok;

	//count number of ok data
	nok = 0;
	for ( i = 0; i < *n; i++ ) {
		if (!isnanorinf(&xreal[i])) {
			nok = nok + 1;
		}
	}
		
	if ((*periodic) & (nok == *n)) {
		//Calculate power spectrum via Fourier transform
		//Only possible when there are no nans, and for periodic analysis
		rfft(n, xreal, c);
		
		for ( i = 0; i < *n; i++ ) {
			power[i] = pow(cabs(c[i]) / *n,2);
		}					
	} else {
		//Calculate power spectrum via autocovariance
		autocovariance(n, xreal, ac, periodic);
		calc_avg(n, xreal, &mu, &nok);

		for ( i = 0; i < *n; i++ ) {
			x2[i] = (ac[i] + pow(mu,2.)) / *n;
		}
		fft(n, x2, acfft);

		for ( i = 0; i < *n; i++ ) {
			power[i] = creal(acfft[i]);
			power[i] = fabs(power[i]);
		}		
	}

	free(c);
	free(acfft);
	free(ac);
	free(x2);
}


void rifft(
	int *n,
	double complex *c,
	double *xreal)
{
	int i;
	
	double complex *x = malloc(*n * sizeof(double complex));
		
	ifft(n, c, x);
	for ( i = 0; i < *n; i++ ) {
		xreal[i] = creal(x[i]);
	}

	free(x);
}

void rfft(
	int		*n,
	double *xreal,
	double complex *c
	)
{
	int i;
	double complex *x = malloc(*n * sizeof(double complex));
	
	for ( i = 0; i < *n; i++ ) {
		x[i] = xreal[i];
	}		
	fft(n, x, c);
	
	free(x);
}


void fft(
	int		*n,
	double complex *x,
	double complex *c)
{
	int i,j;
			
	for ( i = 0; i < *n; i++ ) {
		c[i] = 0.;
		for ( j = 0; j < *n; j++ ) {
			c[i] += 
				x[j] * cexp( I * ((-2.  * M_PI *  i * j) / *n) );
		}		
	}
}


void ifft(
	int		*n,
	double complex *c,
	double complex *x)
{
	int i,j;
	
	for ( i = 0; i < *n; i++ ) {
		x[i] = 0.;
		
		for ( j = 0; j < *n; j++ ) {		
			x[i] = x[i] + 
					c[j] * cexp( (2.0I * M_PI * i * j) / *n);
		}
	}

	for ( j = 0; j < *n; j++ ) {
		x[i] /= *n;
	}
}



void covariance(
	int		*n,
	double 	*x,
	double 	*y,
	double ans)
{
	double mux, muy;
	int i;
	int nok;
	
	nok = 0;
	for ( i = 0; i < *n; i++ ) {
		if (!(isnanorinf(x + i) | isnanorinf(y + i))) {
			nok += 1;
			mux += x[i];
			muy += y[i];
		}
	}
	
	mux /= nok;
	muy /= nok;
	
	for ( i = 0; i < *n; i++ ) {
		if (!(isnanorinf(x + i) | isnanorinf(y + i))) {
			ans =  (x[i] - mux) * (y[i] - muy);
		}
	}	
	ans /= nok;
}

void autocovariance(
	int 		*n,
	double		*x,
	double 		*y,
	int			*periodic)
{
	int i, j, t;
	int *tau = malloc(*n * sizeof(int));
	double mu;
	int nok;
	
	f_tau(n, tau);
		
	//R(t) = <(x(i) - mu)(x(i+t) - mu)>
	//R(t) = <x(i)(x(i+t))> - mu^2

	calc_avg(n, x, &mu, &nok);
	
	if (*periodic) {
		//assume that signal is periodic
		for ( i = 0; i < *n; i++ ) {		
			t = tau[i];
			y[i] = 0.;
			nok = 0;
			if (t > 0) {
				for ( j = 0		; j < (*n - t); j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j + t))) {
						y[i] += x[j] * x[j+t];
						nok++;
					}
				}
				for ( j = *n-t	; j < *n; j++ ) {
					if (!(isnanorinf(x+j) | isnanorinf(x + j + t - *n))) {
						y[i] += x[j] * x[j+t-*n];
						nok++;
					}
				}
			} else {
				for ( j = 0		; j < -t; j++ ) {
					if (!(isnanorinf(x+ j) | isnanorinf(x + j + t + *n))) {
						y[i] += x[j] * x[j+t+*n];
						nok++;
					}
				}
				for ( j = -t	; j < *n; j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j + t))) {
						y[i] += x[j] * x[j+t];
						nok++;
					}
				}
			}
			y[i] = (y[i] / nok) - pow(mu,2);
		}
	} else {
		//assume that signal is non-periodic
		for ( i = 0; i < *n; i++ ) {		
			t = tau[i];
			y[i] = 0.;
			nok = 0;
			if (t > 0) {
				for ( j = 0		; j < (*n - t); j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j+t))) {
						y[i] += (x[j] * x[j+t]);
						nok++;
					}
				}
			} else {
				for ( j = -t	; j < *n; j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j+t))) {
						y[i] += (x[j] * x[j+t]);
						nok++;
					}
				}
			}
			y[i] = (y[i] / nok) - pow(mu,2);
		}
	}

	free(tau);
}


void autocorrelation(
	int 	*n,
	double	*x,
	double	*y,
	int		*periodic)
{
	int i;
	autocovariance(n, x,y, periodic);

	//calculate autocorrelation
	for ( i = 0; i < *n; i++ ) {
		y[i] = y[i] / y[0];
	}	
}


//structure function
void structure_function(
	int			*n,
	double 		*x,
	double 		*y,
	int			*order,
	int			*periodic)
{
	int *tau = malloc(*n * sizeof(int));
	int i, t, j;
	int nok;
	
	f_tau(n, tau);	

	if (*periodic == 1) {
		//assume that signal is periodic
		for ( i = 0; i < *n; i++ ) {		
			t = tau[i];
			y[i] = 0.;
			nok = 0;
			if (t > 0) {
				for ( j = 0		; j < (*n - t); j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j + t))) {
						y[i] += pow(x[j] - x[j+t], *order);
						nok++;
					}
				}
				for ( j = *n-t	; j < *n; j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j + t - *n))) {
						y[i] += pow(x[j] - x[j+t-*n], *order);
						nok++;
					}
				}
			} else {
				for ( j = 0		; j < -t; j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j + t + *n))) {
						y[i] += pow(x[j] - x[j+t+*n], *order);
						nok++;
					}
				}
				for ( j = -t	; j < *n; j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j + t))) {
						y[i] += pow(x[j] - x[j+t], *order);
						nok++;
					}
				}
			}
			y[i] = (y[i] / nok);
		}
	} else {
		//assume that signal is non-periodic
		for ( i = 0; i < *n; i++ ) {		
			t = tau[i];
			y[i] = 0.;
			nok = 0;
			if (t > 0) {
				for ( j = 0		; j < (*n - t); j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j + t))) {
						y[i] += pow(x[j] - x[j+t], *order);
						nok++;
					}
				}
			} else {
				for ( j = -t	; j < *n; j++ ) {
					if (!(isnanorinf(x + j) | isnanorinf(x + j + t))) {
						y[i] += pow(x[j] - x[j+t], *order);
						nok++;
					}
				}
			}
			y[i] = (y[i] / nok);
		}
	}
}
		
		
		

void randompermutation(int n, int **output) 
{
    int i;
    int j;
    int temp;
 
    *output = malloc(n * sizeof(int));
       
    for(i=0;i<n;i++){
        (*output)[i]=i;
    }
    
    for(i=0; i<n;i++){
        j = rand() % n;

		//swap
		temp = (*output)[i];
		(*output)[i] = (*output)[j];
		(*output)[j] = temp;
    }
}

double angleAminB(double *Arad, double *Brad)
{
	return func_modulo(M_PI + *Arad - *Brad, 2. * M_PI) - M_PI;
}

//void array_double_malloc(
void func_dbl_arr_calloc(
	int n,
	double **ptrarr
	)
{
	double *arr;

	if (n == 0) {
		printf("WARNING: double array initialized with size 0");
		*ptrarr = NULL;
	} else {
		arr = calloc(n, sizeof(double));
		if (arr == NULL) {
		  printf("out of memory!\n");
		  exit(1);
		}
		*ptrarr = arr;
	}
}

//void array_int_malloc(
void func_int_arr_malloc(
	int n,
	int **ptrarr
	)
{
	int *arr;

	if (n == 0) {
		printf("WARNING: double array initialized with size 0");
		*ptrarr = NULL;
	} else {
		arr = calloc(n, sizeof(int));
		if (arr == NULL) {
		  printf("out of memory!\n");
		  exit(1);
		}
		*ptrarr = arr;
	}
}

void func_ptr_arr_malloc(
	int n,
	int **ptrarr
	)
{
	int *arr;

	if (n == 0) {
		printf("WARNING: pointer array initialized with size 0");
		*ptrarr = NULL;
	} else {
		arr = calloc(n, sizeof(void*));
		if (arr == NULL) {
		  printf("out of memory!\n");
		  exit(1);
		}
		*ptrarr = arr;
	}
}

/* dot product of an integer array */
void dbldotprod(
	int*	n,
	double*	arr1,
	double*	arr2,
	double*	result)
{
	int i;
	
	*result = 0.;
	for ( i = 0; i < *n; i++ ) {
		*result += arr1[i]*arr2[i];
	}
}


double uniform_random()
{
	double mymax = 1.;
	double mymin = 0.;
	return ((double) rand() / (RAND_MAX+1.)) * (mymax-mymin+1.) + mymin;
}

double read_uniform_random(FILE *fp)
{
	double ans;
	
	fscanf(fp, "%lf", &ans);
	if (feof(fp)) {
		rewind(fp);
		fscanf(fp, "%lf", &ans);
	}
	
	return ans;
}



double func_newton_rapson_solver(
	double (*f)(double, void*),
	double (*df)(double, void*),
	double x0_in,
	double allowed_error,
	int	max_itr,
	void *vp
	)
{
	int itr;
    float h, x0, x1, allerr;
	x0 = x0_in;
	
    for (itr=0; itr < max_itr; itr++)
    {
        h=f(x0, vp)/df(x0, vp);
        
        //printf("\n\n");
        //printf("f = %.2e\n", f(x0, vp));
        //printf("df = %.2e\n", df(x0, vp));
        
        x1=x0-h;
        if (fabs(h) < allowed_error)
        {
			//ok
            return x1;
        }
        x0=x1;
    }
    //maximum iterations where performed but solution not found
	return x1;
}

double func_dB(double val)
{
	if (val == 0.) {
		return -9.999e5;
	} else {
		return 10. * log10(val);
	}
}

double func_dB_inv(double val)
{
	if (val == -9.999e5) {
		return 0.;
	} else {
		return pow(10., val / 10.);
	}
}

double func_modulo(double x, double y)
{
	int i;
	
	if (x < 0) {
		i = abs(x / y);
		return fmod(x + ((1 + i) * y), y);
	} else {
		return fmod(x,y);
	}
}
