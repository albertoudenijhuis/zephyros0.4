/*
Description: 
	EDR retrieval functions

Revision History:
	2014

Functions:
	EDR retrieval functions
	
Author:
	Albert Oude Nijhuis <albertoudenijhuis@gmail.com>

Institute:
	Delft University of Technology
	
Zephyros version:
	0.2

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
#include <string.h>

#include "func.h"
#include "edr.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

void retr_edr_variance(
	int 	*n,				//length of variables v,t,l,v2
	double 	*v,
	double 	*t,
	double 	*l,
	double 	*v2,
	int 	*retr_domain,	//0 = space, 1 = time and U0 from v2
	double 	*retr_C,		//Kolmogorov constant that is used.
	int		*retr_n,		//Number of points used for average and variance calculation
	double	*retr_av,		//output: retrieved average velocity
	double 	*retr_var,		//output: retrieved variance of velocity
	double	*retr_edr,		//output: retrieved edr
	double	*retr_edr13_err	//output: error in retrieved edr^(1/3) (std)
	)
{
	int i, i1, i2, ni;
	int nav, lag;
	int calc_variance;
	
	double *abs_dl = malloc((*n-1) * sizeof(double));
	double *abs_dt = malloc((*n-1) * sizeof(double));
	double *v2_avg = malloc(*n * sizeof(double));
	double *v2_var = malloc(*n * sizeof(double));
	double dlmax, dlmin;
	double dtmax, dtmin;
	double f_u, f_l;
	double g_u, g_l;
	
    for ( i = 0; i < (*n - 1); i++ ) {
		abs_dl[i] = fabs(l[i+1] - l[i]);
		abs_dt[i] = fabs(t[i+1] - t[i]);		
	}

	//calculate running average and running std.
	calc_variance = 1;
	nav 	= fmin(*retr_n, *n);
	lag 	= fmax(1, nav / 2.);
	running_average(n, v, retr_n, retr_av, retr_var, &calc_variance);

	if (*retr_domain == 1) {
		//for time domain, we need the running average and running variance of v2
		calc_variance = 1;
		running_average(n, v2, retr_n, v2_avg, v2_var, &calc_variance);		
	}

    for (i = 0; i < *n; i++ ) {
		i1 = fmax(0, i-lag-1);
		i2 = i1 + (*retr_n -1);

		if (i2 >= *n) {
			i2 = *n -1;
            i1 = fmax(0, i2 - (*retr_n - 1));
		}
		ni = i2 + 1 - i1;
		
		if (*retr_domain == 0) {
			//space domain
			dlmax	= fabs(l[i2] - l[i1]);
			dlmin	= dlmax / (ni - 1);
			f_l		= 2. * M_PI / dlmax;
			f_u		= (((int) (ni/2)) -1.) * f_l;
						
			retr_edr[i] = 	pow(	(3./2.) *
									*retr_C *
									( pow(f_l - (0.5 * f_l),-2./3.)  - pow(f_u + (0.5 * f_l), -2./3.) )
								, -3./2.
								)
								* pow(retr_var[i], 3./2.);

			retr_edr13_err[i] = (1./3.) *
								pow(retr_edr[i] , 1. / 3.) 
								* sqrt(pow(f_l / f_u, 4./3.) +  (9. / (2. * (ni - 1.))));		
		} else {
			//time domain
			dtmax	= fabs(t[i2] - t[i1]);
			dtmin	= dtmax / (ni-1);
			g_l		= 2. * M_PI / dtmax;
			g_u		= (((int) (ni/2)) -1.) * g_l;
						
			retr_edr[i] = 	pow(	(3./2.) *
									*retr_C *
									( pow(g_l - (0.5 * g_l),-2./3.)  - pow(g_u + (0.5 * g_l), -2./3.) )
								, -3./2.
								)
								* pow(retr_var[i], 3./2.) 
								/ fabs(v2_avg[i]);

			retr_edr13_err[i] = (1./3.) *
								pow(retr_edr[i] , 1. / 3.) 
								* sqrt(pow(g_l / g_u, 4./3.) + (retr_var[i] / pow(v2_avg[i],2)) + (9. / (2. * (ni - 1.))));
		}
	}
	
	free(abs_dl);
	free(abs_dt);
	free(v2_avg);
	free(v2_var);
}

void retr_edr_powerspectrum(
	int 	*n,				//length of variables v,t,l,v2
	double 	*v,
	double 	*t,
	double 	*l,
	double 	*v2,
	int 	*retr_domain,		//0 = space, 1 = time and U0 from v2
	double 	*retr_C,			//Kolmogorov constant that is used.
	int		*retr_n,			//Number of points used for average and variance calculation
	int		*periodic,			//1 = periodic, 0 = non-periodic
	int		*retr_nintervals,	//number of intervals used
	double	*retr_edr,			//output: retrieved edr
	double	*retr_edr13_err		//error in retrieved edr^(1/3) (std)
	)
{
	int retr_nfreq;
	int i, i1, i2, ni;
	int ic;
	int j;
	int calc_variance;
	
	double *abs_dl = malloc((*n-1) * sizeof(double));
	double *abs_dt = malloc((*n-1) * sizeof(double));
	double *v2_avg = malloc(*n * sizeof(double));
	double *v2_var = malloc(*n * sizeof(double));

	double *fft_pow = malloc(*n * sizeof(double));
	double *fft_freq = malloc(*n * sizeof(double));

	double dlmax, dlmin;
	double dtmax, dtmin;
	double f_u, f_l;
	double g_u, g_l;
	
    double *lst_eps		= malloc(*n * sizeof(double));
    double *lst_eps13	= malloc(*n * sizeof(double));
	int		lst_n;
	double	freqmin, freqmax;
	double 	thispow;
	int 	k;
	double  avg, var;
    int 	nok;
    
    for ( i = 0; i < (*n - 1); i++ ) {
		abs_dl[i] = fabs(l[i+1] - l[i]);
		abs_dt[i] = fabs(t[i+1] - t[i]);		
	}

	if (*retr_domain == 1) {
		//for time domain, we need the running average and running variance of v2
		calc_variance = 1;
		running_average(n, v2, retr_n, v2_avg, v2_var, &calc_variance);		
	}

        
    for (i1 = 0; i1 < *n; i1 = i1 + *retr_n ) {
		i2 = i1 + (*retr_n -1);

		if (i2 >= *n) {
			i2 = *n -1;
            i1 = fmax(0, i2 - (*retr_n - 1));
		}
		ni = i2 + 1 - i1;
		ic = i1 + ni/2;
		
		retr_nfreq = fmax(1, ceil(((int) (ni / 2)) / (1. * *retr_nintervals)));
				
		if (*retr_domain == 0) {
			//space domain
			dlmax	= fabs(l[i2] - l[i1]);
			dlmin	= dlmax / (ni-1);
			f_l		= 2. * M_PI / dlmax;
			f_u		= (((int) (ni/2)) -1.) * f_l;
			
			rpowspec(&ni, v + i1, fft_pow, periodic);
			f_fftfreq(&ni, &dlmin, fft_freq);
			for (j = 0; j < ni; j++ ) {
				fft_freq[j] = 2. * M_PI * fft_freq[j];
				
				//printf("fft_freq[%i] = %.2e \n", j, fft_freq[j]);
				//printf("fft_pow[%i] 	= %.2e \n", j, fft_pow[j]);
			}
			
			//calculate epsilon for frequency intervals
			lst_n = 0;
			k = 0;
			freqmax = 0.;
			freqmin = 1.e20;
			thispow = 0.;
			for (j = 0; j < ni; j++ ) {
				if (fft_freq[j]>0.) {
					k++;
					thispow = thispow + 2. * fft_pow[j]; //pow of positive and negative frequency
					freqmin = fmin(freqmin, fft_freq[j]);
					freqmax = fmax(freqmax, fft_freq[j]);
				
					if (k == retr_nfreq) {
						lst_eps[lst_n] =  pow(	(3./2.) *
													*retr_C *
													( pow(freqmin - (f_l / 2.),-2./3.)  - pow(freqmax + (f_l / 2.), -2./3.) )
													, -3./2.
												)
												* pow(thispow, 3./2.);
						lst_eps13[lst_n] = pow(lst_eps[lst_n], 1./3.);
						lst_n++;
														
						//initialize next interval
						k = 0;
						freqmax = 0.;
						freqmin = 1.e20;
						thispow = 0.;
					}
				}
			}

			//add last interval
			if (k > 0) {
				lst_eps[lst_n] =  pow(	(3./2.) *
											*retr_C *
											( pow(freqmin - (f_l / 2.),-2./3.)  - pow(freqmax + (f_l / 2.) , -2./3.) )
											, -3./2.
										)
										* pow(thispow, 3./2.);
				lst_eps13[lst_n] = pow(lst_eps[lst_n], 1./3.);
				lst_n++;

			}

			

							
			//calculate statistics of list of edrs
            calc_avg_variance(&lst_n, lst_eps13, &avg, &var, &nok);

			for (j = i1; j <= i2; j++ ) {
				retr_edr[j]			= pow(avg,3);
				retr_edr13_err[j] 	= sqrt(var);
			}
		} else {
			//time domain
			dtmax	= fabs(t[i2] - t[i1]);
			dtmin	= dtmax / (ni-1);
			g_l		= 2. * M_PI / dtmax;
			g_u		= (((int) (ni/2)) -1.) * g_l;

			
			rpowspec(&ni, v + i1, fft_pow, periodic);
			f_fftfreq(&ni, &dtmin, fft_freq);
			for (j = 0; j < ni; j++ ) {
				fft_freq[j] = 2. * M_PI * fft_freq[j];
			}

			//calculate epsilon for frequency intervals
			lst_n = 0;
			k = 0;
			freqmax = 0.;
			freqmin = 1.e20;
			thispow = 0.;

			for (j = 0; j < ni; j++ ) {
				if (fft_freq[j]>0.) {
					k++;
					thispow = thispow + 2. * fft_pow[j]; //pow of positive and negative frequency
					freqmin = fmin(freqmin, fft_freq[j]);
					freqmax = fmax(freqmax, fft_freq[j]);
				
					if (k == retr_nfreq) {
						lst_eps[lst_n] =  pow(	(3./2.) *
													*retr_C *
													( pow(freqmin- (g_l / 2.),-2./3.)  - pow(freqmax + (g_l / 2.), -2./3.) )
													, -3./2.
												)
												* pow(thispow, 3./2.)
												/ fabs(v2_avg[ic]);
						lst_eps13[lst_n] = pow(lst_eps[lst_n], 1./3.);
						lst_n++;
								
						//initialize next
						k = 0;
						freqmax = 0.;
						freqmin = 1.e20;
						thispow = 0.;
					}
				}
			}

			if (k > 0) {
				lst_eps[lst_n] =  pow(	(3./2.) *
											*retr_C *
											( pow(freqmin - (g_l / 2.),-2./3.)  - pow(freqmax + (g_l / 2.), -2./3.) )
											, -3./2.
										)
										* pow(thispow, 3./2.)
										/ fabs(v2_avg[ic]);
				lst_eps13[lst_n] = pow(lst_eps[lst_n], 1./3.);
				lst_n++;
			}
			
            //calculate statistics of list of edrs
            calc_avg_variance(&lst_n, lst_eps13, &avg, &var, &nok);
			for (j = i1; j <= i2; j++ ) {
				retr_edr[j]		= pow(avg,3);
				retr_edr13_err[j] = sqrt(var);
			}                             			
		}
	}
	
	free(abs_dl);
	free(abs_dt);
	free(v2_avg);
	free(v2_var);
	free(fft_pow);
	free(fft_freq);
	free(lst_eps);
	free(lst_eps13);        
}

void retr_edr_2nd_order_structure_function(
	int 	*n,				//length of variables v,t,l,v2
	double 	*v,
	double 	*t,
	double 	*l,
	double 	*v2,
	int 	*retr_domain,		//0 = space, 1 = time and U0 from v2
	double 	*retr_C,			//Kolmogorov constant that is used.
	int		*retr_n,			//Number of points used for average and variance calculation
	int		*periodic,			//1 = periodic, 0 = non-periodic
	double	*retr_edr,			//output: retrieved edr
	double	*retr_edr13_err		//error in retrieved edr^(1/3) (std)
	)
{
	int order = 2;	
	
	int i, i1, i2, ni, ni_2, ni_2m;
	int j;
	int calc_variance;
	
	double *abs_dl 		= malloc((*n-1) * sizeof(double));
	double *abs_dt 		= malloc((*n-1) * sizeof(double));
	double *v2_avg 		= malloc(*n * sizeof(double));
	double *v2_var 		= malloc(*n * sizeof(double));
	double *dll	= malloc(*n * sizeof(double));

	double dlmax, dlmin;
	double dtmax, dtmin;
	
    double *lst_eps		= malloc(*n * sizeof(double));
    double *lst_eps13		= malloc(*n * sizeof(double));
	double  avg, var;
	int 	nok;
       
    for ( i = 0; i < (*n - 1); i++ ) {
		abs_dl[i] = fabs(l[i+1] - l[i]);
		abs_dt[i] = fabs(t[i+1] - t[i]);		
	}

	//calculate running average and running std.
	if (*retr_domain == 1) {
		//for time domain, we need the running average and running variance of v2
		calc_variance = 1;
		running_average(n, v2, retr_n, v2_avg, v2_var, &calc_variance);		
	}
        
    for (i1 = 0; i1 < *n; i1 = i1 + *retr_n ) {
		i2 = i1 + (*retr_n -1);

		if (i2 >= *n) {
			i2 = *n -1;
            i1 = fmax(0, i2 - (*retr_n - 1));
		}
		ni = i2 + 1 - i1;
		ni_2 = ni/2;
		ni_2m = ni_2 - 1;
		
		if (*retr_domain == 0) {
			//space domain
			dlmax	= fabs(l[i2] - l[i1]);
			dlmin	= dlmax / (ni-1.);

			//apply second order structure function
			structure_function(&ni, v + i1, dll, &order, periodic);

			//calculate epsilon for each point
			for (j = 1; j < ni_2; j++ ) {
				lst_eps[j] 		= pow(dll[j] / *retr_C, 3./2.) / (dlmin * j);
				lst_eps13[j]	= pow(lst_eps[j], 1./3.);
			}
						
            //calculate statistics of list of edrs
            calc_avg_variance(&ni_2m, lst_eps13 + 1, &avg, &var, &nok);
			for (j = i1; j <= i2; j++ ) {
				retr_edr[j]			= pow(avg,3.);
				retr_edr13_err[j] 	= sqrt(var);
			}
		} else {
			//time domain
			dtmax	= fabs(t[i2] - t[i1]);
			dtmin	= dtmax / (ni-1.);

			//apply second order structure function
			structure_function(&ni, v + i1, dll, &order, periodic);

			//calculate epsilon for each point
			for (j = 1; j < ni_2; j++ ) {
				lst_eps[j] 		= pow(dll[j] / *retr_C, 3./2.) / (fabs(v2_avg[j]) * dtmin * j);
				lst_eps13[j]	= pow(lst_eps[j], 1./3.);
			}
						
            //calculate statistics of list of edrs
            calc_avg_variance(&ni_2m, lst_eps13 + 1, &avg, &var, &nok);
			for (j = i1; j <= i2; j++ ) {
				retr_edr[j]			= pow(avg,3.);
				retr_edr13_err[j] 	= sqrt(var);
			}                                       			
		}
	}
	
	free(abs_dl);
	free(abs_dt);
	free(v2_avg);
	free(v2_var);
	free(dll);
	free(lst_eps);
	free(lst_eps13);        

}

void retr_edr_3rd_order_structure_function(
	int 	*n,				//length of variables v,t,l,v2
	double 	*v,
	double 	*t,
	double 	*l,
	double 	*v2,
	int 	*retr_domain,		//0 = space, 1 = time and U0 from v2
	double 	*retr_C,			//Kolmogorov constant that is used.
	int		*retr_n,			//Number of points used for average and variance calculation
	int		*periodic,			//1 = periodic, 0 = non-periodic
	double	*retr_edr,			//output: retrieved edr
	double	*retr_edr13_err		//error in retrieved edr^(1/3) (std)
	)
{
	int order = 3;	
	
	int i, i1, i2, ni, ni_2, ni_2m;
	int j;
	int calc_variance;
	
	double *abs_dl 		= malloc((*n-1) * sizeof(double));
	double *abs_dt 		= malloc((*n-1) * sizeof(double));
	double *v2_avg 		= malloc(*n * sizeof(double));
	double *v2_var 		= malloc(*n * sizeof(double));
	double *dlll 	= malloc(*n * sizeof(double));

	double dlmax, dlmin;
	double dtmax, dtmin;
	
    double *lst_eps		= malloc(*n * sizeof(double));
    double *lst_eps13		= malloc(*n * sizeof(double));
	double  avg, var;
	int 	nok;
       
    for ( i = 0; i < (*n - 1); i++ ) {
		abs_dl[i] = fabs(l[i+1] - l[i]);
		abs_dt[i] = fabs(t[i+1] - t[i]);		
	}

	//calculate running average and running std.
	if (*retr_domain == 1) {
		//for time domain, we need the running average and running variance of v2
		calc_variance = 1;
		running_average(n, v2, retr_n, v2_avg, v2_var, &calc_variance);		
	}
        
    for (i1 = 0; i1 < *n; i1 = i1 + *retr_n ) {
		i2 = i1 + (*retr_n -1);

		if (i2 >= *n) {
			i2 = *n -1;
            i1 = fmax(0, i2 - (*retr_n - 1));
		}
		ni = i2 + 1 - i1;
		ni_2 = ni/2.;
		ni_2m = ni_2 - 1;
		
		if (*retr_domain == 0) {
			//space domain
			dlmax	= fabs(l[i2] - l[i1]);
			dlmin	= dlmax / (ni-1);

			//apply third structure function
			structure_function(&ni, v + i1, dlll, &order, periodic);

			//calculate epsilon for each point
			for (j = 1; j < ni_2; j++ ) {
				lst_eps[j] 		= *retr_C * dlll[j] / (dlmin * j);
				lst_eps13[j]	= ((lst_eps[j] > 0.) ? 1. : -1.) *
									pow(fabs(lst_eps[j]), 1./3.);
			}
						
            //calculate statistics of list of edrs
            calc_avg_variance(&ni_2m, lst_eps13 + 1, &avg, &var, &nok);
			for (j = i1; j <= i2; j++ ) {
				retr_edr[j]		= pow(avg,3);
				retr_edr13_err[j] = sqrt(var);
			}
		} else {
			//time domain
			dtmax	= fabs(t[i2] - t[i1]);
			dtmin	= dtmax / (ni - 1);

			//apply third order structure function
			structure_function(&ni, v + i1, dlll, &order, periodic);

			//calculate epsilon for each point
			for (j = 1; j < ni_2; j++ ) {
				lst_eps[j] 		= *retr_C * dlll[j] / (fabs(v2_avg[j]) * dtmin * j);
				lst_eps13[j]	= ((lst_eps[j] > 0.) ? 1. : -1.) *
									pow(fabs(lst_eps[j]), 1./3.);
			}
						
            //calculate statistics of list of edrs
            calc_avg_variance(&ni_2m, lst_eps13 + 1, &avg, &var, &nok);
			for (j = i1; j <= i2; j++ ) {
				retr_edr[j]		= pow(avg,3.);
				retr_edr13_err[j] = sqrt(var);
			}
		}
	}
	
	free(abs_dl);
	free(abs_dt);
	free(v2_avg);
	free(v2_var);
	free(dlll);
	free(lst_eps);
	free(lst_eps13);        

}

//Kolmogorov skewness:
//6.89: S = DLLL(t) / (DLL(t) ** 3/2)
//Kolmogorov constant C:
//6.90: C = (-4 / (5 S)) ** (2/3)
void kolmogorov_skewness(
	int 	*n,
	double 	*x,
	double 	*y,
	int		*periodic
	)
{
	double *dll = malloc(*n * sizeof(double));
	double *dlll = malloc(*n * sizeof(double));
	int order;
	int i;
	
	order = 2;
	structure_function(n, x, dll, &order, periodic);
	order = 3;
	structure_function(n, x, dlll, &order, periodic);
	
	for (i = 0; i < *n; i++ ) {
		y[i] = dlll[i] / pow(dll[i], 3. /2. );
	}
	
	free(dll);
	free(dlll);
}

//retrieve the kolmogorov skewness, kolmogorov constant
void retr_skewness_structure_functions(
	int		*n,
	double 	*v,
	int		*retr_n,
	int		*periodic,
	double	*retr_skewness,
	double	*retr_constant
		)
{
	int i1, i2, nok;
	int ni, ni_2, ni_2m;
	
	double *lst_skewness 	= malloc(*retr_n * sizeof(double));
	double *lst_const 		= malloc(*retr_n * sizeof(double));
	
	int j;
	double avg;
	    
	for (i1 = 0; i1 < *n; i1 = i1 + *retr_n ) {
		i2 = i1 + (*retr_n -1);

		if (i2 >= *n) {
			i2 = *n -1;
			i1 = fmax(0, i2 - (*retr_n - 1));
		}
		ni = i2 + 1 - i1;
		ni_2 = ni/2.;
		ni_2m = ni_2 - 1;
		
		//Kolmogorov skewness:
		//6.89: S = DLLL(t) / (DLL(t) ** 3/2)
		//Kolmogorov constant C:
		//6.90: C = (-4 / (5 S)) ** (2/3)
		kolmogorov_skewness(&ni, v + i1, lst_skewness, periodic);

		for (j = 1; j < ni_2; j++ ) {
			lst_const[j] = pow(fabs(-4. / (5. * lst_skewness[j])) , 2. / 3.);
		}
			
		calc_avg(&ni_2m, v + 1, &avg, &nok);
		for (j = i1; j <= i2; j++ ) {            
			retr_skewness[j] = avg;
		}
		
		calc_avg(&ni_2m, lst_const + 1, &avg, &nok);
		for (j = i1; j <= i2; j++ ) {
			retr_constant[j] = avg;
		}
	}
	
	free(lst_skewness);    
	free(lst_const);    
}


//retrieve the fraction S, kolmogorov constant
void retr_factor_wind_components(
	int		*n,
	double	*u,
	double	*v,
	double 	*w,
	int		*retr_n,
	int		*periodic,
	double	*retr_frachor,	//retrieved factor between components
	double 	*retr_fracver	//retrieved factor between components
	)
{
	int order = 2;
	int i1, i2, nok;
	int ni, ni_2, ni_2m;

	double *u_dir 		= malloc(*n * sizeof(double));

	double *u_lon 		= malloc(*n * sizeof(double));
	double *u_tra1 		= malloc(*n * sizeof(double));
	double *u_tra2 		= malloc(*n * sizeof(double));

	double *dll_lon 	= malloc(*n * sizeof(double));
	double *dll_tra1 	= malloc(*n * sizeof(double));
	double *dll_tra2 	= malloc(*n * sizeof(double));

	double *lst_frac 	= malloc(*n * sizeof(double));
	
	int j;
	double avg;
	    
	//Calculate longitudinal and transverse components of the wind.
	for (j = 0; j <= *n; j++ ) {            
		//wind vector direction w.r.t. east (u-direction)
		u_dir[j] = atan2(v[j], u[j]);
		//now calculate longitudinal component
		u_lon[j] = (cos(u_dir[j]) * u[j]) + (sin(u_dir[j]) * v[j]);

		//and transverse componenents			
		u_tra1[j] = (sin(u_dir[j]) * u[j]) + (sin(u_dir[j]) * v[j]);
		u_tra2[j] = w[j];
	}
            

	for (i1 = 0; i1 < *n; i1 = i1 + *retr_n ) {
		i2 = i1 + (*retr_n -1);

		if (i2 >= *n) {
			i2 = *n -1;
			i1 = fmax(0, i2 - (*retr_n - 1));
		}
		ni = i2 + 1 - i1;
		ni_2 = ni/2.;
		ni_2m = ni_2 - 1;

		//calculate second order structure function.
		structure_function(&ni, u_lon + i1, dll_lon, &order, periodic);
		structure_function(&ni, u_tra1 + i1, dll_tra1, &order, periodic);
		structure_function(&ni, u_tra2 + i1, dll_tra2, &order, periodic);
		
		//using horizontal components
		for (j = 0; j < ni_2; j++ ) {
			lst_frac[j] = dll_lon[j] / dll_tra1[j];  //in theory 4/3
        }
                    
        calc_avg(&ni_2, lst_frac, &avg, &nok);
        for (j = i1; j <= i2; j++ ) {            
			retr_frachor[j] = avg;
		}
        
		//using horizontal + vertical components
		for (j = 0; j < ni_2; j++ ) {
			lst_frac[j] = dll_lon[j] / dll_tra2[j];  //in theory 4/3
        }
        calc_avg(&ni_2, lst_frac, &avg, &nok);
        for (j = i1; j <= i2; j++ ) {            
			retr_fracver[j] = avg;
		}
	}

	free(u_dir);
	
	free(u_lon);
	free(u_tra1);
	free(u_tra2);

	free(dll_lon);
	free(dll_tra1);
	free(dll_tra2);

	free(lst_frac);
}

void kolmogorov_constants(
	char   *choice,
	double *power_c,
	double *struc2_c,
	double *struc3_c
	)
{
	double c, q, c1_div_c2, c2_div_c1;
	
	c 			= 1.5;
	q 			= 2./3.;  //~0.66
	c1_div_c2 	= (1. / M_PI) * lgamma(1.+q) * sin(M_PI * q / 2.);	//~0.25
	c2_div_c1	= 1. / c1_div_c2;  //~4.
		
	if (!strcmp(choice, "full")) {
		*power_c 	= c;
		*struc2_c 	= c2_div_c1 * c;
		*struc3_c 	= NAN;
	}

	if (!strcmp(choice, "longitudinal")) {
		*power_c 	= (18./55.) * c;
		*struc2_c 	= c2_div_c1 * (18./55.) * c;
		*struc3_c 	= -4./5.;
	}

	if (!strcmp(choice, "transverse")) {
		*power_c 	= (4./3.) * (18./55.) * c;
		*struc2_c 	= (4./3.) * c2_div_c1 *  (18./55.) * c;
		*struc3_c 	= -4./5.;		
	}		
}

void radial_kolmogorov_constants(
	double *azimuth_rad,
	double *azimuth0_rad,
	double *elevation_rad,
	double *power_c,
	double *struc2_c,
	double *struc3_c
	)
{
	double delta_azimuth_rad;
	double longi_power_c, longi_struc2_c, longi_struc3_c;
	double trans_power_c, trans_struc2_c, trans_struc3_c;
	
	delta_azimuth_rad = azimuth_rad - azimuth0_rad;
	
	kolmogorov_constants("longitudinal", &longi_power_c, &longi_struc2_c, &longi_struc3_c);
	kolmogorov_constants("transverse",   &trans_power_c, &trans_struc2_c, &trans_struc3_c);

	*power_c 	=  	(pow(cos(*elevation_rad) * cos(delta_azimuth_rad), 2) * longi_power_c)
					+ (pow(cos(*elevation_rad) * sin(delta_azimuth_rad), 2) * trans_power_c)
					+ (pow(sin(*elevation_rad), 2) * trans_power_c);
	*struc2_c 	=   (pow(cos(*elevation_rad) * cos(delta_azimuth_rad), 2) * longi_struc2_c)
					+ (pow(cos(*elevation_rad) * sin(delta_azimuth_rad), 2) * trans_struc2_c)
					+ (pow(sin(*elevation_rad),2) * trans_struc2_c);
	*struc3_c	= -4./5.;
	
}

