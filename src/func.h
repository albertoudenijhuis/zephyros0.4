#ifndef _ZEPHYROS_FUNC
#define _ZEPHYROS_FUNC

#include <complex.h>

void f_tau(
	int 	*n,
	int 	*tau
	);
	
void f_fftfreq(
	int 	*n,
	double 	*d,
	double	*freq
	);
	
int isnanorinf(
	double *x);

double pnan(double *x);
	
void calc_avg(
	int		*n,
	double 	*dat,
	double	*avg,
	int		*nok
	);
	
void calc_avg_ang(
	int		*n,
	double 	*dat,
	double	*avg,
	int		*nok
	);
	
void calc_avg_variance(
	int		*n,
	double 	*dat,
	double	*avg,
	double	*variance,
	int		*nok
	);
	
void calc_avg_variance_ang(
	int		*n,
	double 	*dat,
	double	*avg,
	double	*variance,
	int		*nok
	);

void running_average(
	int		*np,			//number of points
	double 	*dat,
	int		*nav_in,		//number of points to do averaging
	double	*running_avg,
	double 	*running_variance,
	int		*calc_variance
	);
	
void running_average_ang(
	int		*np,			//number of points
	double 	*dat,
	int		*nav_in,		//number of points to do averaging
	double	*running_avg,
	double 	*running_variance,
	int		*calc_variance
	);
	

void rpowspec(
	int 	*n,
	double 	*xreal,
	double 	*pow,
	int		*periodic
	);
	
void rifft(
	int 	*n,
	double complex  *c,
	double *xreal);
	
void rfft(
	int		*n,
	double 	*xreal,
	double complex  *c
	);
	
void fft(
	int		*n,
	double complex *x,
	double complex *c);
	
void ifft(
	int		*n,
	double complex *c,
	double complex *x);
	
void covariance(
	int		*n,
	double 	*x,
	double 	*y,
	double ans);
	
void autocovariance(
	int 	*n,
	double	*x,
	double 	*y,
	int		*periodic);
	
void autocorrelation(
	int 	*n,
	double	*x,
	double	*y,
	int		*periodic);
	
void structure_function(
	int		*n,
	double 	*x,
	double 	*y,
	int		*order,
	int		*periodic);



void randompermutation(int n, int **output);

double angleAminB(double *Arad, double *Brad);

void dbldotprod(
	int*	n,
	double*	arr1,
	double*	arr2,
	double*	result);
	

void func_dbl_arr_malloc(
	int n,
	double **ptrarr
	);

void func_int_arr_malloc(
	int n,
	int **ptrarr
	);

void func_ptr_arr_malloc(
	int n,
	int **ptrarr
	);

double uniform_random();

double read_uniform_random(FILE *fp);

double func_newton_rapson_solver(
	double (*f)(double, void*),
	double (*df)(double, void*),
	double x0_in,
	double allowed_error,
	int	max_itr,
	void *vp
	);
	
double func_dB(double val);
double func_dB_inv(double val);

double func_modulo(double x, double y);

double func_norm(int *n, double *arr);
	
#endif
