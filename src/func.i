%module func

%{
    #define SWIG_FILE_WITH_INIT
	#include "func.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}


%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len101, double* avg101),
						(int len102, double* variance102),
						(int len103, double* nok103)
							};



%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len201, double* avg201),
						(int len202, double* variance202),
						(int len203, double* nok203)
							};


%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len302, double* running_avg302),
						(int len303, double* running_variance303)
							};


%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len402, double* running_avg402),
						(int len403, double* running_variance403)
							};


%apply (int DIM1, double* IN_ARRAY1) {
						(int len100, double* dat100)
						};

%apply (int DIM1, double* IN_ARRAY1) {
						(int len200, double* dat200)
						};


%apply (int DIM1, double* IN_ARRAY1) {
						(int len300, double* dat300),
						(int len301, double* nav_in301)
						};


%apply (int DIM1, double* IN_ARRAY1) {						
						(int len400, double* dat400),					
						(int len401, double* nav_in401)
						};


%rename (calc_avg_variance) my_calc_avg_variance;
%rename (calc_avg_variance_ang) my_calc_avg_variance_ang;
%rename (running_average) my_running_average;
%rename (running_average_ang) my_running_average_ang;

%inline %{
	double my_calc_avg_variance(
		int	len100, double* dat100,
		int	len101, double* avg101,
		int	len102, double* variance102,
		int len103, double* nok103
		) 
	{
		int int_nok;
		calc_avg_variance(&len100, dat100, avg101, variance102, &int_nok);
		*nok103 = int_nok;
	}
		
	double my_calc_avg_variance_ang(
		int	len200,	double* dat200,
		int	len201,	double* avg201,
		int	len202,	double* variance202,
		int len203,	double* nok203
		)
	{
		int int_nok;	
		calc_avg_variance_ang(&len200, dat200, avg201, variance202, &int_nok);
		*nok203 = int_nok;
	}
	
	
	
	double my_running_average(
		int	len300,	double* dat300,
		int	len301,	double* nav_in301,
		int	len302,	double* running_avg302,
		int	len303,	double* running_variance303
		)
	{		
		int nav = (int) *nav_in301;
		int calc_var = 1;
		running_average(&len300, dat300, &nav, running_avg302, running_variance303, &calc_var);
	}
	
		
	double my_running_average_ang(
		int	len400,	double* dat400,
		int	len401,	double* nav_in401,		
		int	len402,	double* running_avg402,
		int	len403,	double* running_variance403)
	{
		int nav = (int) *nav_in401;
		int calc_var = 1;
		running_average_ang(&len400, dat400, &nav,	running_avg402, running_variance403, &calc_var);
	}

%}
