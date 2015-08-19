%module wrapedr

%{
    #define SWIG_FILE_WITH_INIT
	#include "edr.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len200, double *retr_av200),
						(int len201, double *retr_var201),
						(int len202, double *retr_edr202),
						(int len203, double *retr_edr13_err203)
						};
%apply (int DIM1, double* IN_ARRAY1) {
						(int len100, double *v100),
						(int len101, double *t101),
						(int len102, double *l102),
						(int len103, double *v2103),
						(int len104, double *retr_domain104),
						(int len105, double *retr_C105),
						(int len106, double *retr_n106)
						};

%rename (retr_edr_variance) my_retr_edr_variance;

%inline %{
	double my_retr_edr_variance(
		int len100, double *v100,
		int len101, double *t101,
		int len102, double *l102,
		int len103, double *v2103,
		int len104, double *retr_domain104,
		int len105, double *retr_C105,
		int len106, double *retr_n106,
		int len200, double *retr_av200,
		int len201, double *retr_var201,
		int len202, double *retr_edr202,
		int len203, double *retr_edr13_err203)
	{

	int int_retr_domain 	= (int) *retr_domain104;
	int int_retr_n 			= (int) *retr_n106;
	
	retr_edr_variance(	&len100, v100, t101, l102, v2103,
						&int_retr_domain, retr_C105,
						&int_retr_n, retr_av200,
						retr_var201, retr_edr202, retr_edr13_err203);

	}
%}







%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len400, double *retr_edr400),
						(int len401, double *retr_edr13_err401)
						};
%apply (int DIM1, double* IN_ARRAY1) {
						(int len300, double *v300),
						(int len301, double *t301),
						(int len302, double *l302),
						(int len303, double *v2303),
						(int len304, double *retr_domain304),
						(int len305, double *retr_C305),
						(int len306, double *retr_n306),
						(int len307, double *periodic307),
						(int len308, double *retr_nintervals308)
						};

%rename (retr_edr_powerspectrum) my_retr_edr_powerspectrum;

%inline %{
	double my_retr_edr_powerspectrum(
		int len300, double *v300,
		int len301, double *t301,
		int len302, double *l302,
		int len303, double *v2303,
		int len304, double *retr_domain304,
		int len305, double *retr_C305,
		int len306, double *retr_n306,
		int len307, double *periodic307,
		int len308, double *retr_nintervals308,
		int len400, double *retr_edr400,
		int len401, double *retr_edr13_err401
	)
	{
	int int_retr_domain 	= (int) *retr_domain304;
	int int_retr_n 			= (int) *retr_n306;
	int int_periodic 		= (int) *periodic307;
	int int_nintervals 		= (int) *retr_nintervals308;
	
	retr_edr_powerspectrum(	&len300, v300, t301, l302, v2303,
						&int_retr_domain, retr_C305,
						&int_retr_n, &int_periodic, 
						&int_nintervals, retr_edr400, retr_edr13_err401);
	}
%}




%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len600, double *retr_edr600),
						(int len601, double *retr_edr13_err601)
						};
%apply (int DIM1, double* IN_ARRAY1) {
						(int len500, double *v500),
						(int len501, double *t501),
						(int len502, double *l502),
						(int len503, double *v2503),
						(int len504, double *retr_domain504),
						(int len505, double *retr_C505),
						(int len506, double *retr_n506),
						(int len507, double *periodic507),
						(int len508, double *retr_nintervals508)
						};

%rename (retr_edr_2nd_order_structure_function) my_retr_edr_2nd_order_structure_function;

%inline %{
	double my_retr_edr_2nd_order_structure_function(
		int len500, double *v500,
		int len501, double *t501,
		int len502, double *l502,
		int len503, double *v2503,
		int len504, double *retr_domain504,
		int len505, double *retr_C505,
		int len506, double *retr_n506,
		int len507, double *periodic507,
		int len600, double *retr_edr600,
		int len601, double *retr_edr13_err601
	)
	{
	int int_retr_domain 	= (int) *retr_domain504;
	int int_retr_n 			= (int) *retr_n506;
	int int_periodic 		= (int) *periodic507;
	
	retr_edr_2nd_order_structure_function(	&len500, v500, t501, l502, v2503,
						&int_retr_domain, retr_C505,
						&int_retr_n, &int_periodic, 
						retr_edr600, retr_edr13_err601);
	}
%}




%rename (retr_edr_3rd_order_structure_function) my_retr_edr_3rd_order_structure_function;

%inline %{
	double my_retr_edr_3rd_order_structure_function(
		int len500, double *v500,
		int len501, double *t501,
		int len502, double *l502,
		int len503, double *v2503,
		int len504, double *retr_domain504,
		int len505, double *retr_C505,
		int len506, double *retr_n506,
		int len507, double *periodic507,
		int len600, double *retr_edr600,
		int len601, double *retr_edr13_err601
	)
	{
	int int_retr_domain 	= (int) *retr_domain504;
	int int_retr_n 			= (int) *retr_n506;
	int int_periodic 		= (int) *periodic507;
	
	retr_edr_3rd_order_structure_function(	&len500, v500, t501, l502, v2503,
						&int_retr_domain, retr_C505,
						&int_retr_n, &int_periodic, 
						retr_edr600, retr_edr13_err601);
	}
%}
