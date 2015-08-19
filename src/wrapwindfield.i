%module wrapwindfield

%{
    #define SWIG_FILE_WITH_INIT
	#include "wrapwindfield.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len200, double *u),
						(int len201, double *v),
						(int len202, double *w)
								};
								
%apply (int DIM1, double* IN_ARRAY1) {
						(int len100, double *x),
						(int len101, double *y),
						(int len102, double *z),
						(int len103, double *t)
						};

%rename (windfield) my_windfield;

%inline %{
	
	double my_windfield(
		char cfg_filename[8192],
		char additional_output_filename[8192],
		int len100, double *x,
		int len101, double *y,
		int len102, double *z,
		int len103, double *t,
		int len200, double *u,
		int len201, double *v,
		int len202, double *w)
	{

	wrapwindfield(cfg_filename, additional_output_filename, len100, x, y, z, t, u, v, w);
	
	}
%}






