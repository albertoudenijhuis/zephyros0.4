#ifndef _ZEPHYROS_TURBULENCE
#define _ZEPHYROS_TURBULENCE

#include <complex.h>
#include "fields.h"
#include "interpolation.h"

typedef struct t_zephyros_turbulence_widget
{
	int			type;				//(1 = Mann1998, 2 = CTM, 3 = Careta1993, 4 = Pinsky2006, 5 = parametric)
	t_zephyros_field 	*field;
	double 		*grid_edr;
	double 		*grid_edr13;
	t_zephyros_interpolation_bilint_lut	*lut_edr13;

	double 		Lx;							
	double 		Ly;							
	double 		Lz;							
	int			Nx;							
	int			Ny;							
	int			Nz;							
	double		*freqx;						
	double		*freqy;						
	double		*freqz;						
		
	double 		*grid_karman_L;						//mann1998, ctm
	double 		*grid_kolmogorov_constant;			//mann1998, ctm
	t_zephyros_interpolation_bilint_lut	*lut_karman_L;
	t_zephyros_interpolation_bilint_lut	*lut_kolmogorov_constant;
	
	double complex		*random_fourier_cx;		//ctm
	double complex		*random_fourier_cy;		//ctm
	double complex		*random_fourier_cz;		//ctm
	double 				*random_shift_x;		//ctm
	double 				*random_shift_y;		//ctm
	double 				*random_shift_z;		//ctm
		
	double complex		***random_fourier_c0;	//mann1998
	double complex 		***random_fourier_c1;	//mann1998
	double complex 		***random_fourier_c2;	//mann1998
	double complex 		***u_fourier;	//mann1998
	double complex 		***v_fourier;	//mann1998
	double complex 		***w_fourier;	//mann1998
	
	double complex		**random_fourier_beta;	//careta1993
	
	char		pinsky2006_file[8192];					//pinsky2006
	
	double 		minL_div_maxL;						//ctm
	
	int			random_numbers;					//mann1998, ctm
	char		rn_file[8192];					//mann1998, ctm

	double 		lambdax; 						//careta1993
	double 		lambday; 						//careta1993
	
	int			M;								//pinsky2006
	int			K;								//pinsky2006
	double 		**random_a;						//pinsky2006
	double 		**random_b;						//pinsky2006
	double 		*lambdak;						//pinsky2006
	double 		*alphak;						//pinsky2006
	
	double		calibration_factor;
	int			calibration_method;				
	int			calibration_n;
	double		calibration_C;
	int			calibration_periodic;
	int			calibration_nint;
	int			calibration_dir;
	double		calibration_L;
} t_zephyros_turbulence_widget;


typedef struct t_turbulence_karmanspec
{
	double	a;
	double	L;
	double	edr;
} t_turbulence_karmanspec;

void turbulence_prepare_widget(t_zephyros_turbulence_widget *cfg);
void turbulence_free_widget(t_zephyros_turbulence_widget **pcfg);

void turbulence_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double uvw_derivatives[4]);
	
void turbulence_mann1998_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double uvw_derivatives[4]);

void turbulence_mann1998_uvw_particle(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double uvw_derivatives[4],
	double 	*D_maj_mm);

void turbulence_ctm_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double uvw_derivatives[4]);

void turbulence_careta93_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double uvw_derivatives[4]);
	
void turbulence_pinsky2006_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double uvw_derivatives[4]);
	
double turbulence_mann1998_C(
	t_turbulence_karmanspec *sp,
	double kx,
	double ky,
	double kz,
	int i,
	int j,
	double deltax,
	double deltay,
	double deltaz);
	
double turbulence_mann1998_A(
	t_turbulence_karmanspec *sp,
	double kx,
	double ky,
	double kz,
	int i,
	int j);
	
double turbulence_mann1998_Ek(
	t_turbulence_karmanspec *sp,
	double kabs);
	
double turbulence_ctm_C(
	t_turbulence_karmanspec *sp,
	double kabs,
	int 	Nx,
	double Lmax);
	
double turbulence_careta1993_Q(
	t_zephyros_turbulence_widget *cfg,
	double kx,
	double ky);
	
void turbulence_calibrate(t_zephyros_turbulence_widget *cfg);

#endif
