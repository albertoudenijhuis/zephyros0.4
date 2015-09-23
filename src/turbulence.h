#ifndef _ZEPHYROS_TURBULENCE
#define _ZEPHYROS_TURBULENCE

#include <complex.h>
#include "fields.h"
#include "interpolation.h"
#include "util.h"

void turbulence_initialize_widget(t_zephyros_turbulence_widget **pcfg);
void turbulence_prepare_widget(t_zephyros_turbulence_widget *cfg);
void turbulence_free_widget(t_zephyros_turbulence_widget **pcfg);

typedef struct t_turbulence_karmanspec
{
	double	a;
	double	L;
	double	edr;
} t_turbulence_karmanspec;

void turbulence_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);
	
void turbulence_mann1998_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);

void turbulence_mann1998_uvw_particle(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives,
	double 	*D_maj_mm);

void turbulence_ctm_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);

void turbulence_careta93_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);
	
void turbulence_pinsky2006_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);
	
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
