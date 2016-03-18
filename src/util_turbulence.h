#ifndef _ZEPHYROS_TURBULENCE
#define _ZEPHYROS_TURBULENCE

#include <complex.h>
#include "fields.h"
#include "interpolation.h"
#include "util.h"

void util_initialize_turbulence_widget(t_zephyros_turbulence_widget **pcfg);
void util_prepare_turbulence_widget(t_zephyros_turbulence_widget *cfg);
void util_free_turbulence_widget(t_zephyros_turbulence_widget **pcfg);

typedef struct t_util_turbulence_karmanspec
{
	double	a;
	double	L;
	double	edr;
} t_util_turbulence_karmanspec;

void util_turbulence_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);
	
void util_turbulence_mann1998_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);

void util_turbulence_mann1998_uvw_particle(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives,
	double 	*D_maj_mm);

void util_turbulence_ctm_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);

void util_turbulence_careta93_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);
	
void util_turbulence_pinsky2006_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives);
	
double util_turbulence_mann1998_C(
	t_util_turbulence_karmanspec *sp,
	double kx,
	double ky,
	double kz,
	int i,
	int j,
	double deltax,
	double deltay,
	double deltaz);
	
double util_turbulence_mann1998_A(
	t_util_turbulence_karmanspec *sp,
	double kx,
	double ky,
	double kz,
	int i,
	int j);
	
double util_turbulence_mann1998_Ek(
	t_util_turbulence_karmanspec *sp,
	double kabs);
	
double util_turbulence_ctm_C(
	t_util_turbulence_karmanspec *sp,
	double kabs,
	int 	Nx,
	double Lmax);
	
double util_turbulence_careta1993_Q(
	t_zephyros_turbulence_widget *cfg,
	double kx,
	double ky);
	
void util_turbulence_calibrate(t_zephyros_turbulence_widget *cfg);

#endif
