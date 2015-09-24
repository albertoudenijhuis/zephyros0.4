#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <math.h>

#include "util.h"
#include "specialfunctions.h"
#include "func.h"

//uncomment next statement for debug mode
#define _ZEPHYROS_UTIL_DEBUG

void util_field2lut(t_zephyros_field *field, double *variable, int special, t_zephyros_interpolation_bilint_lut **plut)
{
	t_zephyros_interpolation_bilint_lut *lut;
	int i;
	
	interpolation_initialize_lut(&lut);
	
	//LUT interpolation stuff
	lut->n_dim 	= 4;
	lut->n		= field->n;
	lut->shape 	= malloc(lut->n_dim * sizeof(int));	
	lut->mesy_y_allocated = 0;
	lut->periodic = 0;
	
	lut->shape[0] = field->n_x;
	lut->shape[1] = field->n_y;
	lut->shape[2] = field->n_z;
	lut->shape[3] = field->n_t;
				
	lut->ax_values = malloc(lut->n_dim * sizeof(double*));
	lut->ax_values[0] = malloc(lut->shape[0] * sizeof(double));
	lut->ax_values[1] = malloc(lut->shape[1] * sizeof(double));
	lut->ax_values[2] = malloc(lut->shape[2] * sizeof(double));
	lut->ax_values[3] = malloc(lut->shape[3] * sizeof(double));

	for ( i = 0; i < field->n_x; i++ ) {
		lut->ax_values[0][i] = field->vec_x[i];
	}
	for ( i = 0; i < field->n_y; i++ ) {
		lut->ax_values[1][i] = field->vec_y[i];
	}
	for ( i = 0; i < field->n_z; i++ ) {
		lut->ax_values[2][i] = field->vec_z[i];
	}
	for ( i = 0; i < field->n_t; i++ ) {
		lut->ax_values[3][i] = field->vec_t[i];
	}
	
	lut->special = special;
	
	if (lut->special == 2) {
		lut->mesh_y = malloc(lut->n * sizeof(double));
		lut->mesy_y_allocated = 1;
		for ( i = 0; i < lut->n ; i++ ) {
			lut->mesh_y[i] = fabs(variable[i]);
		}		
		lut->special = 0;
	}
	if (lut->special == 3) {
		lut->mesh_y = malloc(lut->n * sizeof(double));
		lut->mesy_y_allocated = 1;
		for ( i = 0; i < lut->n ; i++ ) {
			lut->mesh_y[i] = log(variable[i]);
		}		
		lut->special = 0;		
	} else {
		lut->mesh_y = variable;
	}

	*plut = lut;
}



void util_initialize_windfield(t_zephyros_windfield **pwindfield)
{
	t_zephyros_windfield *windfield = calloc(1, sizeof(t_zephyros_windfield));
	int i;

	*pwindfield = windfield;	
	for ( i = 0; i <= 100; i++ ) {
		windfield->type 	= 0;
		windfield->field[i] 				= NULL;
		
		windfield->grid_hspeed[i] 			= NULL;
		windfield->grid_hdir[i] 			= NULL;
		windfield->grid_u[i] 				= NULL;
		windfield->grid_v[i] 				= NULL;
		windfield->grid_w[i] 				= NULL;
		windfield->vortex[i] 				= NULL;
		windfield->wave[i] 					= NULL;
		windfield->turbulence[i] 			= NULL;
		
		windfield->grid_hspeed_err[i] 		= NULL;
		windfield->grid_hdir_err[i] 		= NULL;
		windfield->grid_u_err[i] 			= NULL;
		windfield->grid_v_err[i] 			= NULL;
		windfield->grid_w_err[i] 			= NULL;
		
		windfield->hdir_ecm[i].mat 			= NULL;
		windfield->hdir_ecm[i].mat_S		= NULL;
		windfield->hdir_ecm[i].mat_N		= NULL;

		windfield->hspeed_ecm[i].mat 		= NULL;
		windfield->hspeed_ecm[i].mat_S		= NULL;
		windfield->hspeed_ecm[i].mat_N		= NULL;

		windfield->u_ecm[i].mat 		= NULL;
		windfield->u_ecm[i].mat_S		= NULL;
		windfield->u_ecm[i].mat_N		= NULL;

		windfield->v_ecm[i].mat 		= NULL;
		windfield->v_ecm[i].mat_S		= NULL;
		windfield->v_ecm[i].mat_N		= NULL;

		windfield->w_ecm[i].mat 		= NULL;
		windfield->w_ecm[i].mat_S		= NULL;
		windfield->w_ecm[i].mat_N		= NULL;
		
		windfield->fit_hspeed[i]		= 0;
		windfield->fit_hdir[i]		= 0;
		windfield->fit_u[i]		= 0;
		windfield->fit_v[i]		= 0;
		windfield->fit_w[i]		= 0;
			
	}
	windfield->nfields = 0;
	windfield->nvortices = 0;
	windfield->nwaves = 0; 
	windfield->nturbulences = 0; 
}


void util_prepare_windfield(t_zephyros_windfield *windfield)
{
	int i;
	double tmp;
	
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_prepare_windfield\n"); fflush(stdout);
	#endif	
	//wind field
	for ( i = 0; i < windfield->nfields; i++ ) {
		if (windfield->field[i] != NULL) {
			fields_prepare(windfield->field[i]);
			if (windfield->grid_hspeed[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_hspeed[i], 0, windfield->lut_hspeed + i);
			if (windfield->grid_hdir[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_hdir[i], 0, windfield->lut_hdir + i);
			
			//TBD: override grid_u and grid_v if hspeed and hdir is set
			if (windfield->grid_u[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_u[i], 0, windfield->lut_u + i);
			if (windfield->grid_v[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_v[i], 0, windfield->lut_v + i);
			if (windfield->grid_w[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_w[i], 0, windfield->lut_w + i);

			if (windfield->grid_hspeed_err[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_hspeed_err[i], 0, windfield->lut_hspeed_err + i);
			if (windfield->grid_hdir_err[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_hdir_err[i], 0, windfield->lut_hdir_err + i);
			if (windfield->grid_u_err[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_u_err[i], 0, windfield->lut_u_err + i);
			if (windfield->grid_v_err[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_v_err[i], 0, windfield->lut_v_err + i);
			if (windfield->grid_w_err[i] != NULL) util_field2lut(windfield->field[i], windfield->grid_w_err[i], 0, windfield->lut_w_err + i);

			if (windfield->fit_hspeed[i]) 
				util_initialize_field_err_cov_matrix(windfield->field[i], windfield->grid_hspeed_err[i], &(windfield->hspeed_ecm[i]));
			if (windfield->fit_hdir[i]) 
				util_initialize_field_err_cov_matrix(windfield->field[i], windfield->grid_hdir_err[i], &(windfield->hdir_ecm[i]));
			if (windfield->fit_u[i]) 
				util_initialize_field_err_cov_matrix(windfield->field[i], windfield->grid_u_err[i], &(windfield->u_ecm[i]));
			if (windfield->fit_v[i]) 
				util_initialize_field_err_cov_matrix(windfield->field[i], windfield->grid_v_err[i], &(windfield->v_ecm[i]));
			if (windfield->fit_w[i]) 
				util_initialize_field_err_cov_matrix(windfield->field[i], windfield->grid_w_err[i], &(windfield->w_ecm[i]));
		}
	}
	for ( i = 0; i < windfield->nvortices; i++ ) {
		if (windfield->vortex[i] != NULL) {
			//Make sure rotation direction is a unit vector.
			tmp = sqrt(
				pow(windfield->vortex[i]->rotation_direction[0], 2.) +
				pow(windfield->vortex[i]->rotation_direction[1], 2.) +
				pow(windfield->vortex[i]->rotation_direction[2], 2.));
			windfield->vortex[i]->rotation_direction[0] /= tmp;
			windfield->vortex[i]->rotation_direction[1] /= tmp;
			windfield->vortex[i]->rotation_direction[2] /= tmp;
		}
	}
	
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("turbulence\n"); fflush(stdout);
	#endif	
	for ( i = 0; i < windfield->nturbulences; i++ ) {
		if (windfield->turbulence[i] != NULL) {
			fields_prepare(windfield->turbulence[i]->field);
			turbulence_prepare_widget(windfield->turbulence[i]);
		}
	}
	
	
}


void util_prepare_post_windfield(t_zephyros_windfield **pdst, t_zephyros_windfield *src)
{
	t_zephyros_windfield *dst;
	int i;
	
	if (src == NULL) {
		dst = NULL;
	} else {
		dst = malloc(sizeof(t_zephyros_windfield));
		 
		memcpy(dst, src, sizeof(t_zephyros_windfield));
		
		for ( i = 0; i < 101; i++ ) {
			fields_copy(&(dst->field[i]), src->field[i]);

			dst->lut_hspeed[i] = NULL;
			dst->lut_hdir[i] = NULL;
			dst->lut_u[i] = NULL;
			dst->lut_v[i] = NULL;
			dst->lut_w[i] = NULL;

			dst->lut_hspeed_err[i] = NULL;
			dst->lut_hdir_err[i] = NULL;
			dst->lut_u_err[i] = NULL;
			dst->lut_v_err[i] = NULL;
			dst->lut_w_err[i] = NULL;
						
			if (src->field[i] != NULL) {
				util_copy_dbl_array(src->field[i]->n, &(dst->grid_hspeed[i]), src->grid_hspeed[i]);
				util_copy_dbl_array(src->field[i]->n, &(dst->grid_hdir[i]), src->grid_hdir[i]);
				util_copy_dbl_array(src->field[i]->n, &(dst->grid_u[i]), src->grid_u[i]);
				util_copy_dbl_array(src->field[i]->n, &(dst->grid_v[i]), src->grid_v[i]);
				util_copy_dbl_array(src->field[i]->n, &(dst->grid_w[i]), src->grid_w[i]);

				if (dst->grid_hspeed[i] != NULL) util_field2lut(dst->field[i], dst->grid_hspeed[i], 0, dst->lut_hspeed + i);
				if (dst->grid_hdir[i] != NULL) util_field2lut(dst->field[i], dst->grid_hdir[i], 0, dst->lut_hdir + i);
				if (dst->grid_u[i] != NULL) util_field2lut(dst->field[i], dst->grid_u[i], 0, dst->lut_u + i);
				if (dst->grid_v[i] != NULL) util_field2lut(dst->field[i], dst->grid_v[i], 0, dst->lut_v + i);
				if (dst->grid_w[i] != NULL) util_field2lut(dst->field[i], dst->grid_w[i], 0, dst->lut_w + i);

				util_copy_dbl_array(src->field[i]->n, &(dst->grid_hspeed_err[i]), src->grid_hspeed_err[i]);
				util_copy_dbl_array(src->field[i]->n, &(dst->grid_hdir_err[i]), src->grid_hdir_err[i]);
				util_copy_dbl_array(src->field[i]->n, &(dst->grid_u_err[i]), src->grid_u_err[i]);
				util_copy_dbl_array(src->field[i]->n, &(dst->grid_v_err[i]), src->grid_v_err[i]);
				util_copy_dbl_array(src->field[i]->n, &(dst->grid_w_err[i]), src->grid_w_err[i]);

				if (dst->grid_hspeed_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_hspeed[i], 0, dst->lut_hspeed_err + i);
				if (dst->grid_hdir_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_hdir[i], 0, dst->lut_hdir_err + i);
				if (dst->grid_u_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_u[i], 0, dst->lut_u_err + i);
				if (dst->grid_v_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_v[i], 0, dst->lut_v_err + i);
				if (dst->grid_w_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_w[i], 0, dst->lut_w_err + i);
			}

			dst->vortex[i] = NULL;
			dst->wave[i] = NULL;
			
			
			if (src->turbulence[i] != NULL) {
				if (src->turbulence[i]->type == 5) {
					turbulence_initialize_widget(&(dst->turbulence[i]));
					dst->turbulence[i]->type = src->turbulence[i]->type;
					fields_copy(&(dst->turbulence[i]->field), src->turbulence[i]->field);

					util_copy_dbl_array(src->turbulence[i]->field->n, &(dst->turbulence[i]->grid_edr), src->turbulence[i]->grid_edr);
					util_copy_dbl_array(src->turbulence[i]->field->n, &(dst->turbulence[i]->grid_edr13), src->turbulence[i]->grid_edr13);
					util_copy_dbl_array(src->turbulence[i]->field->n, &(dst->turbulence[i]->grid_edr13_err), src->turbulence[i]->grid_edr13_err);
					if (dst->turbulence[i]->grid_edr13 != NULL) util_field2lut(dst->turbulence[i]->field, dst->turbulence[i]->grid_edr13, 0, &(dst->turbulence[i]->lut_edr13));
					
					util_copy_dbl_array(src->turbulence[i]->field->n, &(dst->turbulence[i]->grid_kolmogorov_constant), src->turbulence[i]->grid_kolmogorov_constant);
					if (dst->turbulence[i]->grid_kolmogorov_constant != NULL) util_field2lut(dst->turbulence[i]->field, dst->turbulence[i]->grid_kolmogorov_constant, 0, &(dst->turbulence[i]->lut_kolmogorov_constant));
					dst->turbulence[i]->fit_edr13 = src->turbulence[i]->fit_edr13;
					dst->turbulence[i]->edr13_ecm.mat = NULL;
					dst->turbulence[i]->edr13_ecm.mat_N = NULL;
					dst->turbulence[i]->edr13_ecm.mat_S = NULL;
				}
			} else {
				dst->turbulence[i] = NULL;
			}

			//not used in post object
			dst->hdir_ecm[i].mat = NULL;
			dst->hdir_ecm[i].mat_S = NULL;
			dst->hdir_ecm[i].mat_N = NULL;
			dst->hspeed_ecm[i].mat = NULL;
			dst->hspeed_ecm[i].mat_S = NULL;
			dst->hspeed_ecm[i].mat_N = NULL;
			dst->u_ecm[i].mat = NULL;
			dst->u_ecm[i].mat_S = NULL;
			dst->u_ecm[i].mat_N = NULL;
			dst->v_ecm[i].mat = NULL;
			dst->v_ecm[i].mat_S = NULL;
			dst->v_ecm[i].mat_N = NULL;
			dst->w_ecm[i].mat = NULL;
			dst->w_ecm[i].mat_S = NULL;
			dst->w_ecm[i].mat_N = NULL;
		}
	}
	*pdst = dst;
}

void util_windfield_wave_fuvw(t_zephyros_wave *wave, double *xyzt, int i_uvw, double *output, int calcderivatives, double *derivatives)
{
	double phase;
	phase = wave->phi0_rad + (wave->k_minv[0] * xyzt[0]) + (wave->k_minv[1] * xyzt[1]) + (wave->k_minv[2] * xyzt[2]) + (wave->f_sinv * xyzt[3]);
	*output = wave->amplitude_msinv[i_uvw] * sin(phase);
	if (calcderivatives) {
			derivatives[0] = wave->amplitude_msinv[i_uvw] * wave->k_minv[0] * cos(phase);
			derivatives[1] = wave->amplitude_msinv[i_uvw] * wave->k_minv[1] * cos(phase);
			derivatives[2] = wave->amplitude_msinv[i_uvw] * wave->k_minv[2] * cos(phase);
			derivatives[3] = wave->amplitude_msinv[i_uvw] * wave->f_sinv * cos(phase);
	}
}

void util_windfield_vortex_fuvw(t_zephyros_vortex *vortex, double *xyzt, int i_uvw, double *output, int calcderivatives, double *derivatives)
{
	double rv[3], rabs, rc;
	double omega[3], omegaabs;
	double domega_dt[3], domegaabs_dt;
	double tstart;
	double drc_dt;
	double dummy;
	double lamb_oseen_alpha, lamb_oseen_rmax, lamb_oseen_t0;
				
	rv[0] = xyzt[0] - vortex->xyz0_m[0];
	rv[1] = xyzt[1] - vortex->xyz0_m[1];
	rv[2] = xyzt[2] - vortex->xyz0_m[2];
	
	//subtract rotation direction
	dummy = 	  rv[0] * vortex->rotation_direction[0]
				+ rv[1] * vortex->rotation_direction[1]
				+ rv[2] * vortex->rotation_direction[2];
	rv[0] -= dummy * vortex->rotation_direction[0];
	rv[1] -= dummy * vortex->rotation_direction[1];
	rv[2] -= dummy * vortex->rotation_direction[2];
	
	//apply scaling
	rv[0] *= vortex->xyz_scaling_m[0];
	rv[1] *= vortex->xyz_scaling_m[1];
	rv[2] *= vortex->xyz_scaling_m[2];
	 
	rabs = sqrt(pow(rv[0], 2.) + pow(rv[1], 2.) + pow(rv[2], 2.));
	
	if ((rabs == 0.) | (rabs > vortex->maxr_m)) {
		*output = 0;
		if (calcderivatives) {
			derivatives[0] = 0.;
			derivatives[1] = 0.;
			derivatives[2] = 0.;
			derivatives[3] = 0.;
		}
	} else {
		if (vortex->type == 1) {
			//1 = Rankine vortex
			if (rabs < vortex->r_m) {
				omegaabs = (vortex->vmax_msinv * (rabs / vortex->r_m)) / rabs;	
			} else {
				omegaabs = (vortex->vmax_msinv * (vortex->r_m / rabs)) / rabs;	
			}
			domegaabs_dt = 0.;	
		}
		
		if (vortex->type == 2) {
			//2 = Lamb-Oseen vortex
			lamb_oseen_alpha = 1.25643; //Devenport et al. / Wikipedia
			lamb_oseen_t0	 = -1. * pow( vortex->r_m, 2.) / (lamb_oseen_alpha * 4. * vortex->nu);
			
			if ((xyzt[3] - lamb_oseen_t0) > 0.) {
				rc = sqrt(4. * vortex->nu * (xyzt[3] - lamb_oseen_t0));
				lamb_oseen_rmax = sqrt(lamb_oseen_alpha) * rc;
				omegaabs = (vortex->vmax_msinv * (1. + (1. / (2. * lamb_oseen_alpha))) * (lamb_oseen_rmax / rabs ))  * (1. - exp(-1. * lamb_oseen_alpha * pow(rabs / lamb_oseen_rmax,2.))) * (1. / rabs);

				if (calcderivatives) {
					drc_dt 			= 4. * vortex->nu * (1./2.) / rc;
					domegaabs_dt 	= 
						vortex->vmax_msinv * (1. + (1. / (2. * lamb_oseen_alpha))) *
						((drc_dt / rabs)
						+ ((pow(rabs, -2.) + (2. * lamb_oseen_alpha / pow(lamb_oseen_rmax, 2.)))
							 * (exp(-1. * lamb_oseen_alpha * pow(rabs / lamb_oseen_rmax,2.)))));
				}
			} else {
				omegaabs = 0.;
				domegaabs_dt = 0.;
			}
		}

		omega[0] = omegaabs * vortex->rotation_direction[0];
		omega[1] = omegaabs * vortex->rotation_direction[1];
		omega[2] = omegaabs * vortex->rotation_direction[2];
		
		//outer product
		if (i_uvw == 0) {*output = (omega[1] * rv[2]) - (omega[2] * rv[1]);}
		if (i_uvw == 1) {*output = (omega[2] * rv[0]) - (omega[0] * rv[2]);}
		if (i_uvw == 2) {*output = (omega[0] * rv[1]) - (omega[1] * rv[0]);}

		if (calcderivatives) {
			domega_dt[0] = domegaabs_dt * vortex->rotation_direction[0];
			domega_dt[1] = domegaabs_dt * vortex->rotation_direction[1];
			domega_dt[2] = domegaabs_dt * vortex->rotation_direction[2];
			
			if (i_uvw == 0) {
				derivatives[0] = 0.;
				derivatives[1] = vortex->xyz_scaling_m[1] * (-1. * omega[2]);
				derivatives[2] = vortex->xyz_scaling_m[2] * omega[1] ;
				derivatives[3] = (domega_dt[1] * rv[2]) - (domega_dt[2] * rv[1]);
			}
			if (i_uvw == 1) {
				derivatives[0] = vortex->xyz_scaling_m[0] * omega[2];
				derivatives[1] = 0.;
				derivatives[2] = vortex->xyz_scaling_m[2] * (-1. * omega[0]);
				derivatives[3] = (domega_dt[2] * rv[0]) - (domega_dt[0] * rv[2]);
			}
			if (i_uvw == 2) {
				derivatives[0] = vortex->xyz_scaling_m[0] * (-1. * omega[1]);
				derivatives[1] = vortex->xyz_scaling_m[1] * omega[0];
				derivatives[2] = 0.;
				derivatives[3] = (domega_dt[0] * rv[1]) - (domega_dt[1] * rv[0]);
			}
		}
	}
}

void util_windfield_fuvw(t_zephyros_windfield *windfield, double* xyzt, int i_uvw, double *output, int calcderivatives, double *derivatives) 
{
	int bilint_special = 0;
	int i;
	double tmp, tmpderivatives[4];
	double *dummy;
	t_zephyros_interpolation_bilint_lut *my_lut_pointer;
	
	//set zero
	*output = 0.;
	
	if (calcderivatives) {
		derivatives[0] = 0.;
		derivatives[1] = 0.;
		derivatives[2] = 0.;
		derivatives[3] = 0.;
	}
	
	//walk through fields
	for ( i = 0; i < windfield->nfields; i++ ) {
		my_lut_pointer = NULL;
		if (i_uvw == 0) {my_lut_pointer = windfield->lut_u[i];}
		if (i_uvw == 1) {my_lut_pointer = windfield->lut_v[i];}
		if (i_uvw == 2) {my_lut_pointer = windfield->lut_w[i];}

		if (my_lut_pointer != NULL) {
			interpolation_bilint(my_lut_pointer,
					xyzt,
					&tmp,
					0,
					dummy);
			*output = *output + tmp;
			if (calcderivatives) {
				derivatives[0] = derivatives[0] + tmpderivatives[0];
				derivatives[1] = derivatives[1] + tmpderivatives[1];
				derivatives[2] = derivatives[2] + tmpderivatives[2];
				derivatives[3] = derivatives[3] + tmpderivatives[3];
			}
		}
	}
	
	//walk through waves
	for ( i = 0; i < windfield->nwaves; i++ ) {
		if (windfield->wave[i] != NULL) {
			util_windfield_wave_fuvw(windfield->wave[i], xyzt, i_uvw, &tmp, calcderivatives, tmpderivatives);
			*output = *output + tmp;
			if (calcderivatives) {
				derivatives[0] = derivatives[0] + tmpderivatives[0];
				derivatives[1] = derivatives[1] + tmpderivatives[1];
				derivatives[2] = derivatives[2] + tmpderivatives[2];
				derivatives[3] = derivatives[3] + tmpderivatives[3];
			}	
		}
	}
	
	//walk through vortices
	for ( i = 0; i < windfield->nvortices; i++ ) {
		if (windfield->vortex[i] != NULL) {
			util_windfield_vortex_fuvw(windfield->vortex[i], xyzt, i_uvw, &tmp, calcderivatives, tmpderivatives);
			*output = *output + tmp;
			if (calcderivatives) {
				derivatives[0] = derivatives[0] + tmpderivatives[0];
				derivatives[1] = derivatives[1] + tmpderivatives[1];
				derivatives[2] = derivatives[2] + tmpderivatives[2];
				derivatives[3] = derivatives[3] + tmpderivatives[3];
			}	
		}
	}
	
	//walk through turbulences
	for ( i = 0; i < windfield->nturbulences; i++ ) {
		if (windfield->turbulence[i] != NULL) {
			turbulence_uvw((windfield->turbulence[i]), xyzt, i_uvw, &tmp, calcderivatives, tmpderivatives);
			*output = *output + tmp;
			if (calcderivatives) {
				derivatives[0] = derivatives[0] + tmpderivatives[0];
				derivatives[1] = derivatives[1] + tmpderivatives[1];
				derivatives[2] = derivatives[2] + tmpderivatives[2];
				derivatives[3] = derivatives[3] + tmpderivatives[3];
			}
		}
	}
}

void util_free_windfield(t_zephyros_windfield **pwindfield)
{
	t_zephyros_windfield *windfield = *pwindfield;
	int i;
	
	if (windfield != NULL) {
		for ( i = 0; i <= 100; i++ ) {
			fields_free(&windfield->field[i]);
			
			if (windfield->grid_u[i] != NULL) free(windfield->grid_u[i]);
			if (windfield->grid_v[i] != NULL) free(windfield->grid_v[i]);
			if (windfield->grid_w[i] != NULL) free(windfield->grid_w[i]);
			if (windfield->wave[i] != NULL) free(windfield->wave[i]);
			if (windfield->vortex[i] != NULL) free(windfield->vortex[i]); 
			if (windfield->turbulence[i] != NULL) {
				turbulence_free_widget(windfield->turbulence + i);
			}

			if (windfield->grid_u_err[i] != NULL) free(windfield->grid_u_err[i]);
			if (windfield->grid_v_err[i] != NULL) free(windfield->grid_v_err[i]);
			if (windfield->grid_w_err[i] != NULL) free(windfield->grid_w_err[i]);

			if (windfield->hdir_ecm[i].mat != NULL) cs_spfree(windfield->hdir_ecm[i].mat);
			if (windfield->hdir_ecm[i].mat_S != NULL) cs_sfree(windfield->hdir_ecm[i].mat_S);
			if (windfield->hdir_ecm[i].mat_N != NULL) cs_nfree(windfield->hdir_ecm[i].mat_N);

			if (windfield->hspeed_ecm[i].mat != NULL) cs_spfree(windfield->hspeed_ecm[i].mat);
			if (windfield->hspeed_ecm[i].mat_S != NULL) cs_sfree(windfield->hspeed_ecm[i].mat_S);
			if (windfield->hspeed_ecm[i].mat_N != NULL) cs_nfree(windfield->hspeed_ecm[i].mat_N);

			if (windfield->u_ecm[i].mat != NULL) cs_spfree(windfield->u_ecm[i].mat);
			if (windfield->u_ecm[i].mat_S != NULL) cs_sfree(windfield->u_ecm[i].mat_S);
			if (windfield->u_ecm[i].mat_N != NULL) cs_nfree(windfield->u_ecm[i].mat_N);

			if (windfield->v_ecm[i].mat != NULL) cs_spfree(windfield->v_ecm[i].mat);
			if (windfield->v_ecm[i].mat_S != NULL) cs_sfree(windfield->v_ecm[i].mat_S);
			if (windfield->v_ecm[i].mat_N != NULL) cs_nfree(windfield->v_ecm[i].mat_N);

			if (windfield->w_ecm[i].mat != NULL) cs_spfree(windfield->w_ecm[i].mat);
			if (windfield->w_ecm[i].mat_S != NULL) cs_sfree(windfield->w_ecm[i].mat_S);
			if (windfield->w_ecm[i].mat_N != NULL) cs_nfree(windfield->w_ecm[i].mat_N);
		}
		free(windfield);
		windfield = NULL;
	}
}



void util_initialize_scattererfield(t_zephyros_scattererfield **pscattererfield)
{
	t_zephyros_scattererfield *scattererfield = malloc(sizeof(t_zephyros_scattererfield));
	int i;
	
	for ( i = 0; i <= 100; i++ ) {
		scattererfield->psd[i] 	= NULL;
	}
	scattererfield->npsd = 0;

	*pscattererfield = scattererfield;
}

void util_prepare_scattererfield(t_zephyros_scattererfield *scattererfield)
{
	int i;
	
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_prepare_scattererfield\n"); fflush(stdout);
	#endif	
	for ( i = 0; i < scattererfield->npsd; i++ ) {
		if (scattererfield->psd[i] != NULL) {
			util_prepare_psd(scattererfield->psd[i], scattererfield);
		}
	}
}

void util_prepare_post_scattererfield(t_zephyros_scattererfield **pdst, t_zephyros_scattererfield *src)
{
	t_zephyros_scattererfield *dst;
	int i;
	
	if (src == NULL) {
		dst = NULL;
	} else {
		dst = malloc(sizeof(t_zephyros_scattererfield));
		 
		dst->type = 2;
		dst->npsd = src->npsd;
			
		for ( i = 0; i < 101; i++ ) {
			util_prepare_post_psd(&(dst->psd[i]), src->psd[i]);

		}
		
	}
	*pdst = dst;
}

void util_free_scattererfield(t_zephyros_scattererfield **pscattererfield)
{
	t_zephyros_scattererfield *scattererfield = *pscattererfield;
	int i, j;
	
	if (scattererfield != NULL) {
		for ( i = 0; i <= 100; i++ ) {
			if (scattererfield->psd[i] != NULL) {
				fields_free(&scattererfield->psd[i]->field); 
								
				if (scattererfield->psd[i]->grid_lwc_gm3 != NULL) free(scattererfield->psd[i]->grid_lwc_gm3);
				if (scattererfield->psd[i]->grid_dBlwc_err_gm3 != NULL) free(scattererfield->psd[i]->grid_dBlwc_err_gm3);
				if (scattererfield->psd[i]->grid_gammadistribution_mu != NULL) free(scattererfield->psd[i]->grid_gammadistribution_mu);
				if (scattererfield->psd[i]->grid_gammadistribution_D0_mm != NULL) free(scattererfield->psd[i]->grid_gammadistribution_D0_mm);
				
				if (scattererfield->psd[i]->grid_number_density_m3 != NULL) {
					for ( j = 0; j < scattererfield->psd[i]->n_diameters; j++ ) {
						if (scattererfield->psd[i]->grid_number_density_m3[j] != NULL) free(scattererfield->psd[i]->grid_number_density_m3[j]);
					}
					free(scattererfield->psd[i]->grid_number_density_m3);
				}
				if (scattererfield->psd[i]->lut_ln_number_density_m3 != NULL) {
					for ( j = 0; j < scattererfield->psd[i]->n_diameters; j++ ) {
						interpolation_free_lut(&scattererfield->psd[i]->lut_ln_number_density_m3[j]);
					}
					free(scattererfield->psd[i]->lut_ln_number_density_m3);
				}
				
				if (scattererfield->psd[i]->grid_dBnumber_density_err_m3 != NULL) {
					for ( j = 0; j < scattererfield->psd[i]->n_diameters; j++ ) {
						if (scattererfield->psd[i]->grid_dBnumber_density_err_m3[j] != NULL) free(scattererfield->psd[i]->grid_dBnumber_density_err_m3[j]);
					}
					free(scattererfield->psd[i]->grid_dBnumber_density_err_m3);
				}
				if (scattererfield->psd[i]->lut_dBnumber_density_err_m3 != NULL) {
					for ( j = 0; j < scattererfield->psd[i]->n_diameters; j++ ) {
						interpolation_free_lut(&scattererfield->psd[i]->lut_dBnumber_density_err_m3[j]);
					}
					free(scattererfield->psd[i]->lut_dBnumber_density_err_m3);
				}
									
				if (scattererfield->psd[i]->discrete_D_equiv_mm != NULL) {
					free(scattererfield->psd[i]->discrete_D_equiv_mm);
				}
				
				if (scattererfield->psd[i]->grid_gammadistribution_N0 != NULL) free(scattererfield->psd[i]->grid_gammadistribution_N0);
				if (scattererfield->psd[i]->grid_gammadistribution_N0_err != NULL) free(scattererfield->psd[i]->grid_gammadistribution_N0_err);

				if (scattererfield->psd[i]->dBlwc_ecm.mat != NULL) cs_spfree(scattererfield->psd[i]->dBlwc_ecm.mat);
				if (scattererfield->psd[i]->dBlwc_ecm.mat_S != NULL) cs_sfree(scattererfield->psd[i]->dBlwc_ecm.mat_S);
				if (scattererfield->psd[i]->dBlwc_ecm.mat_N != NULL) cs_nfree(scattererfield->psd[i]->dBlwc_ecm.mat_N);

				free(scattererfield->psd[i]);
			}			
		}
		free(scattererfield);
		scattererfield = NULL;
	}
}






/*

//info.xyzt[4] is set.
void scatterers_calc_scattererinfo(
	t_zephyros_scattererfield *scattererfield, 
	t_zephyros_config_simulation_atmosphere	*atmosphere,
	t_zephyros_windfield *windfield,
	t_zephyros_scattererinfo *info)
{
	double tmp;
	double tmpderivatives[4];
	int calcderivatives;
	int bilint_special;
	
	//dBZ (equivalent reflectivity)
	//etc. etc.

	info.dBZ = 0.;
	info.dBZ_der[0] = 0.;
	info.dBZ_der[1] = 0.;
	info.dBZ_der[2] = 0.;
	info.dBZ_der[3] = 0.;
	
	info.attenuation = 0.;
	info.attenuation_der[0] = 0.;
	info.attenuation_der[1] = 0.;
	info.attenuation_der[2] = 0.;
	info.attenuation_der[3] = 0.;
	
	//info.terminal_fall_speed
	
	//walk through particles
	for ( i = 0; i < scattererfield->nparticles; i++ ) {
		//0 = perfect tracer
		if (scattererfield->particle_type[i] == 0) {
			calcderivatives = 0; //todo->calc_point_dBZ_derivatives;
			bilint_special = 0;
			
			bilint(	&scattererfield->field[i]->n_dim,
					&scattererfield->field[i]->n,
					scattererfield->field[i]->int_ax_count,
					scattererfield->field[i]->int_ax_values,
					scattererfield->dBZ[i],
					info.xyzt,
					&tmp,		//memory arithmetic
					&calcderivatives,
					&bilint_special,
					tmpderivatives);
			
			info.dBZ += tmp;
			info.dBZ_der[0] += tmpderivatives[0];
			info.dBZ_der[1] += tmpderivatives[1];
			info.dBZ_der[2] += tmpderivatives[2];
			info.dBZ_der[3] += tmpderivatives[3];

			calcderivatives = 0; //todo->calc_point_dBZ_derivatives;
			bilint_special = 0;


			
			bilint(	&scattererfield->field[i]->n_dim,
					&scattererfield->field[i]->n,
					scattererfield->field[i]->int_ax_count,
					scattererfield->field[i]->int_ax_values,
					scattererfield->attenuation[i],
					info.xyzt,
					&tmp,		//memory arithmetic
					&calcderivatives,
					&bilint_special,
					tmpderivatives);
			
			info.attenuation += tmp;
			info.attenuation_der[0] += tmpderivatives[0];
			info.attenuation_der[1] += tmpderivatives[1];
			info.attenuation_der[2] += tmpderivatives[2];
			info.attenuation_der[3] += tmpderivatives[3];
		}

		//1, or 2 ... ???
		

	}


}
*/



void util_initialize_psd(t_zephyros_psd **pthepsd) 
{
	t_zephyros_psd	*thepsd;
	
	thepsd = malloc(sizeof(t_zephyros_psd));
	
	thepsd->distribution_type = 0;
	thepsd->particle_type = -1;
	fields_initialize(&(thepsd->field));
	strcpy(thepsd->field->name, "scatterer");		
	thepsd->grid_lwc_gm3 = NULL;
	thepsd->grid_dBlwc_err_gm3 = NULL;
	thepsd->discrete_D_equiv_mm = NULL;
	thepsd->grid_number_density_m3 = NULL;
	thepsd->grid_dBnumber_density_err_m3 = NULL;
	thepsd->lut_ln_number_density_m3 = NULL;
	thepsd->lut_dBnumber_density_err_m3 = NULL;
	thepsd->grid_gammadistribution_N0 = NULL;
	thepsd->grid_gammadistribution_N0_err = NULL;
	thepsd->grid_gammadistribution_mu = NULL;
	thepsd->grid_gammadistribution_D0_mm = NULL;
	thepsd->dBlwc_ecm.mat  = NULL;
	thepsd->dBlwc_ecm.mat_N  = NULL;
	thepsd->dBlwc_ecm.mat_S  = NULL;
	thepsd->gammadistribution_mu_err = 5.;
	thepsd->gammadistribution_D0_err_mm = 1.;	
	thepsd->n_diameters = 0;
	
	thepsd->fit_dBlwc = 0;
	thepsd->fit_dBN = 0;

	*pthepsd = thepsd;
}


void util_prepare_psd(t_zephyros_psd *thepsd, t_zephyros_scattererfield *scattererfield) 
{
	int i,j;
	double mincdfP, maxcdfP, difcdfP;
	double *Dl, *Du, *Dcenter;
	double thislwc_gm3;
	
	double dn_dmu;
	double dn_dD0;
	double dn_dN0;
	double mydelta;
	double l_lim, u_lim;
	
	//The following steps are taken
	//1. Prepare field
	//2. Assert grid_number_density is set
	//3. Update the discrete probability density function
	//4. Update grid_number_density with grid_mu grid_D0
	//5. Assert grid_number_density_err_m3 is set
	//6. Update grid_number_density_err_m3
	//7. Update number density with lwc, or update lwc
	//8. Assert grid_dBlwc_err_gm3
	//9. make luts
	//10. make error covariance matrices
	
	//1. Pepare field
	fields_prepare(thepsd->field);

	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("Assert grid_number_density is set\n"); fflush(stdout);
	#endif		

	//2. Assert that grid_number_density is set
	if (thepsd->grid_number_density_m3 == NULL) {
		thepsd->grid_number_density_m3 = malloc(thepsd->n_diameters * sizeof(double*));
		for (i = 0; i < thepsd->n_diameters; i++ )
			if (thepsd->grid_number_density_m3[i] = NULL);
	}
	for (i = 0; i < thepsd->n_diameters; i++ ) {
		if (thepsd->grid_number_density_m3[i] == NULL) {
			thepsd->grid_number_density_m3[i] = malloc(thepsd->field->n * sizeof(double));
			for (j = 0; j < thepsd->field->n; j++ )
				thepsd->grid_number_density_m3[i][j] = 1.;
		}
	}

	//3. Update the discrete probability density function
	//3a. discrete probability density function
	//if (thepsd->distribution_type == 0) {
	//
	//}
	
	//3b. Translate gamma pdf to discrete probability density function
	if (thepsd->distribution_type == 1) {
		#ifdef _ZEPHYROS_UTIL_DEBUG
			printf("translate gamma pdf to discrete probability density function\n"); fflush(stdout);
		#endif	

		Dl = malloc(thepsd->n_diameters * sizeof(double));
		Dcenter = malloc(thepsd->n_diameters * sizeof(double));
		Du = malloc(thepsd->n_diameters * sizeof(double));
		
		//assert thepsd->grid_gammadistribution_N0 is set
		if (thepsd->grid_gammadistribution_N0 == NULL) {
			thepsd->grid_gammadistribution_N0 = malloc(thepsd->field->n * sizeof(double));
			for (j = 0; j < thepsd->field->n; j++ ) {
				thepsd->grid_gammadistribution_N0[j] = 1.;
			}
		}
		
		//assert thepsd->discrete_D_equiv_mm is set
		if (thepsd->discrete_D_equiv_mm != NULL) free(thepsd->discrete_D_equiv_mm);
		thepsd->discrete_D_equiv_mm = malloc(thepsd->n_diameters * sizeof(double));

		mincdfP = 1.e-50;
		maxcdfP = 1. - 1.e-50;
		difcdfP = (maxcdfP - mincdfP) / thepsd->n_diameters;

		for (i = 0; i < thepsd->n_diameters; i++ ) {
			/*
			Dl[i] = util_inverse_gammaDalpha_cdf(thepsd->gammadistribution_mu, thepsd->gammadistribution_D0_mm, thepsd->gammadistribution_dmin_mm , thepsd->gammadistribution_dmax_mm, mincdfP + (i * difcdfP));
			Du[i] = util_inverse_gammaDalpha_cdf(thepsd->gammadistribution_mu, thepsd->gammadistribution_D0_mm, thepsd->gammadistribution_dmin_mm , thepsd->gammadistribution_dmax_mm, mincdfP + ((i + 1.) * difcdfP));
			Dcenter[i] = util_inverse_gammaDalpha_cdf(thepsd->gammadistribution_mu, thepsd->gammadistribution_D0_mm, thepsd->gammadistribution_dmin_mm , thepsd->gammadistribution_dmax_mm, mincdfP + ((i + .5) * difcdfP));
			*/
			
			//for the moment used a simplified model
			Dl[i] = thepsd->gammadistribution_dmin_mm  +
					((thepsd->gammadistribution_dmax_mm - thepsd->gammadistribution_dmin_mm) *
					(mincdfP + (i * difcdfP)));
			Du[i] = thepsd->gammadistribution_dmin_mm  +
					((thepsd->gammadistribution_dmax_mm - thepsd->gammadistribution_dmin_mm) *
					(mincdfP + ((i + 1) * difcdfP)));

			Dcenter[i] = (Du[i] - Dl[i]) / 2.;
			//set diameter size
			thepsd->discrete_D_equiv_mm[i] = Dcenter[i];
			
			//calculate the number density
			for (j = 0; j < thepsd->field->n; j++ ) {
				thepsd->grid_number_density_m3[i][j] = util_gamma_integral(thepsd->grid_gammadistribution_N0[j], thepsd->gammadistribution_mu, thepsd->gammadistribution_D0_mm, Dl[i], Du[i]);
			}						
		}
	}

	
	//4. Update grid_number_density_m3 with grid_gammadistribution_D0_mm grid_gammadistribution_mu
	//NOTE: Dl, Dcenter, Du are reused here.
	if (thepsd->distribution_type == 1) {
		if ((thepsd->grid_gammadistribution_D0_mm != NULL) & (thepsd->grid_gammadistribution_mu != NULL)) {
			for (i = 0; i < thepsd->n_diameters; i++ ) {
				//calculate the number density
				for (j = 0; j < thepsd->field->n; j++ ) {
					thepsd->grid_number_density_m3[i][j] = util_gamma_integral(thepsd->grid_gammadistribution_N0[j], thepsd->grid_gammadistribution_mu[j], thepsd->grid_gammadistribution_D0_mm[j], Dl[i], Du[i]);
				}						
			}
		}
	}

	//5. Assert grid_dBnumber_density_err_m3 is set
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("4. Assert grid_dBnumber_density_err_m3 is set\n"); fflush(stdout);
	#endif	
	
	if (thepsd->grid_dBnumber_density_err_m3 == NULL) {
		thepsd->grid_dBnumber_density_err_m3 = malloc(thepsd->n_diameters * sizeof(double*));
		for (i = 0; i < thepsd->n_diameters; i++ )
			if (thepsd->grid_dBnumber_density_err_m3[i] = NULL);
	}
	
	for (i = 0; i < thepsd->n_diameters; i++ ) {
		if (thepsd->grid_dBnumber_density_err_m3[i] == NULL)
			thepsd->grid_dBnumber_density_err_m3[i] = malloc(thepsd->field->n * sizeof(double));
		for (j = 0; j < thepsd->field->n; j++ )
			thepsd->grid_dBnumber_density_err_m3[i][j] = 40.;
	}

	//6. Update grid_dBnumber_density_err_m3
	if (thepsd->distribution_type == 1) {
		#ifdef _ZEPHYROS_UTIL_DEBUG
			printf("update grid number density error\n"); fflush(stdout);
		#endif	
		for (i = 0; i < thepsd->n_diameters; i++ ) {
			for (j = 0; j < thepsd->field->n; j++ ) {
				mydelta = fmax(1.e-5, thepsd->gammadistribution_mu / 100.);
				l_lim = thepsd->gammadistribution_mu - mydelta/2.; u_lim = l_lim + mydelta;
				if (l_lim < 0.) {u_lim -= l_lim; l_lim = u_lim/1.e4;}					
				dn_dmu = 	(util_gamma_integral(thepsd->grid_gammadistribution_N0[j], l_lim, thepsd->gammadistribution_D0_mm, Dl[i], Du[i]) -
							util_gamma_integral(thepsd->grid_gammadistribution_N0[j], u_lim, thepsd->gammadistribution_D0_mm, Dl[i], Du[i]))
							/ mydelta;
							
				mydelta = fmax(1.e-5, thepsd->gammadistribution_D0_mm / 100.);
				l_lim = thepsd->gammadistribution_D0_mm - mydelta/2.; u_lim = l_lim + mydelta;
				if (l_lim < 0.) {u_lim -= l_lim; l_lim = u_lim/1.e4;}
				dn_dD0 = 	(util_gamma_integral(thepsd->grid_gammadistribution_N0[j], thepsd->gammadistribution_mu, l_lim, Dl[i], Du[i]) -
							util_gamma_integral(thepsd->grid_gammadistribution_N0[j], thepsd->gammadistribution_mu, u_lim, Dl[i], Du[i]))
							/ mydelta;
							
				mydelta = fmax(1.e-10, thepsd->grid_gammadistribution_N0[j] / 100.);
				l_lim = thepsd->grid_gammadistribution_N0[j] - mydelta/2.; u_lim = l_lim + mydelta;
				if (l_lim < 0.) {u_lim -= l_lim; l_lim = u_lim/1.e4;}
				dn_dN0 =	(util_gamma_integral(l_lim, thepsd->gammadistribution_mu, thepsd->gammadistribution_D0_mm, Dl[i], Du[i]) -
							util_gamma_integral(u_lim, thepsd->gammadistribution_mu, thepsd->gammadistribution_D0_mm, Dl[i], Du[i]))
							/ mydelta;
							
				thepsd->grid_dBnumber_density_err_m3[i][j] = 
								pow(dn_dmu * thepsd->gammadistribution_mu_err ,2.) +
								pow(dn_dD0 * thepsd->gammadistribution_D0_err_mm, 2.);
								
				if (thepsd->grid_gammadistribution_N0_err != NULL) {
					thepsd->grid_dBnumber_density_err_m3[i][j] += 
						pow(dn_dN0 * thepsd->grid_gammadistribution_N0_err[j], 2.);
				}

				thepsd->grid_dBnumber_density_err_m3[i][j] = func_dB(sqrt(thepsd->grid_dBnumber_density_err_m3[i][j]));
			}
		}
	}
		
	//7. Update number density with lwc, or update lwc
	if (thepsd->grid_lwc_gm3 != NULL) {
		#ifdef _ZEPHYROS_UTIL_DEBUG
			printf("update number density with lwc\n"); fflush(stdout);
		#endif
		for (j = 0; j < thepsd->field->n; j++ ) {
			thislwc_gm3 = 0.;
			for (i = 0; i < thepsd->n_diameters; i++ ) {
				thislwc_gm3 +=
					thepsd->grid_number_density_m3[i][j] *	//# m-3
					1.e3 *			//g kg-1
					998.2071 * 		//[water density in kg m-3]
					(4./3.) * M_PI * pow((thepsd->discrete_D_equiv_mm[i]/2.) * 1.e-3,3.); //m3 of single droplet
			}
			
			//update N0
			if (thepsd->distribution_type == 1) {
				thepsd->grid_gammadistribution_N0[j] *= thepsd->grid_lwc_gm3[j] / thislwc_gm3;
			}
			
			//update number densities
			for (i = 0; i < thepsd->n_diameters; i++ ) {
				thepsd->grid_number_density_m3[i][j] *= thepsd->grid_lwc_gm3[j] / thislwc_gm3;
				thepsd->grid_dBnumber_density_err_m3[i][j] += func_dB(thepsd->grid_lwc_gm3[j] / thislwc_gm3);
			}
			
		}
	//calculate lwc
	} else {
		thepsd->grid_lwc_gm3 = malloc(thepsd->field->n * sizeof(double));
		for (j = 0; j < thepsd->field->n; j++ ) {
			thislwc_gm3 = 0.;
			for (i = 0; i < thepsd->n_diameters; i++ ) {
				thislwc_gm3 +=
					thepsd->grid_number_density_m3[i][j] *	//# m-3
					1.e3 *			//g kg-1
					998.2071 * 		//[water density in kg m-3]
					(4./3.) * M_PI * pow((thepsd->discrete_D_equiv_mm[i]/2.) * 1.e-3,3.); //m3 of single droplet
			}
			thepsd->grid_lwc_gm3[j] = thislwc_gm3;
		}
	}

#ifdef _ZEPHYROS_UTIL_DEBUG
	printf("update lwc error\n"); fflush(stdout);
#endif	

	//8. Assert grid_dBlwc_err_gm3
	if (thepsd->fit_dBlwc) {
		thepsd->grid_dBlwc_err_gm3 = malloc(thepsd->field->n * sizeof(double));
		for (j = 0; j < thepsd->field->n; j++ ) {	
			//1. calculate lwc_err_gm3, with error propagation
			//2. convert to error in dB
			
			//1.
			thepsd->grid_dBlwc_err_gm3[j] = 0.;
			for (i = 0; i < thepsd->n_diameters; i++ ) {
				thepsd->grid_dBlwc_err_gm3[j] +=
					pow(
					func_dB_inv(thepsd->grid_dBnumber_density_err_m3[i][j]) *	//# m-3
					1.e3 *			//g kg-1
					998.2071 * 		//[water density in kg m-3]
					(4./3.) * M_PI * pow((thepsd->discrete_D_equiv_mm[i]/2.) * 1.e-3,3.)
					, 2.);
			}
			thepsd->grid_dBlwc_err_gm3[j] = sqrt(thepsd->grid_dBlwc_err_gm3[j]);
			//2.
			thepsd->grid_dBlwc_err_gm3[j] = (10. / (thepsd->grid_lwc_gm3[j] * log(10.)) ) * thepsd->grid_dBlwc_err_gm3[j];
		}
	}

	//9. make luts
#ifdef _ZEPHYROS_UTIL_DEBUG
	printf("make luts\n"); fflush(stdout);
#endif 	
	thepsd->lut_ln_number_density_m3 = 
		malloc(thepsd->n_diameters * sizeof(t_zephyros_interpolation_bilint_lut*));
	thepsd->lut_dBnumber_density_err_m3 = 
		malloc(thepsd->n_diameters * sizeof(t_zephyros_interpolation_bilint_lut*));
	for (i = 0; i < thepsd->n_diameters; i++ ) {
		util_field2lut(thepsd->field, thepsd->grid_number_density_m3[i], 3, thepsd->lut_ln_number_density_m3 + i);
		util_field2lut(thepsd->field, thepsd->grid_dBnumber_density_err_m3[i], 0, thepsd->lut_dBnumber_density_err_m3 + i);
	}
	
	//10. make error covariance matrices
	if (thepsd->fit_dBlwc) {
		#ifdef _ZEPHYROS_UTIL_DEBUG
			printf("make error covariance matrices\n"); fflush(stdout);
		#endif	
		util_initialize_field_err_cov_matrix(thepsd->field, thepsd->grid_dBlwc_err_gm3, &thepsd->dBlwc_ecm);
	}
		
	//clean up
	if (thepsd->distribution_type == 1) {
		#ifdef _ZEPHYROS_UTIL_DEBUG
			printf("clean up\n"); fflush(stdout);
		#endif	
		free(Dl);
		free(Dcenter);
		free(Du);
	}
}


void util_prepare_post_psd(t_zephyros_psd **pdst, t_zephyros_psd *src)
{
	t_zephyros_psd *dst;
	int i, j;
	
	if (src == NULL) {
		dst = NULL;
	} else {
		dst = malloc(sizeof(t_zephyros_psd));
		//copy everything
		memcpy(dst, src, sizeof(t_zephyros_psd));

		dst->distribution_type = 0;
			
		//copy field
		fields_copy(&dst->field, src->field);
	   
		dst->grid_lwc_gm3 = malloc(src->field->n * sizeof(double));
		for (j = 0; j < src->field->n; j++ ) 
			dst->grid_lwc_gm3[j] = src->grid_lwc_gm3[j];
		//errors not required in post
		dst->grid_dBlwc_err_gm3 = NULL;
		
		dst->discrete_D_equiv_mm = malloc(src->n_diameters * sizeof(double));
		for (i = 0; i < src->n_diameters; i++ ) 
			dst->discrete_D_equiv_mm[i] = src->discrete_D_equiv_mm[i];
		
		
		dst->grid_number_density_m3 = malloc(src->n_diameters * sizeof(double*));
		dst->grid_dBnumber_density_err_m3 = malloc(src->n_diameters * sizeof(double*));
		for (i = 0; i < src->n_diameters; i++ ) {
			dst->grid_number_density_m3[i] = malloc(src->field->n * sizeof(double));
			dst->grid_dBnumber_density_err_m3[i] = malloc(src->field->n * sizeof(double));
			for (j = 0; j < src->field->n; j++ ) {		
				dst->grid_number_density_m3[i][j] = src->grid_number_density_m3[i][j];
				dst->grid_dBnumber_density_err_m3[i][j] = src->grid_dBnumber_density_err_m3[i][j];
		}}
							
		//make luts
		dst->lut_ln_number_density_m3 = 
			malloc(src->n_diameters * sizeof(t_zephyros_interpolation_bilint_lut*));
		dst->lut_dBnumber_density_err_m3 = 
			malloc(src->n_diameters * sizeof(t_zephyros_interpolation_bilint_lut*));
		for (i = 0; i < src->n_diameters; i++ ) {
			util_field2lut(dst->field, dst->grid_number_density_m3[i], 3, dst->lut_ln_number_density_m3 + i);
			util_field2lut(dst->field, dst->grid_dBnumber_density_err_m3[i], 0, dst->lut_dBnumber_density_err_m3 + i);
		}

		dst->grid_gammadistribution_N0 = NULL;
		dst->grid_gammadistribution_N0_err = NULL;
		dst->grid_gammadistribution_mu = NULL;
		dst->grid_gammadistribution_D0_mm = NULL;
		dst->gammadistribution_mu = 0;
		dst->gammadistribution_mu_err = 0;
		dst->gammadistribution_D0_mm = 0;
		dst->gammadistribution_D0_err_mm = 0;
		dst->gammadistribution_dmin_mm = 0;
		dst->gammadistribution_dmax_mm = 0;
		dst->n_diameters = src->n_diameters;
		dst->fit_dBlwc = src->fit_dBlwc;
		dst->dBlwc_ecm.cl_x_m = 0;
		dst->dBlwc_ecm.cl_y_m = 0;
		dst->dBlwc_ecm.cl_z_m = 0;
		dst->dBlwc_ecm.cl_t_s = 0;
		dst->dBlwc_ecm.c_threshold = 0;
		
		//error covariance matrices
		dst->dBlwc_ecm.mat = NULL;
		dst->dBlwc_ecm.mat_S = NULL;
		dst->dBlwc_ecm.mat_N = NULL;
	}
	
	*pdst = dst;
}


//return D
double util_inverse_gammaDalpha_cdf(double mu, double D0, double Dmin, double Dmax, double cdfP)
{
	double pars[6];
	double err = 1.e-5;
	double maxiter = 25;
	
	double (*f)(double, void*);
	double (*df)(double, void*);
	
	double startD = D0;
	
	pars[0] = mu;
	pars[1] = D0;
	pars[2] = Dmin;
	pars[3] = Dmax;
	pars[4] = 1.;  //D^alpha
	pars[5] = cdfP;
		
	f = &util_inverse_gammaDalpha_cdf_costf;
	df = &util_inverse_gammaDalpha_cdf_costf_der;
		
	//solve with Newton Rapson method
	return fabs(func_newton_rapson_solver(f, df, startD, err, maxiter, (void*) pars));
}

double util_inverse_gammaDalpha_cdf_costf(double D, void *vp)
{
	double cdfP = ((double*)vp)[5];
	return pow(util_gammaDalpha_cdf(fabs(D), vp) - cdfP, 2.);
}

double util_inverse_gammaDalpha_cdf_costf_der(double D, void *vp) 
{
	double deltaD = 1.e-4;
	return (util_inverse_gammaDalpha_cdf_costf(D + deltaD / 2., vp) - util_inverse_gammaDalpha_cdf_costf(D - deltaD / 2., vp)) / deltaD;	
}

double util_gammaDalpha_cdf(double D, void *vp) 
{
	double mu 	= ((double*)vp)[0];
	double D0	= ((double*)vp)[1];
	double Dmin = ((double*)vp)[2];
	double Dmax = ((double*)vp)[3];
	double alpha = ((double*)vp)[4];
	
	return (	(specialfunctions_incompletegamma(mu + 1. + alpha, ((3.67 + mu) / D0) * D)
				- specialfunctions_incompletegamma(mu + 1. + alpha, ((3.67 + mu) / D0) * Dmin)) / 
				(specialfunctions_incompletegamma(mu + 1. + alpha, ((3.67 + mu) / D0) * Dmax)
				- specialfunctions_incompletegamma(mu + 1. + alpha, ((3.67 + mu) / D0) * Dmin) )
			);
}

double util_gamma_integral(double N0, double mu, double D0, double D1, double D2)
{
	double tmp;
	tmp = N0 * (	(specialfunctions_incompletegamma(mu + 1., ((3.67 + mu) / D0) * D2)
				- specialfunctions_incompletegamma(mu + 1., ((3.67 + mu) / D0) * D1)) / 
				pow((3.67 + mu) / D0 , mu+1.));

	if (tmp == 0.) {
		return N0 * 1.e-50;
	} else {
		return tmp;
	}
}








void util_initialize_field_err_cov_matrix(t_zephyros_field *field, double *err, t_zephyros_field_errcovmat *ecm)
{
    cs *errcovmattriplet; 
	double myval;
	int i, j;
	
	ecm->n = field->n;
	
    //Allocate empty matrices    
    errcovmattriplet 		= cs_spalloc (ecm->n, ecm->n, 1, 1, 1) ;

	//Walk through matrix, and save the relevant entries.
	for ( i = 0; i < ecm->n; i++ )
	{
		for ( j = 0; j < ecm->n; j++ )
		{
			myval = exp(-1. *  fmax(0.,
						pow((field->x((void*)field, i) - field->x((void*)field, j))/ ecm->cl_x_m ,2.)
					 + 	pow((field->y((void*)field, i) - field->y((void*)field, j))/ ecm->cl_y_m ,2.)
					 + 	pow((field->z((void*)field, i) - field->z((void*)field, j))/ ecm->cl_z_m ,2.)
					 + 	pow((field->t((void*)field, i) - field->t((void*)field, j))/ ecm->cl_t_s ,2.))
					);
					
			if (myval > ecm->c_threshold) {											
				myval *= err[i] * err[j];
					
				if (!cs_entry (errcovmattriplet, i, j, myval)) {
					printf("Error with matrix allocation");
					exit(1);
				}
			}
		}
	}
	
	//compress matrices
    ecm->mat = cs_compress (errcovmattriplet) ;             

    //inform user
	//printf("errcovmattriplet entries = %i\n" , errcovmattriplet->nzmax);
	
	//clear triplet matrices
    cs_spfree (errcovmattriplet) ; 
	
    ecm->mat_S = cs_sqr( 1, ecm->mat, 1);            	 				
    ecm->mat_N = cs_qr (ecm->mat, ecm->mat_S); 		

    if (!(ecm->mat_S && ecm->mat_N)) {
		printf("Memory error!"); exit(1);
	}
}


void util_safe_free(void **ptr)
{
	if (*ptr != NULL) {
	  free(*ptr);
	  *ptr = NULL;
	}
}

void util_copy_dbl_array(int n, double **pdst, double *src)
{
	double *dst;
	int i;
	
	if ((n == 0) | (src == NULL)) {
		dst = NULL;
	} else {
		func_dbl_arr_calloc(n, &dst);
		for(i=0;i<n;i++){
			dst[i] = src[i];
		}
		*pdst = dst;
	}
}

