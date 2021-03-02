#include <stdlib.h>
#include <string.h> 
#include <execinfo.h>
#include <stdio.h>
#include <math.h>

#include <nlopt.h>

#include "util.h"
#include "util_turbulence.h"
#include "specialfunctions.h"
#include "func.h"


//uncomment next statement for debug mode
//#define _ZEPHYROS_UTIL_DEBUG

void util_initialize_field_errcovmatrix(t_zephyros_field_errcovmat **pecm)
{
	t_zephyros_field_errcovmat	*ecm;
	ecm = calloc(1, sizeof(t_zephyros_field_errcovmat));

	//set standard values
	ecm->cl_x_m 		= 1.e-20;
	ecm->cl_y_m 		= 1.e-20;
	ecm->cl_z_m 		= 1.e-20;
	ecm->cl_t_s 		= 1.e-20;
	ecm->c_threshold 	= 1. - 1.e-5;
	
	ecm->n				= 0;
	ecm->prepared		= 0;
	ecm->mat 			= NULL;
	ecm->mat_S 			= NULL;
	ecm->mat_N 			= NULL;

	strcpy(ecm->name, "Untitled");
	ecm->initialized = 1;
	
	*pecm = ecm;    
}

void util_assert_initialized_errcovmatrix(t_zephyros_field_errcovmat *ecm)
{
	if (ecm == NULL) {
		printf("ECM was NULL pointer. Exiting.\n"); 
		fflush(stdout); exit(0);
	}
	if (ecm->initialized != 1) {
		printf("ECM was not initialized. Exiting.\n");
		fflush(stdout); exit(0);
	}
}

void util_assert_prepared_errcovmatrix(t_zephyros_field_errcovmat *ecm)
{
	util_assert_initialized_errcovmatrix(ecm);
	if (ecm->prepared != 1) {
		printf("ECM `%s' was not prepared. Exiting.\n", ecm->name);
		fflush(stdout); exit(0);
	}
}

void util_free_field_errcovmatrix(t_zephyros_field_errcovmat **pecm)
{
	t_zephyros_field_errcovmat *ecm = *pecm;
	if (ecm != NULL) {
		if (ecm->mat != NULL) 	cs_spfree (ecm->mat) ; 
		if (ecm->mat_S != NULL) cs_sfree (ecm->mat_S) ; 
		if (ecm->mat_N != NULL) cs_nfree (ecm->mat_N) ; 
		
		free(ecm);
		*pecm = NULL;
	}
}

/*
void util_copy_field_errcovmatrix(t_zephyros_field_errcovmat **pdst, t_zephyros_field_errcovmat *src)
{
	int i;
	t_zephyros_field_errcovmat *dst;
		
	if (src == NULL) {
		dst = NULL;
	} else {
		util_assert_initialized_errcovmatrix(src);
		
		dst = calloc(1, sizeof(t_zephyros_field_errcovmat));
		memcpy(dst, src, sizeof(t_zephyros_field_errcovmat));
		
		cs *mat;
		css *mat_S;
		csn *mat_N;      
    
		//not sure how to copy these foreign objects yet ...
				
	}
	*pdst = dst;	
}
*/

void util_prepare_field_errcovmatrix(t_zephyros_field *field, double *err, t_zephyros_field_errcovmat *ecm)
{
    cs *errcovmattriplet; 
	double myval;
	int i, j;
	int nelements;
	
	util_assert_initialized_errcovmatrix(ecm);
	fields_assert_prepared(field);
	if (err == NULL) {
		printf("Given err is NULL pointer\n"); 
		fflush(stdout); exit(0);
	}
	
	ecm->n = field->n;
		
    //Allocate empty matrices    
    errcovmattriplet 		= cs_spalloc (ecm->n, ecm->n, 1, 1, 1) ;
    if (!(errcovmattriplet)) {
		printf("Memory error!"); 
		fflush(stdout); exit(1);
	}

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
					fflush(stdout); exit(1);
				}
			}
		}
	}
	
	//compress matrices
    ecm->mat = cs_compress (errcovmattriplet) ;             

    if (!(ecm->mat)) {
		printf("Memory error!"); 
		fflush(stdout); exit(1);
	}
	
	nelements = errcovmattriplet->nzmax;
	
	//clear triplet matrices
    cs_spfree (errcovmattriplet) ; 
	
    ecm->mat_S = cs_sqr( 1, ecm->mat, 1);            	 				
    ecm->mat_N = cs_qr (ecm->mat, ecm->mat_S); 		

    if (!(ecm->mat_S && ecm->mat_N)) {
		printf("Memory error!"); 
		fflush(stdout); exit(1);
	}

    //inform user
	printf("Error covariance matrix created with %i elements for field `%s' with %i grid points \n" ,
		nelements, field->name, field->n);
		
	ecm->prepared=1;	
}



void util_initialize_radarfilter(t_zephyros_radarfilter **pr, int i_simulation_retrieval)
{
	//simulation: i_simulation_retrieval = 0
	//retrieval: i_simulation_retrieval = 1
	t_zephyros_radarfilter	*r;
	r = calloc(1, sizeof(t_zephyros_radarfilter));	

	r->geostrophic_advection 	= 0;
	r->coriolis_parameter 		= 0.0001148 ;
	r->cross_sections 			= 1;
	
	r->n_beam_range 			= 2;
	r->n_beam_theta 			= 3;
	r->n_beam_phi 				= 3;
	r->n_t 						= 1;

	r->n_parmod_az 				= 8;
	r->n_parmod_el 				= 16;
	
	r->n_spectrum 				= 8;

	r->filter_dBZ_hh 			= 1;
	r->filter_dBZ_hv 			= 0;
	r->filter_dBZ_vh 			= 0;
	r->filter_dBZ_vv 			= 0;
	r->filter_dBZdr 			= 0;
	r->filter_dBLdr 			= 0;
	r->filter_rho_co 			= 0;
	r->filter_rho_cxh 			= 0;
	r->filter_rho_cxv 			= 0;
	r->filter_KDP	 			= 0;
	
	r->filter_Doppler_velocity_hh_ms		= 1;
	r->filter_Doppler_velocity_hv_ms		= 0;
	r->filter_Doppler_velocity_vh_ms		= 0;
	r->filter_Doppler_velocity_vv_ms		= 0;
	r->filter_Doppler_spectralwidth_hh_ms	= 0;
	r->filter_Doppler_spectralwidth_hv_ms	= 0;
	r->filter_Doppler_spectralwidth_vh_ms	= 0;
	r->filter_Doppler_spectralwidth_vv_ms	= 0;
	r->filter_Doppler_spectrum_dBZ_hh		= 0;
	r->filter_Doppler_spectrum_dBZ_hv		= 0;
	r->filter_Doppler_spectrum_dBZ_vh		= 0;
	r->filter_Doppler_spectrum_dBZ_vv		= 0;
	r->filter_specific_dBZdr				= 0;
	r->filter_specific_dBLdr				= 0;
	r->filter_specific_rho_co				= 0;
	r->filter_specific_rho_cxh				= 0;
	r->filter_specific_rho_cxv				= 0;	
	
	r->filter_errors 						= 1;
	r->effective_earth_correction 			= 1;
	r->inertia_effect 						= 0;
	if (i_simulation_retrieval == 0) {
		r->additive_noise 				= 1;
	} else {
		r->additive_noise 				= 0;		
	}
	r->multiplicative_noise 		= 0;
	
	*pr = r;
}



void util_initialize_atmosphere(t_zephyros_atmosphere **pa)
{
	t_zephyros_atmosphere	*a;
	a = calloc(1, sizeof(t_zephyros_atmosphere));

	fields_initialize(&a->field);
	strcpy(a->field->name, "atmosphere");
	a->grid_T_K 				= NULL;
	a->grid_pair_hPa 			= NULL;
	a->grid_pvapor_hPa 			= NULL;
	a->lut_T_K 					= NULL;
	a->lut_ln_grid_pair_hPa 	= NULL;
	a->lut_ln_grid_pvapor_hPa 	= NULL;

	*pa = a;
}
void util_free_atmosphere(t_zephyros_atmosphere **pa)
{
	t_zephyros_atmosphere *a = *pa;
	if (a != NULL) {

		fields_free(&a->field);
		util_safe_free(&a->grid_T_K);
		util_safe_free(&a->grid_pair_hPa);
		util_safe_free(&a->grid_pvapor_hPa);
		interpolation_free_lut(&a->lut_T_K);
		interpolation_free_lut(&a->lut_ln_grid_pair_hPa);
		interpolation_free_lut(&a->lut_ln_grid_pvapor_hPa);

		free(a);
		*pa = NULL;
	}	
}

void util_prepare_atmosphere(t_zephyros_atmosphere *a)
{
	int i;
	
	if (a->grid_T_K != NULL) {
		fields_prepare(a->field);
		
		//add small number, so that interpolation of log values goes ok.
		for ( i = 0; i < a->field->n; i++ ) {		
			a->grid_pair_hPa[i] += 1.e-5;
			a->grid_pvapor_hPa[i] += 1.e-5;
		}
		
		util_field2lut(a->field, a->grid_T_K, 0, &a->lut_T_K);
		util_field2lut(a->field, a->grid_pair_hPa, 3, &a->lut_ln_grid_pair_hPa);
		util_field2lut(a->field, a->grid_pvapor_hPa, 3, &a->lut_ln_grid_pvapor_hPa);
	}
}


//	type;		//0 = normal, 1 = prior, 2 = post
void util_initialize_windfield(t_zephyros_windfield **pwindfield)
{
	t_zephyros_windfield *windfield = calloc(1, sizeof(t_zephyros_windfield));
	int i;
	 
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_initialize_windfield\n"); fflush(stdout);
	#endif	

	//For all types
	windfield->nfields 	= 0;
	windfield->field 	= calloc(101, sizeof(t_zephyros_windfield*));
	
	windfield->grid_u 	= calloc(101, sizeof(double*));
	windfield->grid_v 	= calloc(101, sizeof(double*));
	windfield->grid_w 	= calloc(101, sizeof(double*));
	
	windfield->grid_hspeed 	= calloc(101, sizeof(double*));
	windfield->grid_hdir  	= calloc(101, sizeof(double*));
	
	windfield->lut_u 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	windfield->lut_v 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	windfield->lut_w 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	
	windfield->lut_hspeed 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	windfield->lut_hdir  	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	
	windfield->nturbulences = 0; 
	windfield->turbulence 	= calloc(101, sizeof(t_zephyros_turbulence_widget*));

	windfield->nekmanspirals = 0;
	windfield->nvortices = 0; 
	windfield->nwaves = 0;
		
	windfield->ekmanspiral 	= calloc(101, sizeof(t_zephyros_ekmanspiral*));	
	windfield->vortex 		= calloc(101, sizeof(t_zephyros_vortex*));
	windfield->wave 		= calloc(101, sizeof(t_zephyros_wave*));

	windfield->grid_u_err 	= calloc(101, sizeof(double*));
	windfield->grid_v_err 	= calloc(101, sizeof(double*));
	windfield->grid_w_err 	= calloc(101, sizeof(double*));

	windfield->lut_u_err 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	windfield->lut_v_err 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	windfield->lut_w_err 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));

	windfield->grid_hspeed_err 			= calloc(101, sizeof(double*));
	windfield->grid_hdir_err			= calloc(101, sizeof(double*));

	windfield->lut_hspeed_err 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	windfield->lut_hdir_err 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));

	windfield->fit_u	= calloc(101, sizeof(int));
	windfield->fit_v	= calloc(101, sizeof(int));
	windfield->fit_w	= calloc(101, sizeof(int));
	
	windfield->fit_hspeed_hdir	= calloc(101, sizeof(int));
	
	windfield->fit_u_Knr	= calloc(101, sizeof(int));
	windfield->fit_v_Knr	= calloc(101, sizeof(int));
	windfield->fit_w_Knr	= calloc(101, sizeof(int));
	
	windfield->fit_hspeed_Knr	= calloc(101, sizeof(int));
	windfield->fit_hdir_Knr	= calloc(101, sizeof(int));
	
	windfield->hspeed_ecm	= calloc(101, sizeof(t_zephyros_field_errcovmat*));
	windfield->hdir_ecm		= calloc(101, sizeof(t_zephyros_field_errcovmat*));
	windfield->u_ecm		= calloc(101, sizeof(t_zephyros_field_errcovmat*));
	windfield->v_ecm		= calloc(101, sizeof(t_zephyros_field_errcovmat*));
	windfield->w_ecm		= calloc(101, sizeof(t_zephyros_field_errcovmat*));

	for ( i = 0; i < 101; i++ ) {
		windfield->field[i] 				= NULL;
		
		windfield->grid_u[i] 				= NULL;
		windfield->grid_v[i] 				= NULL;
		windfield->grid_w[i] 				= NULL;
		
		windfield->grid_hspeed[i] 			= NULL;
		windfield->grid_hdir[i] 			= NULL;
		
		windfield->lut_u[i] 				= NULL;
		windfield->lut_v[i] 				= NULL;
		windfield->lut_w[i] 				= NULL;
		
		windfield->lut_hspeed[i] 			= NULL;
		windfield->lut_hdir[i]	 			= NULL;
		
		windfield->turbulence[i] 			= NULL;
		
		windfield->ekmanspiral[i] 			= NULL;
		windfield->vortex[i] 				= NULL;
		windfield->wave[i] 					= NULL;
		
		windfield->grid_u_err[i] 			= NULL;
		windfield->grid_v_err[i] 			= NULL;
		windfield->grid_w_err[i] 			= NULL;
		windfield->lut_u_err[i] 			= NULL;
		windfield->lut_v_err[i] 			= NULL;
		windfield->lut_w_err[i] 			= NULL;

		windfield->grid_hspeed_err[i] 		= NULL;
		windfield->grid_hdir_err[i] 		= NULL;
		windfield->lut_hspeed_err[i] 		= NULL;
		windfield->lut_hdir_err[i] 			= NULL;
		
		windfield->fit_u[i]			= 0;
		windfield->fit_v[i]			= 0;
		windfield->fit_w[i]			= 0;	
		
		windfield->fit_hspeed_hdir[i]			= 0;	
			
		windfield->fit_u_Knr[i]		= -1;
		windfield->fit_v_Knr[i]		= -1;
		windfield->fit_w_Knr[i]		= -1;
		
		windfield->fit_hspeed_Knr[i]	= -1;
		windfield->fit_hdir_Knr[i]		= -1;
		
		windfield->hspeed_ecm[i] 			= NULL;
		windfield->hdir_ecm[i] 				= NULL;
		windfield->u_ecm[i] 				= NULL;
		windfield->v_ecm[i] 				= NULL;
		windfield->w_ecm[i] 				= NULL;
	}

	*pwindfield = windfield;
}

void util_free_windfield(t_zephyros_windfield **pwindfield)
{
	t_zephyros_windfield *windfield = *pwindfield;
	int i;
	
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_free_windfield\n"); fflush(stdout);
	#endif	
	
	if (windfield != NULL) {
		for ( i = 0; i < 101; i++ ) {
			if (windfield->field != NULL) fields_free(&windfield->field[i]);
			
			if (windfield->grid_u != NULL) util_safe_free(&windfield->grid_u[i]);
			if (windfield->grid_v != NULL) util_safe_free(&windfield->grid_v[i]);
			if (windfield->grid_w != NULL) util_safe_free(&windfield->grid_w[i]);
			
			if (windfield->grid_hspeed != NULL) util_safe_free(&windfield->grid_hspeed[i]);
			if (windfield->grid_hdir != NULL) util_safe_free(&windfield->grid_hdir[i]);

			if (windfield->lut_u != NULL) interpolation_free_lut(&windfield->lut_u[i]);
			if (windfield->lut_v != NULL) interpolation_free_lut(&windfield->lut_v[i]);
			if (windfield->lut_w != NULL) interpolation_free_lut(&windfield->lut_w[i]);
			
			if (windfield->lut_hspeed != NULL) interpolation_free_lut(&windfield->lut_hspeed[i]);
			if (windfield->lut_hdir != NULL) interpolation_free_lut(&windfield->lut_hdir[i]);

			if (windfield->ekmanspiral != NULL) 	util_safe_free(&windfield->ekmanspiral[i]);				//TBD: should have dedicated function
			if (windfield->vortex != NULL)			util_safe_free(&windfield->vortex[i]);					//TBD: should have dedicated function
			if (windfield->wave != NULL) 			util_safe_free(&windfield->wave[i]);					//TBD: should have dedicated function
			if (windfield->turbulence != NULL)		util_free_turbulence_widget(windfield->turbulence + i);

			if (windfield->grid_u_err != NULL)		util_safe_free(&windfield->grid_u_err[i]);
			if (windfield->grid_v_err != NULL)		util_safe_free(&windfield->grid_v_err[i]);
			if (windfield->grid_w_err != NULL)		util_safe_free(&windfield->grid_w_err[i]);
			if (windfield->lut_u_err != NULL)		interpolation_free_lut(&windfield->lut_u_err[i]);
			if (windfield->lut_v_err != NULL)		interpolation_free_lut(&windfield->lut_v_err[i]);
			if (windfield->lut_w_err != NULL)		interpolation_free_lut(&windfield->lut_w_err[i]);

			if (windfield->grid_hspeed_err != NULL) 	util_safe_free(&windfield->grid_hspeed_err[i]);
			if (windfield->grid_hdir_err != NULL) 		util_safe_free(&windfield->grid_hdir_err[i]);
			if (windfield->lut_hspeed_err != NULL) 		interpolation_free_lut(&windfield->lut_hspeed_err[i]);
			if (windfield->lut_hdir_err != NULL) 		interpolation_free_lut(&windfield->lut_hdir_err[i]);

			if (windfield->hspeed_ecm != NULL) 			util_free_field_errcovmatrix(&(windfield->hspeed_ecm[i]));
			if (windfield->hdir_ecm != NULL) 			util_free_field_errcovmatrix(&(windfield->hdir_ecm[i]));
			if (windfield->u_ecm != NULL) 				util_free_field_errcovmatrix(&(windfield->u_ecm[i]));
			if (windfield->v_ecm != NULL) 				util_free_field_errcovmatrix(&(windfield->v_ecm[i]));
			if (windfield->w_ecm != NULL) 				util_free_field_errcovmatrix(&(windfield->w_ecm[i]));
		}
		
		util_safe_free(&windfield->field);
		
		util_safe_free(&windfield->grid_u);
		util_safe_free(&windfield->grid_v);
		util_safe_free(&windfield->grid_w);
		
		util_safe_free(&windfield->grid_hspeed);
		util_safe_free(&windfield->grid_hdir);
		
		util_safe_free(&windfield->lut_u);
		util_safe_free(&windfield->lut_v);
		util_safe_free(&windfield->lut_w);
		
		util_safe_free(&windfield->lut_hspeed);
		util_safe_free(&windfield->lut_hdir);
		
		util_safe_free(&windfield->ekmanspiral);
		util_safe_free(&windfield->vortex);
		util_safe_free(&windfield->wave);
		util_safe_free(&windfield->turbulence);
		
		util_safe_free(&windfield->grid_u_err);
		util_safe_free(&windfield->grid_v_err);
		util_safe_free(&windfield->grid_w_err);
		
		util_safe_free(&windfield->lut_u_err);
		util_safe_free(&windfield->lut_v_err);
		util_safe_free(&windfield->lut_w_err);
		
		util_safe_free(&windfield->grid_hspeed_err);
		util_safe_free(&windfield->grid_hdir_err);
		util_safe_free(&windfield->lut_hspeed_err);
		util_safe_free(&windfield->lut_hdir_err);

		util_safe_free(&windfield->fit_u);
		util_safe_free(&windfield->fit_v);
		util_safe_free(&windfield->fit_w);
		util_safe_free(&windfield->fit_hspeed_hdir);
		
		util_safe_free(&windfield->fit_u_Knr);
		util_safe_free(&windfield->fit_v_Knr);
		util_safe_free(&windfield->fit_w_Knr);
		util_safe_free(&windfield->fit_hspeed_Knr);
		util_safe_free(&windfield->fit_hdir_Knr);
		
		util_safe_free(&windfield->hspeed_ecm);
		util_safe_free(&windfield->hdir_ecm);
		util_safe_free(&windfield->u_ecm);
		util_safe_free(&windfield->v_ecm);
		util_safe_free(&windfield->w_ecm);
		
		free(windfield);
		*pwindfield = NULL;
	}
}

void util_prepare_windfield_i(t_zephyros_windfield *windfield, int i)
{

	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_prepare_windfield_i\n"); fflush(stdout);
	#endif	
	
	if (windfield->field[i] == NULL) {
		fields_initialize(&(windfield->field[i]));
		strcpy(windfield->field[i]->name, "windfield");		

		util_initialize_field_errcovmatrix(&(windfield->hspeed_ecm[i]));
		util_initialize_field_errcovmatrix(&(windfield->hdir_ecm[i]));
		util_initialize_field_errcovmatrix(&(windfield->u_ecm[i]));
		util_initialize_field_errcovmatrix(&(windfield->v_ecm[i]));
		util_initialize_field_errcovmatrix(&(windfield->w_ecm[i]));
	}
}

void util_prepare_windfield(t_zephyros_windfield *windfield)
{
	int i, j;
	double tmp;
	double tmp_hspeed;
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_prepare_windfield\n"); fflush(stdout);
	#endif	
	//wind field
	for ( i = 0; i < windfield->nfields; i++ ) {
		if (windfield->field[i] != NULL) {
			fields_prepare(windfield->field[i]);

			if (windfield->fit_hspeed_hdir[i]) {
				//grid_hspeed
				if (windfield->grid_hspeed[i] == NULL) {
					windfield->grid_hspeed[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_hspeed[i][j] = 0.001;
				}							
				if (windfield->grid_hspeed_err[i] == NULL) {
					windfield->grid_hspeed_err[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_hspeed_err[i][j] = 1.e5;
				}

				//grid_hdir
				if (windfield->grid_hdir[i] == NULL) {
					windfield->grid_hdir[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_hdir[i][j] = 0.001;
				}				
				if (windfield->grid_hdir_err[i] == NULL) {
					windfield->grid_hdir_err[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_hdir_err[i][j] = 1.e5;
				}
				
				//update grid_u and grid_v
				//hspeed and hdir are casted to grid_u and grid_v
				util_cast_hspeed_hdir_2_u_v(windfield, i);
			}
			
			if (windfield->fit_u[i]) {
				if (windfield->grid_u[i] == NULL) {
					windfield->grid_u[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_u[i][j] = 0.001;
				}				
				if (windfield->grid_u_err[i] == NULL) {
					windfield->grid_u_err[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_u_err[i][j] =  1.e5;
				}
			}
			
			if (windfield->fit_v[i]) {
				if (windfield->grid_v[i] == NULL) {
					windfield->grid_v[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_v[i][j] = 0.001;
				}			
				if (windfield->grid_v_err[i] == NULL) {
					windfield->grid_v_err[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_v_err[i][j] =  1.e5;
				}
			}
			
			if (windfield->fit_w[i]) {
				if (windfield->grid_w[i] == NULL) {
					windfield->grid_w[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_w[i][j] = 0.001;
				}
				if (windfield->grid_w_err[i] == NULL) {
					windfield->grid_w_err[i] = calloc(windfield->field[i]->n, sizeof(double));
					for (j = 0; j < windfield->field[i]->n; j++) 
						windfield->grid_w_err[i][j] = 1.e5;
				}
			}

			if (windfield->grid_u[i] != NULL) 		util_field2lut(windfield->field[i], windfield->grid_u[i], 0, windfield->lut_u + i);
			if (windfield->grid_v[i] != NULL) 		util_field2lut(windfield->field[i], windfield->grid_v[i], 0, windfield->lut_v + i);
			if (windfield->grid_w[i] != NULL) 		util_field2lut(windfield->field[i], windfield->grid_w[i], 0, windfield->lut_w + i);
			if (windfield->grid_hspeed[i] != NULL) 	util_field2lut(windfield->field[i], windfield->grid_hspeed[i], 0, windfield->lut_hspeed + i);
			if (windfield->grid_hdir[i] != NULL) 	util_field2lut(windfield->field[i], windfield->grid_hdir[i], 1, windfield->lut_hdir + i);

			if (windfield->grid_u_err[i] != NULL) 		util_field2lut(windfield->field[i], windfield->grid_u_err[i], 0, windfield->lut_u_err + i);
			if (windfield->grid_v_err[i] != NULL) 		util_field2lut(windfield->field[i], windfield->grid_v_err[i], 0, windfield->lut_v_err + i);
			if (windfield->grid_w_err[i] != NULL) 		util_field2lut(windfield->field[i], windfield->grid_w_err[i], 0, windfield->lut_w_err + i);
			if (windfield->grid_hspeed_err[i] != NULL) 	util_field2lut(windfield->field[i], windfield->grid_hspeed_err[i], 0, windfield->lut_hspeed_err + i);
			if (windfield->grid_hdir_err[i] != NULL) 	util_field2lut(windfield->field[i], windfield->grid_hdir_err[i], 0, windfield->lut_hdir_err + i);

			if (windfield->grid_u_err[i] != NULL)
				util_prepare_field_errcovmatrix(windfield->field[i], windfield->grid_u_err[i], windfield->u_ecm[i]);
			if (windfield->grid_v_err[i] != NULL)
				util_prepare_field_errcovmatrix(windfield->field[i], windfield->grid_v_err[i], windfield->v_ecm[i]);
			if (windfield->grid_w_err[i] != NULL)
				util_prepare_field_errcovmatrix(windfield->field[i], windfield->grid_w_err[i], windfield->w_ecm[i]);
				
			if (windfield->grid_hspeed_err[i] != NULL)
				util_prepare_field_errcovmatrix(windfield->field[i], windfield->grid_hspeed_err[i], windfield->hspeed_ecm[i]);
			if (windfield->grid_hdir_err[i] != NULL)
				util_prepare_field_errcovmatrix(windfield->field[i], windfield->grid_hdir_err[i], windfield->hdir_ecm[i]);
				
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
			util_prepare_turbulence_widget(windfield->turbulence[i]);
		}
	}
}

void util_prepare_post_windfield(t_zephyros_windfield **pdst, t_zephyros_windfield *src)
{
	t_zephyros_windfield *dst;
	int i;

	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_prepare_post_windfield\n"); fflush(stdout);
	#endif	
	
	dst = calloc(1,sizeof(t_zephyros_windfield));
	memcpy(dst, src, sizeof(t_zephyros_windfield));
	dst->type = 2;

	dst->field 	= calloc(101, sizeof(t_zephyros_windfield*));
	
	dst->grid_u 		= calloc(101, sizeof(double*));
	dst->grid_v 		= calloc(101, sizeof(double*));
	dst->grid_w 		= calloc(101, sizeof(double*));
	dst->grid_hspeed 	= calloc(101, sizeof(double*));
	dst->grid_hdir 		= calloc(101, sizeof(double*));
	
	dst->lut_u 			= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	dst->lut_v 			= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	dst->lut_w 			= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	dst->lut_hspeed 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	dst->lut_hdir 		= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	
	dst->turbulence 	= calloc(101, sizeof(t_zephyros_turbulence_widget*));

	dst->ekmanspiral 	= calloc(101, sizeof(t_zephyros_ekmanspiral*));	
	dst->vortex 		= calloc(101, sizeof(t_zephyros_vortex*));
	dst->wave 			= calloc(101, sizeof(t_zephyros_wave*));

	dst->grid_u_err 	= calloc(101, sizeof(double*));
	dst->grid_v_err 	= calloc(101, sizeof(double*));
	dst->grid_w_err 	= calloc(101, sizeof(double*));

	dst->lut_u_err 		= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	dst->lut_v_err 		= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	dst->lut_w_err 		= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));

	dst->grid_hspeed_err 	= calloc(101, sizeof(double*));
	dst->grid_hdir_err		= calloc(101, sizeof(double*));

	dst->lut_hspeed_err 	= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));
	dst->lut_hdir_err 		= calloc(101, sizeof(t_zephyros_interpolation_bilint_lut*));

	dst->fit_u	= calloc(101, sizeof(int));
	dst->fit_v	= calloc(101, sizeof(int));
	dst->fit_w	= calloc(101, sizeof(int));
	
	dst->fit_hspeed_hdir	= calloc(101, sizeof(int));
	
	dst->fit_u_Knr	= calloc(101, sizeof(int));
	dst->fit_v_Knr	= calloc(101, sizeof(int));
	dst->fit_w_Knr	= calloc(101, sizeof(int));
	
	dst->fit_hspeed_Knr	= calloc(101, sizeof(int));
	dst->fit_hdir_Knr	= calloc(101, sizeof(int));
	
	dst->hspeed_ecm	= calloc(101, sizeof(t_zephyros_field_errcovmat*));
	dst->hdir_ecm	= calloc(101, sizeof(t_zephyros_field_errcovmat*));
	dst->u_ecm		= calloc(101, sizeof(t_zephyros_field_errcovmat*));
	dst->v_ecm		= calloc(101, sizeof(t_zephyros_field_errcovmat*));
	dst->w_ecm		= calloc(101, sizeof(t_zephyros_field_errcovmat*));

	for ( i = 0; i < 101; i++ ) {
		//reset
		dst->field[i] = NULL;
		
		dst->grid_u[i] = NULL;
		dst->grid_v[i] = NULL;
		dst->grid_w[i] = NULL;
		dst->grid_hspeed[i] = NULL;
		dst->grid_hdir[i] = NULL;
		
		dst->lut_u[i] = NULL;
		dst->lut_v[i] = NULL;
		dst->lut_w[i] = NULL;
		dst->lut_hspeed[i] = NULL;
		dst->lut_hdir[i] = NULL;

		dst->turbulence[i] = NULL;
		dst->ekmanspiral[i] = NULL;
		dst->vortex[i] = NULL;
		dst->wave[i] = NULL;					

		dst->grid_u_err[i] = NULL;
		dst->grid_v_err[i] = NULL;
		dst->grid_w_err[i] = NULL;
		dst->lut_u_err[i] = NULL;
		dst->lut_v_err[i] = NULL;
		dst->lut_w_err[i] = NULL;

		dst->grid_hspeed_err[i] = NULL;
		dst->grid_hdir_err[i] = NULL;		
		dst->lut_hspeed_err[i] = NULL;
		dst->lut_hdir_err[i] = NULL;

		dst->fit_u[i]		= 0;
		dst->fit_v[i]		= 0;
		dst->fit_w[i]		= 0;
		dst->fit_hspeed_hdir[i]		= 0;
		
		dst->fit_u_Knr[i]		= -1;
		dst->fit_v_Knr[i]		= -1;
		dst->fit_w_Knr[i]		= -1;
		dst->fit_hspeed_Knr[i]		= -1;
		dst->fit_hdir_Knr[i]		= -1;

		dst->u_ecm[i] = NULL;
		dst->v_ecm[i] = NULL;
		dst->w_ecm[i] = NULL;		
		dst->hspeed_ecm[i] = NULL;
		dst->hdir_ecm[i] = NULL;

		//copy what is needed
		fields_copy(&(dst->field[i]), src->field[i]);
					
		if (src->field[i] != NULL) {
			util_copy_dbl_array(src->field[i]->n, &(dst->grid_u[i]), src->grid_u[i]);
			util_copy_dbl_array(src->field[i]->n, &(dst->grid_v[i]), src->grid_v[i]);
			util_copy_dbl_array(src->field[i]->n, &(dst->grid_w[i]), src->grid_w[i]);
			util_copy_dbl_array(src->field[i]->n, &(dst->grid_hspeed[i]), src->grid_hspeed[i]);
			util_copy_dbl_array(src->field[i]->n, &(dst->grid_hdir[i]), src->grid_hdir[i]);

			//update grid_u and grid_v
			//hspeed and hdir are casted to grid_u and grid_v
			util_cast_hspeed_hdir_2_u_v(dst, i);

			if (dst->grid_u[i] != NULL) util_field2lut(dst->field[i], dst->grid_u[i], 0, dst->lut_u + i);
			if (dst->grid_v[i] != NULL) util_field2lut(dst->field[i], dst->grid_v[i], 0, dst->lut_v + i);
			if (dst->grid_w[i] != NULL) util_field2lut(dst->field[i], dst->grid_w[i], 0, dst->lut_w + i);
			if (dst->grid_hspeed[i] != NULL) util_field2lut(dst->field[i], dst->grid_hspeed[i], 0, dst->lut_hspeed + i);
			if (dst->grid_hdir[i] != NULL) util_field2lut(dst->field[i], dst->grid_hdir[i], 1, dst->lut_hdir + i);

			util_copy_dbl_array(src->field[i]->n, &(dst->grid_u_err[i]), src->grid_u_err[i]);
			util_copy_dbl_array(src->field[i]->n, &(dst->grid_v_err[i]), src->grid_v_err[i]);
			util_copy_dbl_array(src->field[i]->n, &(dst->grid_w_err[i]), src->grid_w_err[i]);

			if (dst->grid_u_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_u[i], 0, dst->lut_u_err + i);
			if (dst->grid_v_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_v[i], 0, dst->lut_v_err + i);
			if (dst->grid_w_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_w[i], 0, dst->lut_w_err + i);

			util_copy_dbl_array(src->field[i]->n, &(dst->grid_hspeed_err[i]), src->grid_hspeed_err[i]);
			util_copy_dbl_array(src->field[i]->n, &(dst->grid_hdir_err[i]), src->grid_hdir_err[i]);

			if (dst->grid_hspeed_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_hspeed_err[i], 0, dst->lut_hspeed_err + i);
			if (dst->grid_hdir_err[i] != NULL) util_field2lut(dst->field[i], dst->grid_hdir_err[i], 0, dst->lut_hdir_err + i);
		}
		
		if (src->turbulence[i] != NULL) {
			if (src->turbulence[i]->type == 5) {
				//TBD: make copy function for turbulence widget ?
				util_initialize_turbulence_widget(&(dst->turbulence[i]));
				dst->turbulence[i]->type = src->turbulence[i]->type;
				fields_copy(&(dst->turbulence[i]->field), src->turbulence[i]->field);

				util_copy_dbl_array(src->turbulence[i]->field->n, &(dst->turbulence[i]->grid_edr), src->turbulence[i]->grid_edr);
				util_copy_dbl_array(src->turbulence[i]->field->n, &(dst->turbulence[i]->grid_edr13), src->turbulence[i]->grid_edr13);
				util_copy_dbl_array(src->turbulence[i]->field->n, &(dst->turbulence[i]->grid_edr13_err), src->turbulence[i]->grid_edr13_err);
				if (dst->turbulence[i]->grid_edr13 != NULL) util_field2lut(dst->turbulence[i]->field, dst->turbulence[i]->grid_edr13, 0, &(dst->turbulence[i]->lut_edr13));
				
				util_copy_dbl_array(src->turbulence[i]->field->n, &(dst->turbulence[i]->grid_kolmogorov_constant), src->turbulence[i]->grid_kolmogorov_constant);
				if (dst->turbulence[i]->grid_kolmogorov_constant != NULL) util_field2lut(dst->turbulence[i]->field, dst->turbulence[i]->grid_kolmogorov_constant, 0, &(dst->turbulence[i]->lut_kolmogorov_constant));
				dst->turbulence[i]->fit_edr13 = src->turbulence[i]->fit_edr13;
				dst->turbulence[i]->edr13_ecm = NULL;
			}
		} 

	}
	
	*pdst = dst;
}


void util_cast_hspeed_hdir_2_u_v(t_zephyros_windfield *windfield, int i) {
	int j;
	
	if (
		(windfield->field[i] != NULL) &&
		(windfield->grid_hspeed[i] != NULL) &&
		(windfield->grid_hdir[i] != NULL)
		) {
			
		//assert that grid_u and grid_v are allocated	
		if (windfield->grid_u[i] == NULL) 
			windfield->grid_u[i] = calloc(windfield->field[i]->n, sizeof(double));
		if (windfield->grid_v[i] == NULL)
			windfield->grid_v[i] = calloc(windfield->field[i]->n, sizeof(double));
			
		for (j = 0; j < windfield->field[i]->n; j++) {
			windfield->grid_u[i][j] = 
				-1. * windfield->grid_hspeed[i][j] * sin(windfield->grid_hdir[i][j]);
			windfield->grid_v[i][j] = 
				-1. * windfield->grid_hspeed[i][j] * cos(windfield->grid_hdir[i][j]);
		}
	}
}

void util_cast_u_v_2_hspeed_hdir(t_zephyros_windfield *windfield, int i) {
	int j;
	
	if (
		(windfield->field[i] != NULL) &&
		(windfield->grid_u[i] != NULL) &&
		(windfield->grid_v[i] != NULL)
		) {
			
		//assert that grid_u and grid_v are allocated	
		if (windfield->grid_hspeed[i] == NULL) 
			windfield->grid_hspeed[i] = calloc(windfield->field[i]->n, sizeof(double));
		if (windfield->grid_hdir[i] == NULL)
			windfield->grid_hdir[i] = calloc(windfield->field[i]->n, sizeof(double));
			
		for (j = 0; j < windfield->field[i]->n; j++) {
			windfield->grid_hspeed[i][j] = 
				sqrt(pow(windfield->grid_u[i][j], 2.) + pow(windfield->grid_v[i][j], 2.));
			windfield->grid_hdir[i][j] = 
				atan2(-1. * windfield->grid_u[i][j], -1. * windfield->grid_v[i][j]);
		}
	}
}






void util_field2lut(t_zephyros_field *field, double *variable, int special, t_zephyros_interpolation_bilint_lut **plut)
{
	t_zephyros_interpolation_bilint_lut *lut;
	int i;
	
	interpolation_initialize_lut(&lut);
	fields_assert_prepared(field);
	
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

	lut->prepared = 1;
	
	*plut = lut;
}


void util_windfield_ekmanspiral_fuvw(t_zephyros_ekmanspiral *ems, double *xyzt, int i_uvw, double *output, int calcderivatives, double *derivatives)
{
	double g;
	g = M_PI / ems->depth_m;
	
	if (i_uvw == 0) {
		*output = 
			(ems->ug_ms * (1. - (exp(-1. * g * xyzt[2]) * cos( g * xyzt[2] ))))
			- (ems->vg_ms * exp(-1. * g * xyzt[2]) * sin( g * xyzt[2] ) );

	}
	if (i_uvw == 1) {
		*output = 
			(ems->ug_ms * (exp(-1. * g * xyzt[2]) * sin(g * xyzt[2]))) +
			(ems->vg_ms * (1. - (exp(-1. * g * xyzt[2]) * cos(g * xyzt[2]))));
	}
 
	if (i_uvw == 2) {
		*output = 0.;
	} 
	/* Not implemented yet ...
	 * 	if (calcderivatives) {
	 * 
	 * }
	 */
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

void util_windfield_fuvw(
	t_zephyros_windfield *windfield, 
	double *xyzt, 
	double *u,
	double *v,
	double *w,
	int calcderivatives,
	double *uderivatives,
	double *vderivatives,
	double *wderivatives,
	int	geostrophic_advection,
	double fcoriolis
	)
{
	int i, j;
	double tmp, *tmpderivatives;
	double tmp2, *tmpderivatives2;
	
	tmpderivatives = malloc(4 * sizeof(double));
	tmpderivatives2 = malloc(4 * sizeof(double));
	
	//set zero
	*u = 0.; *v = 0.; *w = 0.;
	
	if (calcderivatives) {
		for ( i = 0; i < 4; i++ ) {
			uderivatives[i] = 0.;
			vderivatives[i] = 0.;
			wderivatives[i] = 0.;
		}
	}
	
	//walk through fields
	for ( i = 0; i < windfield->nfields; i++ ) {		
		if (windfield->field[i] != NULL) {
			//u
			if (windfield->lut_u[i] != NULL) {
				interpolation_bilint(windfield->lut_u[i],
						xyzt,
						&tmp,
						calcderivatives,
						tmpderivatives);
				*u = *u + tmp;
				if (calcderivatives) {			
					for ( j = 0; j < 4; j++ )
						uderivatives[j] = uderivatives[j] + tmpderivatives[j];
				}
			}
			//v
			if (windfield->lut_v[i] != NULL) {
				interpolation_bilint(windfield->lut_v[i],
						xyzt,
						&tmp,
						calcderivatives,
						tmpderivatives);
				*v = *v + tmp;
				if (calcderivatives) {
					for ( j = 0; j < 4; j++ )
						vderivatives[j] = vderivatives[j] + tmpderivatives[j];
				}
			}
			//w
			if (windfield->lut_w[i] != NULL) {
				interpolation_bilint(windfield->lut_w[i],
						xyzt,
						&tmp,
						calcderivatives,
						tmpderivatives);
				*w = *w + tmp;
				if (calcderivatives) {
					for ( j = 0; j < 4; j++ )
						wderivatives[j] = wderivatives[j] + tmpderivatives[j];
				}			
			}
		}
	}
	
	//walk through ekmanspirals
	for ( i = 0; i < windfield->nekmanspirals; i++ ) {
		if (windfield->ekmanspiral[i] != NULL) {
			util_windfield_ekmanspiral_fuvw(windfield->ekmanspiral[i], xyzt, 0, &tmp, calcderivatives, tmpderivatives);
			*u = *u + tmp;
			/*
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					uderivatives[j] = uderivatives[j] + tmpderivatives[j];
			}
			*/
			util_windfield_ekmanspiral_fuvw(windfield->ekmanspiral[i], xyzt, 1, &tmp, calcderivatives, tmpderivatives);
			*v = *v + tmp;
			/*
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					vderivatives[j] = vderivatives[j] + tmpderivatives[j];
			}
			*/
		}
	}
	
	//walk through waves
	for ( i = 0; i < windfield->nwaves; i++ ) {
		if (windfield->wave[i] != NULL) {
			util_windfield_wave_fuvw(windfield->wave[i], xyzt, 0, &tmp, calcderivatives, tmpderivatives);
			*u = *u + tmp;
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					uderivatives[j] = uderivatives[j] + tmpderivatives[j];
			}
			util_windfield_wave_fuvw(windfield->wave[i], xyzt, 1, &tmp, calcderivatives, tmpderivatives);
			*v = *v + tmp;
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					vderivatives[j] = vderivatives[j] + tmpderivatives[j];
			}
			util_windfield_wave_fuvw(windfield->wave[i], xyzt, 2, &tmp, calcderivatives, tmpderivatives);
			*w = *w + tmp;
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					wderivatives[j] = wderivatives[j] + tmpderivatives[j];
			}
		}
	}
	
	//walk through vortices
	for ( i = 0; i < windfield->nvortices; i++ ) {
		if (windfield->vortex[i] != NULL) {
			util_windfield_vortex_fuvw(windfield->vortex[i], xyzt, 0, &tmp, calcderivatives, tmpderivatives);
			*u = *u + tmp;
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					uderivatives[j] = uderivatives[j] + tmpderivatives[j];
			}
			util_windfield_vortex_fuvw(windfield->vortex[i], xyzt, 1, &tmp, calcderivatives, tmpderivatives);
			*v = *v + tmp;
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					vderivatives[j] = vderivatives[j] + tmpderivatives[j];
			}
			util_windfield_vortex_fuvw(windfield->vortex[i], xyzt, 2, &tmp, calcderivatives, tmpderivatives);
			*w = *w + tmp;
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					wderivatives[j] = wderivatives[j] + tmpderivatives[j];
			}
		}
	}
	
	//walk through turbulences
	for ( i = 0; i < windfield->nturbulences; i++ ) {
		if (windfield->turbulence[i] != NULL) {
			util_turbulence_uvw(windfield->turbulence[i], xyzt, 0, &tmp, calcderivatives, tmpderivatives);
			*u = *u + tmp;
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					uderivatives[j] = uderivatives[j] + tmpderivatives[j];
			}
			util_turbulence_uvw(windfield->turbulence[i], xyzt, 1, &tmp, calcderivatives, tmpderivatives);
			*v = *v + tmp;
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					vderivatives[j] = vderivatives[j] + tmpderivatives[j];
			}
			util_turbulence_uvw(windfield->turbulence[i], xyzt, 2, &tmp, calcderivatives, tmpderivatives);
			*w = *w + tmp;
			if (calcderivatives) {
				for ( j = 0; j < 4; j++ )
					wderivatives[j] = wderivatives[j] + tmpderivatives[j];
			}
		}
	}

	//apply geostrophic correction
	if (geostrophic_advection) {
		*u = (*u * cos(fcoriolis * xyzt[3])) + (*v * sin(fcoriolis * xyzt[3]));
		*v = (-1. * *u * sin(fcoriolis * xyzt[3]))  + (*v * cos(fcoriolis * xyzt[3]));
	}

	free(tmpderivatives);
	free(tmpderivatives2);
	
}




void util_initialize_scattererfield(t_zephyros_scattererfield **pscattererfield)
{
	t_zephyros_scattererfield *scattererfield = calloc(1, sizeof(t_zephyros_scattererfield));
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

	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_prepare_post_scattererfield\n"); fflush(stdout);
	#endif		
	if (src == NULL) {
		dst = NULL;
	} else {
		util_initialize_scattererfield(&dst);
		 
		dst->type = 2;
		dst->npsd = src->npsd;
			
		for ( i = 0; i < 101; i++ ) {
			util_prepare_post_psd(dst->psd + i, src->psd[i]);
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
			util_free_psd(scattererfield->psd + i);
		}
		free(scattererfield);
		*pscattererfield = NULL;
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
	
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_initialize_psd\n"); 
	#endif		
	
	thepsd = calloc(1, sizeof(t_zephyros_psd));
	
	thepsd->distribution_type = 1;
	thepsd->particle_type = -1;
	
	fields_initialize(&(thepsd->field));
	strcpy(thepsd->field->name, "particle size distribution");		
	
	thepsd->grid_lwc_gm3 						= NULL;
	thepsd->grid_dBlwc_err_gm3 					= NULL;
	thepsd->discrete_D_equiv_mm_interval_llimt	= NULL;
	thepsd->discrete_D_equiv_mm					= NULL;
	thepsd->discrete_D_equiv_mm_interval_ulimt	= NULL;
	thepsd->grid_number_density_m3 				= NULL;
	thepsd->grid_dBnumber_density_err_m3 		= NULL;
	thepsd->lut_ln_number_density_m3 			= NULL;
	thepsd->lut_dBnumber_density_err_m3 		= NULL;
	
	thepsd->grid_gammadistribution_N0 		= NULL;
	thepsd->grid_gammadistribution_N0_err 	= NULL;
	thepsd->gammadistribution_mu	 		= 5.;
	thepsd->gammadistribution_mu_err 		= 5.;
	thepsd->gammadistribution_D0_mm		 	= 0.05;	
	thepsd->gammadistribution_D0_err_mm 	= 2.;	
	thepsd->gammadistribution_dmin_mm	 	= 0.01;	
	thepsd->gammadistribution_dmax_mm	 	= 8.;	
	thepsd->n_diameters = 50;

	thepsd->grid_gammadistribution_mu 		= NULL;
	thepsd->grid_gammadistribution_D0_mm 	= NULL;
	
	thepsd->N_constraint 					= 0;
	
	thepsd->fit_dBlwc 						= 0;
	thepsd->fit_dBlwc_Knr 					= -1;
	util_initialize_field_errcovmatrix(&(thepsd->dBlwc_ecm));

	thepsd->fit_dBm 						= 0;
	thepsd->fit_dBm_Knr 					= -1;
	util_initialize_field_errcovmatrix(&(thepsd->dBm_ecm));
	
	thepsd->initialized 					= 1;
	thepsd->prepared	 					= 0;
	
	thepsd->grid_rcs_hh	 				= NULL;
	thepsd->grid_rcs_hv	 				= NULL;
	thepsd->grid_rcs_vv	 				= NULL;
	
	*pthepsd = thepsd;
}

void util_assert_initialized_psd(t_zephyros_psd *thepsd)
{
	if (thepsd == NULL) {
		printf("PSD was NULL pointer. Exiting.\n"); 
		fflush(stdout); exit(0);
	}
	if (thepsd->initialized != 1) {
		printf("PSD was not initialized. Exiting.\n"); 
		fflush(stdout); exit(0);
	}
}

void util_assert_prepared_psd(t_zephyros_psd *thepsd)
{
	util_assert_initialized_psd(thepsd);

	//minimal check
	if (thepsd->discrete_D_equiv_mm_interval_llimt	== NULL) 	thepsd->prepared = 0;
	if (thepsd->discrete_D_equiv_mm	== NULL) 					thepsd->prepared = 0;
	if (thepsd->discrete_D_equiv_mm_interval_ulimt	== NULL) 	thepsd->prepared = 0;
	if (thepsd->grid_number_density_m3	== NULL) 				thepsd->prepared = 0;
		
	if (thepsd->prepared != 1) {
		printf("PSD was not prepared. Exiting.\n");
		fflush(stdout); exit(0);
	}
}

void util_free_psd(t_zephyros_psd **pthepsd)
{
	t_zephyros_psd *thepsd = *pthepsd;
	int j;
	
	if (thepsd != NULL) {
		util_assert_initialized_psd(thepsd);
		
		//walk through deeper structures
		for ( j = 0; j < thepsd->n_diameters; j++ ) {
			if (thepsd->grid_number_density_m3 != NULL) 
				util_safe_free(thepsd->grid_number_density_m3 + j);
			if (thepsd->grid_dBnumber_density_err_m3 != NULL) 
				util_safe_free(thepsd->grid_dBnumber_density_err_m3 + j);
			if (thepsd->lut_ln_number_density_m3 != NULL) 
				interpolation_free_lut(thepsd->lut_ln_number_density_m3 + j);
			if (thepsd->lut_dBnumber_density_err_m3 != NULL) 
				interpolation_free_lut(thepsd->lut_dBnumber_density_err_m3 + j);

			if (thepsd->grid_rcs_hh != NULL) 
				util_safe_free(thepsd->grid_rcs_hh + j);
			if (thepsd->grid_rcs_hv != NULL) 
				util_safe_free(thepsd->grid_rcs_hv + j);
			if (thepsd->grid_rcs_vv != NULL) 
				util_safe_free(thepsd->grid_rcs_vv + j);
		}

		fields_free(&thepsd->field); 
		
		util_safe_free(&thepsd->grid_lwc_gm3);
		util_safe_free(&thepsd->grid_dBlwc_err_gm3);
		util_safe_free(&thepsd->discrete_D_equiv_mm_interval_llimt);
		util_safe_free(&thepsd->discrete_D_equiv_mm);
		util_safe_free(&thepsd->discrete_D_equiv_mm_interval_ulimt);

		util_safe_free(&thepsd->grid_number_density_m3);
		util_safe_free(&thepsd->grid_dBnumber_density_err_m3);
		util_safe_free(&thepsd->lut_ln_number_density_m3);
		util_safe_free(&thepsd->lut_dBnumber_density_err_m3);
		
		util_safe_free(&thepsd->grid_gammadistribution_N0);
		util_safe_free(&thepsd->grid_gammadistribution_N0_err);
		
		util_safe_free(&thepsd->grid_gammadistribution_mu);
		util_safe_free(&thepsd->grid_gammadistribution_D0_mm);

		util_free_field_errcovmatrix(&thepsd->dBlwc_ecm);
		util_free_field_errcovmatrix(&thepsd->dBm_ecm);
		
		util_safe_free(&thepsd->grid_rcs_hh);
		util_safe_free(&thepsd->grid_rcs_hv);
		util_safe_free(&thepsd->grid_rcs_vv);

		free(thepsd);
		*pthepsd = NULL;
	}
}

void util_copy_psd(t_zephyros_psd **pdst, t_zephyros_psd *src)
{
	t_zephyros_psd *dst;
	int i, j;

	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_copy_psd\n");
	#endif	
	
	if (src == NULL) {
		dst = NULL;
	} else {
		dst = calloc(1, sizeof(t_zephyros_psd));
		//copy everything
		memcpy(dst, src, sizeof(t_zephyros_psd));

		//deeper structures
		if (src->grid_number_density_m3 == NULL) {
			dst->grid_number_density_m3 = NULL;
		} else {
			dst->grid_number_density_m3 = calloc(src->n_diameters, sizeof(double*));
		}
		if (src->grid_dBnumber_density_err_m3 == NULL) {
			dst->grid_dBnumber_density_err_m3 = NULL;
		} else {
			dst->grid_dBnumber_density_err_m3 = calloc(src->n_diameters, sizeof(double*));
		}
		if (src->lut_ln_number_density_m3 == NULL) {
			dst->lut_ln_number_density_m3 = NULL;
		} else {
			dst->lut_ln_number_density_m3 = calloc(src->n_diameters, sizeof(t_zephyros_interpolation_bilint_lut*));
		}
		if (src->lut_dBnumber_density_err_m3 == NULL) {
			dst->lut_dBnumber_density_err_m3 = NULL;
		} else {
			dst->lut_dBnumber_density_err_m3 = calloc(src->n_diameters, sizeof(t_zephyros_interpolation_bilint_lut*));
		}
 					
		if (src->grid_rcs_hh == NULL) {
			dst->grid_rcs_hh = NULL;
		} else {
			dst->grid_rcs_hh = calloc(src->n_diameters, sizeof(double*));
		} 					
		if (src->grid_rcs_hv == NULL) {
			dst->grid_rcs_hv = NULL;
		} else {
			dst->grid_rcs_hv = calloc(src->n_diameters, sizeof(double*));
		} 					
		if (src->grid_rcs_vv == NULL) {
			dst->grid_rcs_vv = NULL;
		} else {
			dst->grid_rcs_vv = calloc(src->n_diameters, sizeof(double*));
		} 					
 						
		//copy field
		fields_copy(&dst->field, src->field);
	   
		util_copy_dbl_array(src->field->n, &dst->grid_lwc_gm3, src->grid_lwc_gm3);
		util_copy_dbl_array(src->field->n, &dst->grid_dBlwc_err_gm3, src->grid_dBlwc_err_gm3);
		util_copy_dbl_array(src->n_diameters, &dst->discrete_D_equiv_mm_interval_llimt, src->discrete_D_equiv_mm_interval_llimt);
		util_copy_dbl_array(src->n_diameters, &dst->discrete_D_equiv_mm, src->discrete_D_equiv_mm);
		util_copy_dbl_array(src->n_diameters, &dst->discrete_D_equiv_mm_interval_ulimt, src->discrete_D_equiv_mm_interval_ulimt);
	   
		for (i = 0; i < src->n_diameters; i++ ) {
			if (src->grid_number_density_m3 != NULL)
				util_copy_dbl_array(src->field->n, dst->grid_number_density_m3 + i, src->grid_number_density_m3[i]);
			if (src->grid_dBnumber_density_err_m3 != NULL)
				util_copy_dbl_array(src->field->n, dst->grid_dBnumber_density_err_m3 + i, src->grid_dBnumber_density_err_m3[i]);

			if (src->lut_ln_number_density_m3 != NULL) 
				util_field2lut(dst->field, dst->grid_number_density_m3[i], 3, dst->lut_ln_number_density_m3 + i);
			if (src->lut_dBnumber_density_err_m3 != NULL) 
				util_field2lut(dst->field, dst->grid_dBnumber_density_err_m3[i], 0, dst->lut_dBnumber_density_err_m3 + i);

			if (src->grid_rcs_hh != NULL)
				util_copy_dbl_array(src->field->n, dst->grid_rcs_hh + i, src->grid_rcs_hh[i]);
			if (src->grid_rcs_hv != NULL)
				util_copy_dbl_array(src->field->n, dst->grid_rcs_hv + i, src->grid_rcs_hv[i]);
			if (src->grid_rcs_vv != NULL)
				util_copy_dbl_array(src->field->n, dst->grid_rcs_vv + i, src->grid_rcs_vv[i]);
		}
		
		util_copy_dbl_array(src->field->n, &dst->grid_gammadistribution_N0, src->grid_gammadistribution_N0);
		util_copy_dbl_array(src->field->n, &dst->grid_gammadistribution_N0_err, src->grid_gammadistribution_N0_err);
		util_copy_dbl_array(src->field->n, &dst->grid_gammadistribution_mu, src->grid_gammadistribution_mu);
		util_copy_dbl_array(src->field->n, &dst->grid_gammadistribution_D0_mm, src->grid_gammadistribution_D0_mm);

		//TBD: error covariance matrices are not copied.
		//Typically not necessary when copied from prior to post.
		dst->dBlwc_ecm = NULL;
		dst->dBm_ecm = NULL;
	}
 	
	*pdst = dst;	
}




/*
//The following steps are taken
1. Prepare field
2. Assert that grid_number_density is set
3. Update the discrete probability density function
4. Update grid_number_density with grid_mu grid_D0
5. Assert grid_number_density_err_m3 is set
6. Update grid_number_density_err_m3
7. Update number density with lwc, or update lwc
8. Assert grid_dBlwc_err_gm3
9. make luts
10. make error covariance matrices

*/

void util_prepare_psd(t_zephyros_psd *thepsd, t_zephyros_psd *scattererfield) 
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
	
	//1. Pepare field
	fields_prepare(thepsd->field);

	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("Assert grid_number_density is set\n"); fflush(stdout);
	#endif		

	//2. Assert that grid_number_density is set
	if (thepsd->grid_number_density_m3 == NULL) {
		thepsd->grid_number_density_m3 = calloc(thepsd->n_diameters, sizeof(double*));
		for (i = 0; i < thepsd->n_diameters; i++ )
			if (thepsd->grid_number_density_m3[i] = NULL);
	}
	for (i = 0; i < thepsd->n_diameters; i++ ) {
		if (thepsd->grid_number_density_m3[i] == NULL) {
			thepsd->grid_number_density_m3[i] = calloc(thepsd->field->n, sizeof(double));
			for (j = 0; j < thepsd->field->n; j++ )
				thepsd->grid_number_density_m3[i][j] = 1.;
		}
	}

	//3. Update the discrete probability density function

	//3a. discrete probability density function
	if (thepsd->distribution_type == 0) {

		if (
			(thepsd->discrete_D_equiv_mm_interval_llimt == NULL) |
			(thepsd->discrete_D_equiv_mm_interval_ulimt == NULL)) {

			util_safe_free(&thepsd->discrete_D_equiv_mm_interval_llimt);
			util_safe_free(&thepsd->discrete_D_equiv_mm_interval_ulimt);
			thepsd->discrete_D_equiv_mm_interval_llimt = calloc(thepsd->n_diameters, sizeof(double));
			thepsd->discrete_D_equiv_mm_interval_ulimt = calloc(thepsd->n_diameters, sizeof(double));
		
			
			thepsd->discrete_D_equiv_mm_interval_llimt[0] 	= 
				thepsd->discrete_D_equiv_mm[i] / 100.;
			thepsd->discrete_D_equiv_mm_interval_ulimt[thepsd->n_diameters - 1] = thepsd->gammadistribution_dmax_mm;
			
			for (i = 0; i < thepsd->n_diameters -1 ; i++ ) {
				//set diameter size
				thepsd->discrete_D_equiv_mm_interval_ulimt[i] 
					= thepsd->discrete_D_equiv_mm[i] + ((thepsd->discrete_D_equiv_mm[i + 1] - thepsd->discrete_D_equiv_mm[i]) /2.);
				thepsd->discrete_D_equiv_mm_interval_llimt[i+1] = thepsd->discrete_D_equiv_mm_interval_ulimt[i];
			}
		}
	}
	
	//3b. Translate gamma pdf to discrete probability density function
	if (thepsd->distribution_type == 1) {
		#ifdef _ZEPHYROS_UTIL_DEBUG
			printf("translate gamma pdf to discrete probability density function\n"); fflush(stdout);
		#endif	
	
		util_safe_free(&thepsd->discrete_D_equiv_mm_interval_llimt);
		util_safe_free(&thepsd->discrete_D_equiv_mm_interval_ulimt);

		thepsd->discrete_D_equiv_mm_interval_llimt = calloc(thepsd->n_diameters, sizeof(double));
		thepsd->discrete_D_equiv_mm_interval_ulimt = calloc(thepsd->n_diameters, sizeof(double));


		Dl 		= calloc(thepsd->n_diameters, sizeof(double));
		Dcenter = calloc(thepsd->n_diameters, sizeof(double));
		Du 		= calloc(thepsd->n_diameters, sizeof(double));
		
		//assert thepsd->grid_gammadistribution_N0 is set
		if (thepsd->grid_gammadistribution_N0 == NULL) {
			thepsd->grid_gammadistribution_N0 = malloc(thepsd->field->n * sizeof(double));
			for (j = 0; j < thepsd->field->n; j++ ) {
				thepsd->grid_gammadistribution_N0[j] = 1.;
			}
		}
		
		//assert thepsd->discrete_D_equiv_mm is set
		util_safe_free(&thepsd->discrete_D_equiv_mm);
		thepsd->discrete_D_equiv_mm = calloc(thepsd->n_diameters, sizeof(double));

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

			Dcenter[i] = (Du[i] + Dl[i]) / 2.;
			//set diameter size
			thepsd->discrete_D_equiv_mm[i] = Dcenter[i];
			thepsd->discrete_D_equiv_mm_interval_llimt[i] = Dl[i];
			thepsd->discrete_D_equiv_mm_interval_ulimt[i] = Du[i];
			
			//calculate the number density
			for (j = 0; j < thepsd->field->n; j++ ) {
				thepsd->grid_number_density_m3[i][j] = 
					util_gamma_integral(thepsd->grid_gammadistribution_N0[j],
					 thepsd->gammadistribution_mu, thepsd->gammadistribution_D0_mm, Dl[i], Du[i]);
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
					thepsd->grid_number_density_m3[i][j] = 
						util_gamma_integral(thepsd->grid_gammadistribution_N0[j],
						thepsd->grid_gammadistribution_mu[j], thepsd->grid_gammadistribution_D0_mm[j], Dl[i], Du[i]);
				}						
			}
		}
	}

	//5. Assert grid_dBnumber_density_err_m3 is set
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("4. Assert grid_dBnumber_density_err_m3 is set\n"); fflush(stdout);
	#endif	
	
	if (thepsd->grid_dBnumber_density_err_m3 == NULL) {
		thepsd->grid_dBnumber_density_err_m3 = calloc(thepsd->n_diameters, sizeof(double*));
		for (i = 0; i < thepsd->n_diameters; i++ )
			if (thepsd->grid_dBnumber_density_err_m3[i] = NULL);
	}
	
	for (i = 0; i < thepsd->n_diameters; i++ ) {
		if (thepsd->grid_dBnumber_density_err_m3[i] == NULL)
			thepsd->grid_dBnumber_density_err_m3[i] = calloc(thepsd->field->n, sizeof(double));
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
				thepsd->grid_dBnumber_density_err_m3[i][j] = 1000.;

				/*
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
				*/
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
		thepsd->grid_lwc_gm3 = calloc(thepsd->field->n, sizeof(double));
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
		thepsd->grid_dBlwc_err_gm3 = calloc(thepsd->field->n, sizeof(double));
		for (j = 0; j < thepsd->field->n; j++ ) {
			thepsd->grid_dBlwc_err_gm3[j] = 1000.;
			//When a priori is not set, 
			//essential kill usage of a priori by setting a large error.
		}
	}
			
			/*
			//1. calculate lwc_err_gm3, with error propagation
			//2. convert to error in dB
			
			//1.
			thepsd->grid_dBlwc_err_gm3[j] = 0.;
			for (i = 0; i < thepsd->n_diameters; i++ ) {
				thepsd->grid_dBlwc_err_gm3[j] +=
					pow(
					func_dB_inv(thepsd->grid_dBnumber_density_err_m3[i][j]) *
					1.e3 *			//g kg-1
					998.2071 * 		//[water density in kg m-3]
					(4./3.) * M_PI * pow((thepsd->discrete_D_equiv_mm[i]/2.) * 1.e-3,3.)
					, 2.);
			}
			thepsd->grid_dBlwc_err_gm3[j] = sqrt(thepsd->grid_dBlwc_err_gm3[j]);
			
			//2.
			thepsd->grid_dBlwc_err_gm3[j] = func_dB(thepsd->grid_dBlwc_err_gm3[j] / thepsd->grid_lwc_gm3[j]);
			
			//Error less than 1dB doesn't make sense. Update
			if (thepsd->grid_dBlwc_err_gm3[j] < 1.) {
				thepsd->grid_dBlwc_err_gm3[j] = 1.;
			}
			*/

	//9. make luts
#ifdef _ZEPHYROS_UTIL_DEBUG
	printf("make luts\n"); fflush(stdout);
#endif 	
	thepsd->lut_ln_number_density_m3 = 
		calloc(thepsd->n_diameters, sizeof(t_zephyros_interpolation_bilint_lut*));
	thepsd->lut_dBnumber_density_err_m3 = 
		calloc(thepsd->n_diameters, sizeof(t_zephyros_interpolation_bilint_lut*));
	for (i = 0; i < thepsd->n_diameters; i++ ) {
		util_field2lut(thepsd->field, thepsd->grid_number_density_m3[i], 3, thepsd->lut_ln_number_density_m3 + i);
		util_field2lut(thepsd->field, thepsd->grid_dBnumber_density_err_m3[i], 0, thepsd->lut_dBnumber_density_err_m3 + i);
	}
	
	//10. make error covariance matrices
	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("make error covariance matrices\n"); fflush(stdout);
	#endif	
	if (thepsd->fit_dBlwc) 
		util_prepare_field_errcovmatrix(thepsd->field, thepsd->grid_dBlwc_err_gm3, thepsd->dBlwc_ecm);

	//TBD: should be updated somehow...
	/*
	if (thepsd->fit_dBN) 
		util_prepare_field_errcovmatrix(thepsd->field, thepsd->grid_dBnumber_density_err_m3[0], thepsd->dBN_ecm);
	*/
		
	thepsd->prepared = 1;
	
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

	#ifdef _ZEPHYROS_UTIL_DEBUG
		printf("util_prepare_post_psd\n"); fflush(stdout);
	#endif	
		
	//copy
	util_copy_psd(pdst, src);
	
	//update some settings for post type.
	dst = *pdst;
	if (dst != NULL) dst->distribution_type = 0;	
}

//given number densities, fit the gamma distribution parameters.
//model 1: nothing is fitted. parameters come from the LWC.
//model 2: mu-D0 relation is used. (D0) is fitted.
//model 3: (D0, mu) are fitted.
void util_fit_gammadistribution_parameters(t_zephyros_psd *fit_psd, t_zephyros_psd *ref_psd)
{
	int i,j;
	double *x = NULL;
	double *xu = NULL;
	double *xl = NULL;
	double minf;	
	int inf;
	void 	**somepointers = NULL;
	nlopt_opt opt;

	double sum_a[4];
	double sum_b[4];
	double myintegral;
	int model;
	
	if ((fit_psd == NULL) | (ref_psd == NULL)) {
		printf("Casting not possible.\n");
		fflush(stdout); exit(0);
	}
	somepointers 		= calloc(6, sizeof(void*));
	somepointers[0] 	= (void*)fit_psd;
	somepointers[1] 	= (void*)ref_psd;
	somepointers[2]		= (void*)&model;
	somepointers[3]		= (void*)&j;

	model = fit_psd->N_constraint;
		
	util_assert_prepared_psd(fit_psd);
	util_assert_prepared_psd(ref_psd);
		
	if (fit_psd->field->n != ref_psd->field->n) {
		printf("util_fit_gammadistribution_parameters does not support fitting different field sizes."); 
		fflush(stdout); exit(0);			
	}
		
	if ((model < 1) || (model > 2)) {
		printf("util_fit_gammadistribution_parameters does not support this model"); 
		fflush(stdout); exit(0);
	}
	
	if (model == 1) {
		x 	= calloc(1, sizeof(double));
		xl 	= calloc(1, sizeof(double));
		xu 	= calloc(1, sizeof(double));		
	}

	if (model == 2) {
		x 	= calloc(2, sizeof(double));
		xl 	= calloc(2, sizeof(double));
		xu 	= calloc(2, sizeof(double));		
	}
	
	if (model == 3) {
		x 	= calloc(3, sizeof(double));
		xl 	= calloc(3, sizeof(double));
		xu 	= calloc(3, sizeof(double));		
	}
	
	//assert that gamma distribution parameters are allocated.
	if (fit_psd->grid_gammadistribution_N0 == NULL)
		fit_psd->grid_gammadistribution_N0 =  calloc(fit_psd->field->n, sizeof(double));
	if (fit_psd->grid_gammadistribution_mu == NULL)
		fit_psd->grid_gammadistribution_mu =  calloc(fit_psd->field->n, sizeof(double));
	if (fit_psd->grid_gammadistribution_D0_mm == NULL)
		fit_psd->grid_gammadistribution_D0_mm =  calloc(fit_psd->field->n, sizeof(double));
				
	for (j = 0; j < fit_psd->field->n; j++ ) {
		if (model == 1) {
			if (fit_psd->n_diameters > 0) {			
				//fit D0_mm (which determines N0, mu)
				x[0] = 1.;		xl[0] = 1.e-25; 		xu[0] = 100.;	//D0_mm
						
				//best guess for minimization routine.
				opt = nlopt_create(NLOPT_LN_SBPLX, 1); 

				nlopt_set_min_objective(opt, util_fit_gammadistribution_parameters_costf, (void*)somepointers);
				nlopt_set_lower_bounds(opt, xl);
				nlopt_set_upper_bounds(opt, xu);
				
				//according to manual: nlopt_set_xtol_rel(opt, 1e-4);
				nlopt_set_xtol_rel(opt, 1.e-4);
				nlopt_set_maxtime(opt, 10);

				inf = nlopt_optimize(opt, x, &minf);
				
				nlopt_destroy(opt);
			}
		}


		if (model == 2) {
			if (fit_psd->n_diameters > 1) {
				//first guess
				x[0] = 1.;		xl[0] = 1.e-40; 		xu[0] = 1.e40;	//N0
				x[1] = 1.;		xl[1] = 1.e-25; 		xu[1] = 100.;	//D0_mm
							
				//best guess for minimization routine.
				opt = nlopt_create(NLOPT_LN_SBPLX, 2); 

				nlopt_set_min_objective(opt, util_fit_gammadistribution_parameters_costf, (void*)somepointers);
				nlopt_set_lower_bounds(opt, xl);
				nlopt_set_upper_bounds(opt, xu);
				
				//according to manual: nlopt_set_xtol_rel(opt, 1e-4);
				nlopt_set_xtol_rel(opt, 1.e-10);
				nlopt_set_maxtime(opt, 10);

				inf = nlopt_optimize(opt, x, &minf);

				nlopt_destroy(opt);		
			}		
		}
	
		//~ if (model == 3) {
			//~ if (thepsd->n_diameters > 2) {
				
				//~ x[0] = 1.;		xl[0] = 1.e-20; 		xu[0] = 1.e20;	//N0
				//~ x[1] = 1.;		xl[1] = 1.e-25; 		xu[1] = 100.;	//D0_mm
				//~ x[2] = 5.;		xl[2] = -5.;			xu[2] = 100.;	//mu
							
				//~ //best guess for minimization routine.
				//~ opt = nlopt_create(NLOPT_LN_SBPLX, 3); 

				//~ nlopt_set_min_objective(opt, util_fit_gammadistribution_parameters_costf, (void*)somepointers);
				//~ nlopt_set_lower_bounds(opt, xl);
				//~ nlopt_set_upper_bounds(opt, xu);
				
				//~ //according to manual: nlopt_set_xtol_rel(opt, 1e-4);
				//~ nlopt_set_xtol_rel(opt, 1.e-4);
				//~ nlopt_set_maxtime(opt, 10);

				//~ inf = nlopt_optimize(opt, x, &minf);

				//~ nlopt_destroy(opt);
			//~ }
		//~ }
		
		//Some analysis for debugging 
		//~ sum_a[0] = 0.; sum_b[0] = 0.;
		//~ sum_a[1] = 0.; sum_b[1] = 0.;
		//~ sum_a[2] = 0.; sum_b[2] = 0.;
		//~ sum_a[3] = 0.; sum_b[3] = 0.;
		
		//~ for (i = 0; i < thepsd->n_diameters; i++ ) {			
			//~ myintegral = 
			//~ util_gamma_integral(
			//~ thepsd->grid_gammadistribution_N0[j],		//N0
			//~ thepsd->grid_gammadistribution_mu[j], 		//mu
			//~ thepsd->grid_gammadistribution_D0_mm[j], 	//D0_mm
			//~ thepsd->discrete_D_equiv_mm_interval_llimt[i],
			//~ thepsd->discrete_D_equiv_mm_interval_ulimt[i]);
			
			//~ sum_a[0] += thepsd->grid_number_density_m3[i][j];
			//~ sum_a[1] += thepsd->grid_number_density_m3[i][j] * pow(thepsd->discrete_D_equiv_mm[i], 3.);
			//~ sum_a[2] += thepsd->grid_number_density_m3[i][j] * pow(thepsd->discrete_D_equiv_mm[i], 6.);
			//~ sum_a[3] += thepsd->grid_number_density_m3[i][j] * weights1[i][j];

			//~ sum_b[0] += myintegral;
			//~ sum_b[1] += myintegral * pow(thepsd->discrete_D_equiv_mm[i], 3.);
			//~ sum_b[2] += myintegral * pow(thepsd->discrete_D_equiv_mm[i], 6.);			
			//~ sum_b[3] += myintegral * weights1[i][j];
		//~ }
	
		//~ printf("\n");
		//~ printf("fitted N0 = %.2e\n", thepsd->grid_gammadistribution_N0[j]);
		//~ printf("fitted D0_mm = %.2e\n", thepsd->grid_gammadistribution_D0_mm[j]);
		//~ printf("fitted mu = %.2e\n", thepsd->grid_gammadistribution_mu[j]);
		//~ printf("relative difference in number density: %.2e \%\n", 100. * (sum_b[0] - sum_a[0]) / sum_a[0]);
		//~ printf("relative difference in liqud water content: %.2e \%\n", 100. * (sum_b[1] - sum_a[1]) / sum_a[1]);
		//~ printf("relative difference in RCS (D^6): %.2e \%\n", 100. * (sum_b[2] - sum_a[2]) / sum_a[2]);
		//~ printf("relative difference for fitted weights: %.2e \%\n", 100. * (sum_b[3] - sum_a[3]) / sum_a[3]);
		//~ printf("absoluut difference with fitted weights: %.2e \%\n", (sum_b[3] - sum_a[3]));
		//~ exit(0);
	}
	util_safe_free(&x);
	util_safe_free(&xl);
	util_safe_free(&xu);	
	util_safe_free(&somepointers);
}

double util_fit_gammadistribution_parameters_costf(unsigned n, const double *x, double *grad, void *params)
{
	void 	**somepointers = (void**)params;
	void 	*pointer0		= somepointers[0];
	void 	*pointer1		= somepointers[1];
	void 	*pointer2		= somepointers[2];
	void 	*pointer3		= somepointers[3];
	void 	*pointer4		= somepointers[4];
	void 	*pointer5		= somepointers[5];
	
	double cost;
	int i;
	int ip;
	double sum_a[3];
	double sum_b[3];
	double N0, D0_mm, mu;
	double a_thislwc_gm3;
	double b_thislwc_gm3;
	double *griddep = NULL;
		
	t_zephyros_psd *fit_psd 	= (t_zephyros_psd*) pointer0;
	t_zephyros_psd *ref_psd 	= (t_zephyros_psd*) pointer1;
	int				model 		= (int) *((int*)pointer2);
	int				j 			= (int) *((int*)pointer3);

	double *nd, *nd2;
	double *xyzt = calloc(4, sizeof(double));
	
	double rcs_hh;
	double rcs_vv;
	double *dummy = NULL;
	
	nd 		= calloc(fit_psd->n_diameters, sizeof(double));
	nd2 	= calloc(ref_psd->n_diameters, sizeof(double));

	if (grad) {
		printf("Gradient not supoort when fitting a gamma distribution");
		fflush(stdout); exit(0);
	}
	

	//
	//1. Get the DSD parameters
	//
	
	if (model == 1) {
		//N0 will follow from liquid water content - D0 relation.
		D0_mm	= x[0]; //D0_mm
		mu		= 0.;
		
		//calculate number densities
		for (i = 0; i < fit_psd->n_diameters; i++ ) {
			nd[i] = 
				util_gamma_integral(
				1.,			//N0
				mu, 		//mu
				D0_mm, 		//D0_mm
				fit_psd->discrete_D_equiv_mm_interval_llimt[i],
				fit_psd->discrete_D_equiv_mm_interval_ulimt[i]
				);
		}
		
		//Calculate liquid water content for N0 = 1.
		a_thislwc_gm3 = 0.;
		for (i = 0; i < fit_psd->n_diameters; i++ ) {
			a_thislwc_gm3 +=
				nd[i] *	//# m-3
				1.e3 *			//g kg-1
				998.2071 * 		//[water density in kg m-3]
				(4./3.) * M_PI * pow((fit_psd->discrete_D_equiv_mm[i]/2.) * 1.e-3,3.); //m3 of single droplet
		}
		
		//D0_mm 	= 1.639 * pow(a_thislwc_gm3, 0.25);
		b_thislwc_gm3 	= pow(D0_mm / 1.639, 4.); 
		N0 = b_thislwc_gm3 / a_thislwc_gm3;		
	}
	

	if (model == 2) {
		N0 		= x[0];
		D0_mm	= x[1];
		mu = (6.084 * pow(D0_mm, 2.)) - (29.85 * D0_mm) + 34.64;				
	}

	//~ if (model == 3) {	
		//~ if ((weights1 == NULL) | (weights2 == NULL) | (weights3 == NULL)) {
			//~ printf("fit gamma distribution not correctly initialized...\n");
			//~ fflush(stdout); exit(0);
		//~ }

		//~ N0		= x[0];
		//~ D0_mm	= x[1];
		//~ mu		= x[2];
	//~ }

	//place fitted results back
	fit_psd->grid_gammadistribution_N0[j] 		= N0;
	fit_psd->grid_gammadistribution_D0_mm[j] 	= D0_mm;
	fit_psd->grid_gammadistribution_mu[j] 		= mu;

	//
	//2. cost function.
	//
	sum_a[0]	= 0.;
	sum_a[1]	= 0.;
	sum_b[0]	= 0.;
	sum_b[1]	= 0.;
	
	for (i = 0; i < fit_psd->n_diameters; i++ ) {
		nd[i] = util_gamma_integral(
			fit_psd->grid_gammadistribution_N0[j] ,			//N0
			fit_psd->grid_gammadistribution_mu[j], 			//mu
			fit_psd->grid_gammadistribution_D0_mm[j], 			//D0_mm
			fit_psd->discrete_D_equiv_mm_interval_llimt[i],
			fit_psd->discrete_D_equiv_mm_interval_ulimt[i]
			);				
		sum_a[0] += nd[i] * fit_psd->grid_rcs_hh[i][j];			
		sum_a[1] += nd[i] * fit_psd->grid_rcs_vv[i][j];			
	}

	griddep = calloc(ref_psd->field->n, sizeof(double));
	for (i = 0; i < ref_psd->n_diameters; i++ ) {
		xyzt[0] = fit_psd->field->x((void*)fit_psd->field, j);
		xyzt[1] = fit_psd->field->y((void*)fit_psd->field, j);
		xyzt[2] = fit_psd->field->z((void*)fit_psd->field, j);
		xyzt[3] = fit_psd->field->t((void*)fit_psd->field, j);
		
		//get number density
		interpolation_bilint(
			ref_psd->lut_ln_number_density_m3[i], 
			xyzt,
			&(nd2[i]),
			0,
			dummy);
		nd2[i] = exp(nd2[i]);
		
		//get cross section
		interpolation_bilint_griddep(
			ref_psd->lut_ln_number_density_m3[0],
			xyzt,
			griddep);
		
		rcs_hh = 0.;
		for (ip = 0; ip < ref_psd->field->n; ip++ ) {
			rcs_hh += griddep[ip] * ref_psd->grid_rcs_hh[i][ip];
		}
		
		rcs_vv = 0.;
		for (ip = 0; ip < ref_psd->field->n; ip++ ) {
			rcs_vv += griddep[ip] * ref_psd->grid_rcs_vv[i][ip];
		}

		sum_b[0] += ref_psd->grid_number_density_m3[i][j] * rcs_hh;
		sum_b[1] += ref_psd->grid_number_density_m3[i][j] * rcs_vv;
	}
	util_safe_free(&griddep);

	if (model == 1) {
		cost =  pow(sum_a[0] - sum_b[0], 2.);	
	}

	if (model == 2) {
		cost = 0.;
			
		cost += pow(sum_a[0] - sum_b[0], 2.);	
		cost += pow(sum_a[1] - sum_b[1], 2.);	
	}

	//~ if (model == 3) {
		//~ cost = 0.;
		
		//~ sum_a[0]	= 0.;
		//~ sum_b[0]	= 0.;
		//~ sum_a[1]	= 0.;
		//~ sum_b[1]	= 0.;
		//~ sum_a[2]	= 0.;
		//~ sum_b[2]	= 0.;
		//~ for (i = 0; i < thepsd->n_diameters; i++ ) {
			//~ nd[i] = util_gamma_integral(
				//~ N0,			//N0
				//~ mu, 		//mu
				//~ D0_mm, 		//D0_mm
				//~ thepsd->discrete_D_equiv_mm_interval_llimt[i],
				//~ thepsd->discrete_D_equiv_mm_interval_ulimt[i]
				//~ );

			//~ sum_a[0] += nd[i] * (weights1[i] / weights1_sum);
			//~ sum_b[0] += thepsd->grid_number_density_m3[i][j] * (weights1[i] / weights1_sum);			
			//~ sum_a[1] += nd[i] * (weights2[i] / weights2_sum);
			//~ sum_b[1] += thepsd->grid_number_density_m3[i][j] * (weights2[i] / weights2_sum);
			//~ sum_a[2] += nd[i] * (weights3[i] / weights3_sum);
			//~ sum_b[2] += thepsd->grid_number_density_m3[i][j] * (weights3[i] / weights3_sum);
		//~ }
		
		//~ cost += pow(sum_a[0] - sum_b[0], 2.);	
		//~ cost += pow(sum_a[1] - sum_b[1], 2.);	
		//~ cost += pow(sum_a[2] - sum_b[2], 2.);	
	//~ }
	

	//~ printf("N0 = %.2e\n", N0);
	//~ printf("D0_mm = %.2e\n", D0_mm);
	//~ printf("mu = %.2e\n", mu);

	util_safe_free(&xyzt);
	util_safe_free(&griddep);
	free(nd);
	free(nd2);
	return cost;
}

/*
 	int iflag;
	int n2;
	double result;
	
	n2 = n;
	if (grad) {
		iflag = 2;	//calculate grad
	} else {
		iflag = 0;
	}
	retrieval_fdvar_cost_function(params, &n2, (double*) x, &result, grad, &iflag);
	return result;
*/





//iot = integration over tracjorie.
void util_iot_initialize(t_zephyros_iot **piot)
{	
	t_zephyros_iot *iot = calloc(1, sizeof(t_zephyros_iot));
	
	iot->n 						= 20;
	iot->xyzt 					= NULL;
	iot->u						= NULL;
	iot->v						= NULL;
	iot->w						= NULL;

	iot->store_par			 	= 0;
	iot->wi_u 					= NULL;
	iot->wi_v 					= NULL;
	iot->wi_w 					= NULL;

	iot->prepared 				= 0;
	iot->initialized 			= 1;
	
	*piot = iot;
}

void util_iot_assert_initialized(t_zephyros_iot *iot) 
{
	if (iot == NULL) {
		printf("iot was NULL pointer. Exiting.\n"); 
		fflush(stdout); exit(0);
	}
	if (iot->initialized != 1) {
		printf("iot was not initialized. Exiting.\n"); 
		fflush(stdout); exit(0);
	}
}


void util_iot_assert_prepared(t_zephyros_iot *iot)
{
	util_iot_assert_initialized(iot);
	if (iot->prepared != 1) {
		printf("iot was not prepared. Exiting.\n"); 
		fflush(stdout); exit(0);
	}
}

void util_iot_free(t_zephyros_iot **piot)
{
	t_zephyros_iot *iot = *piot;
	int i;
	
	if (iot != NULL) {
		util_iot_assert_initialized(iot);
		if (iot->u != NULL) {free(iot->u);}
		if (iot->v != NULL) {free(iot->v);}
		if (iot->w != NULL) {free(iot->w);}
		
		if (iot->xyzt != NULL) {
			for ( i = 0; i < iot->n; i++ ) {
				free(iot->xyzt[i]);
			}
			free(iot->xyzt);
		}
		
		if (iot->wi_u != NULL) {free(iot->wi_u);}
		if (iot->wi_v != NULL) {free(iot->wi_v);}
		if (iot->wi_w != NULL) {free(iot->wi_w);}
			
		free(iot);
		*piot = NULL;
	}
}

void util_iot_prepare(t_zephyros_iot *iot)
{
	int i;
	util_iot_assert_initialized(iot);
	
	func_dbl_arr_calloc(iot->n, &iot->u);
	func_dbl_arr_calloc(iot->n, &iot->v);
	func_dbl_arr_calloc(iot->n, &iot->w);

	iot->xyzt = malloc(iot->n * sizeof(double*));
	for ( i = 0; i < iot->n; i++ ) {
		func_dbl_arr_calloc(4, iot->xyzt + i);
	}
	
	if (iot->store_par) {
		func_dbl_arr_calloc(iot->n, &iot->wi_u);
		func_dbl_arr_calloc(iot->n, &iot->wi_v);
		func_dbl_arr_calloc(iot->n, &iot->wi_w);
	}
		
	iot->prepared = 1;	
}

void util_iot_traj(t_zephyros_iot *iot, t_zephyros_windfield *windfield, t_zephyros_particles_widget *scat)
{			
	double traj_d_min;
	double *traj_d;

	double traj_dt;
	int i_traj;
	
	//1. define the trajectorie
	if (scat->particle_inertial_distance_xy < scat->particle_inertial_distance_z_vt_large) {
		traj_d_min 	= 0.25 * scat->particle_inertial_distance_xy;
	} else {
		traj_d_min 	= 0.25 * scat->particle_inertial_distance_z_vt_large;;		
	}

	func_dbl_arr_calloc(iot->n, &traj_d);
	for ( i_traj = 0; i_traj < iot->n; i_traj++ ) {
		traj_d[i_traj] = traj_d_min * pow(2., ((iot->n - 1.) - i_traj)/2.);
	}

	//obtain wind velocity at end point
	i_traj = iot->n - 1;
	util_windfield_fuvw(windfield, iot->xyzt[i_traj], 
		iot->u + i_traj, iot->v + i_traj, iot->w + i_traj,
		0,
		NULL, NULL, NULL,
		0,
		1.
		);
	iot->w[i_traj] -= scat->particle_terminal_fall_speed;
	
	for ( i_traj = iot->n - 2; i_traj >= 0; i_traj-- ) {
		traj_dt = traj_d[i_traj] / sqrt(pow(iot->u[i_traj+1], 2.) + pow(iot->v[i_traj+1], 2.) + pow(iot->w[i_traj+1], 2.));
		iot->xyzt[i_traj][0] = iot->xyzt[i_traj+1][0] - (traj_dt * iot->u[i_traj+1]);
		iot->xyzt[i_traj][1] = iot->xyzt[i_traj+1][1] - (traj_dt * iot->v[i_traj+1]);
		iot->xyzt[i_traj][2] = iot->xyzt[i_traj+1][2] - (traj_dt * iot->w[i_traj+1]);
		iot->xyzt[i_traj][3] = iot->xyzt[i_traj+1][3] - traj_dt;
		util_windfield_fuvw(windfield, iot->xyzt[i_traj], 
			iot->u + i_traj, iot->v + i_traj, iot->w + i_traj,
			0,
			NULL, NULL, NULL,
			0,
			1.
			);
		iot->w[i_traj] -= scat->particle_terminal_fall_speed;
	}
}

void util_iot_calc(t_zephyros_iot *iot, t_zephyros_particles_widget *scat)
{
	double delta_u;
	double delta_v;
	double delta_w;
	int i_traj;
	double tmpdu, tmpdv, tmpdw;
	double sgndu, sgndv, sgndw;
	double traj_dt;
	double tmpvar;
	
	double Ux, Uy, Uz;
	double Q1, Q2x, Q2y, Q3, Q4, Q5;
	
	//zero for the first point
	delta_u = 0.;
	delta_v = 0.;
	delta_w = 0.;

	if (iot->store_par) {
		iot->wi_u[0] = iot->u[0];
		iot->wi_v[0] = iot->v[0];
		iot->wi_w[0] = iot->w[0];		
	}
	
	//integrate to the solution
	for ( i_traj = 0; i_traj < (iot->n - 1); i_traj++ ) {
		//~ tmpdu = 
		//~ tmpdv = 
		//~ tmpdw = iot->w[i_traj + 1] - (delta_w + iot->w[i_traj]);
		//~ sgndu = tmpdu / fabs(tmpdu);
		//~ sgndv = tmpdv / fabs(tmpdv);
		//~ sgndw = tmpdw / fabs(tmpdw);
		//~ if (isnanorinf(&sgndu)) sgndu = 1.;								
		//~ if (isnanorinf(&sgndv)) sgndv = 1.;								
		//~ if (isnanorinf(&sgndw)) sgndw = 1.;	
		
		traj_dt = (iot->xyzt[i_traj + 1][3] - iot->xyzt[i_traj][3]);
		
		Ux = iot->u[i_traj + 1] - (delta_u + iot->u[i_traj]);
		Uy = iot->v[i_traj + 1] - (delta_v + iot->v[i_traj]);
		Uz = iot->w[i_traj + 1] - (delta_w + iot->w[i_traj]);

		Q1 	= (1. / (2. * scat->particle_inertial_eta_xy * traj_dt));
		Q2x = fabs(Ux / (scat->particle_inertial_eta_xy * traj_dt));
		Q2y = fabs(Uy / (scat->particle_inertial_eta_xy * traj_dt));
		
		Q3 	= scat->particle_terminal_fall_speed  + 
				 (1. / (2. * scat->particle_inertial_eta_z * traj_dt));

		Q4 = fabs(Uz / (scat->particle_inertial_eta_z * traj_dt));
		Q5 = Q4 - (2. * pow(scat->particle_terminal_fall_speed, 2.));

		if (Ux == 0.) {
			delta_u = 0.;
		} else if (Ux < 0.) {
			delta_u = (-1. * Q1) + sqrt(pow(Q1, 2.) + Q2x);
		} else { //Ux > 0
			delta_u = Q1 - sqrt(pow(Q1, 2.) + Q2x);
		}

		if (Uy == 0.) {
			delta_v = 0.;
		} else if (Uy < 0.) {
			delta_v = (-1. * Q1) + sqrt(pow(Q1, 2.) + Q2y);
		} else { //Uy > 0
			delta_v = Q1 - sqrt(pow(Q1, 2.) + Q2y);
		}
		
		if (Uz == 0) {
			delta_w = 0.;
		} else if (Uz < 0.) {
			delta_w = (-1. * Q3) + sqrt(pow(Q3, 2.) + Q4);
		} else if ((0 < Uz) & (Uz <= (scat->particle_terminal_fall_speed / 2.))) {
			delta_w = Q3 - sqrt(pow(Q3, 2.) + Q4);
		} else if (Uz > (scat->particle_terminal_fall_speed / 2.)) {
			delta_w = Q3 - sqrt(pow(Q3, 2.) + Q5);
		}
			
		if (iot->store_par) {
			iot->wi_u[i_traj + 1] = iot->u[i_traj + 1] + delta_u;
			iot->wi_v[i_traj + 1] = iot->v[i_traj + 1] + delta_v;
			iot->wi_w[i_traj + 1] = iot->w[i_traj + 1] + delta_w;
		}											
	}

	if (iot->store_par) {
		iot->wi_u[iot->n - 1] = iot->u[iot->n - 1] + delta_u;
		iot->wi_v[iot->n - 1] = iot->v[iot->n - 1] + delta_v;
		iot->wi_w[iot->n - 1] = iot->w[iot->n - 1] + delta_w;
	}	
	
	iot->res_wi_u = iot->u[iot->n - 1] + delta_u;
	iot->res_wi_v = iot->v[iot->n - 1] + delta_v;
	iot->res_wi_w = iot->w[iot->n - 1] + delta_w;
}

void util_iot_write(t_zephyros_iot *iot, char output_fname[8192])
{			
	FILE *fp;
	int i_traj;
	
	fp = fopen(output_fname, "w");
	for ( i_traj = 0; i_traj < iot->n; i_traj++ ) {
		fprintf(fp, "%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e%10.2e \n", 
			iot->xyzt[i_traj][0],
			iot->xyzt[i_traj][1],
			iot->xyzt[i_traj][2],
			iot->xyzt[i_traj][3],
			iot->u[i_traj],
			iot->v[i_traj],
			iot->w[i_traj],
			iot->wi_u[i_traj],
			iot->wi_v[i_traj],
			iot->wi_w[i_traj]
		); fflush(stdout);
	}
	fclose(fp);
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
	double tmp_a;
	double tmp_b;
	double tmp_c;

	tmp_a = specialfunctions_incompletegamma(mu + 1., ((3.67 + mu) / D0) * D2);
	tmp_b = specialfunctions_incompletegamma(mu + 1., ((3.67 + mu) / D0) * D1);
	tmp_c = pow((3.67 + mu) / D0 , mu+1.);
	
	tmp = N0 * ((tmp_a - tmp_b) / tmp_c);
	
	if (tmp == 0.) {
		return N0 * 1.e-50;
	} else {
		return tmp;
	}

}









void util_safe_free_(void **ptr)
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




