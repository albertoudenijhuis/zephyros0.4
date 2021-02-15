/*
Description: 
	4D Var retrieval of wind

Revision History:
	2014

Functions:
	4D Var retrieval of wind
	* 
Author:
	Albert Oude Nijhuis <albertoudenijhuis@gmail.com>

Institute:
	Delft University of Technology
	
Zephyros version:
	0.4

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

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#include <nlopt.h>

#include "retrieval_fdvar.h"
#include "radarfilter.h"
#include "coordinates.h"
#include "interpolation.h"
#include "util.h"
#include "func.h"

//uncomment next statement for debug mode
#define _ZEPHYROS_FDVAR_DEBUG


void retrieval_fdvar_apply(t_fdvar_opc *opc)
{
	prev_costfunction = 0.;
    t_zephyros_psd 					*cfg_psd[101];
	t_zephyros_field 				*cfg_field[101];
    t_zephyros_turbulence_widget 	*cfg_turbulence[101];
	int i, j;
	
	//do requested casts
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_cast_solutions\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_cast_solutions(opc);
	
	//store current post configuration pointers
	for ( i = 0; i <= 100; i++ ) {
		cfg_psd[i]			= opc->zc->retrieval->post_scattererfield->psd[i];
		cfg_field[i]		= opc->zc->retrieval->post_windfield->field[i];
		cfg_turbulence[i]	= opc->zc->retrieval->post_windfield->turbulence[i];
	}
	
	//update active elements
	for ( i = 0; i <= 100; i++ ) {
		opc->zc->retrieval->post_scattererfield->psd[i] = NULL;
		opc->zc->retrieval->post_windfield->turbulence[i] = NULL;
		opc->zc->retrieval->post_windfield->field[i] = NULL;

		for (j = 0; j < opc->c->n_active_scattererfield_nrs; j++) {
			if (i == opc->c->active_scattererfield_nrs[j])
				opc->zc->retrieval->post_scattererfield->psd[i] = cfg_psd[i];
		}
		for (j = 0; j < opc->c->n_active_windfield_grid_nrs; j++) {
			if (i == opc->c->active_windfield_grid_nrs[j])
				opc->zc->retrieval->post_windfield->field[i] = cfg_field[i];			
		}
		for (j = 0; j < opc->c->n_active_windfield_turbulence_nrs; j++) {
			if (i == opc->c->active_windfield_turbulence_nrs[j])
				opc->zc->retrieval->post_windfield->turbulence[i] = cfg_turbulence[i];
		}
	}

	//initialize p
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_initialize_p\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_initialize_p(opc);

	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_minimize_cost_function\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_minimize_cost_function(opc);

	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_estimate_posterior_errors\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_estimate_posterior_errors(opc);

	//restore post configuration pointers
	for ( i = 0; i <= 100; i++ ) {
		opc->zc->retrieval->post_scattererfield->psd[i] = cfg_psd[i];
		opc->zc->retrieval->post_windfield->field[i] = cfg_field[i];
		opc->zc->retrieval->post_windfield->turbulence[i] = cfg_turbulence[i];
	}

	//free p
	retrieval_fdvar_free_p(opc);
}



int retrieval_fdvar_minimize_cost_function(t_fdvar_opc *opc)
{
	clock_t start, end;
	double cpu_time_used;
	t_fdvar_o *o = opc->o;
	t_fdvar_p *p = opc->p;
	t_zephyros_config_retrieval_fdvar_cfg *c = opc->c;
	
	double *x;		
	double minf;	
	int inf;
	double *xu, *xl;
	void *vd_opc = (void*) opc;
	int i;
	
	nlopt_opt opt;
	
	start = clock();

	// Starting point
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_init_x\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_init_x(vd_opc, &x);

	func_dbl_arr_malloc(p->Kn, &xu);
	func_dbl_arr_malloc(p->Kn, &xl);
	
	for (i=0; i < p->Kn; i++) {
		xl[i]	= 0.;
		xu[i]	= 1.;
	}

	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("nlopt_create\n"); fflush(stdout);
	#endif 
	
	//best fit algorithm based on some tests.
	opt = nlopt_create(NLOPT_LD_VAR1, p->Kn);

	//global optimizers, seem not really to work...
	//NLOPT_GN_CRS2_LM
	//NLOPT_G_MLSL_LDS
	//NLOPT_GD_STOGO_RAND
	//NLOPT_GN_ISRES
	//NLOPT_GN_ESCH             ::: seems to work!!!!!
	
	//without derivatives!!
	//opt = nlopt_create(NLOPT_LN_COBYLA, p->Kn); 			//Minimum cost function found at 1,33399e+03
	//opt = nlopt_create(NLOPT_LN_BOBYQA, p->Kn); 			//Minimum cost function found at 5,85130e+02
	//opt = nlopt_create(NLOPT_LN_NEWUOA_BOUND, p->Kn); 	//Minimum cost function found at 7,61164e+02
	//opt = nlopt_create(NLOPT_LN_NELDERMEAD, p->Kn); 		//Minimum cost function found at 1,03584e+03
	//opt = nlopt_create(NLOPT_LN_SBPLX, p->Kn); 			//Minimum cost function found at 4,05836e+01	

	//with derivatives
	//opt = nlopt_create(NLOPT_LD_MMA, p->Kn); 						//Minimum cost function found at 2,44692e-01	<<--
	//opt = nlopt_create(NLOPT_LD_SLSQP, p->Kn); 					//Minimum cost function found at 6,40579e+01
	//opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND_RESTART, p->Kn); 	//Minimum cost function found at 1,05249e+02
	//opt = nlopt_create(NLOPT_LD_VAR1, p->Kn); 					//Minimum cost function found at 1,80310e-01	<<--
	//opt = nlopt_create(NLOPT_LD_VAR2, p->Kn); 					//Minimum cost function found at 1,80323e-01	<<--

	//opt = nlopt_create(NLOPT_LN_PRAXIS, p->Kn); // algorithm and dimensionality			//Minimum cost function found at 1,44567e+01
	
	nlopt_set_min_objective(opt, retrieval_fdvar_cost_function_nlopt, vd_opc);
	nlopt_set_lower_bounds(opt, xl);
	nlopt_set_upper_bounds(opt, xu);
	
	//according to manual: nlopt_set_xtol_rel(opt, 1e-4);
	nlopt_set_xtol_rel(opt, 1.e-4);
	nlopt_set_maxtime(opt, c->maximum_time_s);

	if ((inf = nlopt_optimize(opt, x, &minf)) < 0) {
		printf("nlopt failed!, i = %i\n", inf);
	}
	else {
		printf("Minimum cost function found at %.5e\n", minf);
	}

	//unpack cost final function
	retrieval_fdvar_unpack_x(opc,x);

	nlopt_destroy(opt);
	
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Run time: %f s \n", cpu_time_used);

	free(xu);
	free(xl);
	free(x);
}

double retrieval_fdvar_cost_function_nlopt(unsigned n, const double *x, double *grad, void *params)
{
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
}



void retrieval_fdvar_cost_function(
	void 	*vd_opc,			//optional arguments
	int		*n,				//number of parameters
	double 	*x,				//parameters
	double 	*f,				//function value
	double 	*fd,			//function derivatives
	int 	*iflag			
	)
{
    // if iflag = 0 : calculate the function value, return in f
    // if iflag = 1 : calculate the derivatives, return in fd
    // if iflag = 2 : calculate both

	int calc_costfunction				= ((*iflag == 0) | (*iflag == 2));
	int calc_costfunction_derivatives	= ((*iflag == 1) | (*iflag == 2));

	double f_dBZ_hh;
	double *fd_dBZ_hh = NULL;
	double f_dBZdr;
	double *fd_dBZdr = NULL;
	double f_dBLdr;
	double *fd_dBLdr = NULL;
	double f_Doppler_velocity_hh_ms;
	double *fd_Doppler_velocity_hh_ms = NULL;
	double f_Doppler_spectral_width_hh_ms;
	double *fd_Doppler_spectral_width_hh_ms = NULL;
	
	double weight_dBZ_hh;
	double weight_dBZdr;
	double weight_dBLdr;
	double weight_Doppler_velocity_hh_ms;
	double weight_Doppler_spectral_width_hh_ms;

	double *prior_grid_hspeed, *prior_grid_hdir;
	double *post_grid_hspeed, *post_grid_hdir;
	double *tmp_der_hspeed, *tmp_der_hdir;
	
	double alpha0, alpha, *griddep;
	double *delta, *Sinv_delta, tmp;
	double mytmp, mytmp2, myfactor;
	double wf_u, wf_v, wf_w;
	
	int io, ip, in, i_psd, i_par;
	int i_int;	//i for spectral intervals
	int i_w;
	
	t_fdvar_opc										*opc = ((t_fdvar_opc*)vd_opc);
	t_fdvar_o 										*o = opc->o;
	t_fdvar_p 										*p = opc->p;
	t_zephyros_config_retrieval_fdvar_cfg 			*c = opc->c;
    t_zephyros_config 								*zc = opc->zc;
	t_radarfilter_todolist							*todo;

	double eta_hh, eta_hv, eta_vv;
	t_radarfilter_res_vol 			*res_vol;
	
	//Unpack K
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_unpack_x\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_unpack_x(opc,x);
	
	//set values to zero
	if (calc_costfunction) {
		*f	= 0.;	

		f_dBZ_hh = 0;
		f_dBZdr = 0;
		f_dBLdr = 0;
		f_Doppler_velocity_hh_ms = 0;
		f_Doppler_spectral_width_hh_ms = 0;	
	}
	if (calc_costfunction_derivatives) {
		fd_dBZ_hh = malloc(*n * sizeof(double));
		fd_dBZdr = malloc(*n * sizeof(double));
		fd_dBLdr = malloc(*n * sizeof(double));
		fd_Doppler_velocity_hh_ms = malloc(*n * sizeof(double));
		fd_Doppler_spectral_width_hh_ms = malloc(*n * sizeof(double));
		for (in=0; in < *n; in++) 
		{
			fd[in] =0.;
			fd_dBZ_hh[in] =0.;
			fd_dBZdr[in] =0.;
			fd_dBLdr[in] =0.;
			fd_Doppler_velocity_hh_ms[in] =0.;
			fd_Doppler_spectral_width_hh_ms[in] =0.;			
		}
	}
		
	//initialize
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("radarfilter_initialize_todolist\n"); fflush(stdout);
	#endif
	radarfilter_initialize_todolist(zc, 1, &todo);

	//TBD, i_psd ??
	todo->calc_eta_i_hh = (	(zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) &
							(c->costfunction_dBZ_hh | c->costfunction_dBZdr) & 
							(calc_costfunction_derivatives)
							);
	todo->calc_eta_i_hv = (	(zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) &
							(c->costfunction_dBLdr) & 
							(calc_costfunction_derivatives)
							);
	todo->calc_eta_i_vv = (	(zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) &
							(c->costfunction_dBLdr | c->costfunction_dBZdr) & 
							(calc_costfunction_derivatives)
							);

	todo->calc_spectrum_eta_i_hh 	= (	(zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) &
										(c->costfunction_Doppler_spectrum_dBZ_hh | c->costfunction_specific_dBZdr) & 
										(calc_costfunction_derivatives)
										);
	todo->calc_spectrum_eta_i_hv 	= (	(zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) &
										(c->costfunction_specific_dBLdr) & 
										(calc_costfunction_derivatives)
										);
	todo->calc_spectrum_eta_i_vv 	= (	(zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) &
										(c->costfunction_specific_dBLdr | c->costfunction_specific_dBZdr) & 
										(calc_costfunction_derivatives)
										);
	
	todo->der_edr13 =  (calc_costfunction_derivatives);
	/*
	(	(zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) &
										(c->costfunction_) & 
										(calc_costfunction_derivatives)
										);
	*/
	todo->der_dBZ_hh			= (calc_costfunction_derivatives);
	 							
	//update todo list with logic
	radarfilter_prepare_todolist(zc, 1, todo);

	radarfilter_initialize_resolution_volume(zc, 1, &res_vol, todo);
	radarfilter_exec(zc, 1, todo, o->n, o->model_radarmeasurement, res_vol);
	radarfilter_free_resolution_volume(zc, 1, &res_vol, todo);
	//TBD: reuse res_vol (contains e.g. scattering t-matrix results)
    

    
    
    //***
    //Observation space
    //***
    
	//loop through observations.
	#ifdef _ZEPHYROS_FDVAR_DEBUG
	printf("loop through observations\n"); fflush(stdout);
	#endif
	for (io=0; io < o->n; io++) 
	{
		//observation space of cost function
		if (calc_costfunction) {			
			if (c->costfunction_dBZ_hh) f_dBZ_hh += pow((o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) / o->radarmeasurement[io]->dBZ_hh_err, 2.);
			if (c->costfunction_dBZdr) f_dBZdr += pow((o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) / o->radarmeasurement[io]->dBZdr_err, 2.);
			if (c->costfunction_dBLdr) f_dBLdr += pow((o->model_radarmeasurement[io]->dBLdr - o->radarmeasurement[io]->dBLdr) / o->radarmeasurement[io]->dBLdr_err, 2.);
			if (c->costfunction_Doppler_velocity_hh_ms) f_Doppler_velocity_hh_ms += pow((o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) / o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, 2.);
			if (c->costfunction_Doppler_spectral_width_hh_ms) f_Doppler_spectral_width_hh_ms += pow((o->model_radarmeasurement[io]->Doppler_spectral_width_hh_ms - o->radarmeasurement[io]->Doppler_spectral_width_hh_ms) / o->radarmeasurement[io]->Doppler_spectral_width_hh_ms_err, 2.);

			
			//spectra
			for ( i_int = 0; i_int < o->radarmeasurement[io]->n_spectrum; i_int++ ) {
				//TBD
				if (c->costfunction_Doppler_spectrum_dBZ_hh) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int] - o->radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int]) / (o->radarmeasurement[io]->Doppler_spectrum_dBZ_hh_err[i_int]), 2);
				if (c->costfunction_specific_dBZdr) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->specific_dBZdr[i_int] - o->radarmeasurement[io]->specific_dBZdr[i_int]) / (o->radarmeasurement[io]->specific_dBZdr_err[i_int]), 2);
				if (c->costfunction_specific_dBLdr) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->specific_dBLdr[i_int] - o->radarmeasurement[io]->specific_dBLdr[i_int]) / (o->radarmeasurement[io]->specific_dBLdr_err[i_int]), 2);

			}




			//TBD: all other variables that are in the cost function
		}


		
		if (calc_costfunction_derivatives) {
			//observation space of jacobian

			//loop over scatterers
			for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
				if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {
					//dBlwc
					if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
						griddep = malloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_scattererfield->psd[i_psd]->lut_ln_number_density_m3[0],
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							griddep);

						//dBZ_hh
						if (c->costfunction_dBZ_hh) {
							alpha = 2. * (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) * pow(o->radarmeasurement[io]->dBZ_hh_err , -2.);
							//TBD: find out why we get nans here ... 
							if (isnanorinf(&alpha) == 0) {
								for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
									fd_dBZ_hh[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
										 alpha * griddep[ip];
								}
							} else {
								printf("nan!\n");
								printf("alpha = %.2e\n", alpha);
								printf("pow(o->radarmeasurement[io]->dBZ_hh_err , -2.) = %.2e\n", pow(o->radarmeasurement[io]->dBZ_hh_err , -2.));
							}
						}

						//Doppler_spectrum_dBZ_hh
						//TBD
						for ( i_int = 0; i_int < o->radarmeasurement[io]->n_spectrum; i_int++ ) {
							if (c->costfunction_Doppler_spectrum_dBZ_hh) {
								alpha = (1./ p->cost_no) * 2. * (o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int] - o->radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int]) / pow(o->radarmeasurement[io]->Doppler_spectrum_dBZ_hh_err[i_int], 2.);
								for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
									fd[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
										 alpha * griddep[ip];
								}
							}
						}
						//TBD: all other variables that depend on dBlwc
						
						free(griddep);
					}
					
					//dBN
					if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) {
						griddep = malloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_scattererfield->psd[i_psd]->lut_ln_number_density_m3[0],
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							griddep);

						//dBZ_hh
						if (c->costfunction_dBZ_hh) {
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) * pow(o->radarmeasurement[io]->dBZ_hh_err, -2.);
							eta_hh = func_dB_inv(o->model_radarmeasurement[io]->dBZ_hh) / o->model_radarmeasurement[io]->coef_eta2Z;

							for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
								alpha = alpha0 * (o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] / eta_hh);								
								if (isnanorinf(&alpha) == 0) {
									for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
										fd_dBZ_hh[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
											+ (zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
											+ ip] += 
											 alpha * griddep[ip];
									}
								} else {
									printf("nan!\n");
								}
							}
						}
						
						//dBZdr
						if (c->costfunction_dBZdr) {
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) * pow(o->radarmeasurement[io]->dBZdr_err, -2.);
							eta_hh = func_dB_inv(o->model_radarmeasurement[io]->dBZ_hh) / o->model_radarmeasurement[io]->coef_eta2Z;
							eta_vv = func_dB_inv(o->model_radarmeasurement[io]->dBZ_vv) / o->model_radarmeasurement[io]->coef_eta2Z;

							for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
								alpha = alpha0 * ((o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] / eta_hh) - (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par] / eta_vv));								

								if (isnanorinf(&alpha) == 0) {
									for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
										fd_dBZdr[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
											+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
											+ip] += 
											 alpha * griddep[ip];
									}
								} else {
									printf("nan!\n");
								}
							}
						}
						//dBLdr
						if (c->costfunction_dBLdr) {
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBLdr - o->radarmeasurement[io]->dBLdr) * pow(o->radarmeasurement[io]->dBLdr_err, -2.);
							eta_hv = func_dB_inv(o->model_radarmeasurement[io]->dBZ_hv) / o->model_radarmeasurement[io]->coef_eta2Z;
							eta_vv = func_dB_inv(o->model_radarmeasurement[io]->dBZ_vv) / o->model_radarmeasurement[io]->coef_eta2Z;

							for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
								alpha = alpha0 * ((o->model_radarmeasurement[io]->eta_i_hv[i_psd][i_par] / eta_hv) - (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par] / eta_vv));								
								
								if (isnanorinf(&alpha) == 0) {
									for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
										fd_dBLdr[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
											+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
											+ip] += 
											 alpha * griddep[ip];
									}
								} else {
									printf("nan!\n");
								}
							}
						}
						
						//spectral variables
						//TBD
						for ( i_int = 0; i_int < o->radarmeasurement[io]->n_spectrum; i_int++ ) {
							if (c->costfunction_Doppler_spectrum_dBZ_hh) {							
								alpha0 = (1./ p->cost_no) * 2. * (o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int] - o->radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int]) / pow(o->radarmeasurement[io]->Doppler_spectrum_dBZ_hh_err[i_int], 2.);
								for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
									mytmp 		= func_dB_inv(o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int]) / o->model_radarmeasurement[io]->coef_eta2Z;
									alpha = alpha0 * (o->model_radarmeasurement[io]->spectrum_eta_i_hh[i_psd][i_par][i_int] / mytmp);								
									if (isnanorinf(&alpha) == 0) {
										for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
											fd[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
												+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
												+ip] += alpha * griddep[ip];
										}
									}
								}
							}
							
							if (c->costfunction_specific_dBZdr) {							
								alpha0 = (1./ p->cost_no) * 2. * (o->model_radarmeasurement[io]->specific_dBZdr[i_int] - o->radarmeasurement[io]->specific_dBZdr[i_int]) / pow(o->radarmeasurement[io]->specific_dBZdr_err[i_int], 2.);
								for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
									mytmp 		= func_dB_inv(o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int]) / o->model_radarmeasurement[io]->coef_eta2Z;
									mytmp2 		= func_dB_inv(o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_vv[i_int]) / o->model_radarmeasurement[io]->coef_eta2Z;
									
									if ((mytmp != 0.) & (mytmp2 != 0.)) {
										alpha = alpha0 * ((o->model_radarmeasurement[io]->spectrum_eta_i_hh[i_psd][i_par][i_int] / mytmp) - (o->model_radarmeasurement[io]->spectrum_eta_i_vv[i_psd][i_par][i_int] / mytmp2));								
										
										if (isnanorinf(&alpha) == 0) {
											for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
												fd[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
													+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
													+ip] += alpha * griddep[ip];
											}
										}
									}
								}
							}

							if (c->costfunction_specific_dBLdr) {								
								alpha0 = (1./ p->cost_no) * 2. * (o->model_radarmeasurement[io]->specific_dBLdr[i_int] - o->radarmeasurement[io]->specific_dBLdr[i_int]) / pow(o->radarmeasurement[io]->specific_dBLdr_err[i_int], 2.);
								for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
									mytmp 		= func_dB_inv(o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_hv[i_int]) / o->model_radarmeasurement[io]->coef_eta2Z;
									mytmp2 		= func_dB_inv(o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_vv[i_int]) / o->model_radarmeasurement[io]->coef_eta2Z;
									
									if ((mytmp != 0.) & (mytmp2 != 0.)) {
										alpha = alpha0 * ((o->model_radarmeasurement[io]->spectrum_eta_i_hv[i_psd][i_par][i_int] / mytmp) - (o->model_radarmeasurement[io]->spectrum_eta_i_vv[i_psd][i_par][i_int] / mytmp2));								

										if (isnanorinf(&alpha) == 0) {
											for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
												fd[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
													+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
													+ip] += alpha * griddep[ip];
											}
										}
									}
								}
							}
							
						}
						
						//TBD: other variables

						free(griddep);
					}	
				}
			}
		
			//obtain windfield uvw at central location
			util_windfield_fuvw(zc->retrieval->post_windfield, o->model_radarmeasurement[io]->center_coor->enu_xyzt, 
				&wf_u, &wf_v, &wf_w,
				0,
				NULL, NULL, NULL,
				zc->retrieval->radarfilter->geostrophic_advection,
				zc->retrieval->radarfilter->coriolis_parameter
				);
				
				
			//loop over wind fields
			for ( i_w = 0; i_w <= 100; i_w++ ) {
				if (zc->retrieval->post_windfield->field[i_w] != NULL) {								
					if (zc->retrieval->prior_windfield->fit_u[i_w]) {
						griddep = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_windfield->lut_u[i_w],
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							griddep);

						//grid_u
						if (c->costfunction_dBZ_hh) {
							//assert that todo->der_dBZ_hh is set. 
							alpha = 2. * (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) * pow(o->radarmeasurement[io]->dBZ_hh_err, -2.);							
							alpha = alpha * o->model_radarmeasurement[io]->der_dBZ_hh[0] * (-1. * wf_u);
							for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
								fd_dBZ_hh[zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] += 
									 alpha * griddep[ip];
							}
						}

						//grid_u
						if (c->costfunction_Doppler_velocity_hh_ms) {
							alpha = 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) * pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);
							alpha *= o->model_radarmeasurement[io]->center_coor->radar_enu_dir[0];
							//TBD
							//it would be better to use reflectivity weighted coordinates
							for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
								fd_Doppler_velocity_hh_ms[
									zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] += 
									 alpha * griddep[ip];
							}
						}
						
						free(griddep);
					}
					
					if (zc->retrieval->prior_windfield->fit_v[i_w]) {
						griddep = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_windfield->lut_v[i_w],
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							griddep);

						if (c->costfunction_dBZ_hh) {
							//assert that todo->der_dBZ_hh is set. 
							alpha = 2. * (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) * pow(o->radarmeasurement[io]->dBZ_hh_err, -2.);							
							alpha = alpha * o->model_radarmeasurement[io]->der_dBZ_hh[1] * (-1. * wf_v);
							for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
								fd_dBZ_hh[zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] += 
									 alpha * griddep[ip];
							}
						}
						
						//grid_v
						if (c->costfunction_Doppler_velocity_hh_ms) {
							alpha = 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) * pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);
							alpha *= o->model_radarmeasurement[io]->center_coor->radar_enu_dir[1];
							//TBD
							//it would be better to use reflectivity weighted coordinates
							for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
								fd_Doppler_velocity_hh_ms[
									zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] += 
									 alpha * griddep[ip];
							}
						}
						
						free(griddep);
					}
					
					if (zc->retrieval->prior_windfield->fit_w[i_w]) {
						griddep = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_windfield->lut_w[i_w],
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							griddep);

						if (c->costfunction_dBZ_hh) {
							//assert that todo->der_dBZ_hh is set. 
							alpha = 2. * (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) * pow(o->radarmeasurement[io]->dBZ_hh_err, -2.);							
							alpha = alpha * o->model_radarmeasurement[io]->der_dBZ_hh[2] * (-1. * wf_w);
							for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
								fd_dBZ_hh[zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] += 
									 alpha * griddep[ip];
							}
						}
						
						//grid_w
						if (c->costfunction_Doppler_velocity_hh_ms) {
							alpha = 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) * pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);
							alpha *= o->model_radarmeasurement[io]->center_coor->radar_enu_dir[2];
							//TBD
							//it would be better to use reflectivity weighted coordinates
							for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
								fd_Doppler_velocity_hh_ms[
									zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] += 
									 alpha * griddep[ip];
							}
						}
						
						free(griddep);
					}
				}
			}
			
			//loop over turbulence fields
			for ( i_w = 0; i_w <= 100; i_w++ ) {
				if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
					if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
						griddep = malloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n * sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_windfield->turbulence[i_w]->lut_edr13,
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							griddep);

						//grid_edr13
						if (c->costfunction_Doppler_spectral_width_hh_ms) {
							alpha = 2. * (o->model_radarmeasurement[io]->Doppler_spectral_width_hh_ms - o->radarmeasurement[io]->Doppler_spectral_width_hh_ms) * pow(o->radarmeasurement[io]->Doppler_spectral_width_hh_ms_err, -2.);
							alpha *= o->model_radarmeasurement[io]->der_edr13_Doppler_spectral_width_hh_ms;
							
							for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
								fd_Doppler_spectral_width_hh_ms[
									zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
									 alpha * griddep[ip];
							}
						}
						
						if (c->costfunction_dBZdr) {
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) * pow(o->radarmeasurement[io]->dBZdr_err, -2.);
							alpha = alpha0 * o->model_radarmeasurement[io]->der_edr13_dBZdr;
							
							for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
								fd_dBZdr[
									zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
									 alpha * griddep[ip];
							}

						}
						
						if (c->costfunction_dBLdr) {
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBLdr - o->radarmeasurement[io]->dBLdr) * pow(o->radarmeasurement[io]->dBLdr_err, -2.);
							alpha = alpha0 * o->model_radarmeasurement[io]->der_edr13_dBLdr;

							for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
								fd_dBLdr[
									zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
									 alpha * griddep[ip];
							}
						}

						
						//TBD, other variables that depends on edr13
						free(griddep);						
						
					}
				}
			}


			
			
			//~ if (c->fit_refl) {
				//~ for (ip = 0; ip < p->field->n; ip++) {


	

				
					//~ fd[p->fit_Knr_refl + ip] += (1./p->on) * alpharelf * res_vol->integrated_scatterer_griddep[ip];
				//~ }
				//~ 
				//~ //add derivatives related to advection
				//~ if (c->apply_advection) {
					//~ if (c->fit_hspeed)			drefl_dhspeed	= (-1. * res_vol->integrated_ddBZ_du * sin(res_vol->integrated_hdir)) + (-1. * res_vol->integrated_ddBZ_dv * cos(res_vol->integrated_hdir));
					//~ if (c->fit_hdir)			drefl_dhdir		= (-1. * res_vol->integrated_ddBZ_du * res_vol->integrated_hspeed * cos(res_vol->integrated_hdir)) + (res_vol->integrated_ddBZ_dv * res_vol->integrated_hspeed * sin(res_vol->integrated_hdir));
					//~ if (c->fit_w)				drefl_dw		= res_vol->integrated_ddBZ_dw;
					//~ 
					//~ for (ip = 0; ip < p->field->n; ip++) {
						//~ if (c->fit_hspeed)		fd[p->fit_Knr_hspeed + ip] 	+= (1./p->on) * drefl_dhspeed * alpharelf * res_vol->integrated_scatterer_griddep[ip];
						//~ if (c->fit_hdir)		fd[p->fit_Knr_hdir + ip] 	+= (1./p->on) * drefl_dhdir * alpharelf * res_vol->integrated_scatterer_griddep[ip];
						//~ if (c->fit_w)			fd[p->fit_Knr_w + ip] 		+= (1./p->on) * drefl_dw * alpharelf * res_vol->integrated_scatterer_griddep[ip];
					//~ }
				//~ }
			//~ }
			//~ 
			//~ if (c->fit_vr) {
				//~ if (c->fit_hspeed)			dPvrK_dhspeed 	= 	(-1. * sin(res_vol->integrated_hdir) * res_vol->integrated_beamdir_enu_dx) - (cos(res_vol->integrated_hdir) * res_vol->integrated_beamdir_enu_dy);
				//~ if (c->fit_hdir)			dPvrK_dhdir 	= 	(-1. * res_vol->integrated_hspeed * cos(res_vol->integrated_hdir) * res_vol->integrated_beamdir_enu_dx) + (res_vol->integrated_hspeed * sin(res_vol->integrated_hdir) * res_vol->integrated_beamdir_enu_dy);
//~ 
				//~ if (c->fit_hspeed)			alphahspeed	= dPvrK_dhspeed	 		* 2. * (res_vol->integrated_radial_vel_ms   - o->vr[io] ) 	/ pow(o->svr[io], 2.);
				//~ if (c->fit_hdir)			alphahdir	= dPvrK_dhdir			* 2. * (res_vol->integrated_radial_vel_ms   - o->vr[io] ) 	/ pow(o->svr[io], 2.);
				//~ if (c->fit_w)				alphaw		= res_vol->integrated_beamdir_enu_dz	* 2. * (res_vol->integrated_radial_vel_ms   - o->vr[io] ) 	/ pow(o->svr[io], 2.);
				//~ 
				//~ for (ip = 0; ip < p->field->n; ip++) {
					//~ if (c->fit_hspeed)		fd[p->fit_Knr_hspeed + ip] 	+= (1./p->on) * alphahspeed * res_vol->integrated_wind_griddep[ip];
					//~ if (c->fit_hdir)		fd[p->fit_Knr_hdir + ip] 	+= (1./p->on) * alphahdir * res_vol->integrated_wind_griddep[ip];
					//~ if (c->fit_w)			fd[p->fit_Knr_w + ip] 		+= (1./p->on) * alphaw * res_vol->integrated_wind_griddep[ip];
				//~ }
				//~ 
				//~ //add derivatives related to advection
				//~ if (c->apply_advection) {
					//~ if (c->fit_hspeed)			dradial_vel_ms_dhspeed	= (-1. * res_vol->integrated_dradial_vel_ms_du * sin(res_vol->integrated_hdir)) + (-1. * res_vol->integrated_dradial_vel_ms_dv * cos(res_vol->integrated_hdir));
					//~ if (c->fit_hdir)			dradial_vel_ms_dhdir		= (-1. * res_vol->integrated_dradial_vel_ms_du * res_vol->integrated_hspeed * cos(res_vol->integrated_hdir)) + (res_vol->integrated_dradial_vel_ms_dv * res_vol->integrated_hspeed * sin(res_vol->integrated_hdir));
					//~ if (c->fit_w)				dradial_vel_ms_dw		= res_vol->integrated_dradial_vel_ms_dw;
					//~ 
					//~ for (ip = 0; ip < p->field->n; ip++) {
						//~ if (c->fit_hspeed)		fd[p->fit_Knr_hspeed + ip] 	+= (1./p->on) * dradial_vel_ms_dhspeed * alphahspeed * res_vol->integrated_wind_griddep[ip];
						//~ if (c->fit_hdir)		fd[p->fit_Knr_hdir + ip] 	+= (1./p->on) * dradial_vel_ms_dhdir * alphahdir * res_vol->integrated_wind_griddep[ip];
						//~ if (c->fit_w)			fd[p->fit_Knr_w + ip] 		+= (1./p->on) * dradial_vel_ms_dw * alphaw * res_vol->integrated_wind_griddep[ip];
					//~ }
				//~ }
			//~ }
			
			//calculation of parameter posterior errors, update with observation dependence
			//should be updated
			//observation space of jacobian
			//~ for (ip = 0; ip < p->field->n; ip++) {
				//if ((c->fit_refl)	& (fabs(res_vol->integrated_scatterer_griddep[ip]) > 0.))	p->post_srefl[ip]   += pow(o->svr[io] * fabs(res_vol->integrated_scatterer_griddep[ip]), -2);
				//~ if ((c->fit_hspeed)	& (fabs(res_vol->integrated_wind_griddep[ip]) > 0.))	p->post_shspeed[ip] += pow(o->svr[io] * fabs(res_vol->integrated_wind_griddep[ip]), -2);
				//~ if ((c->fit_hdir)	& (fabs(res_vol->integrated_wind_griddep[ip]) > 0.))	p->post_shdir[ip]  	+= pow(o->svr[io] * fabs(res_vol->integrated_wind_griddep[ip]), -2);
				//~ if ((c->fit_w)		& (fabs(res_vol->integrated_wind_griddep[ip]) > 0.))	p->post_sw[ip] 		+= pow(o->svr[io] * fabs(res_vol->integrated_wind_griddep[ip]), -2);
			//~ }
		}


		//copy resolution volume integration parameters to parameter container.
		//~ o->windvector_u[io]	= res_vol->integrated_u;
		//~ o->windvector_v[io]	= res_vol->integrated_v;
		//~ o->windvector_w[io]	= res_vol->integrated_w;
		//~ 
		//~ if (calc_costfunction_derivatives) {
			//~ if (c->fit_hspeed | c->fit_hdir | c->fit_w) {
				//~ o->windvector_su[io] = 1. / (3. * fabs(res_vol->integrated_beamdir_enu_dx)); //calculation continues...
				//~ o->windvector_sv[io] = 1. / (3. * fabs(res_vol->integrated_beamdir_enu_dy)); //calculation continues...
				//~ o->windvector_sw[io] = 1. / (3. * fabs(res_vol->integrated_beamdir_enu_dz)); //calculation continues...				
			//~ }
		//~ }
	}

    //End of observation space
    //***




	//#ifdef _ZEPHYROS_FDVAR_DEBUG
	//if (calc_costfunction) {printf("current cost value: %.5e (%.2e%)\n"   , *f, 100.*(*f - prev_costfunction)/prev_costfunction);fflush(stdout);}	
	//#endif





    //***
    //Parameter space
    //***
	#ifdef _ZEPHYROS_FDVAR_DEBUG
	printf("Parameter space\n"); fflush(stdout);
	#endif
	
	//loop over scatterers
	for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {		
			//dBlwc
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				delta = malloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * sizeof(double));
				Sinv_delta = malloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
					delta[ip] = 
						func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]) -
							func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]);
				}

				retrieval_fdvar_cs_qrsol(zc->retrieval->prior_scattererfield->psd[i_psd]->dBlwc_ecm, delta, Sinv_delta);							
				
				dbldotprod(&zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_dBZ_hh += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {						
						fd_dBZ_hh[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
							2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);	
			}

			//dBN: TBD
		}
	}

	//loop over windfields
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->field[i_w] != NULL) {
			//hspeed, hdir
			if (zc->retrieval->prior_windfield->fit_u[i_w] & zc->retrieval->prior_windfield->fit_v[i_w] &
					(zc->retrieval->prior_windfield->use_hspeed_hdir_erorrs[i_w] == 1)
				) {

				delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				Sinv_delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				
				//hspeed
				//hdir
				prior_grid_hspeed 	= malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				prior_grid_hdir		= malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				post_grid_hspeed	= malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				post_grid_hdir		= malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				
				tmp_der_hspeed		= malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				tmp_der_hdir		= malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				

				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					prior_grid_hspeed[ip] =
						sqrt(pow(zc->retrieval->prior_windfield->grid_u[i_w][ip],2.) + pow(zc->retrieval->prior_windfield->grid_v[i_w][ip],2.));
					post_grid_hspeed[ip] =
						sqrt(pow(zc->retrieval->post_windfield->grid_u[i_w][ip],2.) + pow(zc->retrieval->post_windfield->grid_v[i_w][ip],2.));
					
					prior_grid_hdir[ip] =
						atan2(-1. * zc->retrieval->prior_windfield->grid_u[i_w][ip], -1. * zc->retrieval->prior_windfield->grid_v[i_w][ip]);
					post_grid_hdir[ip] =
						atan2(-1. * zc->retrieval->post_windfield->grid_u[i_w][ip], -1. * zc->retrieval->post_windfield->grid_v[i_w][ip]);
				}
				

				//***
				//hspeed
				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						post_grid_hspeed[ip]
						-
						prior_grid_hspeed[ip];						
				}

				retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->hspeed_ecm[i_w], delta, Sinv_delta);							
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_Doppler_velocity_hh_ms += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
						tmp_der_hspeed[ip] = 				
							p->Kn * 2. * Sinv_delta[ip];
					}
				}
				
				//***
				//hdir			
				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						angleAminB(post_grid_hdir + ip,
							prior_grid_hdir + ip);
				}

				retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->hdir_ecm[i_w], delta, Sinv_delta);							
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_Doppler_velocity_hh_ms += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {				
						tmp_der_hdir[ip] = 				
						  2. * Sinv_delta[ip];
					}
				}
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {				
						fd_Doppler_velocity_hh_ms[zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] += 
							(
							(tmp_der_hspeed[ip] * (zc->retrieval->post_windfield->grid_u[i_w][ip] / post_grid_hspeed[ip])) +
							(tmp_der_hdir[ip] * (zc->retrieval->post_windfield->grid_v[i_w][ip] / pow(post_grid_hspeed[ip], 2.)))
							);

						fd_Doppler_velocity_hh_ms[zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] += 
							(
							(tmp_der_hspeed[ip] * (zc->retrieval->post_windfield->grid_v[i_w][ip] / post_grid_hspeed[ip])) +
							(tmp_der_hdir[ip] * (-1. * zc->retrieval->post_windfield->grid_u[i_w][ip] / pow(post_grid_hspeed[ip], 2.)) )
							);
					}
				}
				
				free(delta);
				free(Sinv_delta);	

				free(prior_grid_hspeed);
				free(prior_grid_hdir);
				free(post_grid_hspeed);
				free(post_grid_hdir);

				free(tmp_der_hspeed);
				free(tmp_der_hdir);
			}
					
			//grid_u
			if (zc->retrieval->prior_windfield->fit_u[i_w] &
				(zc->retrieval->prior_windfield->use_hspeed_hdir_erorrs[i_w] == 0)
				) {
				delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				Sinv_delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						zc->retrieval->post_windfield->grid_u[i_w][ip]
						-
						zc->retrieval->prior_windfield->grid_u[i_w][ip];						
				}

				retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->u_ecm[i_w], delta, Sinv_delta);							
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_Doppler_velocity_hh_ms += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
						fd_Doppler_velocity_hh_ms[zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] += 
							 2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);	
			}
			
			//grid_v
			if (zc->retrieval->prior_windfield->fit_v[i_w] &
				(zc->retrieval->prior_windfield->use_hspeed_hdir_erorrs[i_w] == 0)
				) {
				delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				Sinv_delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						zc->retrieval->post_windfield->grid_v[i_w][ip]
						-
						zc->retrieval->prior_windfield->grid_v[i_w][ip];						
				}

				retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->v_ecm[i_w], delta, Sinv_delta);							
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_Doppler_velocity_hh_ms += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
						fd_Doppler_velocity_hh_ms[zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] += 
							2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);	
			}

			//grid_w
			if (zc->retrieval->prior_windfield->fit_w[i_w]) {
				delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
				Sinv_delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						zc->retrieval->post_windfield->grid_w[i_w][ip]
						-
						zc->retrieval->prior_windfield->grid_w[i_w][ip];						
				}

				retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->w_ecm[i_w], delta, Sinv_delta);							
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_Doppler_velocity_hh_ms += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
						fd_Doppler_velocity_hh_ms[zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] += 
							2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);	
			}
		}
	}
	

	//turbulence, grid_edr13	
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
			if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
				delta = malloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n * sizeof(double));
				Sinv_delta = malloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n * sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
					delta[ip] = 
						zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[ip]
						-
						zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13[ip];						
				}

				retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->turbulence[i_w]->edr13_ecm, delta, Sinv_delta);							
				dbldotprod(&zc->retrieval->prior_windfield->turbulence[i_w]->field->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_Doppler_spectral_width_hh_ms += tmp;

				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {					
						fd_Doppler_spectral_width_hh_ms[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
							2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);
			}
		}
	}
	
    //End of parameter space
    //***
	
	
	
	
	
	
	
	


	//calculate weights
	weight_dBZ_hh = 0.;
	weight_dBZdr = 0.;
	weight_dBLdr = 0.;
	weight_Doppler_velocity_hh_ms = 0.;
	weight_Doppler_spectral_width_hh_ms = 0.;

	//observation space
	for (io=0; io < o->n; io++) weight_dBZ_hh += pow(o->radarmeasurement[io]->dBZ_hh_err, -2.);
	for (io=0; io < o->n; io++) weight_dBZdr += pow(o->radarmeasurement[io]->dBZdr_err, -2.);
	for (io=0; io < o->n; io++) weight_dBLdr += pow(o->radarmeasurement[io]->dBLdr_err, -2.);
	for (io=0; io < o->n; io++) weight_Doppler_velocity_hh_ms += pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);
	for (io=0; io < o->n; io++) weight_Doppler_spectral_width_hh_ms += pow(o->radarmeasurement[io]->Doppler_spectral_width_hh_ms_err, -2.);

	//parameter space
	//loop over scatterers
	for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {		
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) 
					weight_dBZ_hh += pow(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3[ip], -2.);
			}
		}
	}

	//loop over windfields
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->field[i_w] != NULL) {
			if (zc->retrieval->prior_windfield->fit_u[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) 
					weight_Doppler_velocity_hh_ms += pow(zc->retrieval->prior_windfield->grid_u_err[i_w][ip], -2.);		
			}
			if (zc->retrieval->prior_windfield->fit_v[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) 
					weight_Doppler_velocity_hh_ms += pow(zc->retrieval->prior_windfield->grid_v_err[i_w][ip], -2.);		
			}
			if (zc->retrieval->prior_windfield->fit_w[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) 
					weight_Doppler_velocity_hh_ms += pow(zc->retrieval->prior_windfield->grid_w_err[i_w][ip], -2.);		
			}
		}
	}
	//loop over turbulences
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
			if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {				
				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) 
					weight_Doppler_spectral_width_hh_ms += pow(zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13[ip], -2.);
			}
		}
	}

	//and take inverse
	weight_dBZ_hh = pow(weight_dBZ_hh, -1.);
	weight_dBZdr = pow(weight_dBZdr, -1.);
	weight_dBLdr = pow(weight_dBLdr, -1.);
	weight_Doppler_velocity_hh_ms = pow(weight_Doppler_velocity_hh_ms, -1.);
	weight_Doppler_spectral_width_hh_ms = pow(weight_Doppler_spectral_width_hh_ms, -1.);

	if (calc_costfunction) {
		*f += weight_dBZ_hh * f_dBZ_hh;
		*f += weight_dBZdr * f_dBZdr;
		*f += weight_dBLdr * f_dBLdr;
		*f += weight_Doppler_velocity_hh_ms * f_Doppler_velocity_hh_ms;
		*f += weight_Doppler_spectral_width_hh_ms * f_Doppler_spectral_width_hh_ms;
	}

	if (calc_costfunction_derivatives) {
		for (in=0; in < *n; in++) {
			fd[in] += weight_dBZ_hh * fd_dBZ_hh[in];
			fd[in] += weight_dBZdr * fd_dBZdr[in];
			fd[in] += weight_dBLdr * fd_dBLdr[in];
			fd[in] += weight_Doppler_velocity_hh_ms * fd_Doppler_velocity_hh_ms[in];
			fd[in] += weight_Doppler_spectral_width_hh_ms * fd_Doppler_spectral_width_hh_ms[in];
		}
	}
	
	//account for derivative in cost function to x
	//Hence dK/dx has to be added. dK/dx = (u - l)
	if (calc_costfunction_derivatives) {
		for (in=0; in < *n; in++) {
			fd[in] 	*= (p->Kubound[in] - p->Klbound[in]);
		}
	}
	
	if (calc_costfunction) {printf("current cost value: %.5e (%.2e%)\n"   , *f, 100.*(*f - prev_costfunction)/prev_costfunction);fflush(stdout);}
	#ifdef _ZEPHYROS_FDVAR_DEBUG
	if (calc_costfunction_derivatives) {
		printf("current cost derivatives:\n");
		printf("%i derivatives:\n", *n);
		for (in=0; in < *n; in++) {
			if (fd[in] != 0.) {
				printf("fd[%i] = %.5e \n", in, fd[in]); fflush(stdout);
			}
		}
	}
	#endif
	
	prev_costfunction = *f;

	radarfilter_free_todolist(&todo);
	
	util_safe_free(&fd_dBZ_hh);
	util_safe_free(&fd_dBZdr);
	util_safe_free(&fd_dBLdr);
	util_safe_free(&fd_Doppler_velocity_hh_ms);
	util_safe_free(&fd_Doppler_spectral_width_hh_ms);
}

double retrieval_fdvar_K2x(t_fdvar_p *p, double *K, int Knr)
{
	//rescaling
	return (*K - p->Klbound[Knr]) / (p->Kubound[Knr] - p->Klbound[Knr]);
}

double retrieval_fdvar_x2K(t_fdvar_p *p, double *x, int Knr)
{
	//rescaling
	return p->Klbound[Knr] + (p->Kubound[Knr] - p->Klbound[Knr]) * *x;
}

void retrieval_fdvar_init_x(
	void *vd_opc,
	double **x)	//x = K
{
	int i_psd, i_par;
	int ip, in;
	t_fdvar_o *o = ((t_fdvar_opc*)vd_opc)->o;
	t_fdvar_p *p = ((t_fdvar_opc*)vd_opc)->p;
	t_zephyros_config_retrieval_fdvar_cfg *c = ((t_fdvar_opc*)vd_opc)->c;
	t_zephyros_config *zc = ((t_fdvar_opc*)vd_opc)->zc;
			
	double factor = 4.;
	double tmp;
	double(*K2x)(t_fdvar_p *p, double *x, int Knr);
	K2x = retrieval_fdvar_K2x; 
	int i_w;
	
	func_dbl_arr_malloc(p->Kn, x);

	//set Klbound, Kubound and x value

	//walk through psds
	for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {
			//dBlwc
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
					p->Klbound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] =
						func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]) -
						(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3[ip]);						
					p->Kubound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] =
						func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]) +
						(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3[ip]);
					if (p->Klbound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] < -100.)
						p->Klbound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] = -100.;
					if (p->Kubound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] > 100.)
						p->Kubound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] = 100.;
					tmp = func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]);
					(*x)[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr	+ip] =
						K2x(p, &tmp, zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip);
				}				
			}
			//dBN
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) {
				for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
					for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
						p->Klbound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
									+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									+ip] =
							func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip]) -
							(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBnumber_density_err_m3[i_par][ip]);
						if (p->Klbound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
										+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ip] < -100.)
							p->Klbound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
										+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ip] = -100.;
								
						p->Kubound[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr
									+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									+ip] =
							func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip]) +
							(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBnumber_density_err_m3[i_par][ip]);
						if (p->Kubound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
										+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ip] > 100.)
							p->Kubound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
										+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ip] = 100.;
							
							
						tmp = func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip]);
						(*x)[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr
								+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
								+ip] =
							K2x(p, &tmp, zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr +
							+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par) + ip);
					}
				}
			}
		}
	}
		 
	//walk through wind fields
	for ( i_w = 0; i_w <= 100; i_w++ ) {		
		if (zc->retrieval->post_windfield->field[i_w] != NULL) {

			//grid_u
			if (zc->retrieval->prior_windfield->fit_u[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					p->Klbound[zc->retrieval->prior_windfield->fit_u_Knr[i_w] + ip] =
						zc->retrieval->prior_windfield->grid_u[i_w][ip] -
						(factor * zc->retrieval->prior_windfield->grid_u_err[i_w][ip]);
					p->Kubound[zc->retrieval->prior_windfield->fit_u_Knr[i_w] + ip] =
						zc->retrieval->prior_windfield->grid_u[i_w][ip] +
						(factor * zc->retrieval->prior_windfield->grid_u_err[i_w][ip]);				
					(*x)[zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] =
						K2x(p, &(zc->retrieval->prior_windfield->grid_u[i_w][ip]), zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip);
				}
			}
			
			//grid_v
			if (zc->retrieval->prior_windfield->fit_v[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					p->Klbound[zc->retrieval->prior_windfield->fit_v_Knr[i_w] + ip] =
						zc->retrieval->prior_windfield->grid_v[i_w][ip] -
						(factor * zc->retrieval->prior_windfield->grid_v_err[i_w][ip]);
					p->Kubound[zc->retrieval->prior_windfield->fit_v_Knr[i_w] + ip] =
						zc->retrieval->prior_windfield->grid_v[i_w][ip] +
						(factor * zc->retrieval->prior_windfield->grid_v_err[i_w][ip]);				
					(*x)[zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] =
						K2x(p, &(zc->retrieval->prior_windfield->grid_v[i_w][ip]), zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip);
				}
			}

			//grid_w
			if (zc->retrieval->prior_windfield->fit_w[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					p->Klbound[zc->retrieval->prior_windfield->fit_w_Knr[i_w] + ip] =
						zc->retrieval->prior_windfield->grid_w[i_w][ip] -
						(factor * zc->retrieval->prior_windfield->grid_w_err[i_w][ip]);
					p->Kubound[zc->retrieval->prior_windfield->fit_w_Knr[i_w] + ip] =
						zc->retrieval->prior_windfield->grid_w[i_w][ip] +
						(factor * zc->retrieval->prior_windfield->grid_w_err[i_w][ip]);				
					(*x)[zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] =
						K2x(p, &(zc->retrieval->prior_windfield->grid_w[i_w][ip]), zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip);
				}
			}	
		}
	}
	
	//turbulence, grid_edr13	
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
			if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
				for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
					p->Klbound[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr + ip] =
						zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13[ip] -
						(factor * zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13_err[ip]);
					if (p->Klbound[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr + ip] <= 0)
						p->Klbound[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr + ip] = 
							1.e-10 * zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13_err[ip];
						
					p->Kubound[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr + ip] =
						zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13[ip] +
						(factor * zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13_err[ip]);				
					(*x)[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] =
						K2x(p, &(zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13[ip]), zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip);
				}		
			}
		}
	}	

	//check bounds
	for (in=0; in < p->Kn; in++) {	
		(*x)[in] = fmax((*x)[in], K2x(p, p->Klbound + in, in) );
		(*x)[in] = fmin((*x)[in], K2x(p, p->Kubound + in, in) );
	}
}

void retrieval_fdvar_unpack_x(t_fdvar_opc *opc, double *x)
{
	int i_psd, i_par;
	int ip;
	int i, j;
	int i_w;
	t_fdvar_p *p = opc->p;
	t_zephyros_config_retrieval_fdvar_cfg *c = opc->c;
	t_zephyros_config *zc = opc->zc;

	double(*x2K)(t_fdvar_p *p, double *K, int Knr); 
	x2K = retrieval_fdvar_x2K;
	int in;
	
	//walk through scatterers
	for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {
			//dBlwc
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
					zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[ip] =
						func_dB_inv(x2K(p, x + zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip, zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip));
 				}

 				//update number densities
				for (i = 0; i < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i++ ) {
					for (j = 0; j < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; j++ ) {		
						zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i][j] =
						(zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[j] / zc->retrieval->prior_scattererfield->psd[i_psd]->grid_lwc_gm3[j])
						*  zc->retrieval->prior_scattererfield->psd[i_psd]->grid_number_density_m3[i][j];
				}}
				
 				//update look-up table
				for (i = 0; i < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i++ ) {
					//free lut
					interpolation_free_lut(zc->retrieval->post_scattererfield->psd[i_psd]->lut_ln_number_density_m3 + i);
					util_field2lut(zc->retrieval->post_scattererfield->psd[i_psd]->field, zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i], 3, zc->retrieval->post_scattererfield->psd[i_psd]->lut_ln_number_density_m3 + i);
				}
			}
			//dBN
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) {
				//assert grid_number_density_m3 is set
				if (zc->retrieval->prior_scattererfield->psd[i_psd]->grid_number_density_m3 == NULL) {
					zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3 = malloc(zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters * sizeof(double*));
					for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ )
						zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par] = malloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * sizeof(double));
				}
				
 				//update number densities
				for (i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
					for (ip = 0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++ ) {		
						zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip] =
						func_dB_inv(x2K(p, 
									x + zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr
									+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									+ip,
									zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
									+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									+ip));
				}}
				
 				//update look-up table
				for (i = 0; i < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i++ ) {
					//free lut
					interpolation_free_lut(zc->retrieval->post_scattererfield->psd[i_psd]->lut_ln_number_density_m3 + i);
					util_field2lut(zc->retrieval->post_scattererfield->psd[i_psd]->field, zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i], 3, zc->retrieval->post_scattererfield->psd[i_psd]->lut_ln_number_density_m3 + i);
				}
			}
		}
	}
	
	//walk trough wind fields
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->field[i_w] != NULL) {
			//grid_u
			if (zc->retrieval->prior_windfield->fit_u[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					zc->retrieval->post_windfield->grid_u[i_w][ip] =
					x2K(p, 	x + zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip,
							zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip);
				}			
			}
			
			//grid_v
			if (zc->retrieval->prior_windfield->fit_v[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					zc->retrieval->post_windfield->grid_v[i_w][ip] =
					x2K(p, 	x + zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip,
							zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip);
				}			
			}

			//grid_w
			if (zc->retrieval->prior_windfield->fit_w[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					zc->retrieval->post_windfield->grid_w[i_w][ip] =
						x2K(p, 	x + zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip,
									zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip);
				}
			}
		}
	}
	

	//turbulence, grid_edr13	
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
			if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
				for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
					zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[ip] =
						x2K(p, 	x + zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip,
									zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip);
					
					//update edr
					zc->retrieval->post_windfield->turbulence[i_w]->grid_edr[ip] =
						pow(zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[ip], 3.);
				}
			}
		}
	}
}
	




void retrieval_fdvar_estimate_posterior_errors(t_fdvar_opc *opc)
{
	t_fdvar_o 										*o = opc->o;
	t_fdvar_p 										*p = opc->p;
	t_zephyros_config_retrieval_fdvar_cfg 			*c = opc->c;
    t_zephyros_config 								*zc = opc->zc;
	
	int io, ip, in, i_psd, i_par;
	int i_int;	//i for spectral intervals
	int i_w;
	double myval;
	double myerr;
	double *griddep;
	double weight, weightsum1, weightsum2;
	double model_dBZ_hh_err, model_Doppler_velocity_hh_ms_err;
	double tmp_hspeed;
	
    cs *mat1_triplet; 
    cs *mat1_triplet_a; 
    cs *mat1_triplet_b; 

	//estimate modelling errors
	model_dBZ_hh_err 					= 0.;
	model_Doppler_velocity_hh_ms_err 	= 0.;
	weightsum1 = 0.; weightsum2 = 0.;
	for (io=0; io < o->n; io++) 
	{	
		weight = pow(o->radarmeasurement[io]->dBZ_hh_err, -2.);
		model_dBZ_hh_err 					+= weight *
			fmax(0.,
				pow(o->model_radarmeasurement[io]->dBZ_hh  -o->radarmeasurement[io]->dBZ_hh, 2.)
				-
				pow(o->radarmeasurement[io]->dBZ_hh_err , 2.) );
		weightsum1 += weight;
		
		weight = pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);		
		model_Doppler_velocity_hh_ms_err 	+= weight *
			fmax(0.,
				pow(o->model_radarmeasurement[io]->Doppler_velocity_hh_ms -o->radarmeasurement[io]->Doppler_velocity_hh_ms, 2.)
				-
				pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, 2.) );
		weightsum2 += weight;
	}
	
	model_dBZ_hh_err = sqrt(model_dBZ_hh_err / weightsum1);
	model_Doppler_velocity_hh_ms_err = sqrt(model_Doppler_velocity_hh_ms_err / weightsum2);
	if (isnan(model_dBZ_hh_err)) model_dBZ_hh_err = 0.;
	if (isnan(model_Doppler_velocity_hh_ms_err)) model_Doppler_velocity_hh_ms_err = 0.;

	//loop over scatterers
	for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {
			//dBlwc
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				griddep = malloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * sizeof(double));

				mat1_triplet 		= cs_spalloc (o->n, zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, 1, 1, 1) ;
				for (io=0; io < o->n; io++) 
				{
					interpolation_bilint_griddep(
						zc->retrieval->prior_scattererfield->psd[i_psd]->lut_ln_number_density_m3[0],
						o->model_radarmeasurement[io]->advected_center_enu_xyzt,
						griddep);

					for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
						if (griddep[ip] != 0.) {
							myerr = sqrt(pow(o->radarmeasurement[io]->dBZ_hh_err, 2.) + pow(model_dBZ_hh_err, 2.));
							myval = griddep[ip] / myerr;
							if (!cs_entry (mat1_triplet, io, ip, myval)) {
								printf("Error with matrix allocation");
								exit(1);
							}
						}
					}
				}

				//helper function to solve it
				if (zc->retrieval->post_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3 == NULL) {
					zc->retrieval->post_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3 = 
						calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, sizeof(double));
				}
				retrieval_fdvar_estimate_posterior_errors_solver1(mat1_triplet,
					zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3,
					zc->retrieval->post_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3
				);
				
				//TBD: all other variables that depend on dBlwc
				free(griddep);
			}
			
			//dBN
			//TBD


			//dBZdr
			//TBD
		}
	}


	//loop over wind fields
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->field[i_w] != NULL) {

			//u
			if (zc->retrieval->prior_windfield->fit_u[i_w]) {
				griddep = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

				mat1_triplet 		= cs_spalloc (o->n, zc->retrieval->prior_windfield->field[i_w]->n, 1, 1, 1) ;
				for (io=0; io < o->n; io++) 
				{
					interpolation_bilint_griddep(
						zc->retrieval->prior_windfield->lut_u[i_w],
						o->model_radarmeasurement[io]->advected_center_enu_xyzt,
						griddep);
				
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
						if (griddep[ip] != 0.) {
							myerr = sqrt(pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, 2.) + pow(model_Doppler_velocity_hh_ms_err, 2.));
							myval = griddep[ip] *
									 o->model_radarmeasurement[io]->center_coor->radar_enu_dir[0]
									/ myerr;
							if (!cs_entry (mat1_triplet, io, ip, myval)) {
								printf("Error with matrix allocation");
								exit(1);
							}
						}
					}
				}

				if (zc->retrieval->post_windfield->grid_u_err[i_w] == NULL) {
					zc->retrieval->post_windfield->grid_u_err[i_w] = 
						calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				}
				retrieval_fdvar_estimate_posterior_errors_solver1(mat1_triplet,
					zc->retrieval->prior_windfield->grid_u_err[i_w],
					zc->retrieval->post_windfield->grid_u_err[i_w]
				);
				
				free(griddep);
			}
			//v
			if (zc->retrieval->prior_windfield->fit_v[i_w]) {
				griddep = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

				mat1_triplet 		= cs_spalloc (o->n, zc->retrieval->prior_windfield->field[i_w]->n, 1, 1, 1) ;
				for (io=0; io < o->n; io++) 
				{
					interpolation_bilint_griddep(
						zc->retrieval->prior_windfield->lut_v[i_w],
						o->model_radarmeasurement[io]->advected_center_enu_xyzt,
						griddep);
				
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
						if (griddep[ip] != 0.) {
							myerr = sqrt(pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, 2.) + pow(model_Doppler_velocity_hh_ms_err, 2.));
							myval = griddep[ip] *
									 o->model_radarmeasurement[io]->center_coor->radar_enu_dir[1]
									/ myerr;
							if (!cs_entry (mat1_triplet, io, ip, myval)) {
								printf("Error with matrix allocation");
								exit(1);
							}
						}
					}
				}

				if (zc->retrieval->post_windfield->grid_v_err[i_w] == NULL) {
					zc->retrieval->post_windfield->grid_v_err[i_w] = 
						calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				}
				retrieval_fdvar_estimate_posterior_errors_solver1(mat1_triplet,
					zc->retrieval->prior_windfield->grid_v_err[i_w],
					zc->retrieval->post_windfield->grid_v_err[i_w]
				);
				
				free(griddep);
			}
			//hdir, hspeed
			if (
				zc->retrieval->prior_windfield->fit_u[i_w] &
				zc->retrieval->prior_windfield->fit_v[i_w] &
				(zc->retrieval->prior_windfield->use_hspeed_hdir_erorrs[i_w] == 1)
				) {
				if (zc->retrieval->post_windfield->grid_hspeed_err[i_w] == NULL)
					zc->retrieval->post_windfield->grid_hspeed_err[i_w] = 
						calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				if (zc->retrieval->post_windfield->grid_hdir_err[i_w] == NULL)
					zc->retrieval->post_windfield->grid_hdir_err[i_w] = 
						calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					tmp_hspeed = sqrt(pow(zc->retrieval->post_windfield->grid_u[i_w][ip],2.) + pow(zc->retrieval->post_windfield->grid_v[i_w][ip],2.));				
					zc->retrieval->post_windfield->grid_hspeed_err[i_w][ip] =
					sqrt(
						pow( (zc->retrieval->post_windfield->grid_u[i_w][ip] / tmp_hspeed) * zc->retrieval->post_windfield->grid_u_err[i_w][ip] ,2.) +
						pow( (zc->retrieval->post_windfield->grid_v[i_w][ip] / tmp_hspeed) * zc->retrieval->post_windfield->grid_v_err[i_w][ip] ,2.)
					);
					zc->retrieval->post_windfield->grid_hdir_err[i_w][ip] =
					sqrt(
						pow( (zc->retrieval->post_windfield->grid_v[i_w][ip] / pow(tmp_hspeed, 2.)) * zc->retrieval->post_windfield->grid_u_err[i_w][ip] ,2.) +
						pow( (zc->retrieval->post_windfield->grid_u[i_w][ip] / pow(tmp_hspeed, 2.)) * zc->retrieval->post_windfield->grid_v_err[i_w][ip] ,2.)
					);
				}
			}
			//w
			if (zc->retrieval->prior_windfield->fit_w[i_w]) {
				griddep = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

				mat1_triplet 		= cs_spalloc (o->n, zc->retrieval->prior_windfield->field[i_w]->n, 1, 1, 1) ;
				for (io=0; io < o->n; io++) 
				{
					interpolation_bilint_griddep(
						zc->retrieval->prior_windfield->lut_w[i_w],
						o->model_radarmeasurement[io]->advected_center_enu_xyzt,
						griddep);
				
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
						if (griddep[ip] != 0.) {
							//w
							myerr = sqrt(pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, 2.) + pow(model_Doppler_velocity_hh_ms_err, 2.));
							myval = griddep[ip] *
									 o->model_radarmeasurement[io]->center_coor->radar_enu_dir[2]
									/ myerr;
							if (!cs_entry (mat1_triplet, io, ip, myval)) {
								printf("Error with matrix allocation");
								exit(1);
							}
						}
					}
				}

				if (zc->retrieval->post_windfield->grid_w_err[i_w] == NULL) {
					zc->retrieval->post_windfield->grid_w_err[i_w] = 
						calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				}
				retrieval_fdvar_estimate_posterior_errors_solver1(mat1_triplet,
					zc->retrieval->prior_windfield->grid_w_err[i_w],
					zc->retrieval->post_windfield->grid_w_err[i_w]
				);
				
				free(griddep);
			}
		}
	}

	//TBD: turbulence fields

	
	
	
	
}


void retrieval_fdvar_estimate_posterior_errors_solver1(cs *mat1_triplet, double *prior_errors, double *result_errors)
{
	int m;
	int n;
	int ip;
	
	cs *mat1;
	cs *mat1_T;
	cs *mat2;
	cs *mat3_triplet;
	cs *mat3;
	cs *mat4;
    css *mat4_S;
    csn *mat4_N; 
    
    m = mat1_triplet->m;	//number of observations
    n = mat1_triplet->n;	//number of parameters
    
    //compress matrix
	mat1 = cs_compress (mat1_triplet) ;           //m x n 
	cs_spfree (mat1_triplet) ; 

	mat1_T = cs_transpose (mat1, 1) ;          //n x m
	mat2 = cs_multiply (mat1_T, mat1) ;			//n x n
	cs_spfree(mat1);
	cs_spfree(mat1_T);
	
	mat3_triplet = cs_spalloc (n, n, 1, 1, 1) ;
	for (ip=0; ip < n; ip++) {
		if (!cs_entry (mat3_triplet, ip, ip, pow(prior_errors[ip], -2.))) {
			printf("Error with matrix allocation"); exit(1);
		}
	}
	mat3 = cs_compress (mat3_triplet) ;           //n x n 
	cs_spfree (mat3_triplet) ; 

	mat4 = cs_add (mat2, mat3, 1, 1) ;
	cs_spfree (mat2) ; 
	cs_spfree (mat3) ; 

	for (ip=0; ip < n; ip++) {
		result_errors[ip] = 1.;
	}
	cs_qrsol (1, mat4, result_errors);
	cs_spfree (mat4) ; 

	for (ip=0; ip < n; ip++) {
		result_errors[ip] = sqrt(fabs(result_errors[ip]));
	}	
}





void retrieval_fdvar_cast_solutions(t_fdvar_opc *opc)
{
	t_fdvar_o *o = opc->o;
	t_fdvar_p *p = opc->p;
	t_zephyros_config_retrieval_fdvar_cfg *c = opc->c;
    t_zephyros_config *zc = opc->zc;

	double *xyzt = malloc(4 * sizeof(double));
	double *dummy;
	int i_cast;
	int i_cast_from;
	int i_cast_to;
	int ip;
	double tmpH;
	int i_w;
	
	//wind fields
	for (i_cast = 1; i_cast < c->n_cast_windfield_grid_nrs; i_cast = i_cast + 2) {
		i_cast_from 	= c->cast_windfield_grid_nrs[i_cast - 1];
		i_cast_to 		= c->cast_windfield_grid_nrs[i_cast];
	
		for (ip=0; ip < zc->retrieval->prior_windfield->field[i_cast_to]->n; ip++) {
			xyzt[0] = zc->retrieval->prior_windfield->field[i_cast_to]->x((void*)zc->retrieval->prior_windfield->field[i_cast_to], ip);
			xyzt[1] = zc->retrieval->prior_windfield->field[i_cast_to]->y((void*)zc->retrieval->prior_windfield->field[i_cast_to], ip);
			xyzt[2] = zc->retrieval->prior_windfield->field[i_cast_to]->z((void*)zc->retrieval->prior_windfield->field[i_cast_to], ip);
			xyzt[3] = zc->retrieval->prior_windfield->field[i_cast_to]->t((void*)zc->retrieval->prior_windfield->field[i_cast_to], ip);

			if (zc->retrieval->post_windfield->lut_u != NULL)
				interpolation_bilint(zc->retrieval->post_windfield->lut_u[i_cast_from],
					xyzt,
					zc->retrieval->prior_windfield->grid_u[i_cast_to] + ip,
					0,
					dummy);
			if (zc->retrieval->post_windfield->lut_v != NULL)
				interpolation_bilint(zc->retrieval->post_windfield->lut_v[i_cast_from],
					xyzt,
					zc->retrieval->prior_windfield->grid_v[i_cast_to] + ip,
					0,
					dummy);
			if (zc->retrieval->post_windfield->lut_w != NULL)
				interpolation_bilint(zc->retrieval->post_windfield->lut_w[i_cast_from],
					xyzt,
					zc->retrieval->prior_windfield->grid_w[i_cast_to] + ip,
					0,
					dummy);
			
			/*
			//errors
			if (
				zc->retrieval->post_windfield->fit_u[i_cast_from]
				& zc->retrieval->post_windfield->fit_v[i_cast_from]
				& zc->retrieval->post_windfield->use_hspeed_hdir_erorrs[i_cast_from]) {

				//hspeed, hdir
				if (zc->retrieval->post_windfield->grid_hspeed_err[i_cast_from] != NULL) {
					//make lut if necessary
					if (zc->retrieval->post_windfield->lut_hspeed_err[i_cast_from] == NULL)
						util_field2lut(
							zc->retrieval->post_windfield->field[i_cast_from]
							, zc->retrieval->post_windfield->grid_hspeed_err[i_cast_from], 0, 
							zc->retrieval->post_windfield->lut_hspeed_err + i_cast_from
							);
					interpolation_bilint(zc->retrieval->post_windfield->lut_hspeed_err[i_cast_from],
						xyzt,
						zc->retrieval->prior_windfield->grid_hspeed_err[i_cast_to] + ip,
						0,
						dummy);
				}
				if (zc->retrieval->post_windfield->grid_hdir_err[i_cast_from] != NULL) {
					//make lut if necessary
					if (zc->retrieval->post_windfield->lut_hdir_err[i_cast_from] == NULL)
						util_field2lut(
							zc->retrieval->post_windfield->field[i_cast_from]
							, zc->retrieval->post_windfield->grid_hdir_err[i_cast_from], 0, 
							zc->retrieval->post_windfield->lut_hdir_err + i_cast_from
							);
					interpolation_bilint(zc->retrieval->post_windfield->lut_hdir_err[i_cast_from],
						xyzt,
						zc->retrieval->prior_windfield->grid_hdir_err[i_cast_to] + ip,
						0,
						dummy);
				}
			} else {
				if (zc->retrieval->post_windfield->grid_u_err[i_cast_from] != NULL) {
					//make lut if necessary
					if (zc->retrieval->post_windfield->lut_u_err[i_cast_from] == NULL)
						util_field2lut(
							zc->retrieval->post_windfield->field[i_cast_from]
							, zc->retrieval->post_windfield->grid_u_err[i_cast_from], 0, 
							zc->retrieval->post_windfield->lut_u_err + i_cast_from
							);
					interpolation_bilint(zc->retrieval->post_windfield->lut_u_err[i_cast_from],
						xyzt,
						zc->retrieval->prior_windfield->grid_u_err[i_cast_to] + ip,
						0,
						dummy);
				}	
				if (zc->retrieval->post_windfield->grid_v_err[i_cast_from] != NULL) {
					//make lut if necessary
					if (zc->retrieval->post_windfield->lut_v_err[i_cast_from] == NULL)
						util_field2lut(
							zc->retrieval->post_windfield->field[i_cast_from]
							, zc->retrieval->post_windfield->grid_v_err[i_cast_from], 0, 
							zc->retrieval->post_windfield->lut_v_err + i_cast_from
							);						
					interpolation_bilint(zc->retrieval->post_windfield->lut_v_err[i_cast_from],
						xyzt,
						zc->retrieval->prior_windfield->grid_v_err[i_cast_to] + ip,
						0,
						dummy);
				}
			}
			if (zc->retrieval->post_windfield->grid_w_err[i_cast_from] != NULL) {
				//make lut if necessary
				if (zc->retrieval->post_windfield->lut_w_err[i_cast_from] == NULL)
					util_field2lut(
						zc->retrieval->post_windfield->field[i_cast_from]
						, zc->retrieval->post_windfield->grid_w_err[i_cast_from], 0, 
						zc->retrieval->post_windfield->lut_w_err + i_cast_from
						);
				interpolation_bilint(zc->retrieval->post_windfield->lut_w_err[i_cast_from],
					xyzt,
					zc->retrieval->prior_windfield->grid_w_err[i_cast_to] + ip,
					0,
					dummy);
			}
			*/
			
		}
		printf("windfield casted from %i to %i\n", i_cast_from, i_cast_to);		
	}
	free(xyzt);
	
	
	//update wind field errors
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->prior_windfield->field[i_w] != NULL) {
			if (	zc->retrieval->prior_windfield->fit_u[i_w] &
					zc->retrieval->prior_windfield->fit_v[i_w] &
					zc->retrieval->prior_windfield->use_hspeed_hdir_erorrs[i_w] &
					((c->update_windfield_hspeed_err > 0.) | (c->update_windfield_hdir_err > 0.))
			) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {				
					if (c->update_windfield_hspeed_err > 0.)
						zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip] = c->update_windfield_hspeed_err;
					
					if (c->update_windfield_hdir_err > 0.)
						zc->retrieval->prior_windfield->grid_hdir_err[i_w][ip] = c->update_windfield_hdir_err;
					
					tmpH		= sqrt(	pow(zc->retrieval->prior_windfield->grid_u[i_w][ip], 2.) +
										pow(zc->retrieval->prior_windfield->grid_v[i_w][ip], 2.) );
					
					//update u and v errors	
					zc->retrieval->prior_windfield->grid_u_err[i_w][ip] = sqrt(
						pow((zc->retrieval->prior_windfield->grid_u_err[i_w][ip] / tmpH) * zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip] , 2.) +
						pow(zc->retrieval->prior_windfield->grid_v_err[i_w][ip] * zc->retrieval->prior_windfield->grid_hdir_err[i_w][ip] , 2.)
						);
					
					zc->retrieval->prior_windfield->grid_v_err[i_w][ip] = sqrt(
						pow((zc->retrieval->prior_windfield->grid_v_err[i_w][ip] / tmpH) * zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip] , 2.) +
						pow(zc->retrieval->prior_windfield->grid_u_err[i_w][ip] * zc->retrieval->prior_windfield->grid_hdir_err[i_w][ip] , 2.)
						);					
				}
			}
		
			if (	zc->retrieval->prior_windfield->fit_u[i_w] &
					(c->update_windfield_u_err > 0.)
			) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {				
					zc->retrieval->prior_windfield->grid_u_err[i_w][ip] = c->update_windfield_u_err;
				}
			}		
			
			if (	zc->retrieval->prior_windfield->fit_v[i_w] &
					(c->update_windfield_v_err > 0.)
			) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {				
					zc->retrieval->prior_windfield->grid_v_err[i_w][ip] = c->update_windfield_v_err;
				}
			}	
				
			if (	zc->retrieval->prior_windfield->fit_w[i_w] &
					(c->update_windfield_w_err > 0.)
			) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {				
					zc->retrieval->prior_windfield->grid_w_err[i_w][ip] = c->update_windfield_w_err;
				}
			}		
		}
	}
}



void retrieval_fdvar_additional_output(t_fdvar_opc *opc, FILE *fp)
{
	int ip, i_w, io;
    t_zephyros_config 								*zc = opc->zc;
	t_fdvar_o *o = opc->o;
	double *dummy;
	double *tmpu, *tmpv, *tmpw;
	
	fprintf(fp, "section retrieval\n"); fflush(fp);
	fprintf(fp, "subsection post_scattererfield\n"); fflush(fp);
	zephyros_config_print_scattererfield(zc->retrieval->post_scattererfield, fp); fflush(fp);

	fprintf(fp, "subsection post_windfield\n"); fflush(fp);
	zephyros_config_print_windfield(zc->retrieval->post_windfield, fp); fflush(fp);
	

	radarfilter_write_measurements(zc, 1, o->n, o->model_radarmeasurement, fp);

	//obtain u,v,w
	tmpu = malloc(o->n * sizeof(double));
	tmpv = malloc(o->n * sizeof(double));
	tmpw = malloc(o->n * sizeof(double));
	for ( io = 0; io < o->n; io++ ) {
		util_windfield_fuvw(zc->retrieval->post_windfield, o->model_radarmeasurement[io]->advected_center_enu_xyzt, 
			tmpu + io, tmpv + io, tmpw + io,
			0,
			dummy, dummy, dummy,
			0,
			1.
			);
	}


	//u
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_u", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", tmpu[io]);
	fprintf(fp, " \n");
	//v
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_v", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", tmpv[io]);
	fprintf(fp, " \n");
	//w
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_w", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", tmpw[io]);
	fprintf(fp, " \n");
	
}


void retrieval_fdvar_initialize_o(t_fdvar_o **p_fdvar_o, char measurements_filename[8192])
{
	t_fdvar_o *fdvar_o = malloc(sizeof(t_fdvar_o));
	
	radarfilter_read_measurements(&fdvar_o->n, &fdvar_o->radarmeasurement, measurements_filename);
	radarfilter_prepare_model_radarmeasurement(fdvar_o->n, &fdvar_o->model_radarmeasurement, fdvar_o->radarmeasurement);

	*p_fdvar_o = fdvar_o;
}

void retrieval_fdvar_free_o(t_fdvar_o **p_fdvar_o)
{
	t_fdvar_o *fdvar_o = *p_fdvar_o;

	if (fdvar_o != NULL) { 
		radarfilter_free_radarmeasurement(fdvar_o->n, &fdvar_o->radarmeasurement);	
		radarfilter_free_radarmeasurement(fdvar_o->n, &fdvar_o->model_radarmeasurement);	

		free(fdvar_o);
		*p_fdvar_o = NULL;
	}
}

void retrieval_fdvar_initialize_p(t_fdvar_opc *opc)
{
	int i_psd, i_par;
	
	t_fdvar_o *o = opc->o;
	t_fdvar_p *p = malloc(sizeof(t_fdvar_p));
	t_zephyros_config_retrieval_fdvar_cfg *c = opc->c;
	t_zephyros_config *zc = opc->zc;
	int i_w, io;
	
	opc->p = p; //forward memory allocation of p
    
	p->Kn = 0;	//number of paremeters
	p->cost_no = 0;  //totale number of observations, e.g. 2 x o->n, reflection + vr

	//walk through psds and assign Kn nr's.
	for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {
			//dBlwc
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr = p->Kn;
				p->Kn += zc->retrieval->prior_scattererfield->psd[i_psd]->field->n;
			}
			//dBN
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) {
				zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr = p->Kn;
				p->Kn += zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters;
			}
		}
	}
	
	//walk trough wind fields
	for ( i_w = 0; i_w <= 100; i_w++ ) {			
		if (zc->retrieval->post_windfield->field[i_w] != NULL) {
			//grid_u
			if (zc->retrieval->prior_windfield->fit_u[i_w]) {
				zc->retrieval->prior_windfield->fit_u_Knr[i_w] = p->Kn;
				p->Kn += zc->retrieval->prior_windfield->field[i_w]->n;
			}
			
			//grid_v
			if (zc->retrieval->prior_windfield->fit_v[i_w]) {
				zc->retrieval->prior_windfield->fit_v_Knr[i_w] = p->Kn;
				p->Kn += zc->retrieval->prior_windfield->field[i_w]->n;
			}

			//grid_w
			if (zc->retrieval->prior_windfield->fit_w[i_w]) {
				zc->retrieval->prior_windfield->fit_w_Knr[i_w] = p->Kn;
				p->Kn += zc->retrieval->prior_windfield->field[i_w]->n;
			}
		}
	}
	

	//turbulence, grid_edr13	
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
			if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
				zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr = p->Kn;
				p->Kn += zc->retrieval->prior_windfield->turbulence[i_w]->field->n;
			}
		}
	}
	
	//***
	if (c->costfunction_dBZ_hh) p->cost_no += o->n;
	if (c->costfunction_dBZdr) p->cost_no += o->n;
	if (c->costfunction_dBLdr) p->cost_no += o->n;
	if (c->costfunction_Doppler_velocity_hh_ms) p->cost_no += o->n;
	if (c->costfunction_Doppler_spectral_width_hh_ms) p->cost_no += o->n;

    //spectra
	if (c->costfunction_Doppler_spectrum_dBZ_hh) {
		for (io=0; io < o->n; io++) p->cost_no += o->radarmeasurement[io]->n_spectrum;
	}
	if (c->costfunction_specific_dBZdr) {
		for (io=0; io < o->n; io++) p->cost_no += o->radarmeasurement[io]->n_spectrum;
	}
	if (c->costfunction_specific_dBLdr) {
		for (io=0; io < o->n; io++) p->cost_no += o->radarmeasurement[io]->n_spectrum;
	}

	//TBD: list incomplete
		
		
	//initialize bounds
	func_dbl_arr_malloc( p->Kn, &p->Klbound);
	func_dbl_arr_malloc( p->Kn, &p->Kubound); 	
}

void retrieval_fdvar_free_p(t_fdvar_opc *opc)
{
	t_fdvar_p *p = opc->p;
  
	free(p->Klbound);
	free(p->Kubound);	

	free(p);
}


// last step to solve QR
void retrieval_fdvar_cs_qrsol(t_zephyros_field_errcovmat *ecm, double *b, double *x)
{
	int k;
	double *xtmp;
	
	xtmp = cs_calloc (ecm->mat_S ? ecm->mat_S->m2 : 1, sizeof (CS_ENTRY)) ;    //get workspace 
	cs_ipvec (ecm->mat_S->pinv, b, xtmp, ecm->n) ;   // x(0:m-1) = b(p(0:m-1) 
	for (k = 0 ; k < ecm->n ; k++)       // apply Householder refl. to x
	{
		cs_happly (ecm->mat_N->L, k, ecm->mat_N->B[k], xtmp) ;
	}
	cs_usolve (ecm->mat_N->U, xtmp) ;           // x = R\x 
    //cs_ipvec (S->q, xtmp, b, n) ;      /* b(q(0:n-1)) = x(0:n-1) */

	for (k = 0 ; k < ecm->n ; k++)       // apply Householder refl. to x
		x[k] = xtmp[k];
		
	free(xtmp);
}


/*
cs *cs_qr_obtain_inverse(
	csn *N,
	css *S,
	int n,
	double tol
	)
{
	double* x;
	double* b;
	double  xmaxabs;
	double  tmp;
	int i, j;
	cs *invmat_triplet;
	
    invmat_triplet 		= cs_spalloc (n, n, 1, 1, 1) ;

	func_dbl_arr_malloc(n, &x);
	func_dbl_arr_malloc(n, &b);
	
	for (i = 0 ; i < n ; i++) {
		for (j = 0 ; j < n ; j++) {
			if (i == j) {
				b[j] = 1.;
			} else {
				b[j] = 0.;
			}
		}
			
		cs_qrsol2(N, S, b, x, n);

		//determine max
		xmaxabs = 0.;
		for (j = 0 ; j < n ; j++) {
			tmp = fabs(x[j]);
			if (tmp > xmaxabs) xmaxabs = tmp;
		}
		
		//cast solution
		for (j = 0 ; j < n ; j++) {
			if (fabs(x[j]) > (tol * xmaxabs)) {
					if (!cs_entry (invmat_triplet, j, i, x[j])) {
					printf("Error with matrix allocation");
					exit(1);
				}
			}
		}
	}
	
	free(x);
	free(b);
	
    return cs_compress (invmat_triplet) ;             
}
*/

//calculate b = A x
void cs_matrix_vector (
	cs *A,
	double* x,
	double* y,
	int n
	)
{
	cs *x_matrixtriplet;
	cs *x_matrix;
	cs *sol_matrix;
	int i;
	
	//make matrix from *x
    x_matrixtriplet 		= cs_spalloc (n, 1, 1, 1, 1) ;
	for (i = 0 ; i < n ; i++) {
		if (!cs_entry (x_matrixtriplet, i, 1, x[i])) {
			printf("Error with matrix allocation");
			exit(1);
		}
	}
	
	x_matrix = cs_compress (x_matrixtriplet) ; 
	cs_spfree(x_matrixtriplet);
	
	sol_matrix = cs_multiply(A, x_matrix);
	
	//heuristic
	for (i = 0 ; i < n ; i++) {
		y[i] = sol_matrix->x[i];
	}		
	
}

