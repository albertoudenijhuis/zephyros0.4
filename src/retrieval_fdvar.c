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
	//initialize p
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_initialize_p\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_initialize_p(opc);

	//miminize cost function
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_minimize_cost_function\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_minimize_cost_function(opc);
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

	func_dbl_arr_calloc(p->Kn, &xu);
	func_dbl_arr_calloc(p->Kn, &xl);
	
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

//	double dPvrK_dhspeed, dPvrK_dhdir;

	
	double alpha0, alpha, *griddep;
	double *delta, *Sinv_delta, tmp;
	double mytmp, mytmp2, myfactor;
	
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

	//~ double sigma_vr;
	//~ double drefl_dhspeed, drefl_dhdir, drefl_dw;
	//~ double dradial_vel_ms_dhspeed, dradial_vel_ms_dhdir, dradial_vel_ms_dw;
	//~ 
	
	//Unpack K
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_unpack_x\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_unpack_x(opc,x);

	#ifdef _ZEPHYROS_FDVAR_DEBUG
	/*
	for (in=0; in < *n; in++) 
		printf("K[%i] = %.2e\n", in, p->Klbound[in] + (p->Kubound[in] - p->Klbound[in]) * x[in]);
	for (in=0; in < *n; in++) 
		printf("Klbound[%i] = %.2e\n", in, p->Klbound[in]);
	for (in=0; in < *n; in++) 
		printf("Kubound[%i] = %.2e\n", in, p->Kubound[in]);
	*/
	#endif 
	
	//set values to zero
	if (calc_costfunction) {
		*f	= 0.;	
	}
	if (calc_costfunction_derivatives) {
		for (in=0; in < *n; in++) 
		{
			fd[in] =0.;
		}
	}
	
	//~ if (calc_costfunction_derivatives) {
		//~ //calculation of parameter posterior errors, take inverse of prior estimates
		//~ 
		//~ for (ip=0; ip < p->field->n; ip++) {
			//~ p->post_srefl[ip] 	= 1.e-10; //pow(p->prior_srefl[ip],-2);
			//~ p->post_shspeed[ip] = 1.e-10; //pow(p->prior_shspeed[ip],-2);
			//~ p->post_shdir[ip] 	= 1.e-10; //pow(p->prior_shdir[ip],-2);
			//~ p->post_sw[ip] 		= 1.e-10; //pow(p->prior_sw[ip],-2);
		//~ }
	//~ }
	
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
									
	//update todo list with logic
	radarfilter_prepare_todolist(zc, 1, todo);
	
	radarfilter_exec(zc, 1, todo, o->n, o->model_radarmeasurement);

    
    
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
			if (c->costfunction_dBZ_hh) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) / (o->radarmeasurement[io]->dBZ_hh_err), 2);
			if (c->costfunction_dBZdr) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) / (o->radarmeasurement[io]->dBZdr_err), 2);
			if (c->costfunction_dBLdr) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->dBLdr - o->radarmeasurement[io]->dBLdr) / (o->radarmeasurement[io]->dBLdr_err), 2);
			if (c->costfunction_Doppler_velocity_hh_ms) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) / (o->radarmeasurement[io]->Doppler_velocity_hh_ms_err), 2);
			if (c->costfunction_Doppler_spectral_width_hh_ms) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->Doppler_spectral_width_hh_ms - o->radarmeasurement[io]->Doppler_spectral_width_hh_ms) / (o->radarmeasurement[io]->Doppler_spectral_width_hh_ms_err), 2);

			/*
			#ifdef _ZEPHYROS_FDVAR_DEBUG
			if (c->costfunction_dBZ_hh) printf("o->model_radarmeasurement[%i]->dBZ_hh = %.2e\n", io, o->model_radarmeasurement[io]->dBZ_hh);
			if (c->costfunction_dBZdr) printf("o->model_radarmeasurement[%i]->dBZdr = %.2e\n", io, o->model_radarmeasurement[io]->dBZdr);
			if (c->costfunction_dBLdr) printf("o->model_radarmeasurement[%i]->dBLdr = %.2e\n", io, o->model_radarmeasurement[io]->dBLdr);
			if (c->costfunction_Doppler_velocity_hh_ms) printf("o->model_radarmeasurement[%i]->Doppler_velocity_hh_ms = %.2e\n", io, o->model_radarmeasurement[io]->Doppler_velocity_hh_ms);
			if (c->costfunction_Doppler_spectral_width_hh_ms) printf("o->model_radarmeasurement[%i]->Doppler_spectral_width_hh_ms = %.2e\n", io, o->model_radarmeasurement[io]->Doppler_spectral_width_hh_ms);
			#endif
			*/
			
			//spectra
			for ( i_int = 0; i_int < o->radarmeasurement[io]->n_spectrum; i_int++ ) {
				if (c->costfunction_Doppler_spectrum_dBZ_hh) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int] - o->radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int]) / (o->radarmeasurement[io]->Doppler_spectrum_dBZ_hh_err[i_int]), 2);
				if (c->costfunction_specific_dBZdr) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->specific_dBZdr[i_int] - o->radarmeasurement[io]->specific_dBZdr[i_int]) / (o->radarmeasurement[io]->specific_dBZdr_err[i_int]), 2);
				if (c->costfunction_specific_dBLdr) *f += (1./ p->cost_no) * pow( (o->model_radarmeasurement[io]->specific_dBLdr[i_int] - o->radarmeasurement[io]->specific_dBLdr[i_int]) / (o->radarmeasurement[io]->specific_dBLdr_err[i_int]), 2);

				/*
				#ifdef _ZEPHYROS_FDVAR_DEBUG
				if (c->costfunction_Doppler_spectrum_dBZ_hh) printf("o->model_radarmeasurement[%i]->Doppler_spectrum_dBZ_hh[%i] = %.2e\n", io, i_int, o->model_radarmeasurement[io]->Doppler_spectrum_dBZ_hh[i_int]);
				if (c->costfunction_specific_dBZdr) printf("o->model_radarmeasurement[%i]->specific_dBZdr[%i] = %.2e\n", io, i_int, o->model_radarmeasurement[io]->specific_dBZdr[i_int]);
				if (c->costfunction_specific_dBLdr) printf("o->model_radarmeasurement[%i]->specific_dBLdr[%i] = %.2e\n", io, i_int, o->model_radarmeasurement[io]->specific_dBLdr[i_int]);
				#endif
				*/
			}




			//TBD: all other variables that are in the cost function
		}

		
		//printf("o->model_radarmeasurement[%i]->dBZ_hh = %.2e\n", io, o->model_radarmeasurement[io]->dBZ_hh);
		//~ if (calc_costfunction) {
			//~ //observation space of cost function
			//~ if (c->fit_vr)	*f += (1./ p->on) * pow( (res_vol->integrated_radial_vel_ms  - o->vr[io] ) 	/ o->svr[io], 2);
		//~ }
		
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
							o->model_radarmeasurement[io]->center_coor->enu_xyzt,
							griddep);

						//dBZ_hh
						if (c->costfunction_dBZ_hh) {
							alpha = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) / pow(o->radarmeasurement[io]->dBZ_hh_err, 2);
							for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
								fd[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
									 alpha * griddep[ip];
							}
						}

						//Doppler_spectrum_dBZ_hh
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
							o->model_radarmeasurement[io]->center_coor->enu_xyzt,
							griddep);

						//dBZ_hh
						if (c->costfunction_dBZ_hh) {
							alpha0 = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) / pow(o->radarmeasurement[io]->dBZ_hh_err, 2);
							eta_hh = func_dB_inv(o->model_radarmeasurement[io]->dBZ_hh) / o->model_radarmeasurement[io]->coef_eta2Z;

							for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
								alpha = alpha0 * (o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] / eta_hh);								
								if (isnanorinf(&alpha) == 0) {
									for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
										fd[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
											+ (zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
											+ ip] += 
											 alpha * griddep[ip];
									}
								}
							}
						}
						
						//dBZdr
						if (c->costfunction_dBZdr) {
							alpha0 = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) / pow(o->radarmeasurement[io]->dBZdr_err, 2);
							eta_hh = func_dB_inv(o->model_radarmeasurement[io]->dBZ_hh) / o->model_radarmeasurement[io]->coef_eta2Z;
							eta_vv = func_dB_inv(o->model_radarmeasurement[io]->dBZ_vv) / o->model_radarmeasurement[io]->coef_eta2Z;

							for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
								alpha = alpha0 * ((o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] / eta_hh) - (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par] / eta_vv));								

								if (isnanorinf(&alpha) == 0) {
									for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
										fd[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
											+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
											+ip] += 
											 alpha * griddep[ip];
									}
								}
							}
						}
						//dBLdr
						if (c->costfunction_dBLdr) {
							alpha0 = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->dBLdr - o->radarmeasurement[io]->dBLdr) / pow(o->radarmeasurement[io]->dBLdr_err, 2);
							eta_hv = func_dB_inv(o->model_radarmeasurement[io]->dBZ_hv) / o->model_radarmeasurement[io]->coef_eta2Z;
							eta_vv = func_dB_inv(o->model_radarmeasurement[io]->dBZ_vv) / o->model_radarmeasurement[io]->coef_eta2Z;

							for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
								alpha = alpha0 * ((o->model_radarmeasurement[io]->eta_i_hv[i_psd][i_par] / eta_hv) - (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par] / eta_vv));								
								
								if (isnanorinf(&alpha) == 0) {
									for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
										fd[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
											+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
											+ip] += 
											 alpha * griddep[ip];
									}
								}
							}
						}
						
						//spectral variables
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
		
			//loop over wind fields
			for ( i_w = 0; i_w <= 100; i_w++ ) {
				//TBD
				//hspeed, hdir
				
				if (zc->retrieval->prior_windfield->fit_u[i_w]) {
					griddep = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
					interpolation_bilint_griddep(
						zc->retrieval->prior_windfield->lut_u[i_w],
						o->model_radarmeasurement[io]->center_coor->enu_xyzt,
						griddep);

					//grid_u
					if (c->costfunction_Doppler_velocity_hh_ms) {
						alpha = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) / pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, 2);
						alpha *= o->model_radarmeasurement[io]->center_coor->radar_enu_dir[0];
						//TBD
						//it would be better to use reflectivity weighted coordinates
						for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
							fd[
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
						o->model_radarmeasurement[io]->center_coor->enu_xyzt,
						griddep);

					//grid_v
					if (c->costfunction_Doppler_velocity_hh_ms) {
						alpha = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) / pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, 2);
						alpha *= o->model_radarmeasurement[io]->center_coor->radar_enu_dir[1];
						//TBD
						//it would be better to use reflectivity weighted coordinates
						for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
							fd[
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
						o->model_radarmeasurement[io]->center_coor->enu_xyzt,
						griddep);

					//grid_w
					if (c->costfunction_Doppler_velocity_hh_ms) {
						alpha = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) / pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, 2);
						alpha *= o->model_radarmeasurement[io]->center_coor->radar_enu_dir[2];
						//TBD
						//it would be better to use reflectivity weighted coordinates
						for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
							fd[
								zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] += 
								 alpha * griddep[ip];
						}
					}
					
					free(griddep);
				}
			}
			
			//loop over turbulence fields
			for ( i_w = 0; i_w <= 100; i_w++ ) {
				if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
					if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
						griddep = malloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n * sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_windfield->turbulence[i_w]->lut_edr13,
							o->model_radarmeasurement[io]->center_coor->enu_xyzt,
							griddep);

						//grid_edr13
						if (c->costfunction_Doppler_spectral_width_hh_ms) {
							alpha = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->Doppler_spectral_width_hh_ms - o->radarmeasurement[io]->Doppler_spectral_width_hh_ms) / pow(o->radarmeasurement[io]->Doppler_spectral_width_hh_ms_err, 2);
							alpha *= o->model_radarmeasurement[io]->der_edr13_Doppler_spectral_width_hh_ms;
							
							for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
								fd[
									zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
									 alpha * griddep[ip];
							}
						}

						/*
						if (c->costfunction_dBZ_hh) {
							alpha0 = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) / pow(o->radarmeasurement[io]->dBZ_hh_err, 2);
							alpha = alpha0 * o->model_radarmeasurement[io]->der_edr13_dBZ_hh;
							
							for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
								fd[
									zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
									 alpha * griddep[ip];
							}
						}
						*/
						
						if (c->costfunction_dBZdr) {
							alpha0 = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) / pow(o->radarmeasurement[io]->dBZdr_err, 2);
							alpha = alpha0 * o->model_radarmeasurement[io]->der_edr13_dBZdr;
							
							for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
								fd[
									zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
									 alpha * griddep[ip];
							}

						}
						
						if (c->costfunction_dBLdr) {
							alpha0 = (1./p->cost_no) * 2. * (o->model_radarmeasurement[io]->dBLdr - o->radarmeasurement[io]->dBLdr) / pow(o->radarmeasurement[io]->dBLdr_err, 2);
							alpha = alpha0 * o->model_radarmeasurement[io]->der_edr13_dBLdr;

							for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
								fd[
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

				retrieval_fdvar_cs_qrsol(&zc->retrieval->prior_scattererfield->psd[i_psd]->dBlwc_ecm, delta, Sinv_delta);							
				
				dbldotprod(&zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) *f += (tmp / p->Kn);
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {						
						fd[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
							(1./ p->Kn) * 2. * delta[ip] * Sinv_delta[ip];
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
		//grid_hspeed
		if (zc->retrieval->prior_windfield->fit_hspeed[i_w]) {
			delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
			Sinv_delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

			//delta is post - prior value
			for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
				delta[ip] = 
					zc->retrieval->post_windfield->grid_hspeed[i_w][ip]
					-
					zc->retrieval->prior_windfield->grid_hspeed[i_w][ip];						
			}

			retrieval_fdvar_cs_qrsol(&zc->retrieval->prior_windfield->hspeed_ecm[i_w], delta, Sinv_delta);							
			dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
										
			if (calc_costfunction) *f += (tmp / p->Kn);
			
			if (calc_costfunction_derivatives) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
					fd[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip] += 
						(1./ p->Kn) * 2. * delta[ip] * Sinv_delta[ip];
				}
			}
			
			free(delta);
			free(Sinv_delta);	
		}
		
		//grid_hdir
		if (zc->retrieval->prior_windfield->fit_hdir[i_w]) {
			delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
			Sinv_delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

			//delta is post - prior value
			for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
				delta[ip] = 
					zc->retrieval->post_windfield->grid_hdir[i_w][ip]
					-
					zc->retrieval->prior_windfield->grid_hdir[i_w][ip];						
			}

			retrieval_fdvar_cs_qrsol(&zc->retrieval->prior_windfield->hdir_ecm[i_w], delta, Sinv_delta);							
			dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
										
			if (calc_costfunction) *f += (tmp / p->Kn);
			
			if (calc_costfunction_derivatives) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
					fd[zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip] += 
						(1./ p->Kn) * 2. * delta[ip] * Sinv_delta[ip];
				}
			}
			
			free(delta);
			free(Sinv_delta);	
		}
					
		//grid_u
		if (zc->retrieval->prior_windfield->fit_u[i_w]) {
			delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
			Sinv_delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

			//delta is post - prior value
			for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
				delta[ip] = 
					zc->retrieval->post_windfield->grid_u[i_w][ip]
					-
					zc->retrieval->prior_windfield->grid_u[i_w][ip];						
			}

			retrieval_fdvar_cs_qrsol(&zc->retrieval->prior_windfield->u_ecm[i_w], delta, Sinv_delta);							
			dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
										
			if (calc_costfunction) *f += (tmp / p->Kn);
			
			if (calc_costfunction_derivatives) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
					fd[zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] += 
						(1./ p->Kn) * 2. * delta[ip] * Sinv_delta[ip];
				}
			}
			
			free(delta);
			free(Sinv_delta);	
		}
		
		//grid_v
		if (zc->retrieval->prior_windfield->fit_v[i_w]) {
			delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));
			Sinv_delta = malloc(zc->retrieval->prior_windfield->field[i_w]->n * sizeof(double));

			//delta is post - prior value
			for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
				delta[ip] = 
					zc->retrieval->post_windfield->grid_v[i_w][ip]
					-
					zc->retrieval->prior_windfield->grid_v[i_w][ip];						
			}

			retrieval_fdvar_cs_qrsol(&zc->retrieval->prior_windfield->v_ecm[i_w], delta, Sinv_delta);							
			dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
										
			if (calc_costfunction) *f += (tmp / p->Kn);
			
			if (calc_costfunction_derivatives) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
					fd[zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] += 
						(1./ p->Kn) * 2. * delta[ip] * Sinv_delta[ip];
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

			retrieval_fdvar_cs_qrsol(&zc->retrieval->prior_windfield->w_ecm[i_w], delta, Sinv_delta);							
			dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
										
			if (calc_costfunction) *f += (tmp / p->Kn);
			
			if (calc_costfunction_derivatives) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
					fd[zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] += 
						(1./ p->Kn) * 2. * delta[ip] * Sinv_delta[ip];
				}
			}
			
			free(delta);
			free(Sinv_delta);	
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
					
					//#ifdef _ZEPHYROS_FDVAR_DEBUG
					//printf("zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[%i] = %.2e\n", ip, zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[ip]);
					//#endif

					delta[ip] = 
						zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[ip]
						-
						zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13[ip];						
				}

				retrieval_fdvar_cs_qrsol(&zc->retrieval->prior_windfield->turbulence[i_w]->edr13_ecm, delta, Sinv_delta);							
				dbldotprod(&zc->retrieval->prior_windfield->turbulence[i_w]->field->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) *f += (tmp / p->Kn);

				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {					
						fd[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
							(1./ p->Kn) * 2. * delta[ip] * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);
			}
		}
	}
	
    //End of parameter space
    //***
	
	
	
	
	
	
	
	
	//parameter space of cost function, hspeed, hdir, w
	//~ for (ip = 0; ip < p->field->n; ip++) {
		//~ if (c->fit_refl) 		deltarefl[ip]	= p->post_refl[ip] 			- p->prior_refl[ip];
		//~ if (c->fit_hspeed) 		deltahspeed[ip]	= p->post_hspeed[ip]		- p->prior_hspeed[ip];
		//~ if (c->fit_hdir)		deltahdir[ip]	= angleAminB(p->post_hdir + ip, p->prior_hdir + ip);
		//~ if (c->fit_w)			deltaw[ip]		= p->post_w[ip] 			- p->prior_w[ip];
	//~ }

	//obtain Sreflpinv_deltarefl
	//~ if (c->fit_refl) 		cs_qrsol2(p->refl_errcovmat_N, p->refl_errcovmat_S, deltarefl, Sreflpinv_deltarefl, p->field->n);
	//~ if (c->fit_hspeed)		cs_qrsol2(p->refl_errcovmat_N, p->refl_errcovmat_S, deltahspeed, Shspeedpinv_deltahspeed, p->field->n);
	//~ if (c->fit_hdir)		cs_qrsol2(p->refl_errcovmat_N, p->refl_errcovmat_S, deltahdir, Shdirpinv_deltahdir, p->field->n);
	//~ if (c->fit_w)			cs_qrsol2(p->refl_errcovmat_N, p->refl_errcovmat_S, deltaw, Swpinv_deltaw, p->field->n);
	   
	    
	//~ if (c->fit_refl) 		cs_matrix_vector(p->refl_errcovmatinv, deltarefl,  Sreflpinv_deltarefl, p->field->n);
	//~ if (c->fit_hspeed) 		cs_matrix_vector(p->hdir_errcovmatinv, deltahspeed, Shspeedpinv_deltahspeed, p->field->n);
	//~ if (c->fit_hdir) 		cs_matrix_vector(p->hspeed_errcovmatinv, deltahdir, Shdirpinv_deltahdir, p->field->n);
	//~ if (c->fit_w) 			cs_matrix_vector(p->w_errcovmatinv, deltaw, Swpinv_deltaw, p->field->n);

	    
	//~ if (c->fit_refl)		dbldotprod(&p->field->n, deltarefl, Sreflpinv_deltarefl, &tmp1);
	//~ if (c->fit_hspeed)		dbldotprod(&p->field->n, deltahspeed, Shspeedpinv_deltahspeed, &tmp2);
	//~ if (c->fit_hdir)		dbldotprod(&p->field->n, deltahdir, Shdirpinv_deltahdir, &tmp3);
	//~ if (c->fit_w)			dbldotprod(&p->field->n, deltaw, Swpinv_deltaw, &tmp4);
                   
	//parameter space
    //~ if (c->fit_refl)	tmp1 /= p->Kn;
    //~ if (c->fit_hspeed)	tmp2 /= p->Kn;
    //~ if (c->fit_hdir)	tmp3 /= p->Kn;
    //~ if (c->fit_w)		tmp4 /= p->Kn;
	
	//parameter space
    //~ if (c->fit_refl)	tmp1 /= (c->prior_srefl * c->prior_srefl);
    //~ if (c->fit_hspeed)	tmp2 /= (c->prior_shspeed * c->prior_shspeed);
    //~ if (c->fit_hdir)	tmp3 /= (c->prior_shdir * c->prior_shdir);
    //~ if (c->fit_w)		tmp4 /= (c->prior_sw * c->prior_sw);
	
	//~ if (calc_costfunction) {
		//~ //parameter space
		//~ if (c->fit_refl)	*f += tmp1;
		//~ if (c->fit_hspeed)	*f += tmp2;
		//~ if (c->fit_hdir)	*f += tmp3;
		//~ if (c->fit_w)		*f += tmp4;
		//debugging
		//~ if (c->fit_refl)	printf("tmp1 = %f \n", tmp1);
		//~ if (c->fit_hspeed)	printf("tmp2 = %f \n", tmp2);
		//~ if (c->fit_hdir)	printf("tmp3 = %f \n", tmp3);
		//~ if (c->fit_w)		printf("tmp4 = %f \n", tmp4);
	//~ }
	

	
	
	//~ if (calc_costfunction_derivatives) {
		//~ for (ip=0; ip < p->field->n; ip++) {
			//~ //parameter space
			//~ if (c->fit_refl)	fd[p->fit_Knr_refl  	+ip] += (1./ p->Kn) * 2. * Sreflpinv_deltarefl[ip] * pow(c->prior_srefl, -2.);
			//~ if (c->fit_hspeed)	fd[p->fit_Knr_hspeed	+ip] += (1./ p->Kn) * 2. * Shspeedpinv_deltahspeed[ip] * pow(c->prior_shspeed, -2.);
			//~ if (c->fit_hdir)	fd[p->fit_Knr_hdir		+ip] += (1./ p->Kn) * 2. * Shdirpinv_deltahdir[ip] * pow(c->prior_shdir, -2.);
			//~ if (c->fit_w)		fd[p->fit_Knr_w			+ip] += (1./ p->Kn) * 2. * Swpinv_deltaw[ip] * pow(c->prior_sw, -2.) ;
			//~ 
		//~ }
	//~ }


	//obtain best estimate for the uncertainty in the measurements
	//~ if (calc_costfunction_derivatives) {
		//~ if (c->fit_hspeed | c->fit_hdir | c->fit_w) {
			//~ sigma_vr = 0.;		
			//~ for (io=0; io < o->n; io++) {
				//~ sigma_vr += pow( p->PvrK[io]  - o->vr[io] , 2.);
			//~ }
			//~ sigma_vr = sqrt((1. / p->cost_no) * sigma_vr);		
			//~ //sigma_vr = sqrt((1. / (p->cost_no - DOF)) * sigma_vr); //TBD: update with right degrees of freedom
//~ 
			//~ for (io=0; io < o->n; io++) {
				//~ o->windvector_su[io] *= sigma_vr; //continued calculation
				//~ o->windvector_sv[io] *= sigma_vr; //continued calculation
				//~ o->windvector_sw[io] *= sigma_vr; //continued calculation
			//~ }
		//~ }
	//~ }



	//~ if (calc_costfunction_derivatives) {
		//~ //calculation of parameter posterior errors, take inverse
		//~ //TBD
		//~ p->post_shspeedg		= pow(p->post_shspeedg,-1./2);
		//~ p->post_shdirg			= pow(p->post_shdirg,-1./2);
		//~ p->post_sw0				= pow(p->post_sw0,-1./2)	;
		//~ for (ip=0; ip < p->field->n; ip++) {
			//~ p->post_srefl[ip] 	= pow(p->post_srefl[ip],-1./2);
			//~ p->post_shspeed[ip] = pow(p->post_shspeed[ip],-1./2);
			//~ p->post_shdir[ip] 	= pow(p->post_shdir[ip],-1./2);
			//~ p->post_sw[ip] 		= pow(p->post_sw[ip],-1./2);
		//~ }
	//~ }

	
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
		for (in=0; in < *n; in++) {
			if (fd[in] != 0.) {
				printf("fd[%i] = %.5e \n", in, fd[in]); fflush(stdout);
			}

		}
	}
	#endif
	
	prev_costfunction = *f;
	 
	radarfilter_free_todolist(&todo);
	
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
	
	func_dbl_arr_calloc(p->Kn, x);

	//set Klbound, Kubound and x value
	
	for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		if (zc->retrieval->prior_scattererfield->psd[i_psd] != NULL) {
			//dBlwc
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
					p->Klbound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] =
						func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]) -
						(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3[ip]);
					p->Kubound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] =
						func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]) +
						(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3[ip]);
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
								
						p->Kubound[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr
									+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									+ip] =
							func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip]) +
							(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBnumber_density_err_m3[i_par][ip]);
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
		//grid_hspeed
		if (zc->retrieval->prior_windfield->fit_hspeed[i_w]) {
			for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
				p->Klbound[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] + ip] =
					zc->retrieval->prior_windfield->grid_hspeed[i_w][ip] -
					(factor * zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip]);
					
				if (p->Klbound[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] + ip] <= 0)
					p->Klbound[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] + ip] = 
						1.e-10 * zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip];
					
				p->Kubound[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] + ip] =
					zc->retrieval->prior_windfield->grid_hspeed[i_w][ip] +
					(factor * zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip]);				
				(*x)[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip] =
					K2x(p, &(zc->retrieval->prior_windfield->grid_hspeed[i_w][ip]), zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip);
			}
		}
		
		//grid_hdir
		if (zc->retrieval->prior_windfield->fit_hdir[i_w]) {
			for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
				p->Klbound[zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] + ip] =
					zc->retrieval->prior_windfield->grid_hdir[i_w][ip] -
					(factor * zc->retrieval->prior_windfield->grid_hdir_err[i_w][ip]);
				p->Kubound[zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] + ip] =
					zc->retrieval->prior_windfield->grid_hdir[i_w][ip] +
					(factor * zc->retrieval->prior_windfield->grid_hdir_err[i_w][ip]);				
				(*x)[zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip] =
					K2x(p, &(zc->retrieval->prior_windfield->grid_hdir[i_w][ip]), zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip);
			}
		}
					
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
	
	//turbulence, grid_edr13	
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->prior_windfield->turbulence[i_w] != NULL) {
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
	
	//for (ip=0; ip < p->field->n; ip++) {
	//	p->prior_u[ip] = -1. * p->prior_hspeed[ip] * sin(p->prior_hdir[ip]);
	//	p->prior_v[ip] = -1. * p->prior_hspeed[ip] * cos(p->prior_hdir[ip]);
	//}
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
		if (zc->retrieval->prior_scattererfield->psd[i_psd] != NULL) {
			//dBlwc
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				if (zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3 == NULL) 
					zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3 = malloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * sizeof(double));
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
		//grid_hspeed
		if (zc->retrieval->prior_windfield->fit_hspeed[i_w]) {
			for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
				zc->retrieval->post_windfield->grid_hspeed[i_w][ip] =
				x2K(p, 	x + zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip,
						zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip);
			}
		}
		
		//grid_hdir
		if (zc->retrieval->prior_windfield->fit_hdir[i_w]) {
			for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
				zc->retrieval->post_windfield->grid_hdir[i_w][ip] =
				x2K(p, 	x + zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip,
						zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip);
			}
		}
					
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
		
		//TBD: if hspeed & hdir, then update u and v
		//	if ((c->fit_hspeed) | (c->fit_hdir)) 	p->post_u[ip] = -1. * p->post_hspeed[ip] * sin(p->post_hdir[ip]);
		//	if ((c->fit_hspeed) | (c->fit_hdir)) 	p->post_v[ip] = -1. * p->post_hspeed[ip] * cos(p->post_hdir[ip]);		
		
	}
	

	//turbulence, grid_edr13	
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->prior_windfield->turbulence[i_w] != NULL) {
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
	









void retrieval_fdvar_cast_p(t_fdvar_opc *opc1, t_fdvar_opc *opc2)
{
/*
	int calcderivatives;
	int bilint_special;
	double derivatives[4];
	double outx[4];
	double myval;
	
	int i;
		
	//debugging
	//~ printf("in casting function\n");
	//~ for (i=0; i < p1->field->n; i++) {
		//~ printf("p1->post_refl[%i] = %.2e\n", i, p1->post_refl[i]);
	//~ }
		
	calcderivatives = 0; //no
	for (i=0; i < p2->field->n; i++) {
		outx[0] = p2->field->x[i];
		outx[1] = p2->field->y[i];
		outx[2] = p2->field->z[i];
		outx[3] = p2->field->t[i];
		
		bilint_special	= 0; //normal interpolation
		bilint(&p1->field->n_dim, &p1->field->n, p1->field->int_ax_count, p1->field->int_ax_values, p1->post_refl, outx, p2->prior_refl + i, &calcderivatives, &bilint_special, derivatives);
		bilint(&p1->field->n_dim, &p1->field->n, p1->field->int_ax_count, p1->field->int_ax_values, p1->post_hspeed, outx, p2->prior_hspeed + i, &calcderivatives, &bilint_special, derivatives);
		bilint(&p1->field->n_dim, &p1->field->n, p1->field->int_ax_count, p1->field->int_ax_values, p1->post_w, outx, p2->prior_w + i, &calcderivatives, &bilint_special, derivatives);
		bilint_special	= 1; //angular interpolation
		bilint(&p1->field->n_dim, &p1->field->n, p1->field->int_ax_count, p1->field->int_ax_values, p1->post_hdir, outx, p2->prior_hdir + i, &calcderivatives, &bilint_special, derivatives);

		//debugging
		//~ printf("p1->post_refl[%i] = %.2e\n", i, p1->post_refl[i]);
		//~ printf("p2->prior_refl[%i] = %.2e\n", i, p2->prior_refl[i]);
	}
	
	printf("solution casted\n");
*/
}



void retrieval_fdvar_additional_output(t_fdvar_opc *opc, FILE *fp)
{
    t_zephyros_config 								*zc = opc->zc;
	t_fdvar_o *o = opc->o;

	fprintf(fp, "section retrieval\n"); fflush(fp);
	fprintf(fp, "subsection post_scattererfield\n"); fflush(fp);
	zephyros_config_print_scattererfield(zc->retrieval->post_scattererfield, fp); fflush(fp);

	fprintf(fp, "subsection post_windfield\n"); fflush(fp);
	zephyros_config_print_windfield(zc->retrieval->post_windfield, fp); fflush(fp);
	

	radarfilter_write_measurements(zc, 1, o->n, o->model_radarmeasurement, fp);


/*
	t_fdvar_o *o = opc->o;
	t_fdvar_p *p = opc->p;
	t_zephyros_config_retrieval_fdvar_cfg *c = opc->c;
	
	int io;
	
	//output as !! <varname> <ndim> <dim1> <data>
	if (c->fit_vr) {
		fprintf(fp, "!! %-30s %-15i %-15i", "PvrK", 1, o->n);
		for (io=0; io < o->n; io++) {
			fprintf(fp, " %-15.3e", p->PvrK[io]);
		}
		fprintf(fp, " \n");
	}

	if (c->fit_refl) {	
		fprintf(fp, "!! %-30s %-15i %-15i", "PreflK", 1, o->n);
		for (io=0; io < o->n; io++) {
			fprintf(fp, " %-15.3e", p->PreflK[io]);
		}
		fprintf(fp, " \n");
	}
	
	//vr_residual
	if (c->fit_vr) {
		fprintf(fp, "!! %-30s %-15i %-15i", "vr_residual", 1, o->n);
		for (io=0; io < o->n; io++) {
			fprintf(fp, " %-15.3e", p->PvrK[io] - o->vr[io]);
		
		}
		fprintf(fp, " \n");
	}
	
	//refl_residual
	if (c->fit_refl) {
		fprintf(fp, "!! %-30s %-15i %-15i", "refl_residual", 1, o->n);
		for (io=0; io < o->n; io++) {
			fprintf(fp, " %-15.3e", p->PreflK[io] - o->refl[io]);
		
		}
		fprintf(fp, " \n");
	}
*/
}


void retrieval_fdvar_initialize_o(t_fdvar_o **p_fdvar_o, char measurements_filename[8192])
{
	t_fdvar_o *fdvar_o = malloc(sizeof(t_fdvar_o));
	
	radarfilter_read_measurements(&fdvar_o->n, &fdvar_o->radarmeasurement, measurements_filename);
	//prepare model radar measurements
	radarfilter_prepare_model_radarmeasurement(fdvar_o->n, &fdvar_o->model_radarmeasurement, fdvar_o->radarmeasurement);

	fdvar_o->windvector_u			= malloc(fdvar_o->n * sizeof(double));
	fdvar_o->windvector_v			= malloc(fdvar_o->n * sizeof(double));
	fdvar_o->windvector_w			= malloc(fdvar_o->n * sizeof(double));
	fdvar_o->windvector_u_err		= malloc(fdvar_o->n * sizeof(double));
	fdvar_o->windvector_v_err		= malloc(fdvar_o->n * sizeof(double));
	fdvar_o->windvector_w_err		= malloc(fdvar_o->n * sizeof(double));

	fdvar_o->edr					= malloc(fdvar_o->n * sizeof(double));
	fdvar_o->lwc					= malloc(fdvar_o->n * sizeof(double));


	*p_fdvar_o = fdvar_o;
}

void retrieval_fdvar_free_o(t_fdvar_o **p_fdvar_o)
{
	t_fdvar_o *fdvar_o = *p_fdvar_o;

	if (fdvar_o != NULL) { 
		radarfilter_free_radarmeasurement(fdvar_o->n, &fdvar_o->radarmeasurement);	
		radarfilter_free_radarmeasurement(fdvar_o->n, &fdvar_o->model_radarmeasurement);	

		free(fdvar_o->windvector_u);
		free(fdvar_o->windvector_v);
		free(fdvar_o->windvector_w);
		free(fdvar_o->windvector_u_err);
		free(fdvar_o->windvector_v_err);
		free(fdvar_o->windvector_v_err);
		
		free(fdvar_o->lwc);
		free(fdvar_o->edr);
		
		free(fdvar_o);
		fdvar_o = NULL;
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
		if (zc->retrieval->prior_scattererfield->psd[i_psd] != NULL) {
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
		//grid_hspeed
		if (zc->retrieval->prior_windfield->fit_hspeed[i_w]) {
			zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] = p->Kn;
			p->Kn += zc->retrieval->prior_windfield->field[i_w]->n;
		}
		
		//grid_hdir
		if (zc->retrieval->prior_windfield->fit_hdir[i_w]) {
			zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] = p->Kn;
			p->Kn += zc->retrieval->prior_windfield->field[i_w]->n;
		}
					
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
	

	//turbulence, grid_edr13	
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->prior_windfield->turbulence[i_w] != NULL) {
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
	func_dbl_arr_calloc( p->Kn, &p->Klbound);
	func_dbl_arr_calloc( p->Kn, &p->Kubound);   

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

	func_dbl_arr_calloc(n, &x);
	func_dbl_arr_calloc(n, &b);
	
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

