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
	
	//set to do list for the radar filter
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("radarfilter_initialize_todolist\n"); fflush(stdout);
	#endif
	radarfilter_initialize_todolist(opc->zc, 1, &opc->todo_no_derivatives);
	radarfilter_initialize_todolist(opc->zc, 1, &opc->todo_with_derivatives);

	//TBD, check if dBN is fitted or not. Walk over PSD.
	opc->todo_no_derivatives->calc_eta_i_hh = 		1;
	opc->todo_no_derivatives->calc_eta_i_hv = 		1;
	opc->todo_no_derivatives->calc_eta_i_vv = 		1;
	
	opc->todo_with_derivatives->calc_eta_i_hh = 	1;
	opc->todo_with_derivatives->calc_eta_i_hv = 	1;
	opc->todo_with_derivatives->calc_eta_i_vv =		1;
	
	/*
	TBD ...
	opc->todo_with_derivatives->calc_spectrum_eta_i_hh 	= opc->c->costfunction_Doppler_spectrum_dBZ_hh ;
	opc->todo_with_derivatives->calc_spectrum_eta_i_hv 	= opc->c->costfunction_Doppler_spectrum_dBZ_hv ;
	opc->todo_with_derivatives->calc_spectrum_eta_i_vv 	= opc->c->costfunction_Doppler_spectrum_dBZ_vv ;
	*/
	
	//Ok from here:
	opc->todo_with_derivatives->der_edr13 		= 1;
	opc->todo_no_derivatives->der_edr13 		= 1;
	 		
	opc->todo_no_derivatives->calc_gridaverage_psd_rcs 			= 1;
	opc->todo_with_derivatives->calc_gridaverage_psd_rcs	 	= 1;	 		
	 			
	opc->todo_no_derivatives->der_dBZ_hh	 	= 1;	 		
	opc->todo_with_derivatives->der_dBZ_hh	 	= 1;	 		
	 							
	//update todo list with logic
	radarfilter_prepare_todolist(opc->zc, 1, opc->todo_no_derivatives);
	radarfilter_prepare_todolist(opc->zc, 1, opc->todo_with_derivatives);
	
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar test radar filter\n"); fflush(stdout);
	#endif 	
	
	//test run with all objects, this also calculates the weighted RCS 
	//which is required for casting psd
	radarfilter_initialize_resolution_volume(opc->zc, 1, &opc->res_vol, opc->todo_with_derivatives);
	radarfilter_exec(opc->zc, 1, opc->todo_with_derivatives, opc->o->n, opc->o->model_radarmeasurement, opc->res_vol);
	radarfilter_free_resolution_volume(opc->zc, 1, &opc->res_vol, opc->todo_with_derivatives);
		
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

	printf("active psds = %i\n", opc->c->n_active_psd_nrs);
	printf("active windfields = %i\n", opc->c->n_active_windfield_grid_nrs);
	printf("active turbulences = %i\n", opc->c->n_active_windfield_turbulence_nrs);
	
	//update active elements
	for ( i = 0; i <= 100; i++ ) {
		opc->zc->retrieval->post_scattererfield->psd[i] 	= NULL;
		opc->zc->retrieval->post_windfield->turbulence[i] 	= NULL;
		opc->zc->retrieval->post_windfield->field[i] 		= NULL;

		for (j = 0; j < opc->c->n_active_psd_nrs; j++) {
			if (i == opc->c->active_psd_nrs[j]) {
				opc->zc->retrieval->post_scattererfield->psd[i] = cfg_psd[i];
			}
		}
				
		for (j = 0; j < opc->c->n_active_windfield_grid_nrs; j++) {
			if (i == opc->c->active_windfield_grid_nrs[j]) {			
				opc->zc->retrieval->post_windfield->field[i] = cfg_field[i];			
			}
		}
		for (j = 0; j < opc->c->n_active_windfield_turbulence_nrs; j++) {
			if (i == opc->c->active_windfield_turbulence_nrs[j]) {
				opc->zc->retrieval->post_windfield->turbulence[i] = cfg_turbulence[i];
			}
		}
	}

	//initialize p
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_initialize_p\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_initialize_p(opc);

	//initialize resolution volume
	radarfilter_initialize_resolution_volume(opc->zc, 1, &opc->res_vol, opc->todo_with_derivatives);

	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_minimize_cost_function\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_minimize_cost_function(opc);

	//free reoslution volume
	radarfilter_free_resolution_volume(opc->zc, 1, &opc->res_vol, opc->todo_with_derivatives);

	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_estimate_posterior_errors\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_estimate_posterior_errors(opc);

	//restore post configuration pointers
	for ( i = 0; i <= 100; i++ ) {
		opc->zc->retrieval->post_scattererfield->psd[i] 	= cfg_psd[i];
		opc->zc->retrieval->post_windfield->field[i] 		= cfg_field[i];
		opc->zc->retrieval->post_windfield->turbulence[i] 	= cfg_turbulence[i];
	}

	//one last run to simulate measurements with all objects
	//~ radarfilter_initialize_resolution_volume(opc->zc, 1, &opc->res_vol, opc->todo_with_derivatives);
	//~ radarfilter_exec(opc->zc, 1, opc->todo_with_derivatives, opc->o->n, opc->o->model_radarmeasurement, opc->res_vol);
	//~ radarfilter_free_resolution_volume(opc->zc, 1, &opc->res_vol, opc->todo_with_derivatives);

	//free todo lists
	radarfilter_free_todolist(&opc->todo_no_derivatives);
	radarfilter_free_todolist(&opc->todo_with_derivatives);

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
	double *tol = NULL;
	
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
		xl[i]	= retrieval_fdvar_K2x(p, p->Klbound + i, i);
		xu[i]	= retrieval_fdvar_K2x(p, p->Kubound + i, i);
	}


	
	//What is implemented.
	if (c->use_derivatives) {
		opt = nlopt_create(NLOPT_LD_VAR1, p->Kn);
		printf("nlopt created with derivatives\n"); fflush(stdout);
	} else {
		opt = nlopt_create(NLOPT_LN_SBPLX, p->Kn);
		printf("nlopt created without derivatives\n"); fflush(stdout);
	}
	
	
	//used untill nov 21, 2016.
	//opt = nlopt_create(NLOPT_LD_VAR1, p->Kn);


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
	//opt = nlopt_create(NLOPT_LN_SBPLX, p->Kn); 			//Minimum cost function found at 4,05836e+01 <<--	

	//with derivatives
	//opt = nlopt_create(NLOPT_LD_MMA, p->Kn); 						//Minimum cost function found at 2,44692e-01	
	//opt = nlopt_create(NLOPT_LD_SLSQP, p->Kn); 					//Minimum cost function found at 6,40579e+01
	//opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND_RESTART, p->Kn); 	//Minimum cost function found at 1,05249e+02
	//opt = nlopt_create(NLOPT_LD_VAR1, p->Kn); 					//Minimum cost function found at 1,80310e-01	<<--
	//opt = nlopt_create(NLOPT_LD_VAR2, p->Kn); 					//Minimum cost function found at 1,80323e-01	

	//opt = nlopt_create(NLOPT_LN_PRAXIS, p->Kn); // algorithm and dimensionality			//Minimum cost function found at 1,44567e+01
	

	
	

	
	nlopt_set_min_objective(opt, retrieval_fdvar_cost_function_nlopt, vd_opc);
	nlopt_set_lower_bounds(opt, xl);
	nlopt_set_upper_bounds(opt, xu);
	
	//according to manual: nlopt_set_xtol_rel(opt, 1e-4);
	nlopt_set_xtol_rel(opt, 1.e-10);
	nlopt_set_ftol_rel(opt, 1.e-10);
	nlopt_set_maxtime(opt, c->maximum_time_s);

	//Add constraint function.
	//Tolerance is a best guess
	//~ tol = calloc( p->Kn, sizeof(double));
	//~ for (i=0; i < p->Kn; i++) {
		//~ tol[i] = 1.e-4;
	//~ }
	//~ nlopt_add_inequality_mconstraint(opt, p->Kn, retrieval_fdvar_constraints,
			//~ vd_opc, NULL);	

	//~ //create local optimize for some constraints
	//~ opt_loc = nlopt_create(NLOPT_AUGLAG, p->Kn);
	//~ nlopt_set_min_objective(opt_loc, retrieval_fdvar_cost_function_nlopt, vd_opc);	
	//~ nlopt_set_local_optimizer(opt, opt_loc);
	//~ nlopt_set_xtol_rel(opt_loc, 1.e-5);
	//~ nlopt_set_ftol_rel(opt_loc, 1.e-5);
	


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
	util_safe_free(&tol);
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
	double f_dBZ_hv;
	double *fd_dBZ_hv = NULL;
	double f_dBZ_vh;
	double *fd_dBZ_vh = NULL;
	double f_dBZ_vv;
	double *fd_dBZ_vv = NULL;

	double f_dBZdr;
	double *fd_dBZdr = NULL;
	
	
	double f_Doppler_velocity_hh_ms;
	double *fd_Doppler_velocity_hh_ms = NULL;	
	
	double f_Doppler_spectral_width_hh_ms;
	double *fd_Doppler_spectral_width_hh_ms = NULL;	
	
	double f_dBlwc;
	double *fd_dBlwc = NULL;
	double f_dBm;
	double *fd_dBm = NULL;
	
	double f_u;
	double *fd_u = NULL;
	double f_v;
	double *fd_v = NULL;
	double f_w;
	double *fd_w = NULL;
	double f_edr13;
	double *fd_edr13 = NULL;

	double f_constraints;
	double *fd_constraints = NULL;
	
	//TBD 
	double f_hspeed;
	double *fd_hspeed = NULL;
	double f_hdir;
	double *fd_hdir = NULL;	
	double f_dBLdr;
	double *fd_dBLdr = NULL;
	
	double *prior_grid_hspeed, *prior_grid_hdir;
	double *post_grid_hspeed, *post_grid_hdir;
	double *tmp_der_hspeed, *tmp_der_hdir;
	
	double alpha0, alpha;
	double *griddep = NULL;
	double *delta, *Sinv_delta, tmp;
	double mytmp, mytmp2, myfactor;
	double wf_u, wf_v, wf_w;
	
	double myintegral, mydelta, myscfac;
	
	int io, ip, in, i_psd, i_par;
	int i_int;	//i for spectral intervals
	int i_w;
	
	int total_n;
	int total_n_constraints;
	
	double tmp_weight;
	double tmp_dbz_err;
	double *gridweights = NULL;
	
	t_fdvar_opc										*opc = ((t_fdvar_opc*)vd_opc);
	t_fdvar_o 										*o = opc->o;
	t_fdvar_p 										*p = opc->p;
	t_zephyros_config_retrieval_fdvar_cfg 			*c = opc->c;
    t_zephyros_config 								*zc = opc->zc;
	t_radarfilter_todolist							*todo;

	double eta_hh, eta_hv, eta_vh, eta_vv;
	t_radarfilter_res_vol 			*res_vol;

	double tot_number_density;
	double tmp_diff;
	double tmp_myfactor;
	double *tmp_weight_arr = NULL;
	double *tmp_model_nd = NULL;
	double *dummy = NULL;

	int i_par2;
	double tmp_eta_hh_p, tmp_eta_vv_p, tmp_zdr_p, tmp_delta;
	double tmp_zdr_m, tmp_eta_hh_m, tmp_eta_vv_m, tmp_p, tmp_m;
		
	int parameter_number;
	
	//Unpack K
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		printf("retrieval_fdvar_unpack_x\n"); fflush(stdout);
	#endif 
	retrieval_fdvar_unpack_x(opc,x);
	
	//set values to zero
	if (calc_costfunction) {
		*f	= 0.;	

		//observation space, seperate terms
		f_dBZ_hh 	= 0.;
		f_dBZ_hv 	= 0.;
		f_dBZ_vh 	= 0.;
		f_dBZ_vv 	= 0.;
		
		f_dBZdr 	= 0.;

		//TBD ...
		f_dBLdr = 0.;
				
		f_Doppler_velocity_hh_ms = 0.;
		f_Doppler_spectral_width_hh_ms = 0.;	
		
		//parameter space, seperate terms
		f_dBlwc		= 0.;
		f_dBm		= 0.;
		f_u			= 0.;
		f_v			= 0.;
		f_w			= 0.;
		f_hspeed	= 0.;
		f_hdir		= 0.;
		f_edr13		= 0.;
		
		//~ //constraints
		f_constraints = 0.;
	}
	if (1) {
		//observation space
		fd_dBZ_hh 	= calloc(*n, sizeof(double));
		fd_dBZ_hv 	= calloc(*n, sizeof(double));
		fd_dBZ_vh 	= calloc(*n, sizeof(double));
		fd_dBZ_vv 	= calloc(*n, sizeof(double));
		fd_dBZdr 	= calloc(*n, sizeof(double));
		
		//TBD
		fd_dBLdr = calloc(*n, sizeof(double));
		
		fd_Doppler_velocity_hh_ms = calloc(*n, sizeof(double));
		fd_Doppler_spectral_width_hh_ms = calloc(*n, sizeof(double));
		
		//parameter space
		fd_dBlwc	= calloc(*n, sizeof(double));
		fd_dBm		= calloc(*n, sizeof(double));
		fd_u		= calloc(*n, sizeof(double));
		fd_v		= calloc(*n, sizeof(double));
		fd_w		= calloc(*n, sizeof(double));
		fd_hspeed	= calloc(*n, sizeof(double));
		fd_hdir 	= calloc(*n, sizeof(double));
		fd_edr13	= calloc(*n, sizeof(double));
		
		//constraints
		fd_constraints	= calloc(*n, sizeof(double));	
		
		for (in=0; in < *n; in++) 
		{
			if (calc_costfunction_derivatives) fd[in] = 0.;
			
			//observation space
			fd_dBZ_hh[in]	= 0.;
			fd_dBZ_hv[in]	= 0.;
			fd_dBZ_vh[in]	= 0.;
			fd_dBZ_vv[in]	= 0.;
			
			fd_dBZdr[in] 	= 0.;

			//TBD
			fd_dBLdr[in] 	=	0.;
						
			fd_Doppler_velocity_hh_ms[in] =0.;
			fd_Doppler_spectral_width_hh_ms[in] =0.;			
			
			//parameter space 
			fd_dBlwc[in]	= 0.;
			fd_dBm[in]		= 0.;
			fd_u[in]		= 0.;
			fd_v[in]		= 0.;
			fd_w[in]		= 0.;
			fd_hspeed[in]	= 0.;
			fd_hdir[in]		= 0.;
			fd_edr13[in]	= 0.;
			
			//~ //constraints
			fd_constraints[in] = 0.;
		}
	}
		
	if (calc_costfunction_derivatives) {
		todo = opc->todo_with_derivatives;
	} else {
		todo = opc->todo_no_derivatives;
	}

	//calculate model
	res_vol = opc->res_vol;
	radarfilter_exec(zc, 1, todo, o->n, o->model_radarmeasurement, res_vol);
    
	if (
		(c->costfunction_dBZ_vh) | 
		(c->costfunction_dBLdr) | 
		(c->costfunction_Doppler_spectrum_dBZ_hh) | 
		(c->costfunction_specific_dBZdr) | 
		(c->costfunction_specific_dBLdr)
		)
	{
		printf("cost function for this configuration not implemented yet.\n");
		fflush(stdout); exit(0);
	}
    

    //***
    //Observation space
    //***
    
	#ifdef _ZEPHYROS_FDVAR_DEBUG
	printf("Observation space: dBZ_hh, dBZ_hv, dBZ_vh, dBZ_vv, dBZdr\n"); fflush(stdout);
	#endif		
	for (io=0; io < o->n; io++) 
	{
		//observation space of cost function
		if (calc_costfunction) {			
			if (c->costfunction_dBZ_hh) f_dBZ_hh += pow((o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) / o->radarmeasurement[io]->dBZ_hh_err, 2.);
			if (c->costfunction_dBZ_hv) f_dBZ_hv += pow((o->model_radarmeasurement[io]->dBZ_hv - o->radarmeasurement[io]->dBZ_hv) / o->radarmeasurement[io]->dBZ_hv_err, 2.);
			if (c->costfunction_dBZ_vh) f_dBZ_vh += pow((o->model_radarmeasurement[io]->dBZ_vh - o->radarmeasurement[io]->dBZ_vh) / o->radarmeasurement[io]->dBZ_vh_err, 2.);
			if (c->costfunction_dBZ_vv) f_dBZ_vv += pow((o->model_radarmeasurement[io]->dBZ_vv - o->radarmeasurement[io]->dBZ_vv) / o->radarmeasurement[io]->dBZ_vv_err, 2.);
			if (c->costfunction_dBZdr) 	f_dBZdr  += pow((o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) / o->radarmeasurement[io]->dBZdr_err, 2.);
		}

		if (calc_costfunction_derivatives) {
			//Go through all objects that are fitted. (1. Scatterers, 2. Windfields, 3. Turbulence)			
			
			//1. Scatterers
			for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
				if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {
					//dBlwc
					if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
						griddep = calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_scattererfield->psd[i_psd]->lut_ln_number_density_m3[0],
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							griddep);
						
						//dBZ_hh
						if (c->costfunction_dBZ_hh) {
							alpha = 2. * (o->model_radarmeasurement[io]->dBZ_hh - o->radarmeasurement[io]->dBZ_hh) * pow(o->radarmeasurement[io]->dBZ_hh_err , -2.);
							for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
								//calculate L/Z * sum(sigma_i / alpha_i)
								//alpha = LWC / N
								//first part: sum(sigma_i /alpha_i)
								tmp_myfactor = 0.;
								for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
									tmp_myfactor += 
										(o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip])	 //sigma_i
										/ 
										(
										1.e3 *			//g kg-1
										998.2071 * 		//[water density in kg m-3]
										(4./3.) * M_PI * pow((zc->retrieval->post_scattererfield->psd[i_psd]->discrete_D_equiv_mm[i_par]/2.) * 1.e-3,3.) //m3 of single droplet
										);

								}
								//second part: L/Z
								tmp_myfactor *= (zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[ip] * o->model_radarmeasurement[io]->coef_eta2Z / func_dB_inv(o->model_radarmeasurement[io]->dBZ_hh));

								fd_dBZ_hh[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
									 alpha * griddep[ip] * tmp_myfactor;
							}
						}
						

						//dBZ_hv
						if (c->costfunction_dBZ_hv) {
							alpha = 2. * (o->model_radarmeasurement[io]->dBZ_hv - o->radarmeasurement[io]->dBZ_hv) * pow(o->radarmeasurement[io]->dBZ_hv_err , -2.);
							for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
								tmp_myfactor = 0.;
								for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
									tmp_myfactor += 
										(o->model_radarmeasurement[io]->eta_i_hv[i_psd][i_par] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip])	 //sigma_i
										/ 
										(
										1.e3 *			//g kg-1
										998.2071 * 		//[water density in kg m-3]
										(4./3.) * M_PI * pow((zc->retrieval->post_scattererfield->psd[i_psd]->discrete_D_equiv_mm[i_par]/2.) * 1.e-3,3.) //m3 of single droplet
										);

								}
								//second part: L/Z
								tmp_myfactor *= (zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[ip] * o->model_radarmeasurement[io]->coef_eta2Z / func_dB_inv(o->model_radarmeasurement[io]->dBZ_hv));

								fd_dBZ_hv[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
									 alpha * griddep[ip] * tmp_myfactor;
							}
						}
						
						/*
						//dBZ_vh
						if (c->costfunction_dBZ_vh) {
							alpha = 2. * (o->model_radarmeasurement[io]->dBZ_vh - o->radarmeasurement[io]->dBZ_vh) * pow(o->radarmeasurement[io]->dBZ_vh_err , -2.);
							if (isnanorinf(&alpha) == 0) {
								for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
									fd_dBZ_vh[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
										 alpha * griddep[ip];
								}
							}
						}
						*/
												
						//dBZ_vv
						if (c->costfunction_dBZ_vv) {
							alpha = 2. * (o->model_radarmeasurement[io]->dBZ_vv - o->radarmeasurement[io]->dBZ_vv) * pow(o->radarmeasurement[io]->dBZ_vv_err , -2.);
							for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
								tmp_myfactor = 0.;
								for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
									tmp_myfactor += 
										(o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip])	 //sigma_i
										/ 
										(
										1.e3 *			//g kg-1
										998.2071 * 		//[water density in kg m-3]
										(4./3.) * M_PI * pow((zc->retrieval->post_scattererfield->psd[i_psd]->discrete_D_equiv_mm[i_par]/2.) * 1.e-3,3.) //m3 of single droplet
										);

								}
								//second part: L/Z
								tmp_myfactor *= (zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[ip] * o->model_radarmeasurement[io]->coef_eta2Z / func_dB_inv(o->model_radarmeasurement[io]->dBZ_vv));

								fd_dBZ_vv[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
									 alpha * griddep[ip] * tmp_myfactor;
							}


						}


						util_safe_free(&griddep);
					}
					
					//dBm
					if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm) {
						griddep = calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, sizeof(double));
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
								for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
									fd_dBZ_hh[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
										+ (zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ ip] += 
										 alpha * griddep[ip];
								}
							}
						}
						
						//dBZ_hv
						if (c->costfunction_dBZ_hv) {
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBZ_hv - o->radarmeasurement[io]->dBZ_hv) * pow(o->radarmeasurement[io]->dBZ_hv_err, -2.);
							eta_hv = func_dB_inv(o->model_radarmeasurement[io]->dBZ_hv) / o->model_radarmeasurement[io]->coef_eta2Z;

							for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
								alpha = alpha0 * (o->model_radarmeasurement[io]->eta_i_hv[i_psd][i_par] / eta_hv);								
								for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
									fd_dBZ_hv[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
										+ (zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ ip] += 
										 alpha * griddep[ip];
								}
							}
						}
						
						/*
						//dBZ_vh
						if (c->costfunction_dBZ_vh) {
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBZ_vh - o->radarmeasurement[io]->dBZ_vh) * pow(o->radarmeasurement[io]->dBZ_vh_err, -2.);
							eta_vh = func_dB_inv(o->model_radarmeasurement[io]->dBZ_vh) / o->model_radarmeasurement[io]->coef_eta2Z;

							for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
								alpha = alpha0 * (o->model_radarmeasurement[io]->eta_i_vh[i_psd][i_par] / eta_vh);								
								for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
									fd_dBZ_vh[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
										+ (zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ ip] += 
										 alpha * griddep[ip];
								}
							}
						}
						*/
						
						//dBZ_vv
						if (c->costfunction_dBZ_vv) {
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBZ_vv - o->radarmeasurement[io]->dBZ_vv) * pow(o->radarmeasurement[io]->dBZ_vv_err, -2.);
							eta_vv = func_dB_inv(o->model_radarmeasurement[io]->dBZ_vv) / o->model_radarmeasurement[io]->coef_eta2Z;

							for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
								alpha = alpha0 * (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par] / eta_vv);								
								for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
									fd_dBZ_vv[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
										+ (zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ ip] += 
										 alpha * griddep[ip];
								}
							}
						}
						
						//dBZdr
						if (c->costfunction_dBZdr) {							
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) * pow(o->radarmeasurement[io]->dBZdr_err, -2.);
							//~ eta_hh = func_dB_inv(o->model_radarmeasurement[io]->dBZ_hh) / o->model_radarmeasurement[io]->coef_eta2Z;
							//~ eta_vv = func_dB_inv(o->model_radarmeasurement[io]->dBZ_vv) / o->model_radarmeasurement[io]->coef_eta2Z;
							//dBZdr = dBZ_hh - dBZ_vv

							for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {

								for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
									//~ //calculation by finite differences
									//~ tmp_eta_hh_p 	= 0.;
									//~ tmp_eta_vv_p 	= 0.;
									//~ tmp_eta_hh_m 	= 0.;
									//~ tmp_eta_vv_m 	= 0.;

									//~ tmp_p			= 1.e-40;
									//~ tmp_m			= 1.e-40;
									
									//~ for ( i_par2 = 0; i_par2 < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par2++ ) {
										//~ if (i_par == i_par2) {
											//~ tmp_eta_hh_p += tmp_p * (o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par2] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par2][ip]);
											//~ tmp_eta_vv_p += tmp_p * (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par2] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par2][ip]);	 //sigma_i
											//~ tmp_eta_hh_m += tmp_m * (o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par2] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par2][ip]);
											//~ tmp_eta_vv_m += tmp_m * (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par2] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par2][ip]);	 //sigma_i
										//~ } 
										//~ tmp_eta_hh_p += o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par2];
										//~ tmp_eta_vv_p += o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par2];											
										//~ tmp_eta_hh_m += o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par2];
										//~ tmp_eta_vv_m += o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par2];											
									
										//~ if (i_par == i_par2) {
											//~ tmp_eta_hh_p += tmp_delta * (o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par2] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par2][ip]);	 //sigma_i
											//~ tmp_eta_vv_p += tmp_delta * (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par2] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par2][ip]);	 //sigma_i
										//~ } else {
											//~ tmp_eta_hh_p -= 1.e-5 * (o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip]);	 //sigma_i
											//~ tmp_eta_vv_p -= 1.e-5 * (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par] / zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip]);	 //sigma_i
										//~ }
									//~ }
									//~ tmp_zdr_p = tmp_eta_hh_p / tmp_eta_vv_p;
									//~ tmp_zdr_m = tmp_eta_hh_m / tmp_eta_vv_m;
									//~ alpha = alpha0 * 
										//~ (func_dB(tmp_zdr_p) - func_dB(tmp_zdr_m)) /
										//~ (func_dB(o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] + tmp_p) - func_dB(o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] + tmp_m) );



									alpha = alpha0 * 
										((o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] / eta_hh)
										- (o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par] / eta_vv));
									fd_dBZdr[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
										+ (zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ ip] += 
										 alpha * griddep[ip];
								}
							}	
						}

						
						util_safe_free(&griddep);
					}
				}
			}
		
			//2. Wind fields				
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
							alpha = alpha * o->model_radarmeasurement[io]->der_dBZ_hh[0] * (-1. * o->model_radarmeasurement[io]->center_air_u);
							for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
								fd_dBZ_hh[zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] += 
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
							alpha = alpha * o->model_radarmeasurement[io]->der_dBZ_hh[1] * (-1. * o->model_radarmeasurement[io]->center_air_v);
							for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
								fd_dBZ_hh[zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] += 
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
							alpha = alpha * o->model_radarmeasurement[io]->der_dBZ_hh[2] * (-1. * o->model_radarmeasurement[io]->center_air_w);
							for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
								fd_dBZ_hh[zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] += 
									 alpha * griddep[ip];
							}
						}
					}
				}
			}
			
			//loop over turbulence fields
			/*
			not implemented yet ...
			for ( i_w = 0; i_w <= 100; i_w++ ) {
				if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
					if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
						griddep = malloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n * sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_windfield->turbulence[i_w]->lut_edr13,
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							griddep);


						//also dBZ_hh here.
						if (c->costfunction_dBZdr) {
							alpha0 = 2. * (o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) * pow(o->radarmeasurement[io]->dBZdr_err, -2.);
							alpha = alpha0 * o->model_radarmeasurement[io]->der_edr13_dBZdr;
							
							for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
								fd_dBZdr[
									zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
									 alpha * griddep[ip];
							}

						}


						
						//TBD, other variables that depends on edr13
						free(griddep);						
						
					}
				}
			}
			*/
		}		
	}


	#ifdef _ZEPHYROS_FDVAR_DEBUG
	printf("Observation space: Doppler_velocity_hh_ms\n"); fflush(stdout);
	#endif		
	for (io=0; io < o->n; io++) 
	{
		//observation space of cost function
		if (calc_costfunction) {	
			if (c->costfunction_Doppler_velocity_hh_ms) f_Doppler_velocity_hh_ms += pow((o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) / o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, 2.);
		}

		if (calc_costfunction_derivatives) {
			//Go through all objects that are fitted. (1. Scatterers, 2. Windfields, 3. Turbulence)			
			//1. Scatterers
			//Will not influence Doppler velocities
			
			//2. Windfields
			//TBD: it would be better to use reflectivity weighted coordinates
			
			if (c->costfunction_Doppler_velocity_hh_ms) {			
				for ( i_w = 0; i_w <= 100; i_w++ ) {
					if (zc->retrieval->post_windfield->field[i_w] != NULL) {							
						//grid_u
						if (zc->retrieval->prior_windfield->fit_u[i_w]) {
							griddep = calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
							interpolation_bilint_griddep(
								zc->retrieval->prior_windfield->lut_u[i_w],
								o->model_radarmeasurement[io]->advected_center_enu_xyzt,
								griddep);

							if (c->costfunction_Doppler_velocity_hh_ms) {
								alpha = 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) * pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);
								alpha *= o->model_radarmeasurement[io]->center_coor->radar_enu_dir[0];
								for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
									fd_Doppler_velocity_hh_ms[
										zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] += 
										 alpha * griddep[ip];
								}
							}
							util_safe_free(&griddep);
						}
						
						//grid_v					
						if (zc->retrieval->prior_windfield->fit_v[i_w]) {
							griddep = calloc(zc->retrieval->prior_windfield->field[i_w]->n , sizeof(double));
							interpolation_bilint_griddep(
								zc->retrieval->prior_windfield->lut_v[i_w],
								o->model_radarmeasurement[io]->advected_center_enu_xyzt,
								griddep);

							if (c->costfunction_Doppler_velocity_hh_ms) {
								alpha = 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) * pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);
								alpha *= o->model_radarmeasurement[io]->center_coor->radar_enu_dir[1];
								for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
									fd_Doppler_velocity_hh_ms[
										zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] += 
										 alpha * griddep[ip];
								}
							}
							util_safe_free(&griddep);
						}
						
						//grid_w					
						if (zc->retrieval->prior_windfield->fit_w[i_w]) {
							griddep = calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
							interpolation_bilint_griddep(
								zc->retrieval->prior_windfield->lut_w[i_w],
								o->model_radarmeasurement[io]->advected_center_enu_xyzt,
								griddep);

							if (c->costfunction_Doppler_velocity_hh_ms) {
								alpha = 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) * pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);
								alpha *= o->model_radarmeasurement[io]->center_coor->radar_enu_dir[2];
								for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
									fd_Doppler_velocity_hh_ms[
										zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] += 
										 alpha * griddep[ip];
								}
							}
							util_safe_free(&griddep);
						}
						
						if (zc->retrieval->prior_windfield->fit_hspeed_hdir[i_w]) {
							griddep = calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
							interpolation_bilint_griddep(
								zc->retrieval->prior_windfield->lut_u[i_w],
								o->model_radarmeasurement[io]->advected_center_enu_xyzt,
								griddep);

							//grid_hspeed
							if (c->costfunction_Doppler_velocity_hh_ms) {
								alpha0 = 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) * pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);
								for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
									alpha = alpha0 * (
										(o->model_radarmeasurement[io]->center_coor->radar_enu_dir[0] * -1. * sin(zc->retrieval->post_windfield->grid_hdir[i_w][ip])) +
										(o->model_radarmeasurement[io]->center_coor->radar_enu_dir[1] * -1. * cos(zc->retrieval->post_windfield->grid_hdir[i_w][ip]))
										);
									fd_Doppler_velocity_hh_ms[
										zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip] += 
										 alpha * griddep[ip];
								}
							}
							
							//grid_hdir
							if (c->costfunction_Doppler_velocity_hh_ms) {
								alpha0 = 2. * (o->model_radarmeasurement[io]->Doppler_velocity_hh_ms - o->radarmeasurement[io]->Doppler_velocity_hh_ms) * pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err, -2.);
								for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
									alpha = alpha0 * (
										(o->model_radarmeasurement[io]->center_coor->radar_enu_dir[0] * -1. * zc->retrieval->post_windfield->grid_hspeed[i_w][ip] * cos(zc->retrieval->post_windfield->grid_hdir[i_w][ip])) +
										(o->model_radarmeasurement[io]->center_coor->radar_enu_dir[1] * zc->retrieval->post_windfield->grid_hspeed[i_w][ip] * sin(zc->retrieval->post_windfield->grid_hdir[i_w][ip]))
										);
									fd_Doppler_velocity_hh_ms[
										zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip] += 
										 alpha * griddep[ip];
								}
							}
							
							util_safe_free(&griddep);
						}
					}
				}
			}
		}
	}

	#ifdef _ZEPHYROS_FDVAR_DEBUG
	printf("Observation space: Doppler_spectral_width_hh_ms\n"); fflush(stdout);
	#endif	
	for (io=0; io < o->n; io++) 
	{
		//observation space of cost function
		if (calc_costfunction) {			
			if (c->costfunction_Doppler_spectral_width_hh_ms) f_Doppler_spectral_width_hh_ms += pow((o->model_radarmeasurement[io]->Doppler_spectral_width_hh_ms - o->radarmeasurement[io]->Doppler_spectral_width_hh_ms) / o->radarmeasurement[io]->Doppler_spectral_width_hh_ms_err, 2.);
		}

		if (calc_costfunction_derivatives) {
			if (c->costfunction_Doppler_spectral_width_hh_ms) {
				//loop over turbulence fields
				for ( i_w = 0; i_w <= 100; i_w++ ) {
					if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
						if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {						
							griddep = calloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n, sizeof(double));
							interpolation_bilint_griddep(
								zc->retrieval->prior_windfield->turbulence[i_w]->lut_edr13,
								o->model_radarmeasurement[io]->advected_center_enu_xyzt,
								griddep);

							//grid_edr13
							if (c->costfunction_Doppler_spectral_width_hh_ms) {
								alpha = 2. * (o->model_radarmeasurement[io]->Doppler_spectral_width_hh_ms - o->radarmeasurement[io]->Doppler_spectral_width_hh_ms) * pow(o->radarmeasurement[io]->Doppler_spectral_width_hh_ms_err, -2.);							
								for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
									fd_Doppler_spectral_width_hh_ms[
										zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
										 o->model_radarmeasurement[io]->der_edr13_Doppler_spectral_width_hh_ms *
										 alpha * griddep[ip];
								}
							}
							
							util_safe_free(&griddep);													
						}
					}
				}
			}
		}
	}
	






	/*
	{
		//observation space of cost function
		if (calc_costfunction) {			
			if (c->costfunction_dBZdr) f_dBZdr += pow((o->model_radarmeasurement[io]->dBZdr - o->radarmeasurement[io]->dBZdr) / o->radarmeasurement[io]->dBZdr_err, 2.);
			if (c->costfunction_dBLdr) f_dBLdr += pow((o->model_radarmeasurement[io]->dBLdr - o->radarmeasurement[io]->dBLdr) / o->radarmeasurement[io]->dBLdr_err, 2.);
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
				

			
			//loop over turbulence fields
			for ( i_w = 0; i_w <= 100; i_w++ ) {
				if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
					if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
						griddep = malloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n * sizeof(double));
						interpolation_bilint_griddep(
							zc->retrieval->prior_windfield->turbulence[i_w]->lut_edr13,
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							griddep);


						
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

	*/
			
			
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
		//}


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
	//}

    //End of observation space
    //***




	//#ifdef _ZEPHYROS_FDVAR_DEBUG
	//if (calc_costfunction) {printf("current cost value: %.5e (%.2e%)\n"   , *f, 100.*(*f - prev_costfunction)/prev_costfunction);fflush(stdout);}	
	//#endif





    //***
    //Parameter space
    //***
	//~ #ifdef _ZEPHYROS_FDVAR_DEBUG
	//~ printf("Parameter space: dBlwc, dBN\n"); fflush(stdout);
	//~ #endif	
	//~ for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		//~ if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {		
			//~ //dBlwc
			//~ if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				//~ delta 		= calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, sizeof(double));
				//~ Sinv_delta 	= calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, sizeof(double));

				//~ //delta is post - prior value
				//~ for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
					//~ delta[ip] = 
						//~ func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]) -
							//~ func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]);
				//~ }

				//~ retrieval_fdvar_cs_qrsol(zc->retrieval->prior_scattererfield->psd[i_psd]->dBlwc_ecm, delta, Sinv_delta);							
				
				//~ dbldotprod(&zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, delta, Sinv_delta, &tmp);
																	
				//~ if (calc_costfunction) f_dBlwc += tmp;
				
				//~ if (calc_costfunction_derivatives) {
					//~ for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {						
						//~ fd_dBlwc[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] += 
							//~ 2. * Sinv_delta[ip];
					//~ }
				//~ }
				
				//~ util_safe_free(&delta);
				//~ util_safe_free(&Sinv_delta);	
			//~ }

			//~ //dBN
			//~ if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) {			
				//~ delta 		= calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, sizeof(double));
				//~ Sinv_delta 	= calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, sizeof(double));

				//~ for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
					//~ //delta is post - prior value
					//~ for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
						//~ delta[ip] = 
							//~ func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip]) -
								//~ func_dB(zc->retrieval->prior_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip]);
					//~ }

					//~ for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
						//~ Sinv_delta[ip] = delta[ip] * 
							//~ pow(zc->retrieval->post_scattererfield->psd[i_psd]->grid_dBnumber_density_err_m3[i_par][ip], -2.);
					//~ }
					
					//~ dbldotprod(&zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, delta, Sinv_delta, &tmp);
					
					//~ if (calc_costfunction) f_dBN += tmp;
					
					//~ if (calc_costfunction_derivatives) {
						//~ for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {						
							//~ fd_dBN[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
									//~ + (zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									//~ + ip] =
								//~ 2. * Sinv_delta[ip];
						//~ }
					//~ }	
					
					//~ //TBD: work with error covariance matrices.
					//~ //retrieval_fdvar_cs_qrsol(zc->retrieval->prior_scattererfield->psd[i_psd]->dBlwc_ecm, delta, Sinv_delta);							
					//~ //dbldotprod(&zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, delta, Sinv_delta, &tmp);
				
				//~ }
				//~ util_safe_free(&delta);
				//~ util_safe_free(&Sinv_delta);

			//~ }		
		//~ }
	//~ }

	#ifdef _ZEPHYROS_FDVAR_DEBUG
	printf("Parameter space: grid_u, grid_v, grid_w, grid_hspeed, grid_hdir\n"); fflush(stdout);
	#endif
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->field[i_w] != NULL) {
			//grid_u
			if (zc->retrieval->prior_windfield->fit_u[i_w]) {
				delta 		= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				Sinv_delta 	= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						zc->retrieval->post_windfield->grid_u[i_w][ip]
						-
						zc->retrieval->prior_windfield->grid_u[i_w][ip];						
				}

				//use the ecm
				//retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->u_ecm[i_w], delta, Sinv_delta);							
				
				//do not use the ecm
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					Sinv_delta[ip] = pow(zc->retrieval->prior_windfield->grid_u_err[i_w][ip], -2.) * delta[ip];
				}
				
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_u += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {		
						fd_u[zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] += 
							 2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);	
			}
			
			//grid_v
			if (zc->retrieval->prior_windfield->fit_v[i_w]) {
				delta 		= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				Sinv_delta 	= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						zc->retrieval->post_windfield->grid_v[i_w][ip]
						-
						zc->retrieval->prior_windfield->grid_v[i_w][ip];						
				}

				//use the ecm
				//retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->v_ecm[i_w], delta, Sinv_delta);							
				
				//do not use the ecm
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					Sinv_delta[ip] = pow(zc->retrieval->prior_windfield->grid_v_err[i_w][ip], -2.) * delta[ip];
				}
								
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_v += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
						fd_v[zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] += 
							2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);	
			}

			//grid_w
			if (zc->retrieval->prior_windfield->fit_w[i_w]) {
				delta 		= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				Sinv_delta 	= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						zc->retrieval->post_windfield->grid_w[i_w][ip]
						-
						zc->retrieval->prior_windfield->grid_w[i_w][ip];						
				}

				//use the ecm
				//retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->w_ecm[i_w], delta, Sinv_delta);							
				
				//do not use the ecm
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					Sinv_delta[ip] = pow(zc->retrieval->prior_windfield->grid_w_err[i_w][ip], -2.) * delta[ip];
				}
				
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_w += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
						fd_w[zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] += 
							2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);	
			}
			
			//grid_hspeed
			if (zc->retrieval->prior_windfield->fit_hspeed_hdir[i_w]) {
				delta 		= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				Sinv_delta 	= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						zc->retrieval->post_windfield->grid_hspeed[i_w][ip]
						-
						zc->retrieval->prior_windfield->grid_hspeed[i_w][ip];						
				}

				//use the ecm
				//retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->hspeed_ecm[i_w], delta, Sinv_delta);							
				
				//do not use the ecm
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					Sinv_delta[ip] = pow(zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip], -2.) * delta[ip];
				}				
				
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);								
											
				if (calc_costfunction) f_hspeed += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
						fd_hspeed[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip] += 
							2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);	
			}
			
			//grid_hdir
			if (zc->retrieval->prior_windfield->fit_hspeed_hdir[i_w]) {
				delta 		= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				Sinv_delta 	= calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));

				//delta is post - prior value
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					delta[ip] = 
						angleAminB(zc->retrieval->post_windfield->grid_hdir[i_w] + ip,
							zc->retrieval->prior_windfield->grid_hdir[i_w] + ip);					
				}

				//use the ecm
				//retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->hdir_ecm[i_w], delta, Sinv_delta);							

				//do not use the ecm
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					Sinv_delta[ip] = pow(zc->retrieval->prior_windfield->grid_hdir_err[i_w][ip], -2.) * delta[ip];
				}
				
				dbldotprod(&zc->retrieval->prior_windfield->field[i_w]->n, delta, Sinv_delta, &tmp);
											
				if (calc_costfunction) f_hdir += tmp;
				
				if (calc_costfunction_derivatives) {
					for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {					
						fd_hdir[zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip] += 
							2. * Sinv_delta[ip];
					}
				}
				
				free(delta);
				free(Sinv_delta);	
			}
			
		}
	}
			
			
 	
	
	//~ #ifdef _ZEPHYROS_FDVAR_DEBUG
	//~ printf("Parameter space: grid_edr13 (turbulence)\n"); fflush(stdout);
	//~ #endif
	//~ for ( i_w = 0; i_w <= 100; i_w++ ) {
		//~ if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
			//~ if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
				//~ delta 		= calloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n, sizeof(double));
				//~ Sinv_delta 	= calloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n, sizeof(double));

				//~ //delta is post - prior value
				//~ for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
					//~ delta[ip] = 
						//~ zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[ip]
						//~ -
						//~ zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13[ip];						
				//~ }

				//~ retrieval_fdvar_cs_qrsol(zc->retrieval->prior_windfield->turbulence[i_w]->edr13_ecm, delta, Sinv_delta);							
				//~ dbldotprod(&zc->retrieval->prior_windfield->turbulence[i_w]->field->n, delta, Sinv_delta, &tmp);
											
				//~ if (calc_costfunction) f_edr13 += tmp;

				//~ if (calc_costfunction_derivatives) {
					//~ for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {					
						//~ fd_edr13[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] += 
							//~ 2. * Sinv_delta[ip];
					//~ }
				//~ }
				
				//~ free(delta);
				//~ free(Sinv_delta);
			//~ }
		//~ }
	//~ }

	
    //End of parameter space
    //***
	
	

	



	//***
    //Some additional constraints
    //***
	
	/*
	#ifdef _ZEPHYROS_FDVAR_DEBUG
	printf("Some additional constraints: fit gamma distribution parameters\n"); fflush(stdout);
	#endif
	total_n_constraints = 0.;
	for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN) {
			if (
				(	(zc->retrieval->prior_scattererfield->psd[i_psd]->N_constraint == 1)
					&& (zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters > 0) ) 
				|
				(	(zc->retrieval->prior_scattererfield->psd[i_psd]->N_constraint == 2)
					&& (zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters > 1) ) 
				|
				(	(zc->retrieval->prior_scattererfield->psd[i_psd]->N_constraint == 3)
					&& (zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters > 2) )
				) {
							
				//Initialize
				tmp_model_nd = calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters , sizeof(double));
				tmp_weight_arr = calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters, sizeof(double));
				tmp_sigma_hh = calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters, sizeof(double));
				tmp_sigma_vv = calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters, sizeof(double));

				total_n_constraints += 
					zc->retrieval->post_scattererfield->psd[i_psd]->field->n
					* zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters;

				//Calculate some average tmp_sigma_hh and tmp_sigma_vv for particles
				//TBD: make averaged sigma grid dependend. Probably not very relevant ... ?			
				for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
					tmp_sigma_hh[i_par] = 0.;
					tmp_sigma_vv[i_par] = 0.;
					for (io=0; io < o->n; io++) {
						interpolation_bilint(
							zc->retrieval->prior_scattererfield->psd[i_psd]->lut_ln_number_density_m3[i_par], 
							o->model_radarmeasurement[io]->advected_center_enu_xyzt,
							&tmp,
							0,
							dummy);
						tmp = exp(tmp); // number density
						
						if (o->model_radarmeasurement[io]->eta_i_hh != NULL)
							tmp_sigma_hh[i_par] += o->model_radarmeasurement[io]->eta_i_hh[i_psd][i_par] / tmp;
						if (o->model_radarmeasurement[io]->eta_i_vv != NULL)
							tmp_sigma_vv[i_par] += o->model_radarmeasurement[io]->eta_i_vv[i_psd][i_par] / tmp;
					}
					tmp_sigma_hh[i_par] /= o->n;
					tmp_sigma_vv[i_par] /= o->n;
				}

				//set weighting function in the psd based on rcs
				//~ for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
					//~ tmp_weight_arr[i_par] = 
							//~ pow(zc->retrieval->prior_scattererfield->psd[i_psd]->discrete_D_equiv_mm[i_par], 6.);

					//~ //TBD: maybe better to use measurements somehow.

				//~ }	
																			
				//model 1: D0 fitted for the same Z_hh.
				//model 2: mu-D0 relation is used. (N0, D0) is fitted for same (Z_hh, Z_vv).
				//model 3: (D0, mu) are fitted. Not sure yet what to do...
				if (zc->retrieval->prior_scattererfield->psd[i_psd]->N_constraint) {
					util_fit_gammadistribution_parameters(zc->retrieval->post_scattererfield->psd[i_psd], 
						zc->retrieval->prior_scattererfield->psd[i_psd]->N_constraint,
						tmp_sigma_hh, tmp_sigma_vv, NULL);
				 }
							
				//add
				for (ip=0; ip < zc->retrieval->post_scattererfield->psd[i_psd]->field->n; ip++) {				
					printf("fitted_N0 = %.2e\n",
						zc->retrieval->post_scattererfield->psd[i_psd]->grid_gammadistribution_N0[ip]);
					printf("fitted D0_mm = %.2e\n",
						zc->retrieval->post_scattererfield->psd[i_psd]->grid_gammadistribution_D0_mm[ip]);
					printf("fitted mu = %.2e\n",
						zc->retrieval->post_scattererfield->psd[i_psd]->grid_gammadistribution_mu[ip]);

					//~ tot_number_density	= 0.;
					//~ for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
						//~ tot_number_density += zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip];
					//~ }
									 
					for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
						tmp_model_nd[i_par] = 
							util_gamma_integral(
							zc->retrieval->post_scattererfield->psd[i_psd]->grid_gammadistribution_N0[ip],		//N0
							zc->retrieval->post_scattererfield->psd[i_psd]->grid_gammadistribution_mu[ip], 		//mu
							zc->retrieval->post_scattererfield->psd[i_psd]->grid_gammadistribution_D0_mm[ip], 	//D0_mm
							zc->retrieval->post_scattererfield->psd[i_psd]->discrete_D_equiv_mm_interval_llimt[i_par],
							zc->retrieval->post_scattererfield->psd[i_psd]->discrete_D_equiv_mm_interval_ulimt[i_par]);
					}
					
					for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {

						//constraint in terms of (Z - Z_mod) / Z

						//constraint in terms of (n_fit - n_model) / n_tot

						//~ mydelta 	= pow(tmp_diff, 2.) * pow(tot_number_density, -2.) * tmp_weight_arr[i_par];
						//~ alpha	 	= 2. * (tmp_diff) * pow(tot_number_density, -2.) * tmp_weight_arr[i_par];
						//~ alpha		*= (log(10.) / 10.) * zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip];



						//constraint in terms of (db(n_fit) - db(n_model) * (n_model / n_tot)
						//As db is fitted, this in linear and therefore preferrable.
						//~ tmp_diff = (
							//~ zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip]
							//~ -
							//~ tmp_model_nd[i_par]);	

						tmp_dbz_err = .1;
						tmp_diff = (
							func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip])
							-
							func_dB(tmp_model_nd[i_par])) ;	

						mydelta 	= pow(tmp_diff, 2.) * pow(tmp_dbz_err, -2.);
						alpha	 	= 2. * tmp_diff * pow(tmp_dbz_err, -2.);
								
						//add to cost function
						f_constraints += mydelta;
						if (calc_costfunction_derivatives) {
							parameter_number = zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBN_Knr 
													+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
													+ip;						
							fd_constraints[parameter_number] += alpha;
						}
							
						//}
						// debugging
						//~ printf("D0 = %.2e; Dl = %.2e, Du = %.2e, f1 = %.2e, f2 = %.2e \n",
						//~ zc->retrieval->post_scattererfield->psd[i_psd]->discrete_D_equiv_mm[i_par],
						//~ zc->retrieval->post_scattererfield->psd[i_psd]->discrete_D_equiv_mm_interval_llimt[i_par],
						//~ zc->retrieval->post_scattererfield->psd[i_psd]->discrete_D_equiv_mm_interval_ulimt[i_par],
						//~ zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip] / tot_number_density,
						//~ myintegral / tot_number_density);
						//~ printf("mydelta = %.2e ; alpha = %.2e\n", mydelta, alpha);
						
					}
				}	
				util_safe_free(&tmp_model_nd);
				util_safe_free(&tmp_weight_arr);
				util_safe_free(&tmp_sigma_hh);
				util_safe_free(&tmp_sigma_vv);
			}
		}
	}
	*/
	//End of some additional constraints
    //***




	//For analysis
	//calcualte grid zetaI
	//assert that it is set.
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->post_windfield->turbulence[i_w] != NULL) {
			if (zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13) {
				//assert that zeta_I is set.
				if (zc->retrieval->post_windfield->turbulence[i_w]->grid_zetaI == NULL)
					zc->retrieval->post_windfield->turbulence[i_w]->grid_zetaI = calloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n, sizeof(double));

				gridweights = calloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n, sizeof(double));						
				for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
					zc->retrieval->post_windfield->turbulence[i_w]->grid_zetaI[ip] = 0.;
					gridweights[ip] = 0.;
				}
				
				for (io=0; io < o->n; io++) {				
					griddep = calloc(zc->retrieval->prior_windfield->turbulence[i_w]->field->n, sizeof(double));
					interpolation_bilint_griddep(
						zc->retrieval->prior_windfield->turbulence[i_w]->lut_edr13,
						o->model_radarmeasurement[io]->advected_center_enu_xyzt,
						griddep);

					for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
						zc->retrieval->post_windfield->turbulence[i_w]->grid_zetaI[ip] +=
							griddep[ip] * o->model_radarmeasurement[io]->zetaI;
						gridweights[ip] += griddep[ip];
					}

					util_safe_free(&griddep);
				}
				
				for (ip=0; ip < zc->retrieval->prior_windfield->turbulence[i_w]->field->n; ip++) {
					zc->retrieval->post_windfield->turbulence[i_w]->grid_zetaI[ip] /= gridweights[ip];
					if (gridweights[ip] == 0.) {
						zc->retrieval->post_windfield->turbulence[i_w]->grid_zetaI[ip] = -999.9;
					}
				}				
				util_safe_free(&gridweights);				
			}
		}
	}

	
	
	
		
	//Calculate total number of measurements used in the cost function
	total_n = 0;
	if (c->costfunction_dBZ_hh) total_n += o->n;
	if (c->costfunction_dBZ_hv) total_n += o->n;
	if (c->costfunction_dBZ_vh) total_n += o->n;
	if (c->costfunction_dBZ_vv) total_n += o->n;
	if (c->costfunction_dBZdr) 	total_n += o->n;
	if (c->costfunction_Doppler_velocity_hh_ms) total_n += o->n;
	if (c->costfunction_Doppler_spectral_width_hh_ms) total_n += o->n;
	
	//do some scaling
	if (calc_costfunction) {
		//scale observation space by number of measurements
		if (c->costfunction_dBZ_hh) f_dBZ_hh *= (1. / total_n);
		if (c->costfunction_dBZ_hv) f_dBZ_hv *= (1. / total_n);
		if (c->costfunction_dBZ_vh) f_dBZ_vh *= (1. / total_n);
		if (c->costfunction_dBZ_vv) f_dBZ_vv *= (1. / total_n);
		if (c->costfunction_dBZdr) 	f_dBZdr	 *= (1. / total_n);
		
		if (c->costfunction_Doppler_velocity_hh_ms) 		f_Doppler_velocity_hh_ms *= (1. / total_n);
		if (c->costfunction_Doppler_spectral_width_hh_ms) 	f_Doppler_spectral_width_hh_ms *= (1. / total_n);
		
		//scale parameter space by number of parameters
		f_dBlwc 	*= (1. / p->Kn);
		f_dBm 		*= (1. / p->Kn);
		f_u 		*= (1. / p->Kn);
		f_v 		*= (1. / p->Kn);
		f_w 		*= (1. / p->Kn);
		f_hspeed	*= (1. / p->Kn);
		f_hdir 		*= (1. / p->Kn);
		f_edr13 	*= (1. / p->Kn);
		
		//~ //scale constraints by factor 
		f_constraints *= (1. / total_n_constraints);
	}

	if (calc_costfunction_derivatives) {
		for (in=0; in < *n; in++) {
			//observation space
			if (c->costfunction_dBZ_hh) fd_dBZ_hh[in] 	*= (1. / total_n);
			if (c->costfunction_dBZ_hv) fd_dBZ_hv[in] 	*= (1. / total_n);
			if (c->costfunction_dBZ_vh) fd_dBZ_vh[in] 	*= (1. / total_n);
			if (c->costfunction_dBZ_vv) fd_dBZ_vv[in] 	*= (1. / total_n);
			if (c->costfunction_dBZdr) 	fd_dBZdr[in] 	*= (1. / total_n);
			
			if (c->costfunction_Doppler_velocity_hh_ms) fd_Doppler_velocity_hh_ms[in] 				*= (1. / total_n);
			if (c->costfunction_Doppler_spectral_width_hh_ms) fd_Doppler_spectral_width_hh_ms[in] 	*= (1. / total_n);
			
			//parameter space
			fd_dBlwc[in] 	*= (1. / p->Kn);
			fd_dBm[in]		*= (1. / p->Kn);
			fd_u[in] 		*= (1. / p->Kn);
			fd_v[in] 		*= (1. / p->Kn);
			fd_w[in] 		*= (1. / p->Kn);
			fd_hspeed[in] 	*= (1. / p->Kn);
			fd_hdir[in] 	*= (1. / p->Kn);
			fd_edr13[in] 	*= (1. / p->Kn);	
			
			//~ //constraints
			fd_constraints[in] *= (1. / total_n_constraints);					
		}
	}





	if (calc_costfunction) {
		//observation space
		if (c->costfunction_dBZ_hh) *f += f_dBZ_hh;
		if (c->costfunction_dBZ_hv) *f += f_dBZ_hv;
		if (c->costfunction_dBZ_vh) *f += f_dBZ_vh;
		if (c->costfunction_dBZ_vv) *f += f_dBZ_vv;
		if (c->costfunction_dBZdr) 	*f += f_dBZdr;
		
		if (c->costfunction_Doppler_velocity_hh_ms) 		*f += f_Doppler_velocity_hh_ms;
		if (c->costfunction_Doppler_spectral_width_hh_ms) 	*f += f_Doppler_spectral_width_hh_ms;
		
		//parameter space
		*f += f_dBlwc;
		*f += f_dBm;
		*f += f_u;
		*f += f_v;
		*f += f_w;
		*f += f_hspeed;
		*f += f_hdir;
		*f += f_edr13;
		
		//~ //constraints
		if (total_n_constraints != 0)
			*f += f_constraints;
		
		
		/*
		*f += weight_dBZdr * f_dBZdr;
		*f += weight_dBLdr * f_dBLdr;
		*f += weight_Doppler_velocity_hh_ms * f_Doppler_velocity_hh_ms;
		*f += weight_Doppler_spectral_width_hh_ms * f_Doppler_spectral_width_hh_ms;
		*/
	}
	
	if (calc_costfunction_derivatives) {
		for (in=0; in < *n; in++) {
			
			if (c->costfunction_dBZ_hh) 						fd[in] += fd_dBZ_hh[in];
			if (c->costfunction_dBZ_hv) 						fd[in] += fd_dBZ_hv[in];
			if (c->costfunction_dBZ_vh) 						fd[in] += fd_dBZ_vh[in];
			if (c->costfunction_dBZ_vv) 						fd[in] += fd_dBZ_vv[in];
			if (c->costfunction_dBZdr)	 						fd[in] += fd_dBZdr[in];
			if (c->costfunction_Doppler_velocity_hh_ms) 		fd[in] += fd_Doppler_velocity_hh_ms[in];
			if (c->costfunction_Doppler_spectral_width_hh_ms) 	fd[in] += fd_Doppler_spectral_width_hh_ms[in];
		
			//parameter space
			fd[in] += fd_dBlwc[in];
			fd[in] += fd_dBm[in];
			fd[in] += fd_u[in];
			fd[in] += fd_v[in];
			fd[in] += fd_w[in];
			fd[in] += fd_hspeed[in];
			fd[in] += fd_hdir[in];
			fd[in] += fd_edr13[in];

			//~ //constraints
			if (total_n_constraints != 0)
				fd[in] += fd_constraints[in];
			
			/*
			fd[in] += weight_dBZdr * fd_dBZdr[in];
			fd[in] += weight_dBLdr * fd_dBLdr[in];
			fd[in] += weight_Doppler_velocity_hh_ms * fd_Doppler_velocity_hh_ms[in];
			fd[in] += weight_Doppler_spectral_width_hh_ms * fd_Doppler_spectral_width_hh_ms[in];
			*/
		}
	}

	if (calc_costfunction_derivatives) {
		printf("norm fd = %.2e\n", func_norm(n, fd));
		printf("info: %s %i \n", __FILE__, __LINE__); fflush(stdout);		
	}

	/*
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
	*/
	

	
	
	#ifdef _ZEPHYROS_FDVAR_DEBUG
		if (calc_costfunction_derivatives) {

			printf("current cost function parts:\n");
			printf("f_dBZ_hh = %.2e							, norm der = %.2e\n", f_dBZ_hh, func_norm(n, fd_dBZ_hh));
			//printf("f_dBZ_hv = %.2e							, norm der = %.2e\n", f_dBZ_hv, func_norm(n, fd_dBZ_hv));
			//printf("f_dBZ_vh = %.2e							, norm der = %.2e\n", f_dBZ_vh, func_norm(n, fd_dBZ_vh));
			printf("f_dBZ_vv = %.2e							, norm der = %.2e\n", f_dBZ_vv, func_norm(n, fd_dBZ_vv));
			printf("f_dBZdr = %.2e							, norm der = %.2e\n", f_dBZdr, func_norm(n, fd_dBZdr));
			printf("f_Doppler_velocity_hh_ms = %.2e			, norm der = %.2e\n", f_Doppler_velocity_hh_ms, func_norm(n, fd_Doppler_velocity_hh_ms));
			printf("f_Doppler_spectral_width_hh_ms = %.2e	, norm der = %.2e\n", f_Doppler_spectral_width_hh_ms, func_norm(n, fd_Doppler_spectral_width_hh_ms));
			
			printf("\n");
			//printf("f_dBlwc = %.2e						, norm der = %.2e\n", f_dBlwc, func_norm(n, fd_dBlwc));
			printf("f_dBm = %.2e							, norm der = %.2e\n", f_dBm, func_norm(n, fd_dBm));
			printf("f_u = %.2e							, norm der = %.2e\n", f_u, func_norm(n, fd_u));
			printf("f_v = %.2e							, norm der = %.2e\n", f_v, func_norm(n, fd_v));
			printf("f_w = %.2e							, norm der = %.2e\n", f_w, func_norm(n, fd_w));
			printf("f_hspeed = %.2e							, norm der = %.2e\n", f_hspeed, func_norm(n, fd_hspeed));
			printf("f_hdir = %.2e							, norm der = %.2e\n", f_hdir, func_norm(n, fd_hdir));
			printf("f_edr13 = %.2e						, norm der = %.2e\n", f_edr13, func_norm(n, fd_edr13));
			
			printf("f_constraints = %.2e						, norm der = %.2e\n", f_constraints, func_norm(n, fd_constraints));
		}

		//~ for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
			//~ if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {	
				//~ for ( i_par = 0; i_par < zc->retrieval->post_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {

					//~ printf("\nD_mm = %.2e\n",
						//~ zc->retrieval->post_scattererfield->psd[i_psd]->discrete_D_equiv_mm[i_par]);
					//~ printf("grid_rcs_hh[i_par]/grid_rcs_vv[i_par] = %.2e (= %.2e dB)\n",
						//~ zc->retrieval->post_scattererfield->psd[i_psd]->grid_rcs_hh[i_par][0]
						//~ / zc->retrieval->post_scattererfield->psd[i_psd]->grid_rcs_vv[i_par][0],
						//~ func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_rcs_hh[i_par][0]
						//~ / zc->retrieval->post_scattererfield->psd[i_psd]->grid_rcs_vv[i_par][0])
						//~ );
				//~ }
			//~ }
		//~ }

		/*
		for (io=0; io < o->n; io++) 
			printf("o->model_radarmeasurement[%i]->der_edr13_Doppler_spectral_width_hh_ms = %.2e\n",
			io, o->model_radarmeasurement[io]->der_edr13_Doppler_spectral_width_hh_ms);
		*/
		
		//~ if (calc_costfunction_derivatives) {
			//~ printf("current cost derivatives:\n");
			//~ printf("%i derivatives:\n", *n);
			//~ for (in=0; in < *n; in++) {
				//~ if (fd[in] != 0.) {
					//~ printf("fd[%i] = %.5e \n", in, fd[in]); fflush(stdout);
				//~ }
			//~ }		
		//~ }


		//~ if (calc_costfunction_derivatives) {
			//~ printf("current cost derivatives:\n");
			//~ printf("%i derivatives:\n", *n);
			//~ for (in=0; in < *n; in++) {
				//~ if (fd[in] != 0.) {
					//~ printf("fd_dBZdr[%i] = %.5e \n", in, fd_dBZdr[in]); fflush(stdout);
				//~ }
			//~ }		
		//~ }

	#endif
	
	
	
	
	
	
	
	
	
	//~ //account for derivative in cost function to x
	//~ //Hence dK/dx has to be added. dK/dx = (u - l)
	//~ if (calc_costfunction_derivatives) {
		//~ for (in=0; in < *n; in++) {
			//~ fd[in] 	*= (p->Kubound[in] - p->Klbound[in]);
		//~ }
	//~ }
	
	if (calc_costfunction) {printf("current cost value: %.5e (%.2e%)\n"   , *f, 100.*(*f - prev_costfunction)/prev_costfunction);fflush(stdout);}
	

	
	prev_costfunction = *f;
	
	util_safe_free(&fd_dBZ_hh);
	util_safe_free(&fd_dBZ_hv);
	util_safe_free(&fd_dBZ_vh);
	util_safe_free(&fd_dBZ_vv);
	util_safe_free(&fd_dBZdr);
	util_safe_free(&fd_Doppler_velocity_hh_ms);
	util_safe_free(&fd_Doppler_spectral_width_hh_ms);
	
	util_safe_free(&fd_dBlwc);
	util_safe_free(&fd_dBm);
	util_safe_free(&fd_u);
	util_safe_free(&fd_v);
	util_safe_free(&fd_w);
	util_safe_free(&fd_hspeed);
	util_safe_free(&fd_hdir);
	util_safe_free(&fd_edr13);
	
	util_safe_free(&fd_dBLdr);
	util_safe_free(&tmp_weight_arr);
}




	






double retrieval_fdvar_K2x(t_fdvar_p *p, double *K, int Knr)
{
	//rescaling
	//return (*K - p->Klbound[Knr]) / (p->Kubound[Knr] - p->Klbound[Knr]);

	//no rescaling
	return *K;
}

double retrieval_fdvar_x2K(t_fdvar_p *p, double *x, int Knr)
{
	//~ //rescaling
	//~ return p->Klbound[Knr] + (p->Kubound[Knr] - p->Klbound[Knr]) * *x;

	//no rescaling
	return *x;
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

	//walk through psds
	for ( i_psd = 0; i_psd < zc->retrieval->prior_scattererfield->npsd; i_psd++ ) {
		if (zc->retrieval->post_scattererfield->psd[i_psd] != NULL) {
			//dBlwc
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc) {
				for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
					p->Klbound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] =
						func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]) -
						(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3[ip]);						
					p->Kubound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] =
						func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]) +
						(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBlwc_err_gm3[ip]);
					if (p->Klbound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] < -100.)
						p->Klbound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] = -100.;
					if (p->Kubound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] > 100.)
						p->Kubound[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip] = 100.;
					tmp = func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_lwc_gm3[ip]);
					(*x)[zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr	+ip] =
						K2x(p, &tmp, zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBlwc_Knr +ip);
				}				
			}
			//dBm
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm) {
				for ( i_par = 0; i_par < zc->retrieval->prior_scattererfield->psd[i_psd]->n_diameters; i_par++ ) {
					for (ip=0; ip < zc->retrieval->prior_scattererfield->psd[i_psd]->field->n; ip++) {
						p->Klbound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
									+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									+ip] =
							(func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip] * zc->retrieval->post_scattererfield->psd[i_psd]->grid_rcs_hh[i_par][ip]) -
							(factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBnumber_density_err_m3[i_par][ip]));
						if (p->Klbound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
										+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ip] < -100.)
							p->Klbound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
										+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ip] = -100.;
								
						p->Kubound[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr
									+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									+ip] =
							func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip] * zc->retrieval->post_scattererfield->psd[i_psd]->grid_rcs_hh[i_par][ip]) +
							factor * zc->retrieval->prior_scattererfield->psd[i_psd]->grid_dBnumber_density_err_m3[i_par][ip];
						if (p->Kubound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
										+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ip] > 100.)
							p->Kubound[ zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
										+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
										+ip] = 100.;
							
							
						tmp = func_dB(zc->retrieval->post_scattererfield->psd[i_psd]->grid_number_density_m3[i_par][ip] * zc->retrieval->post_scattererfield->psd[i_psd]->grid_rcs_hh[i_par][ip]);
						(*x)[	zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr
								+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
								+ip] =
							K2x(p, &tmp, zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr +
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
						zc->retrieval->post_windfield->grid_u[i_w][ip] -
						(factor * zc->retrieval->prior_windfield->grid_u_err[i_w][ip]);
					p->Kubound[zc->retrieval->prior_windfield->fit_u_Knr[i_w] + ip] =
						zc->retrieval->post_windfield->grid_u[i_w][ip] +
						(factor * zc->retrieval->prior_windfield->grid_u_err[i_w][ip]);				
					(*x)[zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip] =
						K2x(p, &(zc->retrieval->post_windfield->grid_u[i_w][ip]), zc->retrieval->prior_windfield->fit_u_Knr[i_w] +ip);
				}
			}
			
			//grid_v
			if (zc->retrieval->prior_windfield->fit_v[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					p->Klbound[zc->retrieval->prior_windfield->fit_v_Knr[i_w] + ip] =
						zc->retrieval->post_windfield->grid_v[i_w][ip] -
						(factor * zc->retrieval->prior_windfield->grid_v_err[i_w][ip]);
					p->Kubound[zc->retrieval->prior_windfield->fit_v_Knr[i_w] + ip] =
						zc->retrieval->post_windfield->grid_v[i_w][ip] +
						(factor * zc->retrieval->prior_windfield->grid_v_err[i_w][ip]);				
					(*x)[zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip] =
						K2x(p, &(zc->retrieval->post_windfield->grid_v[i_w][ip]), zc->retrieval->prior_windfield->fit_v_Knr[i_w] +ip);
				}
			}

			//grid_w
			if (zc->retrieval->prior_windfield->fit_w[i_w]) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					p->Klbound[zc->retrieval->prior_windfield->fit_w_Knr[i_w] + ip] =
						zc->retrieval->post_windfield->grid_w[i_w][ip] -
						(factor * zc->retrieval->prior_windfield->grid_w_err[i_w][ip]);
					p->Kubound[zc->retrieval->prior_windfield->fit_w_Knr[i_w] + ip] =
						zc->retrieval->post_windfield->grid_w[i_w][ip] +
						(factor * zc->retrieval->prior_windfield->grid_w_err[i_w][ip]);				
					(*x)[zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip] =
						K2x(p, &(zc->retrieval->post_windfield->grid_w[i_w][ip]), zc->retrieval->prior_windfield->fit_w_Knr[i_w] +ip);
				}
			}
			
			if (zc->retrieval->prior_windfield->fit_hspeed_hdir[i_w]) {
				//grid_hspeed
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					p->Klbound[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] + ip] =
						zc->retrieval->post_windfield->grid_hspeed[i_w][ip] -
						(factor * zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip]);
					if (p->Klbound[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] + ip] < 0.) {
						p->Klbound[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] + ip] = 1.e-10;
					}
					p->Kubound[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] + ip] =
						zc->retrieval->post_windfield->grid_hspeed[i_w][ip] +
						(factor * zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip]);				
					(*x)[zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip] =
						K2x(p, &(zc->retrieval->post_windfield->grid_hspeed[i_w][ip]), zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip);
				}
				
				//grid_hdir
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					p->Klbound[zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] + ip] =
						zc->retrieval->post_windfield->grid_hdir[i_w][ip] -
						(factor * zc->retrieval->prior_windfield->grid_hdir_err[i_w][ip]);
					p->Kubound[zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] + ip] =
						zc->retrieval->post_windfield->grid_hdir[i_w][ip] +
						(factor * zc->retrieval->prior_windfield->grid_hdir_err[i_w][ip]);				
					(*x)[zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip] =
						K2x(p, &(zc->retrieval->post_windfield->grid_hdir[i_w][ip]), zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip);
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
						zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[ip] -
						(factor * zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13_err[ip]);	
					if (p->Klbound[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr + ip] <= 0)
						p->Klbound[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr + ip] = 
							1.e-10 * zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13_err[ip];
						
					p->Kubound[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr + ip] =
						zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[ip] +
						(factor * zc->retrieval->prior_windfield->turbulence[i_w]->grid_edr13_err[ip]);				
					(*x)[zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip] =
						K2x(p, &(zc->retrieval->post_windfield->turbulence[i_w]->grid_edr13[ip]), zc->retrieval->prior_windfield->turbulence[i_w]->fit_edr13_Knr +ip);
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
			//dBm
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm) {
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
									x + zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr
									+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									+ip,
									zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr 
									+(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n * i_par)
									+ip))
									/ zc->retrieval->post_scattererfield->psd[i_psd]->grid_rcs_hh[i_par][ip];
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
			
			if (zc->retrieval->prior_windfield->fit_hspeed_hdir[i_w]) {
				//hspeed
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					zc->retrieval->post_windfield->grid_hspeed[i_w][ip] =
						x2K(p, 	x + zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip,
									zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] +ip);
				}

				//hdir
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					zc->retrieval->post_windfield->grid_hdir[i_w][ip] =
						x2K(p, 	x + zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip,
									zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] +ip);
				}

				//update u and v
				util_cast_hspeed_hdir_2_u_v(zc->retrieval->post_windfield, i_w);
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
				griddep = calloc(zc->retrieval->prior_scattererfield->psd[i_psd]->field->n, sizeof(double));

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
								fflush(stdout); exit(1);
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
								fflush(stdout); exit(1);
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
								fflush(stdout); exit(1);
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
			//TBD
			//~ if (
				//~ zc->retrieval->prior_windfield->fit_u[i_w] &
				//~ zc->retrieval->prior_windfield->fit_v[i_w] &
				//~ (zc->retrieval->prior_windfield->use_hspeed_hdir_errors[i_w] == 1)
				//~ ) {
				//~ if (zc->retrieval->post_windfield->grid_hspeed_err[i_w] == NULL)
					//~ zc->retrieval->post_windfield->grid_hspeed_err[i_w] = 
						//~ calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				//~ if (zc->retrieval->post_windfield->grid_hdir_err[i_w] == NULL)
					//~ zc->retrieval->post_windfield->grid_hdir_err[i_w] = 
						//~ calloc(zc->retrieval->prior_windfield->field[i_w]->n, sizeof(double));
				//~ for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {
					//~ tmp_hspeed = sqrt(pow(zc->retrieval->post_windfield->grid_u[i_w][ip],2.) + pow(zc->retrieval->post_windfield->grid_v[i_w][ip],2.));				
					//~ zc->retrieval->post_windfield->grid_hspeed_err[i_w][ip] =
					//~ sqrt(
						//~ pow( (zc->retrieval->post_windfield->grid_u[i_w][ip] / tmp_hspeed) * zc->retrieval->post_windfield->grid_u_err[i_w][ip] ,2.) +
						//~ pow( (zc->retrieval->post_windfield->grid_v[i_w][ip] / tmp_hspeed) * zc->retrieval->post_windfield->grid_v_err[i_w][ip] ,2.)
					//~ );
					//~ zc->retrieval->post_windfield->grid_hdir_err[i_w][ip] =
					//~ sqrt(
						//~ pow( (zc->retrieval->post_windfield->grid_v[i_w][ip] / pow(tmp_hspeed, 2.)) * zc->retrieval->post_windfield->grid_u_err[i_w][ip] ,2.) +
						//~ pow( (zc->retrieval->post_windfield->grid_u[i_w][ip] / pow(tmp_hspeed, 2.)) * zc->retrieval->post_windfield->grid_v_err[i_w][ip] ,2.)
					//~ );
				//~ }
			//~ }
			
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
								fflush(stdout); exit(1);
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
			printf("Error with matrix allocation"); 
			fflush(stdout); exit(1);
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

	double *dummy;
	int i_cast;
	int i_cast_from;
	int i_cast_to;
	int ip;
	double tmpH;
	int i_w;
	int i_par;
	int j;
	double tmp;
	double *xyzt = calloc(4, sizeof(double));
	
	t_zephyros_interpolation_bilint_lut	*tmp_lut_N0 	= NULL;
	t_zephyros_interpolation_bilint_lut	*tmp_lut_D0_mm	= NULL;
	t_zephyros_interpolation_bilint_lut	*tmp_lut_mu 	= NULL;

	t_zephyros_psd *psd_src 		= NULL;
	t_zephyros_psd *psd_dst 		= NULL;
	
	
	
	//wind fields
	for (i_cast = 1; i_cast < c->n_cast_windfield_grid_nrs; i_cast = i_cast + 2) {
		i_cast_from 	= c->cast_windfield_grid_nrs[i_cast - 1];
		i_cast_to 		= c->cast_windfield_grid_nrs[i_cast];
	
		for (ip=0; ip < zc->retrieval->prior_windfield->field[i_cast_to]->n; ip++) {
			xyzt[0] = zc->retrieval->prior_windfield->field[i_cast_to]->x((void*)zc->retrieval->prior_windfield->field[i_cast_to], ip);
			xyzt[1] = zc->retrieval->prior_windfield->field[i_cast_to]->y((void*)zc->retrieval->prior_windfield->field[i_cast_to], ip);
			xyzt[2] = zc->retrieval->prior_windfield->field[i_cast_to]->z((void*)zc->retrieval->prior_windfield->field[i_cast_to], ip);
			xyzt[3] = zc->retrieval->prior_windfield->field[i_cast_to]->t((void*)zc->retrieval->prior_windfield->field[i_cast_to], ip);

			if ((zc->retrieval->post_windfield->lut_u[i_cast_from] != NULL)
				& (zc->retrieval->prior_windfield->grid_u[i_cast_to] != NULL)) {
				interpolation_bilint(
					zc->retrieval->post_windfield->lut_u[i_cast_from],
					xyzt,
					zc->retrieval->post_windfield->grid_u[i_cast_to] + ip,
					0,
					dummy);
				zc->retrieval->prior_windfield->grid_u[i_cast_to][ip] = zc->retrieval->post_windfield->grid_u[i_cast_to][ip];
			}
					
			if ((zc->retrieval->post_windfield->lut_v[i_cast_from] != NULL)
				& (zc->retrieval->prior_windfield->grid_v[i_cast_to] != NULL)) {
				interpolation_bilint(
					zc->retrieval->post_windfield->lut_v[i_cast_from],
					xyzt,
					zc->retrieval->post_windfield->grid_v[i_cast_to] + ip,
					0,
					dummy);
				zc->retrieval->prior_windfield->grid_v[i_cast_to][ip] = zc->retrieval->post_windfield->grid_v[i_cast_to][ip];					
			}

			if ((zc->retrieval->post_windfield->lut_w[i_cast_from] != NULL)
				& (zc->retrieval->prior_windfield->grid_w[i_cast_to] != NULL)) {
				interpolation_bilint(
					zc->retrieval->post_windfield->lut_w[i_cast_from],
					xyzt,
					zc->retrieval->post_windfield->grid_w[i_cast_to] + ip,
					0,
					dummy);
				zc->retrieval->prior_windfield->grid_w[i_cast_to][ip] = zc->retrieval->post_windfield->grid_w[i_cast_to][ip];					
			}
			
			
			if (zc->retrieval->prior_windfield->fit_hspeed_hdir[i_cast_to]) {
				util_cast_u_v_2_hspeed_hdir(zc->retrieval->prior_windfield, i_cast_to);
				util_cast_u_v_2_hspeed_hdir(zc->retrieval->post_windfield, i_cast_to);
			}
			
			/*
			//errors
			if (
				zc->retrieval->post_windfield->fit_u[i_cast_from]
				& zc->retrieval->post_windfield->fit_v[i_cast_from]
				& zc->retrieval->post_windfield->use_hspeed_hdir_errors[i_cast_from]) {

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
		printf("windfield casted from %i to %i\n", i_cast_from, i_cast_to);	fflush(stdout);	
	}
	
	
	//update wind field errors
	for ( i_w = 0; i_w <= 100; i_w++ ) {
		if (zc->retrieval->prior_windfield->field[i_w] != NULL) {		
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
			
			if (	zc->retrieval->prior_windfield->fit_hspeed_hdir[i_w] &
					(c->update_windfield_hspeed_err > 0.)
			) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {				
					zc->retrieval->prior_windfield->grid_hspeed_err[i_w][ip] = c->update_windfield_hspeed_err;
				}
			}		
			
			if (	zc->retrieval->prior_windfield->fit_hspeed_hdir[i_w] &
					(c->update_windfield_hdir_err > 0.)
			) {
				for (ip=0; ip < zc->retrieval->prior_windfield->field[i_w]->n; ip++) {				
					zc->retrieval->prior_windfield->grid_hdir_err[i_w][ip] = c->update_windfield_hdir_err;
				}
			}				
		}
	}

	
	//psd from post to post and prior
	for (i_cast = 1; i_cast < c->n_cast_psd_nrs; i_cast = i_cast + 2) {	
		i_cast_from 	= c->cast_psd_nrs[i_cast - 1];
		i_cast_to 		= c->cast_psd_nrs[i_cast];
		psd_src = zc->retrieval->post_scattererfield->psd[i_cast_from];
		psd_dst = zc->retrieval->post_scattererfield->psd[i_cast_to];
	
		if (
			(psd_src->field->n == psd_dst->field->n) &&
			(psd_src->n_diameters == psd_dst->n_diameters)) {
			//simply copy
			for (ip=0; ip < psd_dst->field->n; ip++) {				
				for ( i_par = 0; i_par < psd_dst->n_diameters; i_par++ ) {
					psd_dst->grid_number_density_m3[i_par][ip] =
						psd_src->grid_number_density_m3[i_par][ip];
				}
			}
		} else {		
			//some fitting
			util_fit_gammadistribution_parameters(psd_dst, psd_src);
								
			//Update dst number densities.
			for (ip=0; ip < psd_dst->field->n; ip++) {				
				for ( i_par = 0; i_par < psd_dst->n_diameters; i_par++ ) {
						psd_dst->grid_number_density_m3[i_par][ip] =
						util_gamma_integral(
						psd_dst->grid_gammadistribution_N0[ip],		//N0
						psd_dst->grid_gammadistribution_mu[ip], 	//mu
						psd_dst->grid_gammadistribution_D0_mm[ip], 	//D0_mm
						psd_dst->discrete_D_equiv_mm_interval_llimt[i_par],
						psd_dst->discrete_D_equiv_mm_interval_ulimt[i_par]);
						
						//also update the prior
						zc->retrieval->prior_scattererfield->psd[i_cast_to]->grid_number_density_m3[i_par][ip] =
							psd_dst->grid_number_density_m3[i_par][ip];
				}
			}
		}
			
		//Update the LUTs
		for ( i_par = 0; i_par < psd_dst->n_diameters; i_par++ ) {
			//prior
			interpolation_free_lut(psd_dst->lut_ln_number_density_m3 + i_par);
			util_field2lut(psd_dst->field, psd_dst->grid_number_density_m3[i_par],
				3, psd_dst->lut_ln_number_density_m3 + i_par);
	
			//post
			interpolation_free_lut(zc->retrieval->post_scattererfield->psd[i_cast_to]->lut_ln_number_density_m3 + i_par);
			util_field2lut(zc->retrieval->post_scattererfield->psd[i_cast_to]->field, psd_dst->grid_number_density_m3[i_par],
				3, zc->retrieval->post_scattererfield->psd[i_cast_to]->lut_ln_number_density_m3 + i_par);			
		}

		printf("PSD casted from %i to %i\n", i_cast_from, i_cast_to);	fflush(stdout);	
	}


	//free
	util_safe_free(&xyzt);
}



void retrieval_fdvar_additional_output(t_fdvar_opc *opc, FILE *fp)
{
	int ip, i_w, io;
    t_zephyros_config 								*zc = opc->zc;
	t_fdvar_o *o = opc->o;
	
	fprintf(fp, "section retrieval\n"); fflush(fp);
	fprintf(fp, "subsection post_scattererfield\n"); fflush(fp);
	zephyros_config_print_scattererfield(zc->retrieval->post_scattererfield, fp); fflush(fp);

	fprintf(fp, "subsection post_windfield\n"); fflush(fp);
	zephyros_config_print_windfield(zc->retrieval->post_windfield, fp); fflush(fp);
	
	radarfilter_write_measurements(zc, 1, o->n, o->model_radarmeasurement, fp);

	//u
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_u", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->model_radarmeasurement[io]->center_air_u);
	fprintf(fp, " \n");
	//v
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_v", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->model_radarmeasurement[io]->center_air_v);
	fprintf(fp, " \n");
	//w
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_w", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->model_radarmeasurement[io]->center_air_w);
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
			if (zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm) {			
				zc->retrieval->prior_scattererfield->psd[i_psd]->fit_dBm_Knr = p->Kn;
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
			
			//grid_hspeed
			//grid_hdir
			if (zc->retrieval->prior_windfield->fit_hspeed_hdir[i_w]) {								
				zc->retrieval->prior_windfield->fit_hspeed_Knr[i_w] = p->Kn;
				p->Kn += zc->retrieval->prior_windfield->field[i_w]->n;

				zc->retrieval->prior_windfield->fit_hdir_Knr[i_w] = p->Kn;
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
	if (c->costfunction_dBZ_hh) 						p->cost_no += o->n;
	if (c->costfunction_dBZ_hv) 						p->cost_no += o->n;
	if (c->costfunction_dBZ_vh) 						p->cost_no += o->n;
	if (c->costfunction_dBZ_vv) 						p->cost_no += o->n;
	if (c->costfunction_dBZdr) 							p->cost_no += o->n;
	if (c->costfunction_dBLdr) 							p->cost_no += o->n;
	if (c->costfunction_Doppler_velocity_hh_ms) 		p->cost_no += o->n;
	if (c->costfunction_Doppler_spectral_width_hh_ms) 	p->cost_no += o->n;

    //spectra
	if (c->costfunction_Doppler_spectrum_dBZ_hh) {
		for (io=0; io < o->n; io++) p->cost_no += o->radarmeasurement[io]->n_spectrum;
	}
	if (c->costfunction_Doppler_spectrum_dBZ_hv) {
		for (io=0; io < o->n; io++) p->cost_no += o->radarmeasurement[io]->n_spectrum;
	}
	if (c->costfunction_Doppler_spectrum_dBZ_vh) {
		for (io=0; io < o->n; io++) p->cost_no += o->radarmeasurement[io]->n_spectrum;
	}
	if (c->costfunction_Doppler_spectrum_dBZ_vv) {
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
					fflush(stdout); exit(1);
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
			fflush(stdout); exit(1);
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

