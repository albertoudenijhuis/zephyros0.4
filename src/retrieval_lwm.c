/*
Description: 
	Linear wind model, retrieval of wind

Revision History:
	2014

Functions:
	Linear wind model, retrieval of wind

Author:
	Albert Oude Nijhuis <albertoudenijhuis@gmail.com>

Institute:
	Delft University of Technology
	
Zephyros version:
	0.3

Project:
	EU FP7 program, the UFO project

Dissemination:
	Confidential, only for members of the UFO project. Potentially public in the future.

Acknowledgement and citation:
	Whenever this code used for publication of scientific results,
	the code writer should be informed, acknowledged and referenced.

Note:
	If you have any suggestions for improvements or amendments, please inform the author of this code.

//to do:
//option to select fitting algorithm
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>


#include "retrieval_lwm.h"

#include "radarfilter.h"
#include "util.h"
#include "ltqnorm.h"


#include <stdlib.h>

#include <nlopt.h>

//uncomment next statement for debug mode
#define _ZEPHYROS_LWM_DEBUG

void retrieval_lwm_apply(t_lwm_opc *opc)
{
	t_lwm_p *p;
	t_lwm_o *o = opc->o;

	#ifdef _ZEPHYROS_LWM_DEBUG
		printf("retrieval_lwm_initialize_p\n"); fflush(stdout);
	#endif 
	retrieval_lwm_initialize_p(opc); p= opc->p;
	#ifdef _ZEPHYROS_LWM_DEBUG
		printf("retrieval_lwm_calculate_xyz\n"); fflush(stdout);
	#endif 
	retrieval_lwm_calculate_xyz(opc);
	
	//walk through grid
	for (p->i=0; p->i < p->field->n; p->i++) {
		//calculate coefficients
		retrieval_lwm_calculate_coefficients(opc);

		//fit
		if ((o->in_analysis_volume_n > 0) & (o->used_for_analysis_n > p->Kn)) {					
			printf("retrieval_lwm_apply, step %i/%i\n", p->i, p->field->n);
			retrieval_lwm_minimize_chisquared(opc);

			//store results
			retrieval_lwm_store_results(opc);
		}
	}
	#ifdef _ZEPHYROS_LWM_DEBUG
		printf("Clean up memory.\n"); fflush(stdout);
	#endif 	
	retrieval_lwm_free_p(opc);

	#ifdef _ZEPHYROS_LWM_DEBUG
		printf("Linear wind model retrieval done.\n"); fflush(stdout);
	#endif 
}

int retrieval_lwm_minimize_chisquared(t_lwm_opc *opc) {
	clock_t start, end;
	double cpu_time_used;
	t_lwm_o *o = opc->o;
	t_lwm_p *p = opc->p;
	t_zephyros_config_retrieval_lwm_cfg *c = opc->c;
	
	double *K;		
	double *Klbound;		
	double *Kubound;		
	int inf;
	double minf;
	void *vd_opc = (void*) opc;
	nlopt_opt opt;
	int i;
	
	start = clock();

	// Starting point
	func_dbl_arr_calloc(p->Kn, &K);
	func_dbl_arr_calloc(p->Kn, &Klbound);
	func_dbl_arr_calloc(p->Kn, &Kubound);
	
	for (i=0; i < p->Kn; i++) 
	{
		K[i] =1.e-10;
		Klbound[i] 	= -1.e5;	
		Kubound[i] 	= 1.e5;	
	}
	
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
		
	nlopt_set_lower_bounds(opt, Klbound);
	nlopt_set_upper_bounds(opt, Kubound);
	nlopt_set_min_objective(opt, retrieval_lwm_chisquared_nlopt, vd_opc);


	//according to manual: nlopt_set_xtol_rel(opt, 1e-4);
	nlopt_set_xtol_rel(opt, 1.e-4);
	nlopt_set_maxtime(opt, c->maximum_time_s);	

	if ((inf = nlopt_optimize(opt, K, &minf)) < 0) {
		printf("nlopt failed!, i = %i\n", inf);
	}
	else {
		printf("Minimum chisquared function found at %.5e\n", minf);
	}

	//for (i=0; i < p->Kn; i++) 
	//{
	//	printf("K[%i] = %.2e\n", i, K[i]);
	//}

	retrieval_lwm_unpack_K(opc,K);
	nlopt_destroy(opt);
	
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Run time: %f s \n", cpu_time_used);
	
	free(K);
	free(Klbound);
	free(Kubound);
}


double retrieval_lwm_chisquared_nlopt(
	unsigned n,
	const double *x,
	double *grad,
	void *params)
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
	retrieval_lwm_chisquared(params, &n2, (double*) x, &result, grad, &iflag);
	return result;
}


void retrieval_lwm_chisquared(
	void 	*vd_opc,		//optional arguments
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

	t_lwm_o *o = ((t_lwm_opc*)vd_opc)->o;
	t_lwm_p *p = ((t_lwm_opc*)vd_opc)->p;
	t_zephyros_config_retrieval_lwm_cfg *c = ((t_lwm_opc*)vd_opc)->c;
    
    double weight, v_r;
    double myval, myval2;
    int io, iK, ir;
        
	//set values to zero
	if (calc_costfunction) {
		*f	= 0.;	
	}
	if (calc_costfunction_derivatives) {
		for (iK=0; iK < p->Kn; iK++) 
		{
			fd[iK] =0.;
		}
	}
    
	//loop through observations.
	for (io=0; io < o->n; io++) {
		if (o->used_for_analysis[io]) {
			if (c->apply_weights) {
				weight = pow(o->radarmeasurement[io]->Doppler_velocity_hh_ms_err,-2.);
			} else {
				weight = 1.;
			}

			//Calculate project Doppler velocity
			o->lwm_vr[io] = 0.;
			for (ir=0; ir < p->n_sv; ir++) {
				if (c->fit_u0)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_u0[io][ir] * 		x[p->fit_Knr_u0];
				if (c->fit_u_x)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_u_x[io][ir] * 	x[p->fit_Knr_u_x];
				if (c->fit_u_z)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_u_z[io][ir] * 	x[p->fit_Knr_u_z];
				if (c->fit_v0)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_v0[io][ir] * 		x[p->fit_Knr_v0];
				if (c->fit_v_y)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_v_y[io][ir] * 	x[p->fit_Knr_v_y];
				if (c->fit_v_z)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_v_z[io][ir] * 	x[p->fit_Knr_v_z];
				if (c->fit_u_y_plus_v_x)			o->lwm_vr[io] += (1. / p->n_sv) * p->coef_u_y_plus_v_x[io][ir] * x[p->fit_Knr_u_y_plus_v_x];
				if (c->fit_w0)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_w0[io][ir] * 		x[p->fit_Knr_w0];
				if (c->fit_w_x)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_w_x[io][ir] * 	x[p->fit_Knr_w_x];
				if (c->fit_w_y)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_w_y[io][ir] * 	x[p->fit_Knr_w_y];
				if (c->fit_w_z)						o->lwm_vr[io] += (1. / p->n_sv) * p->coef_w_z[io][ir] * 	x[p->fit_Knr_w_z];
				if (c->fit_u_t_plus_v_t_plus_w_t)	o->lwm_vr[io] += (1. / p->n_sv) * p->coef_u_t_plus_v_t_plus_w_t[io][ir] * 	x[p->fit_Knr_u_t_plus_v_t_plus_w_t];
			}
			
			myval = weight * pow(o->lwm_vr[io] - o->radarmeasurement[io]->Doppler_velocity_hh_ms, 2.);
			myval2 = 2. * weight * (o->lwm_vr[io] - o->radarmeasurement[io]->Doppler_velocity_hh_ms);
			
			//debug
			//~ printf("o->vr[%i] = %.2e\n", io, o->vr[io]);
			//~ printf("o->lwm_vr[%i] = %.2e\n", io, o->lwm_vr[io]);

			if (calc_costfunction) {
				*f	+= myval;	
			}
			
			if (calc_costfunction_derivatives) {
				//loop over parameters
				for (iK = 0; iK < p->Kn; iK++) {
					for (ir=0; ir < p->n_sv; ir++) {
						if ((c->fit_u0) & (iK == p->fit_Knr_u0)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_u0[io][ir];}
						if ((c->fit_u_x) & (iK == p->fit_Knr_u_x)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_u_x[io][ir];}
						if ((c->fit_u_z) & (iK == p->fit_Knr_u_z)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_u_z[io][ir];}
						if ((c->fit_v0) & (iK == p->fit_Knr_v0)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_v0[io][ir];}
						if ((c->fit_v_y) & (iK == p->fit_Knr_v_y)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_v_y[io][ir];}
						if ((c->fit_v_z) & (iK == p->fit_Knr_v_z)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_v_z[io][ir];}
						if ((c->fit_u_y_plus_v_x) & (iK == p->fit_Knr_u_y_plus_v_x)) 	{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_u_y_plus_v_x[io][ir];}
						if ((c->fit_w0) & (iK == p->fit_Knr_w0)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_w0[io][ir];}
						if ((c->fit_w_x) & (iK == p->fit_Knr_w_x)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_w_x[io][ir];}
						if ((c->fit_w_y) & (iK == p->fit_Knr_w_y)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_w_y[io][ir];}
						if ((c->fit_w_z) & (iK == p->fit_Knr_w_z)) 						{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_w_z[io][ir];}
						if ((c->fit_u_t_plus_v_t_plus_w_t) & (iK == p->fit_Knr_u_t_plus_v_t_plus_w_t)) 	{fd[iK] += (1. / p->n_sv) * myval2 * p->coef_u_t_plus_v_t_plus_w_t[io][ir];}
					}
				}
			}
		}
	}
	
	//debug
	//~ if (calc_costfunction) {
		//~ printf("*f = %.2e\n", *f);
//~ 
		//~ if (c->fit_u0)						printf("u0 = %.2e\n", x[p->fit_Knr_u0]);
		//~ if (c->fit_u_x)						printf("u_x = %.2e\n", x[p->fit_Knr_u_x]);
		//~ if (c->fit_u_z)						printf("u_z = %.2e\n", x[p->fit_Knr_u_z]);
		//~ if (c->fit_v0)						printf("v0 = %.2e\n", x[p->fit_Knr_v0]);
		//~ if (c->fit_v_y)						printf("v_y = %.2e\n", x[p->fit_Knr_v_y]);
		//~ if (c->fit_v_z)						printf("v_z = %.2e\n", x[p->fit_Knr_v_z]);
		//~ if (c->fit_u_y_plus_v_x)			printf("u_y_plus_v_x = %.2e\n", x[p->fit_Knr_u_y_plus_v_x]);
		//~ if (c->fit_w0)						printf("w0 = %.2e\n", x[p->fit_Knr_w0]);
		//~ if (c->fit_w_x)						printf("w_x = %.2e\n", x[p->fit_Knr_w_x]);
		//~ if (c->fit_w_y)						printf("w_y = %.2e\n", x[p->fit_Knr_w_y]);
		//~ if (c->fit_w_z)						printf("w_z = %.2e\n", x[p->fit_Knr_w_z]);
		//~ if (c->fit_u_t_plus_v_t_plus_w_t)	printf("u_t_plus_v_t_plus_w_t = %.2e\n", x[p->fit_Knr_u_t_plus_v_t_plus_w_t]);
		//~ printf("\n");
	//~ }
	//~ 
	//~ if (calc_costfunction_derivatives) {
		//~ for (iK=0; iK < p->Kn; iK++) 
		//~ {
			//~ printf("fd[%i] = %.2e\n", iK, fd[iK]);
		//~ }
	//~ }

	//~ if (calc_costfunction) {printf("current chisquared value: %.5e \n"   , *f);}
}

void retrieval_lwm_calculate_xyz(t_lwm_opc *opc)
{
	t_lwm_o *o = opc->o;
	t_lwm_p *p = opc->p;
	t_zephyros_config_retrieval_lwm_cfg *c = opc->c;
    t_zephyros_config *zcfg = (t_zephyros_config*)c->vd_zcfg;

	//calculate center coordinates
	int io, ir;
	double myr1, myr2;
	
	int i_beam_range, i_beam_theta, i_beam_phi, i_t, remainder;
	double cdfP, FWHM_rad, mincdfP, maxcdfP, difcdfP, sigma_gauss;
	
	for (io=0; io < o->n; io++) {
		myr1	= o->radarmeasurement[io]->azel_r1_m;
		myr2	= o->radarmeasurement[io]->azel_r2_m;
		if (myr1 == 0.) {
			myr1 = myr2 / 10.;
		}
		o->radarmeasurement[io]->center_coor->radar_range 		= pow(pow(myr1 , -1.) - (0.5 * (pow(myr1 , -1.) - pow(myr2, -1.))), -1.);		
		o->radarmeasurement[io]->center_coor->radar_azel_alpha 	= o->radarmeasurement[io]->azel_alpha_rad;
		o->radarmeasurement[io]->center_coor->radar_azel_gamma 	= o->radarmeasurement[io]->azel_gamma_rad;

		coordinates_radar_azel2enu(o->radarmeasurement[io]->center_coor);
		coordinates_radar_azelrangedir2enu(o->radarmeasurement[io]->center_coor);

		//*****
		//calculate subvolume coordinates
		for ( ir = 0; ir < p->n_sv; ir++ ) {
			//take over radar location
			p->subvolume_coor[io][ir]->enu_radar_location_xyzt[0] = o->radarmeasurement[io]->center_coor->enu_radar_location_xyzt[0];
			p->subvolume_coor[io][ir]->enu_radar_location_xyzt[1] = o->radarmeasurement[io]->center_coor->enu_radar_location_xyzt[1];
			p->subvolume_coor[io][ir]->enu_radar_location_xyzt[2] = o->radarmeasurement[io]->center_coor->enu_radar_location_xyzt[2];
			p->subvolume_coor[io][ir]->enu_radar_location_xyzt[3] = o->radarmeasurement[io]->center_coor->enu_radar_location_xyzt[3];

			//get indices for range, phi, theta and time
			//i_beam_range 	0 ... n_beam_range -1
			//i_beam_theta 	0 ... n_beam_theta -1
			//i_beam_phi 	0 ... n_beam_phi -1 
			//i_t		 	0 ... n_t -1

			i_beam_range	= ir % zcfg->retrieval->radarfilter->n_beam_range;
			remainder 		= (ir - i_beam_range) / zcfg->retrieval->radarfilter->n_beam_range;
			i_beam_theta	= remainder % zcfg->retrieval->radarfilter->n_beam_theta;
			remainder 		= (remainder - i_beam_theta) / zcfg->retrieval->radarfilter->n_beam_theta;
			i_beam_phi		= remainder % zcfg->retrieval->radarfilter->n_beam_phi;
			remainder 		= (remainder - i_beam_phi) / zcfg->retrieval->radarfilter->n_beam_phi;
			i_t				= remainder % zcfg->retrieval->radarfilter->n_t;
			
			//beam range is chosen such that equal power is there.
			//for volumetric distributed scatterers the weighting factor is (r^-2)
			cdfP = (i_beam_range + 0.5) / zcfg->retrieval->radarfilter->n_beam_range; 
			myr1	= o->radarmeasurement[io]->azel_r1_m;
			myr2	= o->radarmeasurement[io]->azel_r2_m;
			if (myr1 == 0.) {
				myr1 = myr2 / 10.;
			}
			p->subvolume_coor[io][ir]->radar_range = 
				pow(  pow(myr1 , -1.) - (cdfP * (pow(myr1 , -1.) - pow(myr2, -1.))), -1.);		
		
			//calculate point_beam_phi_rad
			p->subvolume_coor[io][ir]->radar_beam_phi =
					(2. * M_PI * i_beam_phi ) / zcfg->retrieval->radarfilter->n_beam_phi;

			//assuming elips form for the beam:
			//given phi_rad, calculate FWHM
			//note:
			//phi is defined s.t.:
			//phi = 0, pi			: beam distortion in azimutal direction
			//phi = pi/2, -pi/2		: beam distortion in elevation direction
			//hence that: 
			//FWHM0: is FWHM in azimutal direction
			//FWHM1: is FWHM in elevation direction
			FWHM_rad = 	pow(
						pow(cos(p->subvolume_coor[io][ir]->radar_beam_phi) * o->radarmeasurement[io]->beam_FWHM0_rad,2) +
						pow(sin(p->subvolume_coor[io][ir]->radar_beam_phi) * o->radarmeasurement[io]->beam_FWHM1_rad,2)
						, 1./2.);

			//beam theta is chosen such that equal power is there.
			if ((zcfg->retrieval->radarfilter->n_beam_theta == 1) & (zcfg->retrieval->radarfilter->n_beam_phi == 1))
			{			
				//special case, n_beam_theta = 1, n_beam_phi = 1
				p->subvolume_coor[io][ir]->radar_beam_theta = 0.;
			} else {					
				//examples:
				//n_beam_theta = 1, cdfP = [0.75]
				//n_beam_theta = 2, cdfP = [0.625, 875]
				//n_beam_theta = 4, cdfP = [0.5625, 0.6875, 0.8125, 0.9375]

				mincdfP = 0.5	;
				maxcdfP = 1.	;
				difcdfP = (maxcdfP - mincdfP) / zcfg->retrieval->radarfilter->n_beam_theta;
				cdfP = mincdfP + (i_beam_theta + 0.5) * difcdfP; 

				sigma_gauss = FWHM_rad / sqrt(64. * log(2.));
				p->subvolume_coor[io][ir]->radar_beam_theta = sigma_gauss * ltqnorm(cdfP);
			}

			//beam to azel
			coordinates_radar_beam2azel(p->subvolume_coor[io][ir]);
			
			//add central values of the azel sytem
			p->subvolume_coor[io][ir]->radar_azel_alpha 		+= o->radarmeasurement[io]->azel_alpha_rad;
			p->subvolume_coor[io][ir]->radar_azel_gamma 		+= o->radarmeasurement[io]->center_coor->radar_azel_gamma;
			p->subvolume_coor[io][ir]->radar_azel_gamma_cor 	+= o->radarmeasurement[io]->center_coor->radar_azel_gamma_cor;
			
			//azu -> enu
			coordinates_radar_azel2enu(p->subvolume_coor[io][ir]);
			coordinates_radar_azelrangedir2enu(p->subvolume_coor[io][ir]);
		}
	}
}	

void retrieval_lwm_calculate_coefficients(t_lwm_opc *opc)
{
	t_lwm_o *o = opc->o;
	t_lwm_p *p = opc->p;
	t_zephyros_config_retrieval_lwm_cfg *c = opc->c;
	
	int io, ir;
	int i;
	
	int *indices_extra_points = malloc(o->n * sizeof(int));
	int n_ep;
	int n_ep2;
	int *my_permutation = NULL;

	o->in_analysis_volume_n = 0;
	//loop through observations.
	for (io=0; io < o->n; io++) {
		//check whether the observation is in the analysis volume
		o->in_analysis_volume[io] = 
			((p->field->x(p->field, p->i) - (c->xvec_m[2]/2.)) <= o->radarmeasurement[io]->center_coor->enu_xyzt[0]) &
				(o->radarmeasurement[io]->center_coor->enu_xyzt[0] < (p->field->x(p->field, p->i) + (c->xvec_m[2]/2.))) &
			((p->field->y(p->field, p->i) - (c->yvec_m[2]/2.)) <= o->radarmeasurement[io]->center_coor->enu_xyzt[1]) &
				(o->radarmeasurement[io]->center_coor->enu_xyzt[1] < (p->field->y(p->field, p->i) + (c->yvec_m[2]/2.))) &
			((p->field->z(p->field, p->i) - (c->zvec_m[2]/2.)) <= o->radarmeasurement[io]->center_coor->enu_xyzt[2]) &
				(o->radarmeasurement[io]->center_coor->enu_xyzt[2] < (p->field->z(p->field, p->i) + (c->zvec_m[2]/2.))) &
			((p->field->t(p->field, p->i) - (c->tvec_s[2]/2.)) <= o->radarmeasurement[io]->center_coor->enu_xyzt[3]) &
				(o->radarmeasurement[io]->center_coor->enu_xyzt[3] < (p->field->t(p->field, p->i) + (c->tvec_s[2]/2.)));
		o->used_for_analysis[io] = o->in_analysis_volume[io];
		if (o->in_analysis_volume[io]) {o->in_analysis_volume_n++;}
	}
	o->used_for_analysis_n = o->in_analysis_volume_n;
			
	//add random extra points
	n_ep = 0;
	for (io=0; io < o->n; io++) {
		//check whether the observation is in the extra analysis volume
		if 	(
			((p->field->x(p->field, p->i) - (c->extra_points_dx/2.)) <= o->radarmeasurement[io]->center_coor->enu_xyzt[0]) & (o->radarmeasurement[io]->center_coor->enu_xyzt[0] < (p->field->x(p->field, p->i) + (c->extra_points_dx/2.))) &
			((p->field->y(p->field, p->i) - (c->extra_points_dy/2.)) <= o->radarmeasurement[io]->center_coor->enu_xyzt[1]) & (o->radarmeasurement[io]->center_coor->enu_xyzt[1] < (p->field->y(p->field, p->i) + (c->extra_points_dy/2.))) &
			((p->field->z(p->field, p->i) - (c->extra_points_dz/2.)) <= o->radarmeasurement[io]->center_coor->enu_xyzt[2]) & (o->radarmeasurement[io]->center_coor->enu_xyzt[2] < (p->field->z(p->field, p->i) + (c->extra_points_dz/2.))) &
			((p->field->t(p->field, p->i) - (c->extra_points_dt/2.)) <= o->radarmeasurement[io]->center_coor->enu_xyzt[3]) & (o->radarmeasurement[io]->center_coor->enu_xyzt[3] < (p->field->t(p->field, p->i) + (c->extra_points_dt/2.))) &
			(o->in_analysis_volume[io] == 0)
			)
			{
				indices_extra_points[n_ep] = io;
				n_ep++;
			}
	}

	if (n_ep > 0) {
		//make random perturbation
		randompermutation(n_ep, &my_permutation);

		if (c->extra_points_n < n_ep) {
			n_ep2 = c->extra_points_n;
		} else {
			n_ep2 = n_ep;
		}
				
		//add points
		for (i=0; i < n_ep2; i++) {
			io = indices_extra_points[my_permutation[i]];
			o->used_for_analysis[io] = 1;
			o->used_for_analysis_n++;
		}
	}


	//loop through observations.
	for (io=0; io < o->n; io++) {
		if (o->used_for_analysis[io] == 1) {
			
			//debug
			//~ printf("o->used_for_analysis[%i] =%i\n", io, o->used_for_analysis[io]);
			//~ printf("o->in_analysis_volume[%i] = %i\n", io, o->in_analysis_volume[io]);
			//~ printf("p->o_x[%i] = %.2e\n", io, p->o_x[io]);
			//~ printf("p->o_y[%i] = %.2e\n", io, p->o_y[io]);
			//~ printf("p->o_z[%i] = %.2e\n", io, p->o_z[io]);
			//~ printf("\n");

			//coefficients
			if (c->fit_u0) 						{p->center_coef_u0[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[0];}
			if (c->fit_u_x) 					{p->center_coef_u_x[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[0] * (o->radarmeasurement[io]->center_coor->enu_xyzt[0] - p->field->x(p->field, p->i));}
			if (c->fit_u_z) 					{p->center_coef_u_z[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[0] * (o->radarmeasurement[io]->center_coor->enu_xyzt[2] - p->field->z(p->field, p->i));}
			if (c->fit_v0) 						{p->center_coef_v0[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[1];}
			if (c->fit_v_y) 					{p->center_coef_v_y[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[1] * (o->radarmeasurement[io]->center_coor->enu_xyzt[1] - p->field->y(p->field, p->i));}
			if (c->fit_v_z) 					{p->center_coef_v_z[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[1] * (o->radarmeasurement[io]->center_coor->enu_xyzt[2] - p->field->z(p->field, p->i));}
			if (c->fit_u_y_plus_v_x) 			{p->center_coef_u_y_plus_v_x[io] = 
													o->radarmeasurement[io]->center_coor->radar_enu_dir[0] * (o->radarmeasurement[io]->center_coor->enu_xyzt[1] - p->field->y(p->field, p->i)) +
													o->radarmeasurement[io]->center_coor->radar_enu_dir[1] * (o->radarmeasurement[io]->center_coor->enu_xyzt[0] - p->field->x(p->field, p->i));}
			if (c->fit_w0) 						{p->center_coef_w0[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[2];}
			if (c->fit_w_x) 					{p->center_coef_w_x[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[2] * (o->radarmeasurement[io]->center_coor->enu_xyzt[0] - p->field->x(p->field, p->i));}
			if (c->fit_w_y) 					{p->center_coef_w_y[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[2] * (o->radarmeasurement[io]->center_coor->enu_xyzt[1] - p->field->y(p->field, p->i));}
			if (c->fit_w_z) 					{p->center_coef_w_z[io]		= o->radarmeasurement[io]->center_coor->radar_enu_dir[2] * (o->radarmeasurement[io]->center_coor->enu_xyzt[2] - p->field->z(p->field, p->i));}
			if (c->fit_u_t_plus_v_t_plus_w_t) 	{p->center_coef_u_t_plus_v_t_plus_w_t[io]	= o->radarmeasurement[io]->center_coor->enu_xyzt[3] - p->field->t(p->field, p->i);}

			for (ir=0; ir < p->n_sv; ir++) {
				if (c->fit_u0) 						{p->coef_u0[io][ir]		= p->subvolume_coor[io][ir]->radar_enu_dir[0];}
				if (c->fit_u_x) 					{p->coef_u_x[io][ir]	= p->subvolume_coor[io][ir]->radar_enu_dir[0] * (p->subvolume_coor[io][ir]->enu_xyzt[0] - p->field->x(p->field, p->i));}
				if (c->fit_u_z) 					{p->coef_u_z[io][ir]	= p->subvolume_coor[io][ir]->radar_enu_dir[0] * (p->subvolume_coor[io][ir]->enu_xyzt[2] - p->field->z(p->field, p->i));}
				if (c->fit_v0) 						{p->coef_v0[io][ir]		= p->subvolume_coor[io][ir]->radar_enu_dir[1];}
				if (c->fit_v_y) 					{p->coef_v_y[io][ir]	= p->subvolume_coor[io][ir]->radar_enu_dir[1] * (p->subvolume_coor[io][ir]->enu_xyzt[1] - p->field->y(p->field, p->i));}
				if (c->fit_v_z) 					{p->coef_v_z[io][ir]	= p->subvolume_coor[io][ir]->radar_enu_dir[1] * (p->subvolume_coor[io][ir]->enu_xyzt[2] - p->field->z(p->field, p->i));}
				if (c->fit_u_y_plus_v_x) 			{p->coef_u_y_plus_v_x[io][ir] = 
														p->subvolume_coor[io][ir]->radar_enu_dir[0] * (p->subvolume_coor[io][ir]->enu_xyzt[1] - p->field->y(p->field, p->i)) +
														p->subvolume_coor[io][ir]->radar_enu_dir[1] * (p->subvolume_coor[io][ir]->enu_xyzt[0] - p->field->x(p->field, p->i));}
				if (c->fit_w0) 						{p->coef_w0[io][ir]		= p->subvolume_coor[io][ir]->radar_enu_dir[2];}
				if (c->fit_w_x) 					{p->coef_w_x[io][ir]	= p->subvolume_coor[io][ir]->radar_enu_dir[2] * (p->subvolume_coor[io][ir]->enu_xyzt[0] - p->field->x(p->field, p->i));}
				if (c->fit_w_y) 					{p->coef_w_y[io][ir]	= p->subvolume_coor[io][ir]->radar_enu_dir[2] * (p->subvolume_coor[io][ir]->enu_xyzt[1] - p->field->y(p->field, p->i));}
				if (c->fit_w_z) 					{p->coef_w_z[io][ir]	= p->subvolume_coor[io][ir]->radar_enu_dir[2] * (p->subvolume_coor[io][ir]->enu_xyzt[2] - p->field->z(p->field, p->i));}
				if (c->fit_u_t_plus_v_t_plus_w_t) 	{p->coef_u_t_plus_v_t_plus_w_t[io][ir]	= p->subvolume_coor[io][ir]->enu_xyzt[3] - p->field->t(p->field, p->i);}
			}
		}
	}

	free(indices_extra_points);
	if (my_permutation != NULL) {
		free(my_permutation);
	}
}

void retrieval_lwm_store_results(t_lwm_opc *opc)
{
	t_lwm_o *o = opc->o;
	t_lwm_p *p = opc->p;
	t_zephyros_config_retrieval_lwm_cfg *c = opc->c;
	
	int ip;
	int io;

	double sigma_vr;
	double sigma_u0, sigma_u_x, sigma_u_z, sigma_v0, sigma_v_y, sigma_v_z, sigma_u_y_plus_v_x, sigma_w0, sigma_w_x, sigma_w_y, sigma_w_z, sigma_u_t_plus_v_t_plus_w_t;
	
	//loop through observations.
	for (io=0; io < o->n; io++) {
		if (o->in_analysis_volume[io]) {
												o->windvector_u[io] 	= 	0.;
			if (c->fit_u0) 						o->windvector_u[io] 	+= 	p->u0[p->i];
			if (c->fit_u_x)						o->windvector_u[io] 	+= 	p->u_x[p->i] * (o->radarmeasurement[io]->center_coor->enu_xyzt[0] - p->field->x(p->field, p->i));
			if (c->fit_u_y_plus_v_x) 			o->windvector_u[io] 	+= 	p->u_y_plus_v_x[p->i] * .5 * (o->radarmeasurement[io]->center_coor->enu_xyzt[1] - p->field->y(p->field, p->i));
			if (c->fit_u_z) 					o->windvector_u[io] 	+= 	p->u_z[p->i] * (o->radarmeasurement[io]->center_coor->enu_xyzt[2] - p->field->z(p->field, p->i));
			if (c->fit_u_t_plus_v_t_plus_w_t)	o->windvector_u[io] 	+= 	p->u_t_plus_v_t_plus_w_t[p->i] * (1./3.) * (o->radarmeasurement[io]->center_coor->enu_xyzt[3] - p->field->t(p->field, p->i));
			
												o->windvector_v[io] 	= 	0.;
			if (c->fit_v0)						o->windvector_v[io] 	+=	p->v0[p->i];
			if (c->fit_u_y_plus_v_x)			o->windvector_v[io] 	+=	p->u_y_plus_v_x[p->i] * .5 * (o->radarmeasurement[io]->center_coor->enu_xyzt[0] - p->field->x(p->field, p->i));
			if (c->fit_v_y)						o->windvector_v[io] 	+=	p->v_y[p->i] * (o->radarmeasurement[io]->center_coor->enu_xyzt[1] - p->field->y(p->field, p->i));
			if (c->fit_v_z) 					o->windvector_v[io] 	+=	p->v_z[p->i] * (o->radarmeasurement[io]->center_coor->enu_xyzt[2] - p->field->z(p->field, p->i));
			if (c->fit_u_t_plus_v_t_plus_w_t)	o->windvector_v[io] 	+= 	p->u_t_plus_v_t_plus_w_t[p->i] * (1./3.) * (o->radarmeasurement[io]->center_coor->enu_xyzt[3] - p->field->t(p->field, p->i));
							
												o->windvector_w[io] 	= 	0.;
			if (c->fit_w0) 						o->windvector_w[io] 	+=	p->w0[p->i];
			if (c->fit_w_x)						o->windvector_w[io] 	+=	p->w_x[p->i] * (o->radarmeasurement[io]->center_coor->enu_xyzt[0] - p->field->x(p->field, p->i));
			if (c->fit_w_y)						o->windvector_w[io] 	+=	p->w_y[p->i] * (o->radarmeasurement[io]->center_coor->enu_xyzt[1] - p->field->y(p->field, p->i));
			if (c->fit_w_z) 					o->windvector_w[io] 	+=	p->w_z[p->i] * (o->radarmeasurement[io]->center_coor->enu_xyzt[2] - p->field->z(p->field, p->i));
			if (c->fit_u_t_plus_v_t_plus_w_t)	o->windvector_w[io] 	+= 	p->u_t_plus_v_t_plus_w_t[p->i] * (1./3.) * (o->radarmeasurement[io]->center_coor->enu_xyzt[3] - p->field->t(p->field, p->i));
		}          
	}
	
	//error calculations
	//obtain best estimate for the uncertainty in the measurements
	sigma_vr = 0;
	for (io=0; io < o->n; io++) {
		if (o->in_analysis_volume[io]) {	
			sigma_vr += pow(o->lwm_vr[io] - o->radarmeasurement[io]->Doppler_velocity_hh_ms, 2.);
		}
	}
	sigma_vr 	= sqrt(   (1. / (o->used_for_analysis_n - p->Kn)) * sigma_vr);
	
	for (io=0; io < o->n; io++) {
		if (o->in_analysis_volume[io]) {
			if (c->fit_u0)						sigma_u0	= ((p->center_coef_u0[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_u0[io])) * sigma_vr));
			if (c->fit_u_x)						sigma_u_x	= ((p->center_coef_u_x[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_u_x[io])) * sigma_vr));
			if (c->fit_u_z)						sigma_u_z	= ((p->center_coef_u_z[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_u_z[io])) * sigma_vr));
			if (c->fit_v0)						sigma_v0	= ((p->center_coef_v0[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_v0[io])) * sigma_vr));
			if (c->fit_v_y)						sigma_v_y	= ((p->center_coef_v_y[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_v_y[io])) * sigma_vr));
			if (c->fit_v_z)						sigma_v_z	= ((p->center_coef_v_z[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_v_z[io])) * sigma_vr));
			if (c->fit_u_y_plus_v_x)			sigma_u_y_plus_v_x	= ((p->center_coef_u_y_plus_v_x[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_u_y_plus_v_x[io])) * sigma_vr));
			if (c->fit_w0)						sigma_w0	= ((p->center_coef_w0[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_w0[io])) * sigma_vr));
			if (c->fit_w_x)						sigma_w_x	= ((p->center_coef_w_x[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_w_x[io])) * sigma_vr));
			if (c->fit_w_y)						sigma_w_y	= ((p->center_coef_w_y[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_w_y[io])) * sigma_vr));
			if (c->fit_w_z)						sigma_w_z	= ((p->center_coef_w_z[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_w_z[io])) * sigma_vr));
			if (c->fit_u_t_plus_v_t_plus_w_t)	sigma_u_t_plus_v_t_plus_w_t	= ((p->center_coef_u_t_plus_v_t_plus_w_t[io] < 1.e-10)?0.:(fabs( 1. / (p->center_coef_u_t_plus_v_t_plus_w_t[io])) * sigma_vr));
			
												o->windvector_u_err[io] 	= 	0.;
			if (c->fit_u0) 						o->windvector_u_err[io] 	+=  pow(sigma_u0,2.);
			if (c->fit_u_x)						o->windvector_u_err[io] 	+= 	pow(sigma_u_x * (o->radarmeasurement[io]->center_coor->enu_xyzt[0] - p->field->x(p->field, p->i)),2.);
			if (c->fit_u_y_plus_v_x) 			o->windvector_u_err[io] 	+= 	pow(sigma_u_y_plus_v_x * .5 * (o->radarmeasurement[io]->center_coor->enu_xyzt[1] - p->field->y(p->field, p->i)),2.);
			if (c->fit_u_z) 					o->windvector_u_err[io] 	+= 	pow(sigma_u_z * (o->radarmeasurement[io]->center_coor->enu_xyzt[2] - p->field->z(p->field, p->i)),2.);
			if (c->fit_u_t_plus_v_t_plus_w_t)	o->windvector_u_err[io] 	+= 	pow(sigma_u_t_plus_v_t_plus_w_t * (1./3.) * (o->radarmeasurement[io]->center_coor->enu_xyzt[3] - p->field->t(p->field, p->i)),2.);
												o->windvector_u_err[io] 	= 	sqrt(o->windvector_u_err[io]);

												o->windvector_v_err[io] 	= 	0.;
			if (c->fit_v0)						o->windvector_v_err[io] 	+=	pow(sigma_v0,2.);
			if (c->fit_u_y_plus_v_x)			o->windvector_v_err[io] 	+=	pow(sigma_u_y_plus_v_x * .5 * (o->radarmeasurement[io]->center_coor->enu_xyzt[0] - p->field->x(p->field, p->i)),2.);
			if (c->fit_v_y)						o->windvector_v_err[io] 	+=	pow(sigma_v_y * (o->radarmeasurement[io]->center_coor->enu_xyzt[1] - p->field->y(p->field, p->i)),2.);
			if (c->fit_v_z) 					o->windvector_v_err[io] 	+=	pow(sigma_v_z * (o->radarmeasurement[io]->center_coor->enu_xyzt[2] - p->field->z(p->field, p->i)),2.);
			if (c->fit_u_t_plus_v_t_plus_w_t)	o->windvector_v_err[io] 	+= 	pow(sigma_u_t_plus_v_t_plus_w_t * (1./3.) * (o->radarmeasurement[io]->center_coor->enu_xyzt[3] - p->field->t(p->field, p->i)),2.);
												o->windvector_v_err[io] 	= 	sqrt(o->windvector_v_err[io]);

												o->windvector_w_err[io] 	= 	0.;
			if (c->fit_w0) 						o->windvector_w_err[io] 	+=	pow(sigma_w0,2.);
			if (c->fit_w_x)						o->windvector_w_err[io] 	+=	pow(sigma_w_x * (o->radarmeasurement[io]->center_coor->enu_xyzt[0] - p->field->x(p->field, p->i)),2.);
			if (c->fit_w_y)						o->windvector_w_err[io] 	+=	pow(sigma_w_y * (o->radarmeasurement[io]->center_coor->enu_xyzt[1] - p->field->y(p->field, p->i)),2.);
			if (c->fit_w_z) 					o->windvector_w_err[io] 	+=	pow(sigma_w_z * (o->radarmeasurement[io]->center_coor->enu_xyzt[2] - p->field->z(p->field, p->i)),2.);
			if (c->fit_u_t_plus_v_t_plus_w_t)	o->windvector_w_err[io] 	+= 	pow(sigma_u_t_plus_v_t_plus_w_t * (1./3.) * (o->radarmeasurement[io]->center_coor->enu_xyzt[3] - p->field->t(p->field, p->i)),2.);
												o->windvector_w_err[io] 	= 	sqrt(o->windvector_w_err[io]);
		}
	}											
}

/*
    #shear along the beam, dvr_dr(x,y,z)
    p_vr_sh = {}
    p_vr_sh['u_x']	        = d['cos_theta_e_prime'] * d['sin_phi'] * d['dx_dr']
    p_vr_sh['u_z']	        = d['cos_theta_e_prime'] * d['sin_phi'] * d['dz_dr']
    p_vr_sh['v_y']  	    = d['cos_theta_e_prime'] * d['cos_phi'] * d['dy_dr']
    p_vr_sh['v_z']	        = d['cos_theta_e_prime'] * d['cos_phi'] * d['dz_dr']
    p_vr_sh['u_y_plus_v_x'] = d['cos_theta_e_prime'] * d['cos_theta_e_prime'] * d['sin_phi'] * d['cos_phi']
    p_vr_sh['w_x']	        = d['sin_theta_e_prime'] * d['dx_dr']
    p_vr_sh['w_y']	        = d['sin_theta_e_prime'] * d['dy_dr']
    p_vr_sh['w_z']          = d['sin_theta_e_prime'] * d['dz_dr']

    #horizontal divergence
    p_div = {}
    p_div['u_x'] = 1.
    p_div['v_y'] = 1.

    #horizontal shear
    p_sh = {}
    p_sh['u_y_plus_v_x'] = 1.

    #horizontal stretching
    p_str = {}
    p_str['u_x'] = 1.
    p_str['v_y'] = -1.
*/


void retrieval_lwm_additional_output(t_lwm_opc *opc, FILE *fp)
{
	t_lwm_o *o = opc->o;
	t_lwm_p *p = opc->p;
	t_zephyros_config_retrieval_lwm_cfg *c = opc->c;
    t_zephyros_config *zcfg = (t_zephyros_config*)c->vd_zcfg;

	int io;
	
	//print coordinates only
	radarfilter_write_measurements(zcfg, 2, o->n, o->radarmeasurement, fp);
	
	//output as !! <varname> <ndim> <dim1> <data>
	//o->lwm_vr
	fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_velocity_hh_ms", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->lwm_vr[io]);
	fprintf(fp, " \n");
	
	//residual
	//fprintf(fp, "!! %-30s %-15i %-15i", "residual", 1, o->n);
	//for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->lwm_vr[io] -o->radarmeasurement[io]->Doppler_velocity_hh_ms);
	//fprintf(fp, " \n");
	
	//windvector_u
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_u", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->windvector_u[io]);
	fprintf(fp, " \n");
	//windvector_u_err
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_u_err", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->windvector_u_err[io]);
	fprintf(fp, " \n");
	//windvector_v
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_v", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->windvector_v[io]);
	fprintf(fp, " \n");
	//windvector_v_err
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_v_err", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->windvector_v_err[io]);
	fprintf(fp, " \n");
	//windvector_w
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_w", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->windvector_w[io]);
	fprintf(fp, " \n");
	//windvector_w_err
	fprintf(fp, "!! %-30s %-15i %-15i", "retrieved_w_err", 1, o->n);
	for (io=0; io < o->n; io++) fprintf(fp, " %-15.3e", o->windvector_w_err[io]);
	fprintf(fp, " \n");
}


void retrieval_lwm_unpack_K(
	t_lwm_opc *opc,
	double *K)
{
	t_lwm_p *p = opc->p;
	t_zephyros_config_retrieval_lwm_cfg *c = opc->c;
	
	if (c->fit_u0) 						{p->u0[p->i]					= K[p->fit_Knr_u0];}
	if (c->fit_u_x) 					{p->u_x[p->i] 					= K[p->fit_Knr_u_x];}
	if (c->fit_u_z) 					{p->u_z[p->i]					= K[p->fit_Knr_u_z];}
	if (c->fit_v0) 						{p->v0[p->i] 					= K[p->fit_Knr_v0];}
	if (c->fit_v_y) 					{p->v_y[p->i] 					= K[p->fit_Knr_v_y];}
	if (c->fit_v_z) 					{p->v_z[p->i] 					= K[p->fit_Knr_v_z];}
	if (c->fit_u_y_plus_v_x) 			{p->u_y_plus_v_x[p->i]			= K[p->fit_Knr_u_y_plus_v_x];}
	if (c->fit_w0) 						{p->w0[p->i] 					= K[p->fit_Knr_w0];}
	if (c->fit_w_x) 					{p->w_x[p->i] 					= K[p->fit_Knr_w_x];}
	if (c->fit_w_y) 					{p->w_y[p->i] 					= K[p->fit_Knr_w_y];}
	if (c->fit_w_z) 					{p->w_z[p->i]					= K[p->fit_Knr_w_z];}
	if (c->fit_u_t_plus_v_t_plus_w_t) 	{p->u_t_plus_v_t_plus_w_t[p->i]	= K[p->fit_Knr_u_t_plus_v_t_plus_w_t];}
}

void retrieval_lwm_initialize_o(t_lwm_o **p_lwm_o, char measurements_filename[8192])
{
	t_lwm_o *o = malloc(sizeof(t_lwm_o));
	
	radarfilter_read_measurements(&o->n, &o->radarmeasurement, measurements_filename);
	
	func_dbl_arr_calloc(o->n, &o->windvector_u );	
	func_dbl_arr_calloc(o->n, &o->windvector_v );	
	func_dbl_arr_calloc(o->n, &o->windvector_w );	
	func_dbl_arr_calloc(o->n, &o->windvector_u_err );	
	func_dbl_arr_calloc(o->n, &o->windvector_v_err );	
	func_dbl_arr_calloc(o->n, &o->windvector_w_err );	
	
	func_dbl_arr_calloc(o->n, &o->lwm_vr );	
	func_int_arr_malloc(o->n, &o->used_for_analysis);	
	func_int_arr_malloc(o->n, &o->in_analysis_volume);	

	*p_lwm_o = o;
}

void retrieval_lwm_free_o(t_lwm_o **p_lwm_o)
{
	t_lwm_o *o = *p_lwm_o;

	if (o != NULL) {
		radarfilter_free_radarmeasurement(o->n, &o->radarmeasurement);	

		free(o->windvector_u);			
		free(o->windvector_v);
		free(o->windvector_w);

		free(o->windvector_u_err);
		free(o->windvector_v_err);
		free(o->windvector_w_err);
		
		free(o->lwm_vr);
		free(o->used_for_analysis);
		free(o->in_analysis_volume);

		free(o);
		*p_lwm_o = NULL;
	}
}

void retrieval_lwm_initialize_p(t_lwm_opc *opc)
{
	t_lwm_o *o = opc->o;
	t_lwm_p *p = malloc(sizeof(t_lwm_p));
	t_zephyros_config_retrieval_lwm_cfg *c = opc->c;
    t_zephyros_config *zcfg = (t_zephyros_config*)c->vd_zcfg;

	int io;
	int ir;
	
	opc->p = p; //forward memory allocation of p

	//allocate field
	fields_prepare2(&p->field,
		c->xvec_m[0], c->xvec_m[1], c->xvec_m[2],
		c->yvec_m[0], c->yvec_m[1], c->yvec_m[2],
		c->zvec_m[0], c->zvec_m[1], c->zvec_m[2],
		c->tvec_s[0], c->tvec_s[1], c->tvec_s[2]);

	p->Kn = 0;	//number of paremeters
	if (c->fit_u0) 						{p->fit_Knr_u0 	= p->Kn; p->Kn++;}
	if (c->fit_u_x) 					{p->fit_Knr_u_x	= p->Kn; p->Kn++;}
	if (c->fit_u_z) 					{p->fit_Knr_u_z	= p->Kn; p->Kn++;}
	if (c->fit_v0) 						{p->fit_Knr_v0	= p->Kn; p->Kn++;}
	if (c->fit_v_y) 					{p->fit_Knr_v_y	= p->Kn; p->Kn++;}
	if (c->fit_v_z) 					{p->fit_Knr_v_z	= p->Kn; p->Kn++;}
	if (c->fit_u_y_plus_v_x) 			{p->fit_Knr_u_y_plus_v_x = p->Kn; p->Kn++;}
	if (c->fit_w0) 						{p->fit_Knr_w0 	= p->Kn; p->Kn++;}
	if (c->fit_w_x) 					{p->fit_Knr_w_x = p->Kn; p->Kn++;}
	if (c->fit_w_y) 					{p->fit_Knr_w_y = p->Kn; p->Kn++;}
	if (c->fit_w_z) 					{p->fit_Knr_w_z = p->Kn; p->Kn++;}
	if (c->fit_u_t_plus_v_t_plus_w_t) 	{p->fit_Knr_u_t_plus_v_t_plus_w_t = p->Kn; p->Kn++;}
	
	func_dbl_arr_calloc(p->field->n, &p->u0 );	
	func_dbl_arr_calloc(p->field->n, &p->u_x );	
	func_dbl_arr_calloc(p->field->n, &p->u_z );	
	func_dbl_arr_calloc(p->field->n, &p->v0 );	
	func_dbl_arr_calloc(p->field->n, &p->v_y );	
	func_dbl_arr_calloc(p->field->n, &p->v_z );	
	func_dbl_arr_calloc(p->field->n, &p->u_y_plus_v_x );	
	func_dbl_arr_calloc(p->field->n, &p->w0 );	
	func_dbl_arr_calloc(p->field->n, &p->w_x );	
	func_dbl_arr_calloc(p->field->n, &p->w_y );	
	func_dbl_arr_calloc(p->field->n, &p->w_z );	
	func_dbl_arr_calloc(p->field->n, &p->u_t_plus_v_t_plus_w_t );

	p->n_sv = zcfg->retrieval->radarfilter->n_beam_range * 
				zcfg->retrieval->radarfilter->n_beam_theta *
				zcfg->retrieval->radarfilter->n_beam_phi *
				zcfg->retrieval->radarfilter->n_t;
	
	p->coef_u0 				= (double**)malloc(o->n * sizeof(double*));
	p->coef_u_x 			= (double**)malloc(o->n * sizeof(double*));
	p->coef_u_z 			= (double**)malloc(o->n * sizeof(double*));
	p->coef_v0 				= (double**)malloc(o->n * sizeof(double*));
	p->coef_v_y 			= (double**)malloc(o->n * sizeof(double*));
	p->coef_v_z 			= (double**)malloc(o->n * sizeof(double*));
	p->coef_u_y_plus_v_x	= (double**)malloc(o->n * sizeof(double*));
	p->coef_w0				= (double**)malloc(o->n * sizeof(double*));
	p->coef_w_x				= (double**)malloc(o->n * sizeof(double*));
	p->coef_w_y				= (double**)malloc(o->n * sizeof(double*));
	p->coef_w_z				= (double**)malloc(o->n * sizeof(double*));
	p->coef_u_t_plus_v_t_plus_w_t	= (double**)malloc(o->n * sizeof(double*));

	for (io = 0; io < o->n; io++ ) {
		func_dbl_arr_calloc(p->n_sv, p->coef_u0 + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_u_x + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_u_z + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_v0 + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_v_y + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_v_z + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_u_y_plus_v_x + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_w0 + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_w_x + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_w_y + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_w_z + io );
		func_dbl_arr_calloc(p->n_sv, p->coef_u_t_plus_v_t_plus_w_t + io );
	}

	func_dbl_arr_calloc(o->n, &p->center_coef_u0 );
	func_dbl_arr_calloc(o->n, &p->center_coef_u_x );
	func_dbl_arr_calloc(o->n, &p->center_coef_u_z );
	func_dbl_arr_calloc(o->n, &p->center_coef_v0 );
	func_dbl_arr_calloc(o->n, &p->center_coef_v_y );
	func_dbl_arr_calloc(o->n, &p->center_coef_v_z );
	func_dbl_arr_calloc(o->n, &p->center_coef_u_y_plus_v_x );
	func_dbl_arr_calloc(o->n, &p->center_coef_w0 );
	func_dbl_arr_calloc(o->n, &p->center_coef_w_x );
	func_dbl_arr_calloc(o->n, &p->center_coef_w_y );
	func_dbl_arr_calloc(o->n, &p->center_coef_w_z );
	func_dbl_arr_calloc(o->n, &p->center_coef_u_t_plus_v_t_plus_w_t );

	p->subvolume_coor			= malloc(o->n * sizeof(t_zephyros_coordinates**));
	for (io = 0; io < o->n; io++ ) {
		p->subvolume_coor[io] = malloc(p->n_sv * sizeof(t_zephyros_coordinates*));
		for (ir = 0; ir < p->n_sv; ir++ ) {
			coordinates_initialize_coor(&(p->subvolume_coor[io][ir]));
		}
	}
}
	
void retrieval_lwm_free_p(t_lwm_opc *opc)
{	
	t_lwm_o *o = opc->o;
	t_lwm_p *p = opc->p;
	
	int io;
	int ir;
	
	fields_free(&p->field);
	
	free(p->u0);
	free(p->u_x);
	free(p->u_z);
	free(p->v0);
	free(p->v_y);
	free(p->v_z);
	free(p->u_y_plus_v_x);
	free(p->w0);
	free(p->w_x);
	free(p->w_y);
	free(p->w_z);
	free(p->u_t_plus_v_t_plus_w_t);
	
	for (io = 0; io < o->n; io++ ) {
		free(p->coef_u0[io]);
		free(p->coef_u_x[io]);
		free(p->coef_u_z[io]);
		free(p->coef_v0[io]);
		free(p->coef_v_y[io]);
		free(p->coef_v_z[io]);
		free(p->coef_u_y_plus_v_x[io]);
		free(p->coef_w0[io]);
		free(p->coef_w_x[io]);
		free(p->coef_w_y[io]);
		free(p->coef_w_z[io]);
		free(p->coef_u_t_plus_v_t_plus_w_t[io]);
	}

	free(p->coef_u0);
	free(p->coef_u_x);
	free(p->coef_u_z);
	free(p->coef_v0);
	free(p->coef_v_y);
	free(p->coef_v_z);
	free(p->coef_u_y_plus_v_x);
	free(p->coef_w0);
	free(p->coef_w_x);
	free(p->coef_w_y);
	free(p->coef_w_z);
	free(p->coef_u_t_plus_v_t_plus_w_t);

	free(p->center_coef_u0);
	free(p->center_coef_u_x);
	free(p->center_coef_u_z);
	free(p->center_coef_v0);
	free(p->center_coef_v_y);
	free(p->center_coef_v_z);
	free(p->center_coef_u_y_plus_v_x);
	free(p->center_coef_w0);
	free(p->center_coef_w_x);
	free(p->center_coef_w_y);
	free(p->center_coef_w_z);
	free(p->center_coef_u_t_plus_v_t_plus_w_t);

	for (io = 0; io < o->n; io++ ) {
		for (ir = 0; ir < p->n_sv; ir++ ) {
			coordinates_free_coor(&(p->subvolume_coor[io][ir]));
		}
		free(p->subvolume_coor[io]);
	}
	free(p->subvolume_coor);
	
	free(p);
}


