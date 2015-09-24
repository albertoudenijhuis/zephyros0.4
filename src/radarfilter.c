//next:
//- add spectrum_eta_i


/*
Description: 
	RADAR filter functions

Revision History:
	2014

Functions:
	RADAR filter functions
	
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "radarfilter.h"

#include "interpolation.h"
#include "ltqnorm.h"
#include "util.h"
#include "func.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

//uncomment next statement for debug mode
//#define _ZEPHYROS_RADARFILTER_DEBUG

void radarfilter_exec(
	t_zephyros_config 				*cfg,
	int								i_mode,	//0 = simulation mode, 1 = retrieval mode
	t_radarfilter_todolist			*todo,	
	int								n_measurements,
	t_radarmeasurement				**radarmeasurement
)
{
	int i;		//used now and then for loops
	int i_m; 	//i measurement
	int i_psd;	//i particle side distribution
	int i_par; 	//i particle
	int i_res;	//i resolution subvolume
	int i_parmod;  //i parametric model

	t_radarfilter_res_vol 		*res_vol;
	t_zephyros_radarfilter		*rcfg;
	t_zephyros_windfield		*mywindfield;
	t_zephyros_scattererfield	*myscattererfield;

	//needed for calculations
	double eta2Z;
	
	double eta_hh;
	double eta_hv;
	double eta_vh;
	double eta_vv;

	double eta_ShhSvvc;
	double eta_ShhShvc;
	double eta_SvvSvhc;
	
	double n_Re_Shh_min_Svv;

	double *Doppler_spectrum_ShhSvvc = NULL;
	double *Doppler_spectrum_ShhShvc = NULL;
	double *Doppler_spectrum_SvvSvhc = NULL;

	int 		remainder;
	int			i_beam_range, i_beam_theta, i_beam_phi, i_t;
	double 		myr1, myr2;
	double 		FWHM_rad;
	double 		mincdfP, maxcdfP, difcdfP, cdfP, sigma_gauss;

	double tmpvar;
	double tmp;
	double *dummy;

	double tmp_turbulence_u;
	double tmp_turbulence_v;
	double tmp_turbulence_w;
	
	double tmp_total_u;
	double tmp_total_v;
	double tmp_total_w;

	double tmpa, tmpb, tmpwindspeed, tmpL;
	double tmpedr, tmpedrplus;
	double tmpIsqrt;
	double tmpkolmogorovconstant;
	double tmparr[10];
	int tmpi;
	int i_parmod_az, i_parmod_el;
	double parmod_az, parmod_el;

	int i_int;	//i for spectral intervals

	double epsilon_r = creal(cfg->derived_quantities->radar_water_refractive_index * cfg->derived_quantities->radar_water_refractive_index);
	double kwsq = pow((epsilon_r - 1.) / (epsilon_r - 2.), 2.);
	//translates eta to Z

	double lbound, ubound, center;

	//related to parametric model of turbulence
	double parametric_turbulence_sigmaT;
	double parametric_turbulence_sigmaT_plus; //for calculation of derivatives via finite difference
	
	int traj_n, i_traj;
	double traj_dxy, traj_dz, traj_d, traj_dt;
	double traj_u[100], traj_v[100], traj_w[100];
	double traj_xyzt[100][4];
	double traj_delta_u, traj_delta_v, traj_delta_w;

	double tmpdu, tmpdv, tmpdw, sgndu, sgndv, sgndw;
	
	//needed to calculate derivative to edr^(1/3)
	double		eta_hh_plus;
	double		eta_hv_plus;
	double		eta_vh_plus;
	double		eta_vv_plus;
	
	double		eta_ShhSvvc_plus;
	double		eta_ShhShvc_plus;
	double		eta_SvvSvhc_plus;

	double		Doppler_velocity_hh_ms_plus;
	double		Doppler_spectral_width_hh_ms_plus;
	
	
	
	//int debug = 0;
	//clock_t start, diff;
	
	//point radar filter configuration in the right way
	if (i_mode == 0) {
		rcfg = cfg->simulation->radarfilter;
		mywindfield = cfg->simulation->windfield;
		myscattererfield = cfg->simulation->scattererfield;
	}
	if (i_mode == 1) {
		rcfg = cfg->retrieval->radarfilter;
		mywindfield = cfg->retrieval->post_windfield;
		myscattererfield = cfg->retrieval->post_scattererfield;		
	}

	//initialize resolution volume
	#ifdef _ZEPHYROS_RADARFILTER_DEBUG
		printf("radarfilter_initialize_resolution_volume\n"); fflush(stdout);
	#endif 
	radarfilter_initialize_resolution_volume(cfg, i_mode, &res_vol, todo);

	//*****
	//main loop over the measurements
	printf("calculate measurement %5i/%5i\n", 0, 0); fflush(stdout);
	for (i_m = 0; i_m < n_measurements; i_m++ ) {
		
		radarmeasurement[i_m]->coef_eta2Z = 
			1.e18 * pow(cfg->derived_quantities->central_wavelength_m, 4.) / (pow(M_PI, 5.) * kwsq);
		eta2Z = radarmeasurement[i_m]->coef_eta2Z;
		
		tmpi = (n_measurements / 10); if (tmpi == 0) tmpi = 1;
		if (	(i_m < 10) |
				( (i_m % tmpi) == 0)) {
			for (i = 0;  i < 35;  i++, printf("\b")); printf("calculate measurement %5i/%5i\n", i_m + 1, n_measurements); fflush(stdout);
		}
		
		//****
		//take over psd structure from configuration
		radarmeasurement[i_m]->n_psd = res_vol->n_psd;
		util_safe_free(&(radarmeasurement[i_m]->n_diameters));		
		radarmeasurement[i_m]->n_diameters = malloc(res_vol->n_psd * sizeof(int));
		for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
			radarmeasurement[i_m]->n_diameters[i_psd] = res_vol->n_diameters[i_psd];
		}
		
		//set spectrum size
		if (i_mode == 0) radarmeasurement[i_m]->n_spectrum = rcfg->n_spectrum;
		
		//*****
		//calculation of radar center resolution volume coordinates
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("calculation of radar center resolution volume coordinates\n"); fflush(stdout);
		#endif 
		if (todo->calc_center_coordinates) {
			myr1	= radarmeasurement[i_m]->azel_r1_m;
			myr2	= radarmeasurement[i_m]->azel_r2_m;
			if (myr1 == 0.) {
				myr1 = myr2 / res_vol->n_beam_range;
			}
			radarmeasurement[i_m]->center_coor->radar_range = pow(pow(myr1 , -1.) - (0.5 * (pow(myr1 , -1.) - pow(myr2, -1.))), -1.);		
			radarmeasurement[i_m]->center_coor->radar_azel_alpha = radarmeasurement[i_m]->azel_alpha_rad;
			radarmeasurement[i_m]->center_coor->radar_azel_gamma = radarmeasurement[i_m]->azel_gamma_rad;

			coordinates_radar_azel2enu(radarmeasurement[i_m]->center_coor);
			coordinates_radar_azelrangedir2enu(radarmeasurement[i_m]->center_coor);
			coordinates_radar_pol_dir(radarmeasurement[i_m]->center_coor);
		}

		//*****
		//free and allocate memory for spectra
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("free and allocate memory for spectra\n"); fflush(stdout);
		#endif 
		util_safe_free(&(radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh));
		util_safe_free(&(radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh_err));
		util_safe_free(&(radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv));
		util_safe_free(&(radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv_err));
		util_safe_free(&(radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh));
		util_safe_free(&(radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh_err));
		util_safe_free(&(radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv));
		util_safe_free(&(radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv_err));
		util_safe_free(&(radarmeasurement[i_m]->specific_dBZdr));
		util_safe_free(&(radarmeasurement[i_m]->specific_dBZdr_err));
		util_safe_free(&(radarmeasurement[i_m]->specific_dBLdr));
		util_safe_free(&(radarmeasurement[i_m]->specific_dBLdr_err));
		util_safe_free(&(radarmeasurement[i_m]->specific_rho_co));
		util_safe_free(&(radarmeasurement[i_m]->specific_rho_co_err));
		util_safe_free(&(radarmeasurement[i_m]->specific_rho_cxh));
		util_safe_free(&(radarmeasurement[i_m]->specific_rho_cxh_err));
		util_safe_free(&(radarmeasurement[i_m]->specific_rho_cxv));
		util_safe_free(&(radarmeasurement[i_m]->specific_rho_cxv_err));

		util_safe_free(&(radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hh));
		util_safe_free(&(radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hv));
		util_safe_free(&(radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vh));
		util_safe_free(&(radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vv));
					
		if (rcfg->filter_Doppler_spectrum_dBZ_hh) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh);
		if ((rcfg->filter_Doppler_spectrum_dBZ_hh) & (rcfg->filter_errors)) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh_err); 		
		if (rcfg->filter_Doppler_spectrum_dBZ_hv) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv);
		if ((rcfg->filter_Doppler_spectrum_dBZ_hv) & (rcfg->filter_errors)) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv_err);
		if (rcfg->filter_Doppler_spectrum_dBZ_vh) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh);
		if ((rcfg->filter_Doppler_spectrum_dBZ_vh) & (rcfg->filter_errors)) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh_err);
		if (rcfg->filter_Doppler_spectrum_dBZ_vv) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv);
		if ((rcfg->filter_Doppler_spectrum_dBZ_vv) & (rcfg->filter_errors)) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv_err);
		if (todo->calc_Doppler_spectrum_ShhSvvc) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &Doppler_spectrum_ShhSvvc);
		if (todo->calc_Doppler_spectrum_ShhShvc) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &Doppler_spectrum_ShhShvc);
		if (todo->calc_Doppler_spectrum_SvvSvhc) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &Doppler_spectrum_SvvSvhc);		
		if (rcfg->filter_specific_dBZdr) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_dBZdr);
		if ((rcfg->filter_specific_dBZdr)  & (rcfg->filter_errors)) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_dBZdr_err);
		if (rcfg->filter_specific_dBLdr) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_dBLdr);
		if ((rcfg->filter_specific_dBLdr)  & (rcfg->filter_errors)) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_dBLdr_err);
		if (rcfg->filter_specific_rho_co) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_rho_co);
		if ((rcfg->filter_specific_rho_co)  & (rcfg->filter_errors)) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_rho_co_err);
		if (rcfg->filter_specific_rho_cxh) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_rho_cxh);
		if ((rcfg->filter_specific_rho_cxh)  & (rcfg->filter_errors)) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_rho_cxh_err);
		if (rcfg->filter_specific_rho_cxv) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_rho_cxv);
		if ((rcfg->filter_specific_rho_cxv)  & (rcfg->filter_errors)) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->specific_rho_cxv_err);
		
		if (todo->der_edr13) {
			if (rcfg->filter_Doppler_spectrum_dBZ_hh) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hh);
			if (rcfg->filter_Doppler_spectrum_dBZ_hv) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hv);
			if (rcfg->filter_Doppler_spectrum_dBZ_vh) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vh);
			if (rcfg->filter_Doppler_spectrum_dBZ_vv) func_dbl_arr_calloc(radarmeasurement[i_m]->n_spectrum, &radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vv);
		}
		
		//free and allocate memory for eta_i
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("free and allocate memory for eta_i\n"); fflush(stdout);
		#endif 
		if (radarmeasurement[i_m]->eta_i_hh != NULL) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				util_safe_free(&(radarmeasurement[i_m]->eta_i_hh[i_psd]));
			}
			util_safe_free(&(radarmeasurement[i_m]->eta_i_hh));
		}
		if (radarmeasurement[i_m]->eta_i_hv != NULL) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				util_safe_free(&(radarmeasurement[i_m]->eta_i_hv[i_psd]));
			}
			util_safe_free(&(radarmeasurement[i_m]->eta_i_hv));
		}
		if (radarmeasurement[i_m]->eta_i_vh != NULL) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				util_safe_free(&(radarmeasurement[i_m]->eta_i_vh[i_psd]));
			}
			util_safe_free(&(radarmeasurement[i_m]->eta_i_vh));
		}
		if (radarmeasurement[i_m]->eta_i_vv != NULL) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				util_safe_free(&(radarmeasurement[i_m]->eta_i_vv[i_psd]));
			}
			util_safe_free(&(radarmeasurement[i_m]->eta_i_vv));
		}

		if (todo->calc_eta_i_hh) {
			radarmeasurement[i_m]->eta_i_hh = malloc(radarmeasurement[i_m]->n_psd * sizeof(double*));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) 
				radarmeasurement[i_m]->eta_i_hh[i_psd] = malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double));
		}
		if (todo->calc_eta_i_hv) {
			radarmeasurement[i_m]->eta_i_hv = malloc(radarmeasurement[i_m]->n_psd * sizeof(double*));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) 
				radarmeasurement[i_m]->eta_i_hv[i_psd] = malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double));
		}
		if (todo->calc_eta_i_vh) {
			radarmeasurement[i_m]->eta_i_vh = malloc(radarmeasurement[i_m]->n_psd * sizeof(double*));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) 
				radarmeasurement[i_m]->eta_i_vh[i_psd] = malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double));
		}
		if (todo->calc_eta_i_vv) {
			radarmeasurement[i_m]->eta_i_vv = malloc(radarmeasurement[i_m]->n_psd * sizeof(double*));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) 
				radarmeasurement[i_m]->eta_i_vv[i_psd] = malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double));
		}
		
		//free and allocate memory for spectrum_eta_i
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("free and allocate memory for spectrum_eta_i\n"); fflush(stdout);
		#endif 
		if (radarmeasurement[i_m]->spectrum_eta_i_hh != NULL) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ )
					util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_hh[i_psd][i_par]));
				util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_hh[i_psd]));
			}
			util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_hh));
		}
		if (radarmeasurement[i_m]->spectrum_eta_i_hv != NULL) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ )
					util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_hv[i_psd][i_par]));
				util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_hv[i_psd]));
			}
			util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_hv));
		}
		if (radarmeasurement[i_m]->spectrum_eta_i_vh != NULL) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ )
					util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_vh[i_psd][i_par]));
				util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_vh[i_psd]));
			}
			util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_vh));
		}
		if (radarmeasurement[i_m]->spectrum_eta_i_vv != NULL) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ )
					util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_vv[i_psd][i_par]));
				util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_vv[i_psd]));
			}
			util_safe_free(&(radarmeasurement[i_m]->spectrum_eta_i_vv));
		}
		
		if (todo->calc_spectrum_eta_i_hh) {
			radarmeasurement[i_m]->spectrum_eta_i_hh = malloc(radarmeasurement[i_m]->n_psd * sizeof(double**));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				radarmeasurement[i_m]->spectrum_eta_i_hh[i_psd] = malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double*));
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ )
					radarmeasurement[i_m]->spectrum_eta_i_hh[i_psd][i_par] = 
						malloc(radarmeasurement[i_m]->n_spectrum * sizeof(double));
			}
		}
		if (todo->calc_spectrum_eta_i_hv) {
			radarmeasurement[i_m]->spectrum_eta_i_hv = malloc(radarmeasurement[i_m]->n_psd * sizeof(double**));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				radarmeasurement[i_m]->spectrum_eta_i_hv[i_psd] = malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double*));
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ )
					radarmeasurement[i_m]->spectrum_eta_i_hv[i_psd][i_par] = 
						malloc(radarmeasurement[i_m]->n_spectrum * sizeof(double));
			}
		}
		if (todo->calc_spectrum_eta_i_vh) {
			radarmeasurement[i_m]->spectrum_eta_i_vh = malloc(radarmeasurement[i_m]->n_psd * sizeof(double**));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				radarmeasurement[i_m]->spectrum_eta_i_vh[i_psd] = malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double*));
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ )
					radarmeasurement[i_m]->spectrum_eta_i_vh[i_psd][i_par] = 
						malloc(radarmeasurement[i_m]->n_spectrum * sizeof(double));
			}
		}
		if (todo->calc_spectrum_eta_i_vv) {
			radarmeasurement[i_m]->spectrum_eta_i_vv = malloc(radarmeasurement[i_m]->n_psd * sizeof(double**));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				radarmeasurement[i_m]->spectrum_eta_i_vv[i_psd] = malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double*));
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ )
					radarmeasurement[i_m]->spectrum_eta_i_vv[i_psd][i_par] = 
						malloc(radarmeasurement[i_m]->n_spectrum * sizeof(double));
			}
		}
		
		//*****
		//calculate subvolume coordinates
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("calculate subvolume coordinates\n"); fflush(stdout);
		#endif 
		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			//take over radar location
			res_vol->subvolume_coor[i_res]->enu_radar_location_xyzt[0] = radarmeasurement[i_m]->center_coor->enu_radar_location_xyzt[0];
			res_vol->subvolume_coor[i_res]->enu_radar_location_xyzt[1] = radarmeasurement[i_m]->center_coor->enu_radar_location_xyzt[1];
			res_vol->subvolume_coor[i_res]->enu_radar_location_xyzt[2] = radarmeasurement[i_m]->center_coor->enu_radar_location_xyzt[2];
			res_vol->subvolume_coor[i_res]->enu_radar_location_xyzt[3] = radarmeasurement[i_m]->center_coor->enu_radar_location_xyzt[3];
			//TBD: account for variation in time

			//get indices for range, phi, theta and time
			//i_beam_range 	0 ... n_beam_range -1
			//i_beam_theta 	0 ... n_beam_theta -1
			//i_beam_phi 	0 ... n_beam_phi -1 
			//i_t		 	0 ... n_t -1

			i_beam_range	= i_res % res_vol->n_beam_range;
			remainder 		= (i_res - i_beam_range) / res_vol->n_beam_range;
			i_beam_theta	= remainder % res_vol->n_beam_theta;
			remainder 		= (remainder - i_beam_theta) / res_vol->n_beam_theta;
			i_beam_phi		= remainder % res_vol->n_beam_phi;
			remainder 		= (remainder - i_beam_phi) / res_vol->n_beam_phi;
			i_t				= remainder % res_vol->n_t;
			
			//beam range is chosen such that equal power is there.
			//for volumetric distributed scatterers the weighting factor is (r^-2)
			cdfP = (i_beam_range + 0.5) / res_vol->n_beam_range; 
			myr1	= radarmeasurement[i_m]->azel_r1_m;
			myr2	= radarmeasurement[i_m]->azel_r2_m;
			if (myr1 == 0.) {
				myr1 = myr2 / res_vol->n_beam_range;
			}
			res_vol->subvolume_coor[i_res]->radar_range = 
				pow(  pow(myr1 , -1.) - (cdfP * (pow(myr1 , -1.) - pow(myr2, -1.))), -1.);		
			
			//calculate point_beam_phi_rad
			res_vol->subvolume_coor[i_res]->radar_beam_phi =
					(2. * M_PI * i_beam_phi ) / res_vol->n_beam_phi;

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
						pow(cos(res_vol->subvolume_coor[i_res]->radar_beam_phi) * radarmeasurement[i_m]->beam_FWHM0_rad,2) +
						pow(sin(res_vol->subvolume_coor[i_res]->radar_beam_phi) * radarmeasurement[i_m]->beam_FWHM1_rad,2)
						, 1./2.);

			//beam theta is chosen such that equal power is there.
			if ((res_vol->n_beam_theta == 1) & (res_vol->n_beam_phi == 1))
			{			
				//special case, n_beam_theta = 1, n_beam_phi = 1
				res_vol->subvolume_coor[i_res]->radar_beam_theta = 0.;
			} else {					
				//examples:
				//n_beam_theta = 1, cdfP = [0.75]
				//n_beam_theta = 2, cdfP = [0.625, 875]
				//n_beam_theta = 4, cdfP = [0.5625, 0.6875, 0.8125, 0.9375]

				mincdfP = 0.5	;
				maxcdfP = 1.	;
				difcdfP = (maxcdfP - mincdfP) / res_vol->n_beam_theta;
				cdfP = mincdfP + (i_beam_theta + 0.5) * difcdfP; 

				sigma_gauss = FWHM_rad / sqrt(64. * log(2.));
				res_vol->subvolume_coor[i_res]->radar_beam_theta = sigma_gauss * ltqnorm(cdfP);
			}

			//beam to azel
			coordinates_radar_beam2azel(res_vol->subvolume_coor[i_res]);
			
			//add central values of the azel sytem
			res_vol->subvolume_coor[i_res]->radar_azel_alpha 		+= radarmeasurement[i_m]->azel_alpha_rad;
			res_vol->subvolume_coor[i_res]->radar_azel_gamma 		+= radarmeasurement[i_m]->center_coor->radar_azel_gamma;
			res_vol->subvolume_coor[i_res]->radar_azel_gamma_cor 	+= radarmeasurement[i_m]->center_coor->radar_azel_gamma_cor;
			
			//azel -> enu
			coordinates_radar_azel2enu(res_vol->subvolume_coor[i_res]);
			coordinates_radar_azelrangedir2enu(res_vol->subvolume_coor[i_res]);
			coordinates_radar_pol_dir(res_vol->subvolume_coor[i_res]);
		}
		
		//to calculate the cross sections, we need the orientation of the particle first
		//1. calculate air velocity
		//2. calculate terminal fall speeds
		//3. calculate particle velocities
		//4. calculate particle orientations
		//5. calculate cross sections

		//*****
		//calculate air velocity
		if (todo->calc_air_velocity == 1) {
			//for whole resolution volume
			util_windfield_fuvw(mywindfield, radarmeasurement[i_m]->center_coor->enu_xyzt, 0, &(res_vol->air_u), todo->calc_air_velocity_der, res_vol->air_u_der);
			util_windfield_fuvw(mywindfield, radarmeasurement[i_m]->center_coor->enu_xyzt, 1, &(res_vol->air_v), todo->calc_air_velocity_der, res_vol->air_v_der);
			util_windfield_fuvw(mywindfield, radarmeasurement[i_m]->center_coor->enu_xyzt, 2, &(res_vol->air_w), todo->calc_air_velocity_der, res_vol->air_w_der);
			
			//for subvolume
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				util_windfield_fuvw(mywindfield, res_vol->subvolume_coor[i_res]->enu_xyzt, 0, res_vol->subvolume_air_u + i_res, todo->calc_air_velocity_der, res_vol->subvolume_air_u_der[i_res]);
				util_windfield_fuvw(mywindfield, res_vol->subvolume_coor[i_res]->enu_xyzt, 1, res_vol->subvolume_air_v + i_res, todo->calc_air_velocity_der, res_vol->subvolume_air_u_der[i_res]);
				util_windfield_fuvw(mywindfield, res_vol->subvolume_coor[i_res]->enu_xyzt, 2, res_vol->subvolume_air_w + i_res, todo->calc_air_velocity_der, res_vol->subvolume_air_u_der[i_res]);
			}
		}
				
		//calculate terminal fall speeds (& related parameters)
		if (todo->calc_terminal_fall_speeds == 1) {
			//obtain T_K
			interpolation_bilint(cfg->general->atmosphere->lut_T_K, radarmeasurement[i_m]->center_coor->enu_xyzt, &tmp, 0, dummy);
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					res_vol->subvolume_scat[i_psd][i_par]->air_temperature_K = tmp;
				}
			}
			
			//obtain air_drypressure_hPa
			interpolation_bilint(cfg->general->atmosphere->lut_ln_grid_pair_hPa, radarmeasurement[i_m]->center_coor->enu_xyzt, &tmp, 0, dummy);
			tmp = exp(tmp);
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					res_vol->subvolume_scat[i_psd][i_par]->air_drypressure_hPa = tmp;
				}
			}

			//obtain air_vaporpressure_hPa
			interpolation_bilint(cfg->general->atmosphere->lut_ln_grid_pvapor_hPa, radarmeasurement[i_m]->center_coor->enu_xyzt, &tmp, 0, dummy);
			tmp = exp(tmp);
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					res_vol->subvolume_scat[i_psd][i_par]->air_vaporpressure_hPa = tmp;
				}
			}

			//terminal fall speed (and related parameters)
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					res_vol->subvolume_scat[i_psd][i_par]->radar_wavelength_m = cfg->derived_quantities->central_wavelength_m;
					
					particles_air_parameters(res_vol->subvolume_scat[i_psd][i_par], radarmeasurement[i_m]->center_coor);

					res_vol->subvolume_scat[i_psd][i_par]->particle_type		= myscattererfield->psd[i_psd]->particle_type;
					res_vol->subvolume_scat[i_psd][i_par]->particle_D_eqvol_mm	= myscattererfield->psd[i_psd]->discrete_D_equiv_mm[i_par];
					particles_spheroid_geometry_beard1987(res_vol->subvolume_scat[i_psd][i_par]);
					particles_terminal_fall_speed_khvorostyanov2005(res_vol->subvolume_scat[i_psd][i_par]);
					//debug particle_print_widget(res_vol->subvolume_scat[i_psd][i_par]);
				}
			}
		}

		//Prepare parametric turbulence model
		//Calculate parametric_turbulence_sigmaT
		if (todo->calc_particle_velocity) {
			if (res_vol->parametric_turbulence == 1) {
				parametric_turbulence_sigmaT = 0.;
				if (todo->der_edr13) parametric_turbulence_sigmaT_plus = 0.;
				for ( i = 0; i < mywindfield->nturbulences; i++ ) {
					if ((mywindfield->turbulence[i] != NULL) & (mywindfield->turbulence[i]->type == 5)) {						
						//interpolate EDR
						interpolation_bilint(mywindfield->turbulence[i]->lut_edr13,
								radarmeasurement[i_m]->center_coor->enu_xyzt,
								&tmp,
								0, //no derivatives
								dummy);
						if (tmp < 0.) tmp = 0.;
						tmpedr 		= pow(tmp, 3.);
						tmpedrplus  = pow(tmp * 1.001, 3.);
						interpolation_bilint(mywindfield->turbulence[i]->lut_kolmogorov_constant,
								radarmeasurement[i_m]->center_coor->enu_xyzt,
								&tmpkolmogorovconstant,
								0, //no derivatives
								dummy);						
						
						//interpolate white integral
						tmpa = radarmeasurement[i_m]->center_coor->radar_range * 
							sqrt( (radarmeasurement[i_m]->beam_FWHM0_rad * radarmeasurement[i_m]->beam_FWHM1_rad)
							/ (64. * log(2.)));
						tmpb = fabs(radarmeasurement[i_m]->azel_r2_m - radarmeasurement[i_m]->azel_r1_m) / 2.;
						tmpwindspeed = sqrt(pow(res_vol->air_u, 2.) + pow(res_vol->air_v, 2.) + pow(res_vol->air_w, 2.));
						tmpL = tmpwindspeed * radarmeasurement[i_m]->dt;
						
						tmparr[0] = log(fmax(1.e-2, tmpa));
						tmparr[1] = log(fmax(1.e-2, tmpb));
						tmparr[2] = log(fmax(1.e-2, tmpL));
						interpolation_bilint(cfg->general->white1999_integral->lut_integral_sqrt,
								tmparr,
								&tmpIsqrt,
								0, //no derivatives
								dummy);						
						
						//calculate sigmaT, and add up
						parametric_turbulence_sigmaT = sqrt(
							pow(parametric_turbulence_sigmaT, 2.) + 
							(((tmpkolmogorovconstant * pow(tmpedr, 2./3.)) / (4. * M_PI)) * pow(tmpIsqrt,2.)));
													
						if (todo->der_edr13) {
							parametric_turbulence_sigmaT_plus = sqrt(
							pow(parametric_turbulence_sigmaT, 2.) + 
							(((tmpkolmogorovconstant * pow(tmpedrplus, 2./3.)) / (4. * M_PI)) * pow(tmpIsqrt,2.)));
						}
							
					}
				}
			}
		}
				
		//calculate particle velocities
		if (todo->calc_particle_velocity) {
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
					for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
						if (res_vol->parametric_turbulence == 0) {
							//no parametric turbulence
							//assumption is that v_par = v_air + v_terminal
							res_vol->subvolume_particle_u[i_res][i_psd][i_par][0] = res_vol->subvolume_air_u[i_res];
							res_vol->subvolume_particle_v[i_res][i_psd][i_par][0] = res_vol->subvolume_air_v[i_res];
							res_vol->subvolume_particle_w[i_res][i_psd][i_par][0] = res_vol->subvolume_air_w[i_res] - res_vol->subvolume_scat[i_psd][i_par]->particle_terminal_fall_speed;

							if (rcfg->inertia_effect) {
								//account for inertia effect
								//assumption is that v_par = v_air + v_terminal + v_par^prime
								//where v_par^prime is solved via a small trajectorie.
							
								//solve for atrajectorie
								//1. define the trajectorie
								traj_n 		= 20;
								traj_dxy 	= 0.25 * res_vol->subvolume_scat[i_psd][i_par]->particle_inertial_distance_xy;
								traj_dz 	= 0.25 * res_vol->subvolume_scat[i_psd][i_par]->particle_inertial_distance_z;
								traj_d 		= sqrt(2. * pow(traj_dxy,2.) + pow(traj_dz, 2.));

								i_traj = traj_n - 1;
								traj_u[i_traj] = res_vol->subvolume_air_u[i_res];
								traj_v[i_traj] = res_vol->subvolume_air_v[i_res];
								traj_w[i_traj] = res_vol->subvolume_air_w[i_res] - res_vol->subvolume_scat[i_psd][i_par]->particle_terminal_fall_speed;
								traj_xyzt[i_traj][0] = res_vol->subvolume_coor[i_res]->enu_xyzt[0];
								traj_xyzt[i_traj][1] = res_vol->subvolume_coor[i_res]->enu_xyzt[1];
								traj_xyzt[i_traj][2] = res_vol->subvolume_coor[i_res]->enu_xyzt[2];
								traj_xyzt[i_traj][3] = res_vol->subvolume_coor[i_res]->enu_xyzt[3];
								for ( i_traj = traj_n - 2; i_traj >= 0; i_traj-- ) {
									traj_dt = traj_d / sqrt(pow(traj_u[i_traj+1], 2.) + pow(traj_v[i_traj+1], 2.) + pow(traj_w[i_traj+1], 2.));
									traj_xyzt[i_traj][0] = traj_xyzt[i_traj+1][0] - (traj_dt * traj_u[i_traj+1]);
									traj_xyzt[i_traj][1] = traj_xyzt[i_traj+1][1] - (traj_dt * traj_v[i_traj+1]);
									traj_xyzt[i_traj][2] = traj_xyzt[i_traj+1][2] - (traj_dt * traj_w[i_traj+1]);
									traj_xyzt[i_traj][3] = traj_xyzt[i_traj+1][3] - traj_dt;
									util_windfield_fuvw(mywindfield, traj_xyzt[i_traj], 0, traj_u + i_traj, 0, dummy);
									util_windfield_fuvw(mywindfield, traj_xyzt[i_traj], 1, traj_v + i_traj, 0, dummy);
									util_windfield_fuvw(mywindfield, traj_xyzt[i_traj], 2, traj_w + i_traj, 0, dummy);
									traj_w[i_traj] -= res_vol->subvolume_scat[i_psd][i_par]->particle_terminal_fall_speed;
								}

								//2. integrate to the solution
								traj_delta_u = 0.;
								traj_delta_v = 0.;
								traj_delta_w = 0.;									
								for ( i_traj = 0; i_traj < traj_n - 1; i_traj++ ) {
									tmpdu = traj_u[i_traj+1] - (traj_delta_u + traj_u[i_traj]);
									tmpdv = traj_v[i_traj+1] - (traj_delta_v + traj_v[i_traj]);
									tmpdw = traj_w[i_traj+1] - (traj_delta_w + traj_w[i_traj]);
									sgndu = tmpdu / fabs(tmpdu);
									sgndv = tmpdv / fabs(tmpdv);
									sgndw = tmpdw / fabs(tmpdw);
									traj_dt = (traj_xyzt[i_traj + 1][3] - traj_xyzt[i_traj][3]);
									

									tmpvar = (-1. / (2. * res_vol->subvolume_scat[i_psd][i_par]->particle_inertial_eta_xy * traj_dt));
									traj_delta_u = 
										sgndu * 
										( tmpvar +
											sqrt(pow(tmpvar,2.) + fabs(tmpdu) / (res_vol->subvolume_scat[i_psd][i_par]->particle_inertial_eta_xy * traj_dt))
										);
									traj_delta_v = 
										sgndv * 
										( tmpvar +
											sqrt(pow(tmpvar,2.) + fabs(tmpdv) / (res_vol->subvolume_scat[i_psd][i_par]->particle_inertial_eta_xy * traj_dt))
										);
								
									tmpvar = (-1. * res_vol->subvolume_scat[i_psd][i_par]->particle_terminal_fall_speed)
											+ (-1. / (2. * res_vol->subvolume_scat[i_psd][i_par]->particle_inertial_eta_z * traj_dt));
									traj_delta_w = 
										sgndw * 
										( tmpvar + 
											sqrt(pow(tmpvar,2.) + (fabs(tmpdw) / (res_vol->subvolume_scat[i_psd][i_par]->particle_inertial_eta_z * traj_dt)))
										);									
								}

								//fix to zeros.
								if (isnanorinf(&traj_delta_u)) traj_delta_u = 0.;								
								if (isnanorinf(&traj_delta_v)) traj_delta_v = 0.;								
								if (isnanorinf(&traj_delta_w)) traj_delta_w = 0.;	
															
								//add to the solution
								res_vol->subvolume_particle_u[i_res][i_psd][i_par][0] += traj_delta_u;
								res_vol->subvolume_particle_v[i_res][i_psd][i_par][0] += traj_delta_v;
								res_vol->subvolume_particle_w[i_res][i_psd][i_par][0] += traj_delta_w;		

								//debugging
								/*
								if (i_par == (res_vol->n_diameters[i_psd] - 1)) {
									particle_print_widget(res_vol->subvolume_scat[i_psd][i_par]);
									for ( i_traj = 0; i_traj < traj_n; i_traj++ ) {
										printf("traj_xyzt[%i][0] = %.2e\n",i_traj, traj_xyzt[i_traj][0]);
										printf("traj_xyzt[%i][1] = %.2e\n",i_traj, traj_xyzt[i_traj][1]);
										printf("traj_xyzt[%i][2] = %.2e\n",i_traj, traj_xyzt[i_traj][2]);
										printf("traj_xyzt[%i][3] = %.2e\n",i_traj, traj_xyzt[i_traj][3]);
										printf("traj_u[%i] = %.2e\n",i_traj, traj_u[i_traj]);
										printf("traj_v[%i] = %.2e\n",i_traj, traj_v[i_traj]);
										printf("traj_w[%i] = %.2e\n",i_traj, traj_w[i_traj]);
									}								
									printf("traj_delta_u = %.2e\n", traj_delta_u);
									printf("traj_delta_v = %.2e\n", traj_delta_v);
									printf("traj_delta_w = %.2e\n", traj_delta_w);
									printf("\n\n\n");
									exit(0);
								}
								*/

							}
						}
						
						if (res_vol->parametric_turbulence == 1) {
							//parametric turbulence
							//assumption is that v_par = v_air + v_terminal + v_turb
							//particle orientation for the parametric turbulence model
							for ( i_parmod = 0; i_parmod < res_vol->n_parmod; i_parmod++ ) {
								i_parmod_az		= i_parmod % res_vol->n_parmod_az;
								remainder 		= (i_parmod - i_parmod_az) / res_vol->n_parmod_az;
								i_parmod_el		= remainder % res_vol->n_parmod_el;
													
								parmod_az = (2. * M_PI * i_parmod_az ) / res_vol->n_parmod_az;
								tmp = (i_parmod_el + 0.5) / res_vol->n_parmod_el;
								parmod_el = asin((2. * x) - 1.);

								//tmp = -1.0 + (2. * (i_parmod_el + 0.5) / res_vol->n_parmod_el);
								//parmod_el = ((tmp > 0.)?1.:-1.) * acos(fabs(tmp));

								//calculate turbulence vector
								tmp_turbulence_u 	= parametric_turbulence_sigmaT * sqrt(3.) * cos(parmod_az) * sin(parmod_el);
								tmp_turbulence_v 	= parametric_turbulence_sigmaT * sqrt(3.) * sin(parmod_az) * sin(parmod_el);
								tmp_turbulence_w 	= parametric_turbulence_sigmaT * sqrt(3.) * cos(parmod_el);
							
								res_vol->subvolume_particle_u[i_res][i_psd][i_par][i_parmod] = res_vol->subvolume_air_u[i_res] + tmp_turbulence_u;
								res_vol->subvolume_particle_v[i_res][i_psd][i_par][i_parmod] = res_vol->subvolume_air_v[i_res] + tmp_turbulence_v;
								res_vol->subvolume_particle_w[i_res][i_psd][i_par][i_parmod] = res_vol->subvolume_air_w[i_res] - res_vol->subvolume_scat[i_psd][i_par]->particle_terminal_fall_speed + tmp_turbulence_w;
								
								//calculate turbulence vector plus (for calculation of derivatives via finite difference)
								if (todo->der_edr13) {
									tmp_turbulence_u 	= parametric_turbulence_sigmaT_plus * sqrt(3.) * cos(parmod_az) * sin(parmod_el);
									tmp_turbulence_v 	= parametric_turbulence_sigmaT_plus * sqrt(3.) * sin(parmod_az) * sin(parmod_el);
									tmp_turbulence_w 	= parametric_turbulence_sigmaT_plus * sqrt(3.) * cos(parmod_el);
								
									res_vol->subvolume_particle_u[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod] = res_vol->subvolume_air_u[i_res] + tmp_turbulence_u;
									res_vol->subvolume_particle_v[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod] = res_vol->subvolume_air_v[i_res] + tmp_turbulence_v;
									res_vol->subvolume_particle_w[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod] = res_vol->subvolume_air_w[i_res] - res_vol->subvolume_scat[i_psd][i_par]->particle_terminal_fall_speed + tmp_turbulence_w;
								}
							}
						}
					}
				}
			}
		}

		//calculate particle direction
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("calculate particle direction\n"); fflush(stdout);
		#endif 
		//assumption is minor axis is parallel to the motion of the droplet
		if (todo->calc_particle_direction) {
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
					for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
						for ( i_parmod = 0; i_parmod < res_vol->n_parmod_plus; i_parmod++ ) {
							//particle orientation of the minor axis
							tmp = sqrt(
									pow(res_vol->subvolume_particle_u[i_res][i_psd][i_par][i_parmod], 2.) +
									pow(res_vol->subvolume_particle_v[i_res][i_psd][i_par][i_parmod], 2.) +
									pow(res_vol->subvolume_particle_w[i_res][i_psd][i_par][i_parmod], 2.)
									);
							if (tmp == 0.) {
								res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][0] = 0.;
								res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][1] = 0.;
								res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][2] = -1.;
								res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][3] = 0.;
							} else {
								res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][0] = res_vol->subvolume_particle_u[i_res][i_psd][i_par][i_parmod] / tmp;
								res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][1] = res_vol->subvolume_particle_v[i_res][i_psd][i_par][i_parmod] / tmp;
								res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][2] = res_vol->subvolume_particle_w[i_res][i_psd][i_par][i_parmod] / tmp;
								res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][3] = 0.;
							}
						}
					}
				}
			}
		}

#ifdef _ZEPHYROS_RADARFILTER_DEBUG
	printf("initialize cross sections\n"); fflush(stdout);
#endif 
			
		//*****
		//initialize cross sections
		if (todo->calc_cross_sections == 1) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					//set water refractive index
					res_vol->subvolume_scat[i_psd][i_par]->particle_refractive_index = cfg->derived_quantities->radar_water_refractive_index;
					
					//calculate cross sections
					//particles_cross_sections_dewolf1990(res_vol->subvolume_scat[i_psd][i_par], radarmeasurement[i_m]->center_coor);			
					
					particles_cross_sections_mischenko2000(res_vol->subvolume_scat[i_psd][i_par], radarmeasurement[i_m]->center_coor, 0);			
				}
			}			
		}
			
#ifdef _ZEPHYROS_RADARFILTER_DEBUG
	printf("calculate particle number densities\n"); fflush(stdout);
#endif 
								
		//*****
		//calculate particle number densities
		if (todo->calc_number_density == 1) {
			//loop over particles
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
					for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
						//interpolate number_density_m3
						interpolation_bilint(
							myscattererfield->psd[i_psd]->lut_ln_number_density_m3[i_par], 
							res_vol->subvolume_coor[i_res]->enu_xyzt,
							&(res_vol->subvolume_number_density_m3[i_res][i_psd][i_par]),
							todo->calc_number_density_der,
							res_vol->subvolume_ln_number_density_m3_der[i_res][i_psd][i_par]);
						res_vol->subvolume_number_density_m3[i_res][i_psd][i_par] = exp(res_vol->subvolume_number_density_m3[i_res][i_psd][i_par]);
					}
				}
			}
		}
		
 		
		
		//*****		
		//initialize everything that needs the cross sections
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("initialize everything that needs the cross sections\n"); fflush(stdout);
		#endif
		if (todo->calc_eta_hh) eta_hh = 0.;
		if (todo->calc_eta_hv) eta_hv = 0.;
		if (todo->calc_eta_vh) eta_vh = 0.;
		if (todo->calc_eta_vv) eta_vv = 0.;
		if (todo->calc_eta_ShhSvvc) eta_ShhSvvc = 0.;
		if (todo->calc_eta_ShhShvc) eta_ShhShvc = 0.;
		if (todo->calc_eta_SvvSvhc) eta_SvvSvhc = 0.;

		if (rcfg->filter_KDP) n_Re_Shh_min_Svv = 0.;

		if (rcfg->filter_Doppler_velocity_hh_ms) radarmeasurement[i_m]->Doppler_velocity_hh_ms = 0.;
		if (rcfg->filter_Doppler_velocity_hv_ms) radarmeasurement[i_m]->Doppler_velocity_hv_ms = 0.;
		if (rcfg->filter_Doppler_velocity_vh_ms) radarmeasurement[i_m]->Doppler_velocity_vh_ms = 0.;
		if (rcfg->filter_Doppler_velocity_vv_ms) radarmeasurement[i_m]->Doppler_velocity_vv_ms = 0.;

		if (rcfg->filter_Doppler_spectralwidth_hh_ms) radarmeasurement[i_m]->Doppler_spectral_width_hh_ms = 0.;
		if (rcfg->filter_Doppler_spectralwidth_hv_ms) radarmeasurement[i_m]->Doppler_spectral_width_hv_ms = 0.;
		if (rcfg->filter_Doppler_spectralwidth_vh_ms) radarmeasurement[i_m]->Doppler_spectral_width_vh_ms = 0.;
		if (rcfg->filter_Doppler_spectralwidth_vv_ms) radarmeasurement[i_m]->Doppler_spectral_width_vv_ms = 0.;

		if (todo->der_edr13) {
			if (todo->calc_eta_hh) eta_hh_plus = 0.;
			if (todo->calc_eta_hv) eta_hv_plus = 0.;
			if (todo->calc_eta_vh) eta_vh_plus = 0.;
			if (todo->calc_eta_vv) eta_vv_plus = 0.;
			if (todo->calc_eta_ShhSvvc) eta_ShhSvvc_plus = 0.;
			if (todo->calc_eta_ShhShvc) eta_ShhShvc_plus = 0.;
			if (todo->calc_eta_SvvSvhc) eta_SvvSvhc_plus = 0.;
			
			if (rcfg->filter_Doppler_velocity_hh_ms) Doppler_velocity_hh_ms_plus = 0.;
			if (rcfg->filter_Doppler_spectralwidth_hh_ms) Doppler_spectral_width_hh_ms_plus = 0.;
		}

				
		//*****
		//calculate eta
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("calculate eta\n"); fflush(stdout);
		#endif 	
		if (todo->calc_cross_sections == 1) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					//initialize further
					if (todo->calc_eta_i_hh) radarmeasurement[i_m]->eta_i_hh[i_psd][i_par] = 0.;
					if (todo->calc_eta_i_hv) radarmeasurement[i_m]->eta_i_hv[i_psd][i_par] = 0.;
					if (todo->calc_eta_i_vh) radarmeasurement[i_m]->eta_i_vh[i_psd][i_par] = 0.;
					if (todo->calc_eta_i_vv) radarmeasurement[i_m]->eta_i_vv[i_psd][i_par] = 0.;
					
					for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
						for ( i_parmod = 0; i_parmod < res_vol->n_parmod_plus; i_parmod++ ) {
							//forward particle direction to scattering widget
							res_vol->subvolume_scat[i_psd][i_par]->particle_dir[0] = res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][0];
							res_vol->subvolume_scat[i_psd][i_par]->particle_dir[1] = res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][1];
							res_vol->subvolume_scat[i_psd][i_par]->particle_dir[2] = res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][2];
							//update cross sections with particle orientation
							//particles_cross_sections_dewolf1990_update_coor(res_vol->subvolume_scat[i_psd][i_par], res_vol->subvolume_coor[i_res]);
							//printf("debugging. using other cross sections...");
	
							particles_cross_sections_mischenko2000_update_coor(res_vol->subvolume_scat[i_psd][i_par], res_vol->subvolume_coor[i_res], 0);

							if (todo->calc_eta_hh) res_vol->subvolume_eta_hh[i_res][i_psd][i_par][i_parmod] = 
								res_vol->subvolume_number_density_m3[i_res][i_psd][i_par] * res_vol->subvolume_scat[i_psd][i_par]->particle_sigma_hh / res_vol->n_parmod;
							if (todo->calc_eta_hv) res_vol->subvolume_eta_hv[i_res][i_psd][i_par][i_parmod] = 
								res_vol->subvolume_number_density_m3[i_res][i_psd][i_par] * res_vol->subvolume_scat[i_psd][i_par]->particle_sigma_hv / res_vol->n_parmod;
							if (todo->calc_eta_vh) res_vol->subvolume_eta_vh[i_res][i_psd][i_par][i_parmod] = 
								res_vol->subvolume_number_density_m3[i_res][i_psd][i_par] * res_vol->subvolume_scat[i_psd][i_par]->particle_sigma_vh / res_vol->n_parmod;
							if (todo->calc_eta_vv) res_vol->subvolume_eta_vv[i_res][i_psd][i_par][i_parmod] = 
								res_vol->subvolume_number_density_m3[i_res][i_psd][i_par] * res_vol->subvolume_scat[i_psd][i_par]->particle_sigma_vv / res_vol->n_parmod;
							if (todo->calc_eta_ShhSvvc) res_vol->subvolume_eta_ShhSvvc[i_res][i_psd][i_par][i_parmod] = 
								res_vol->subvolume_number_density_m3[i_res][i_psd][i_par] * res_vol->subvolume_scat[i_psd][i_par]->particle_sigma_ShhSvvc / res_vol->n_parmod;
							if (todo->calc_eta_ShhShvc) res_vol->subvolume_eta_ShhShvc[i_res][i_psd][i_par][i_parmod] = 
								res_vol->subvolume_number_density_m3[i_res][i_psd][i_par] * res_vol->subvolume_scat[i_psd][i_par]->particle_sigma_ShhShvc / res_vol->n_parmod;
							if (todo->calc_eta_SvvSvhc) res_vol->subvolume_eta_SvvSvhc[i_res][i_psd][i_par][i_parmod] = 
								res_vol->subvolume_number_density_m3[i_res][i_psd][i_par] * res_vol->subvolume_scat[i_psd][i_par]->particle_sigma_SvvSvhc / res_vol->n_parmod;


							if (todo->calc_Doppler_mean_velocity == 1) {
								res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod] = 
									(res_vol->subvolume_particle_u[i_res][i_psd][i_par][i_parmod] * res_vol->subvolume_coor[i_res]->radar_enu_dir[0]) +
									(res_vol->subvolume_particle_v[i_res][i_psd][i_par][i_parmod] * res_vol->subvolume_coor[i_res]->radar_enu_dir[1]) +
									(res_vol->subvolume_particle_w[i_res][i_psd][i_par][i_parmod] * res_vol->subvolume_coor[i_res]->radar_enu_dir[2]);
							}
							
							if (i_parmod <  res_vol->n_parmod) {
								if (todo->calc_eta_hh) eta_hh += res_vol->subvolume_eta_hh[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_hv) eta_hv += res_vol->subvolume_eta_hv[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_vh) eta_vh += res_vol->subvolume_eta_vh[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_vv) eta_vv += res_vol->subvolume_eta_vv[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_ShhSvvc) eta_ShhSvvc += res_vol->subvolume_eta_ShhSvvc[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_ShhShvc) eta_ShhShvc += res_vol->subvolume_eta_ShhShvc[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_SvvSvhc) eta_SvvSvhc += res_vol->subvolume_eta_SvvSvhc[i_res][i_psd][i_par][i_parmod];
		
								if (todo->calc_eta_i_hh) radarmeasurement[i_m]->eta_i_hh[i_psd][i_par] += res_vol->subvolume_eta_hh[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_i_hv) radarmeasurement[i_m]->eta_i_hv[i_psd][i_par] += res_vol->subvolume_eta_hv[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_i_vh) radarmeasurement[i_m]->eta_i_vh[i_psd][i_par] += res_vol->subvolume_eta_vh[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_i_vv) radarmeasurement[i_m]->eta_i_vv[i_psd][i_par] += res_vol->subvolume_eta_vv[i_res][i_psd][i_par][i_parmod];
							} else {
								if (todo->calc_eta_hh) eta_hh_plus += res_vol->subvolume_eta_hh[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_hv) eta_hv_plus += res_vol->subvolume_eta_hv[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_vh) eta_vh_plus += res_vol->subvolume_eta_vh[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_vv) eta_vv_plus += res_vol->subvolume_eta_vv[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_ShhSvvc) eta_ShhSvvc_plus += res_vol->subvolume_eta_ShhSvvc[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_ShhShvc) eta_ShhShvc_plus += res_vol->subvolume_eta_ShhShvc[i_res][i_psd][i_par][i_parmod];
								if (todo->calc_eta_SvvSvhc) eta_SvvSvhc_plus += res_vol->subvolume_eta_SvvSvhc[i_res][i_psd][i_par][i_parmod];								
							}

							//for calculation of KdP, specific differential phase
							//calculate forward scattering amplitude
							if (rcfg->filter_KDP) {
								particles_cross_sections_mischenko2000_update_coor(res_vol->subvolume_scat[i_psd][i_par], res_vol->subvolume_coor[i_res], 1);
								tmp = res_vol->subvolume_number_density_m3[i_res][i_psd][i_par] * res_vol->subvolume_scat[i_psd][i_par]->particle_Re_Shh_min_Svv / res_vol->n_parmod;
								n_Re_Shh_min_Svv += tmp;
							}
						}
					}

					if (todo->calc_eta_i_hh) radarmeasurement[i_m]->eta_i_hh[i_psd][i_par] /= res_vol->n;
					if (todo->calc_eta_i_hv) radarmeasurement[i_m]->eta_i_hv[i_psd][i_par] /= res_vol->n;
					if (todo->calc_eta_i_vh) radarmeasurement[i_m]->eta_i_vh[i_psd][i_par] /= res_vol->n;
					if (todo->calc_eta_i_vv) radarmeasurement[i_m]->eta_i_vv[i_psd][i_par] /= res_vol->n;
					
				}
			}
			
			if (todo->calc_eta_hh) eta_hh /= res_vol->n;
			if (todo->calc_eta_hv) eta_hv /= res_vol->n;
			if (todo->calc_eta_vh) eta_vh /= res_vol->n;
			if (todo->calc_eta_vv) eta_vv /= res_vol->n;
			if (todo->calc_eta_ShhSvvc) eta_ShhSvvc /= res_vol->n;
			if (todo->calc_eta_ShhShvc) eta_ShhShvc /= res_vol->n;
			if (todo->calc_eta_SvvSvhc) eta_SvvSvhc /= res_vol->n;

			if (todo->der_edr13) {
				if (todo->calc_eta_hh) eta_hh_plus /= res_vol->n;
				if (todo->calc_eta_hv) eta_hv_plus /= res_vol->n;
				if (todo->calc_eta_vh) eta_vh_plus /= res_vol->n;
				if (todo->calc_eta_vv) eta_vv_plus /= res_vol->n;
				if (todo->calc_eta_ShhSvvc) eta_ShhSvvc_plus /= res_vol->n;
				if (todo->calc_eta_ShhShvc) eta_ShhShvc_plus /= res_vol->n;
				if (todo->calc_eta_SvvSvhc) eta_SvvSvhc_plus /= res_vol->n;
			}

			if (rcfg->filter_KDP) n_Re_Shh_min_Svv /= res_vol->n;
		}
		
		//*****
		//calculate particle Doppler mean velocity
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("calculate particle Doppler mean velocity\n"); fflush(stdout);
		#endif 	
		if (todo->calc_Doppler_mean_velocity == 1) {
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
					for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
						for ( i_parmod = 0; i_parmod < res_vol->n_parmod; i_parmod++ ) {
							if (rcfg->filter_Doppler_velocity_hh_ms) radarmeasurement[i_m]->Doppler_velocity_hh_ms += 
								(res_vol->subvolume_eta_hh[i_res][i_psd][i_par][i_parmod] *
								res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod]);
							if (rcfg->filter_Doppler_velocity_hv_ms) radarmeasurement[i_m]->Doppler_velocity_hv_ms += 
								(res_vol->subvolume_eta_hv[i_res][i_psd][i_par][i_parmod] *
								res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod]);
							if (rcfg->filter_Doppler_velocity_vh_ms) radarmeasurement[i_m]->Doppler_velocity_vh_ms += 
								(res_vol->subvolume_eta_vh[i_res][i_psd][i_par][i_parmod] *
								res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod]);
							if (rcfg->filter_Doppler_velocity_vv_ms) radarmeasurement[i_m]->Doppler_velocity_vv_ms += 
								(res_vol->subvolume_eta_vv[i_res][i_psd][i_par][i_parmod] *
								res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod]);

							if (todo->der_edr13) {
								if (rcfg->filter_Doppler_velocity_hh_ms)
									Doppler_velocity_hh_ms_plus += 
									(res_vol->subvolume_eta_hh[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod] *
									res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod]);									
							}
						}	
					}	
				}
			}
												
			if (rcfg->filter_Doppler_velocity_hh_ms) radarmeasurement[i_m]->Doppler_velocity_hh_ms /= (res_vol->n * eta_hh);
			if (rcfg->filter_Doppler_velocity_hv_ms) radarmeasurement[i_m]->Doppler_velocity_hv_ms /= (res_vol->n * eta_hv);
			if (rcfg->filter_Doppler_velocity_vh_ms) radarmeasurement[i_m]->Doppler_velocity_vh_ms /= (res_vol->n * eta_vh);
			if (rcfg->filter_Doppler_velocity_vv_ms) radarmeasurement[i_m]->Doppler_velocity_vv_ms /= (res_vol->n * eta_vv);		

			if (todo->der_edr13) {
				if (rcfg->filter_Doppler_velocity_hh_ms) Doppler_velocity_hh_ms_plus /= (res_vol->n * eta_hh_plus);
				radarmeasurement[i_m]->der_edr13_Doppler_velocity_hh_ms =
					(Doppler_velocity_hh_ms_plus - radarmeasurement[i_m]->Doppler_velocity_hh_ms) / 0.001;
				
			}
		}

		//*****
		//Calculate Doppler spectral widths
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("calculate particle Doppler spectral widths\n"); fflush(stdout);
		#endif 
		if (todo->calc_Doppler_spectrum_widths == 1) {
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
					for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
						for ( i_parmod = 0; i_parmod < res_vol->n_parmod; i_parmod++ ) {
							if (rcfg->filter_Doppler_spectralwidth_hh_ms) radarmeasurement[i_m]->Doppler_spectral_width_hh_ms += 
								(res_vol->subvolume_eta_hh[i_res][i_psd][i_par][i_parmod] *
								pow(res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod] - radarmeasurement[i_m]->Doppler_velocity_hh_ms, 2.));
							if (rcfg->filter_Doppler_spectralwidth_hv_ms) radarmeasurement[i_m]->Doppler_spectral_width_hv_ms += 
								(res_vol->subvolume_eta_hv[i_res][i_psd][i_par][i_parmod] *
								pow(res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod] - radarmeasurement[i_m]->Doppler_velocity_hv_ms, 2.));
							if (rcfg->filter_Doppler_spectralwidth_vh_ms) radarmeasurement[i_m]->Doppler_spectral_width_vh_ms += 
								(res_vol->subvolume_eta_vh[i_res][i_psd][i_par][i_parmod] *
								pow(res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod] - radarmeasurement[i_m]->Doppler_velocity_vh_ms, 2.));
							if (rcfg->filter_Doppler_spectralwidth_vv_ms) radarmeasurement[i_m]->Doppler_spectral_width_vv_ms += 
								(res_vol->subvolume_eta_vv[i_res][i_psd][i_par][i_parmod] *
								pow(res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod] - radarmeasurement[i_m]->Doppler_velocity_vv_ms, 2.));

							if (todo->der_edr13) {
								if (rcfg->filter_Doppler_spectralwidth_hh_ms)
									Doppler_spectral_width_hh_ms_plus += 
									(res_vol->subvolume_eta_hh[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod] *
									pow(res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod] - Doppler_velocity_hh_ms_plus, 2.));
							}
						}
					}
				}
			}
			
			if (rcfg->filter_Doppler_spectralwidth_hh_ms) radarmeasurement[i_m]->Doppler_spectral_width_hh_ms /= (res_vol->n * eta_hh);
			if (rcfg->filter_Doppler_spectralwidth_hv_ms) radarmeasurement[i_m]->Doppler_spectral_width_hv_ms /= (res_vol->n * eta_hv);
			if (rcfg->filter_Doppler_spectralwidth_vh_ms) radarmeasurement[i_m]->Doppler_spectral_width_vh_ms /= (res_vol->n * eta_vh);		
			if (rcfg->filter_Doppler_spectralwidth_vv_ms) radarmeasurement[i_m]->Doppler_spectral_width_vv_ms /= (res_vol->n * eta_vv);		
			
			if (rcfg->filter_Doppler_spectralwidth_hh_ms) radarmeasurement[i_m]->Doppler_spectral_width_hh_ms = sqrt(radarmeasurement[i_m]->Doppler_spectral_width_hh_ms);
			if (rcfg->filter_Doppler_spectralwidth_hv_ms) radarmeasurement[i_m]->Doppler_spectral_width_hv_ms = sqrt(radarmeasurement[i_m]->Doppler_spectral_width_hv_ms);
			if (rcfg->filter_Doppler_spectralwidth_vh_ms) radarmeasurement[i_m]->Doppler_spectral_width_vh_ms = sqrt(radarmeasurement[i_m]->Doppler_spectral_width_vh_ms);
			if (rcfg->filter_Doppler_spectralwidth_hh_ms) radarmeasurement[i_m]->Doppler_spectral_width_vv_ms = sqrt(radarmeasurement[i_m]->Doppler_spectral_width_vv_ms);

			if (todo->der_edr13) {
				if (rcfg->filter_Doppler_spectralwidth_hh_ms) Doppler_spectral_width_hh_ms_plus /= (res_vol->n * eta_hh_plus);
				if (rcfg->filter_Doppler_spectralwidth_hh_ms) Doppler_spectral_width_hh_ms_plus = sqrt(Doppler_spectral_width_hh_ms_plus);
				radarmeasurement[i_m]->der_edr13_Doppler_spectral_width_hh_ms =
					(Doppler_spectral_width_hh_ms_plus - radarmeasurement[i_m]->Doppler_spectral_width_hh_ms) / 0.001;
			}
		}
		
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("calculate Doppler spectrum\n"); fflush(stdout);
		#endif 
		if (todo->calc_Doppler_spectrum) {
			if (i_mode == 0) {
				//in case of simulation, set lbound, center and ubound
				util_safe_free(&(radarmeasurement[i_m]->spectrum_lbound));
				util_safe_free(&(radarmeasurement[i_m]->spectrum_ubound));
				util_safe_free(&(radarmeasurement[i_m]->spectrum_center));
				
				if (todo->calc_Doppler_spectrum) radarmeasurement[i_m]->spectrum_lbound = malloc(radarmeasurement[i_m]->n_spectrum * sizeof(double));
				if (todo->calc_Doppler_spectrum) radarmeasurement[i_m]->spectrum_ubound = malloc(radarmeasurement[i_m]->n_spectrum * sizeof(double));
				if (todo->calc_Doppler_spectrum) radarmeasurement[i_m]->spectrum_center = malloc(radarmeasurement[i_m]->n_spectrum * sizeof(double));

				//Calculate spectral shape
				//mincdfP = 1.e-50;
				//maxcdfP = 1. - 1.e-50;
				if (radarmeasurement[i_m]->n_spectrum <= 8) {
					//98%
					mincdfP = .01;
					maxcdfP = .99;
				} else {
					//99.9%
					mincdfP = .001;
					maxcdfP = .999;				
				}
				difcdfP = (maxcdfP - mincdfP) / radarmeasurement[i_m]->n_spectrum;
				
				//walk through intervals
				for ( i_int = 0; i_int < radarmeasurement[i_m]->n_spectrum; i_int++ ) {
					//unscaled bounds of this interval
					lbound = ltqnorm(mincdfP + (i_int * difcdfP));
					ubound = ltqnorm(mincdfP + ((i_int + 1.)* difcdfP));
					center = ltqnorm(mincdfP + ((i_int + .5)* difcdfP));
					
					//scaled bounds based on hh
					radarmeasurement[i_m]->spectrum_lbound[i_int] = (radarmeasurement[i_m]->Doppler_spectral_width_hh_ms * lbound) + radarmeasurement[i_m]->Doppler_velocity_hh_ms;
					radarmeasurement[i_m]->spectrum_ubound[i_int] = (radarmeasurement[i_m]->Doppler_spectral_width_hh_ms * ubound) + radarmeasurement[i_m]->Doppler_velocity_hh_ms;
					radarmeasurement[i_m]->spectrum_center[i_int] = (radarmeasurement[i_m]->Doppler_spectral_width_hh_ms * center) + radarmeasurement[i_m]->Doppler_velocity_hh_ms;
				}
			}
							
			//walk through intervals
			for ( i_int = 0; i_int < radarmeasurement[i_m]->n_spectrum; i_int++ ) {
				//TBD: think about Doppler folding here ...

				//reset
				if (rcfg->filter_Doppler_spectrum_dBZ_hh) radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh[i_int] = 0.;
				if (rcfg->filter_Doppler_spectrum_dBZ_hv) radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv[i_int] = 0.;
				if (rcfg->filter_Doppler_spectrum_dBZ_vh) radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh[i_int] = 0.;
				if (rcfg->filter_Doppler_spectrum_dBZ_vv) radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int] = 0.;
				if (todo->calc_Doppler_spectrum_ShhSvvc) Doppler_spectrum_ShhSvvc[i_int] = 0.;
				if (todo->calc_Doppler_spectrum_ShhShvc) Doppler_spectrum_ShhShvc[i_int] = 0.;
				if (todo->calc_Doppler_spectrum_SvvSvhc) Doppler_spectrum_SvvSvhc[i_int] = 0.;
				
				if (todo->der_edr13) {
					if (rcfg->filter_Doppler_spectrum_dBZ_hh) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hh[i_int] = 0.;
					if (rcfg->filter_Doppler_spectrum_dBZ_hv) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hv[i_int] = 0.;
					if (rcfg->filter_Doppler_spectrum_dBZ_vh) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vh[i_int] = 0.;
					if (rcfg->filter_Doppler_spectrum_dBZ_vv) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vv[i_int] = 0.;					
				}
					
				for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
					for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
						if (todo->calc_spectrum_eta_i_hh) radarmeasurement[i_m]->spectrum_eta_i_hh[i_psd][i_par][i_int] = 0.;
						if (todo->calc_spectrum_eta_i_hv) radarmeasurement[i_m]->spectrum_eta_i_hv[i_psd][i_par][i_int] = 0.;
						if (todo->calc_spectrum_eta_i_vh) radarmeasurement[i_m]->spectrum_eta_i_vh[i_psd][i_par][i_int] = 0.;
						if (todo->calc_spectrum_eta_i_vv) radarmeasurement[i_m]->spectrum_eta_i_vv[i_psd][i_par][i_int] = 0.;

						for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
							for ( i_parmod = 0; i_parmod < res_vol->n_parmod; i_parmod++ ) {
								if ((radarmeasurement[i_m]->spectrum_lbound[i_int] <= res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod]) & (res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par][i_parmod] < radarmeasurement[i_m]->spectrum_ubound[i_int])) {
									if (rcfg->filter_Doppler_spectrum_dBZ_hh) radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh[i_int] += res_vol->subvolume_eta_hh[i_res][i_psd][i_par][i_parmod];
									if (rcfg->filter_Doppler_spectrum_dBZ_hv) radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv[i_int] += res_vol->subvolume_eta_hv[i_res][i_psd][i_par][i_parmod];
									if (rcfg->filter_Doppler_spectrum_dBZ_vh) radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh[i_int] += res_vol->subvolume_eta_vh[i_res][i_psd][i_par][i_parmod];
									if (rcfg->filter_Doppler_spectrum_dBZ_vv) radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int] += res_vol->subvolume_eta_vv[i_res][i_psd][i_par][i_parmod];

									if (todo->der_edr13) {
										if (rcfg->filter_Doppler_spectrum_dBZ_hh) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hh[i_int] += res_vol->subvolume_eta_hh[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod];
										if (rcfg->filter_Doppler_spectrum_dBZ_hv) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hv[i_int] += res_vol->subvolume_eta_hv[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod];
										if (rcfg->filter_Doppler_spectrum_dBZ_vh) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vh[i_int] += res_vol->subvolume_eta_vh[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod];
										if (rcfg->filter_Doppler_spectrum_dBZ_vv) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vv[i_int] += res_vol->subvolume_eta_vv[i_res][i_psd][i_par][res_vol->n_parmod + i_parmod];
									}

									if (todo->calc_spectrum_eta_i_hh) radarmeasurement[i_m]->spectrum_eta_i_hh[i_psd][i_par][i_int] += res_vol->subvolume_eta_hh[i_res][i_psd][i_par][i_parmod];
									if (todo->calc_spectrum_eta_i_hv) radarmeasurement[i_m]->spectrum_eta_i_hv[i_psd][i_par][i_int] += res_vol->subvolume_eta_hv[i_res][i_psd][i_par][i_parmod];
									if (todo->calc_spectrum_eta_i_vh) radarmeasurement[i_m]->spectrum_eta_i_vh[i_psd][i_par][i_int] += res_vol->subvolume_eta_vh[i_res][i_psd][i_par][i_parmod];
									if (todo->calc_spectrum_eta_i_vv) radarmeasurement[i_m]->spectrum_eta_i_vv[i_psd][i_par][i_int] += res_vol->subvolume_eta_vv[i_res][i_psd][i_par][i_parmod];

									if (todo->calc_Doppler_spectrum_ShhSvvc) Doppler_spectrum_ShhSvvc[i_int] += res_vol->subvolume_eta_ShhSvvc[i_res][i_psd][i_par][i_parmod];
									if (todo->calc_Doppler_spectrum_ShhShvc) Doppler_spectrum_ShhShvc[i_int] += res_vol->subvolume_eta_ShhShvc[i_res][i_psd][i_par][i_parmod];
									if (todo->calc_Doppler_spectrum_SvvSvhc) Doppler_spectrum_SvvSvhc[i_int] += res_vol->subvolume_eta_SvvSvhc[i_res][i_psd][i_par][i_parmod];
								}
							}
						}
						if (todo->calc_spectrum_eta_i_hh) radarmeasurement[i_m]->spectrum_eta_i_hh[i_psd][i_par][i_int] *= (1.  / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));
						if (todo->calc_spectrum_eta_i_hv) radarmeasurement[i_m]->spectrum_eta_i_hv[i_psd][i_par][i_int] *= (1.  / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));
						if (todo->calc_spectrum_eta_i_vh) radarmeasurement[i_m]->spectrum_eta_i_vh[i_psd][i_par][i_int] *= (1.  / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));
						if (todo->calc_spectrum_eta_i_vv) radarmeasurement[i_m]->spectrum_eta_i_vv[i_psd][i_par][i_int] *= (1.  / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));
					}
				}
				
				//1.e18: conversion of m^3 to mm^6 m^-3
				if (rcfg->filter_Doppler_spectrum_dBZ_hh) radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh[i_int] = func_dB(eta2Z * radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh[i_int] / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));
				if (rcfg->filter_Doppler_spectrum_dBZ_hv) radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv[i_int] = func_dB(eta2Z * radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv[i_int] / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));
				if (rcfg->filter_Doppler_spectrum_dBZ_vh) radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh[i_int] = func_dB(eta2Z * radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh[i_int]/ ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));				
				if (rcfg->filter_Doppler_spectrum_dBZ_vv) radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int] = func_dB(eta2Z * radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int]/ ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));				

				if (todo->der_edr13) {
					//plus term is calculated
					if (rcfg->filter_Doppler_spectrum_dBZ_hh) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hh[i_int] = func_dB(eta2Z * radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh[i_int] / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));
					if (rcfg->filter_Doppler_spectrum_dBZ_hv) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hv[i_int] = func_dB(eta2Z * radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv[i_int] / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));
					if (rcfg->filter_Doppler_spectrum_dBZ_vh) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vh[i_int] = func_dB(eta2Z * radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh[i_int]/ ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));				
					if (rcfg->filter_Doppler_spectrum_dBZ_vv) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vv[i_int] = func_dB(eta2Z * radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int]/ ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n));				
					//derivative to edr is calculated
					if (rcfg->filter_Doppler_spectrum_dBZ_hh) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hh[i_int] = (radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hh[i_int] - radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh[i_int]) / 0.001;
					if (rcfg->filter_Doppler_spectrum_dBZ_hv) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hv[i_int] = (radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_hv[i_int] - radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv[i_int]) / 0.001;
					if (rcfg->filter_Doppler_spectrum_dBZ_vh) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vh[i_int] = (radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vh[i_int] - radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh[i_int]) / 0.001;
					if (rcfg->filter_Doppler_spectrum_dBZ_vv) radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vv[i_int] = (radarmeasurement[i_m]->der_edr13_Doppler_spectrum_dBZ_vv[i_int] - radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int]) / 0.001;
					
				}
				
				if (todo->calc_Doppler_spectrum_ShhSvvc) Doppler_spectrum_ShhSvvc[i_int] *= eta2Z / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n);				
				if (todo->calc_Doppler_spectrum_ShhShvc) Doppler_spectrum_ShhShvc[i_int] *= eta2Z / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n);				
				if (todo->calc_Doppler_spectrum_SvvSvhc) Doppler_spectrum_SvvSvhc[i_int] *= eta2Z / ((radarmeasurement[i_m]->spectrum_ubound[i_int] - radarmeasurement[i_m]->spectrum_lbound[i_int]) * res_vol->n);				
			}
		}


		//post calculations
		//*****
		#ifdef _ZEPHYROS_RADARFILTER_DEBUG
			printf("post calculations\n"); fflush(stdout);
		#endif
		
		//1.e18: conversion of m^3 to mm^6 m^-3
		if (rcfg->filter_dBZ_hh) radarmeasurement[i_m]->dBZ_hh = func_dB(eta2Z * eta_hh);
		if (rcfg->filter_dBZ_hv) radarmeasurement[i_m]->dBZ_hv = func_dB(eta2Z * eta_hv);
		if (rcfg->filter_dBZ_vh) radarmeasurement[i_m]->dBZ_vh = func_dB(eta2Z * eta_vh);
		if (rcfg->filter_dBZ_vv) radarmeasurement[i_m]->dBZ_vv = func_dB(eta2Z * eta_vv);

		//*****
		if (rcfg->filter_dBZdr) radarmeasurement[i_m]->dBZdr = func_dB(eta_hh / eta_vv);
		if (rcfg->filter_dBLdr) radarmeasurement[i_m]->dBLdr = func_dB(eta_hv / eta_vv);

		//*****
		if (rcfg->filter_specific_dBZdr) {
			for ( i_int = 0; i_int < radarmeasurement[i_m]->n_spectrum; i_int++ ) {
				radarmeasurement[i_m]->specific_dBZdr[i_int] = radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh[i_int] - radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int];
			}
		}
		if (rcfg->filter_specific_dBLdr) {
			for ( i_int = 0; i_int < radarmeasurement[i_m]->n_spectrum; i_int++ ) {
				radarmeasurement[i_m]->specific_dBLdr[i_int] = radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv[i_int] - radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int];
			}
		}

		//*****
		if (rcfg->filter_rho_co) 	radarmeasurement[i_m]->rho_co 	= eta_ShhSvvc / sqrt(eta_hh * eta_vv);
		if (rcfg->filter_rho_cxh) 	radarmeasurement[i_m]->rho_cxh 	= eta_ShhShvc / sqrt(eta_hh * eta_hv);
		if (rcfg->filter_rho_cxv) 	radarmeasurement[i_m]->rho_cxv 	= eta_SvvSvhc / sqrt(eta_vv * eta_vh);

		//*****
		if (rcfg->filter_specific_rho_co) {
			for ( i_int = 0; i_int < radarmeasurement[i_m]->n_spectrum; i_int++ ) {
				radarmeasurement[i_m]->specific_rho_co[i_int] = Doppler_spectrum_ShhSvvc[i_int] / sqrt(func_dB_inv(radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh[i_int]) * func_dB_inv(radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int]));
			}
		}
		if (rcfg->filter_specific_rho_cxh) {
			for ( i_int = 0; i_int < radarmeasurement[i_m]->n_spectrum; i_int++ ) {
				radarmeasurement[i_m]->specific_rho_cxh[i_int] = Doppler_spectrum_ShhShvc[i_int] / sqrt(func_dB_inv(radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh[i_int]) * func_dB_inv(radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv[i_int]));
			}
		}
		if (rcfg->filter_specific_rho_cxv) {
			for ( i_int = 0; i_int < radarmeasurement[i_m]->n_spectrum; i_int++ ) {
				radarmeasurement[i_m]->specific_rho_cxv[i_int] = Doppler_spectrum_SvvSvhc[i_int] / sqrt(func_dB_inv(radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv[i_int]) * func_dB_inv(radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh[i_int]));
			}
		}

		if (todo->der_edr13) {
			if (rcfg->filter_dBZ_hh) radarmeasurement[i_m]->der_edr13_dBZ_hh = 
				func_dB(eta_hh_plus/eta_hh) / 0.001;
			if (rcfg->filter_dBZdr) radarmeasurement[i_m]->der_edr13_dBZdr =
				(func_dB(eta_hh_plus / eta_vv_plus) -
				func_dB(eta_hh / eta_vv)) / 0.001;
			if (rcfg->filter_dBLdr) radarmeasurement[i_m]->der_edr13_dBLdr =
				(func_dB(eta_hv_plus / eta_vv_plus) -
				func_dB(eta_hv / eta_vv)) / 0.001;
		}


		//KDP
		if (rcfg->filter_KDP) radarmeasurement[i_m]->KDP = 1.e3 * cfg->derived_quantities->central_wavelength_m * n_Re_Shh_min_Svv;
				
		if (rcfg->filter_errors) {
			//TBD
			//simply put some numbers in
			if (rcfg->filter_dBZ_hh) radarmeasurement[i_m]->dBZ_hh_err = 1.;
			if (rcfg->filter_dBZ_hv) radarmeasurement[i_m]->dBZ_hv_err = 1.;
			if (rcfg->filter_dBZ_vh) radarmeasurement[i_m]->dBZ_vh_err = 1.;
			if (rcfg->filter_dBZ_vv) radarmeasurement[i_m]->dBZ_vv_err = 1.;
			if (rcfg->filter_dBZdr) radarmeasurement[i_m]->dBZdr_err = 1.;
			if (rcfg->filter_dBLdr) radarmeasurement[i_m]->dBLdr_err = 1.;
			if (rcfg->filter_rho_co) radarmeasurement[i_m]->rho_co_err = 1.;
			if (rcfg->filter_rho_cxh) radarmeasurement[i_m]->rho_cxh_err = 1.;
			if (rcfg->filter_rho_cxv) radarmeasurement[i_m]->rho_cxv_err = 1.;
			if (rcfg->filter_Doppler_velocity_hh_ms) radarmeasurement[i_m]->Doppler_velocity_hh_ms_err = 1.;
			if (rcfg->filter_Doppler_velocity_hv_ms) radarmeasurement[i_m]->Doppler_velocity_hv_ms_err = 1.;
			if (rcfg->filter_Doppler_velocity_vh_ms) radarmeasurement[i_m]->Doppler_velocity_vh_ms_err = 1.;
			if (rcfg->filter_Doppler_velocity_vv_ms) radarmeasurement[i_m]->Doppler_velocity_vv_ms_err = 1.;
			if (rcfg->filter_Doppler_spectralwidth_hh_ms) radarmeasurement[i_m]->Doppler_spectral_width_hh_ms_err = 1.;
			if (rcfg->filter_Doppler_spectralwidth_hv_ms) radarmeasurement[i_m]->Doppler_spectral_width_hv_ms_err = 1.;
			if (rcfg->filter_Doppler_spectralwidth_vh_ms) radarmeasurement[i_m]->Doppler_spectral_width_vh_ms_err = 1.;
			if (rcfg->filter_Doppler_spectralwidth_vv_ms) radarmeasurement[i_m]->Doppler_spectral_width_vv_ms_err = 1.;
			
			for ( i_int = 0; i_int < radarmeasurement[i_m]->n_spectrum; i_int++ ) {
				if (rcfg->filter_Doppler_spectrum_dBZ_hh) radarmeasurement[i_m]->Doppler_spectrum_dBZ_hh_err[i_int] = 1.;
				if (rcfg->filter_Doppler_spectrum_dBZ_hv) radarmeasurement[i_m]->Doppler_spectrum_dBZ_hv_err[i_int] = 1.;
				if (rcfg->filter_Doppler_spectrum_dBZ_vh) radarmeasurement[i_m]->Doppler_spectrum_dBZ_vh_err[i_int] = 1.;
				if (rcfg->filter_Doppler_spectrum_dBZ_vv) radarmeasurement[i_m]->Doppler_spectrum_dBZ_vv_err[i_int] = 1.;
				if (rcfg->filter_specific_dBZdr) radarmeasurement[i_m]->specific_dBZdr_err[i_int] = 1.;
				if (rcfg->filter_specific_dBLdr) radarmeasurement[i_m]->specific_dBLdr_err[i_int] = 1.;
				if (rcfg->filter_specific_rho_co) radarmeasurement[i_m]->specific_rho_co_err[i_int] = 1.;
				if (rcfg->filter_specific_rho_cxh) radarmeasurement[i_m]->specific_rho_cxh_err[i_int] = 1.;
				if (rcfg->filter_specific_rho_cxv) radarmeasurement[i_m]->specific_rho_cxv_err[i_int] = 1.;
			}
			
			if (rcfg->filter_KDP) radarmeasurement[i_m]->KDP_err = 1.;

		}
				
		if ((i_mode == 0) & (cfg->general->additional_output->print_detailed_analysis)) {
			//analysis of canting angles
			//initalize
			radarmeasurement[i_m]->analysis->discrete_D_equiv_mm 								= malloc(res_vol->n_psd * sizeof(double*));
			radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_mean 		= malloc(res_vol->n_psd * sizeof(double*));
			radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_variance 	= malloc(res_vol->n_psd * sizeof(double*));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				radarmeasurement[i_m]->analysis->discrete_D_equiv_mm[i_psd] 							= malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double));
				radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_mean[i_psd] 		= malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double));
				radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_variance[i_psd] 	= malloc(radarmeasurement[i_m]->n_diameters[i_psd] * sizeof(double));				
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					radarmeasurement[i_m]->analysis->discrete_D_equiv_mm[i_psd][i_par] = res_vol->subvolume_scat[i_psd][i_par]->particle_D_eqvol_mm;
					radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_mean[i_psd][i_par]	 	= 0.;
					radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_variance[i_psd][i_par] 	= 0.;
				}
			}

			//mean
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
						for ( i_parmod = 0; i_parmod < res_vol->n_parmod; i_parmod++ ) {
							//calculate canting angle w.r.t. vertical
							//cos theta = b . a , b = vertical, a = minor axis
							radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_mean[i_psd][i_par] += 
								acos(fabs(res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][2]));
						}
					}
					radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_mean[i_psd][i_par] /=
						(1. * res_vol->n * res_vol->n_parmod);
				}
			}

			//variance
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
						for ( i_parmod = 0; i_parmod < res_vol->n_parmod; i_parmod++ ) {
							radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_variance[i_psd][i_par] += 
								pow(acos(fabs(res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod][2])) - 
								radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_mean[i_psd][i_par], 2.);
						}
					}
					radarmeasurement[i_m]->analysis->unweighted_canting_angle_wrt_vertical_variance[i_psd][i_par] /=
						(1. * res_vol->n * res_vol->n_parmod);
				}
			}
		}
		
		//free variables
		util_safe_free(&(Doppler_spectrum_ShhSvvc));
		util_safe_free(&(Doppler_spectrum_ShhShvc));
		util_safe_free(&(Doppler_spectrum_SvvSvhc));

	}

	//free resolution volume
	#ifdef _ZEPHYROS_RADARFILTER_DEBUG
		printf("free resolution volume\n"); fflush(stdout);
	#endif 
	radarfilter_free_resolution_volume(cfg, i_mode, &res_vol, todo);
	
	#ifdef _ZEPHYROS_RADARFILTER_DEBUG
		printf("end of radarfilter_exec\n"); fflush(stdout);
	#endif
}



/*





	if (debug) start = clock();	

	if (debug) diff = clock() - start;
	if (debug) printf("calculation of coordinates: %.6e seconds\n", 1000. * diff / CLOCKS_PER_SEC);


	//calculate cross sections ...
	// .... TBD ....




		
		
		/*

	if (todo->calc_integrated_radial_vel_ms_derivatives == 1) {
		res_vol->integrated_dradial_vel_ms_dx	= 0.;
		res_vol->integrated_dradial_vel_ms_dy	= 0.;
		res_vol->integrated_dradial_vel_ms_dz	= 0.;
		res_vol->integrated_dradial_vel_ms_dt	= 0.;

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->integrated_dradial_vel_ms_dx += 
				res_vol->point_Z[i_res] * res_vol->point_dradial_vel_ms_dx[i_res];
			res_vol->integrated_dradial_vel_ms_dy += 
				res_vol->point_Z[i_res] * res_vol->point_dradial_vel_ms_dy[i_res];
			res_vol->integrated_dradial_vel_ms_dz += 
				res_vol->point_Z[i_res] * res_vol->point_dradial_vel_ms_dz[i_res];
			res_vol->integrated_dradial_vel_ms_dt += 
				res_vol->point_Z[i_res] * res_vol->point_dradial_vel_ms_dt[i_res];
		}
		
		res_vol->integrated_dradial_vel_ms_dx	/= res_vol->integrated_Z;
		res_vol->integrated_dradial_vel_ms_dy	/= res_vol->integrated_Z;
		res_vol->integrated_dradial_vel_ms_dz	/= res_vol->integrated_Z;
		res_vol->integrated_dradial_vel_ms_dt	/= res_vol->integrated_Z;
	}

		if (todo->apply_advection == 1) {
			//advection
			//xnew = f_adv_x(xold, ...)
			//ynew = f_adv_y(yold, ...)
			//znew = f_adv_z(zold, ...)

			//advection derivatives
			res_vol->point_df_adv_x_dx[i_res] = 1.;
			res_vol->point_df_adv_y_dy[i_res] = 1.;
			res_vol->point_df_adv_z_dz[i_res] = 1.;

			dt = outx[3];
			//res_vol->advection_integration_steps = 5;
			for ( i = 0; i < res_vol->advection_integration_steps; i++ ) {
				calcderivatives = 1;
				
				//obtain advection: uadv, vadv, wadv
				windfield_fuvw(windfield, outx, 0, &u_adv, calcderivatives, derivatives);
				du_adv_dx = derivatives[0];
				
				windfield_fuvw(windfield, outx, 1, &v_adv, calcderivatives, derivatives);
				dv_adv_dy = derivatives[1];
						
				windfield_fuvw(windfield, outx, 2, &w_adv, calcderivatives, derivatives);
				dw_adv_dz = derivatives[2];
				
				if (todo->apply_geostrophic_correction == 1) {
					ug = sqrt(pow(u_adv,2.) + pow(v_adv, 2.));
					tg = M_PI + atan2(u_adv, v_adv);
					u_adv = ug * sin((res_vol->fcoriolis * (-1. * dt / res_vol->advection_integration_steps)) - tg);
					v_adv = ug * cos((res_vol->fcoriolis * (-1. * dt / res_vol->advection_integration_steps)) - tg + M_PI);
				}
				
				outx[0] = outx[0] + (u_adv * (-dt / res_vol->advection_integration_steps));
				outx[1] = outx[1] + (v_adv * (-dt / res_vol->advection_integration_steps));
				outx[2] = outx[2] + (w_adv * (-dt / res_vol->advection_integration_steps));
				outx[3] = outx[3] + (-dt / res_vol->advection_integration_steps);
				
				//derivatives
				res_vol->point_df_adv_x_dx[i_res] += (du_adv_dx * (-dt / res_vol->advection_integration_steps));
				res_vol->point_df_adv_y_dy[i_res] += (dv_adv_dy * (-dt / res_vol->advection_integration_steps));
				res_vol->point_df_adv_z_dz[i_res] += (dw_adv_dz * (-dt / res_vol->advection_integration_steps));


			}
			
			//restore original time
			outx[3] = dt;



			if (todo->apply_advection == 1) {
				res_vol->point_ddBZ_dx[i_res] *= res_vol->point_df_adv_x_dx[i_res];	
				res_vol->point_ddBZ_dy[i_res] *= res_vol->point_df_adv_y_dy[i_res];
				res_vol->point_ddBZ_dz[i_res] *= res_vol->point_df_adv_z_dz[i_res];
			}


			if (todo->apply_advection == 1) {
				res_vol->point_du_dx[i_res] *= res_vol->point_df_adv_x_dx[i_res];	
				res_vol->point_du_dy[i_res] *= res_vol->point_df_adv_y_dy[i_res];
				res_vol->point_du_dz[i_res] *= res_vol->point_df_adv_z_dz[i_res];
			}

		}

		if (todo->apply_geostrophic_correction == 1) {
			if ((todo->calc_point_u == 1) & (todo->calc_point_v == 1)) {
				point_u 		= res_vol->point_u[i_res];
				point_v			= res_vol->point_v[i_res];				
				ug = sqrt(pow(point_u,2.) + pow(point_v, 2.));
				tg = M_PI + atan2(point_u, point_v);

				//if derivatives
				if ((todo->calc_point_u_derivatives == 1) & (todo->calc_point_v_derivatives == 1)) {
					point_du_dx		= res_vol->point_du_dx[i_res];
					point_du_dy		= res_vol->point_du_dy[i_res];
					point_dv_dx		= res_vol->point_dv_dx[i_res];
					point_dv_dy		= res_vol->point_dv_dy[i_res];
					
					dug_dx = (1. / ug) * (point_u * point_du_dx + point_v * point_dv_dx);
					dug_dy = (1. / ug) * (point_u * point_du_dx + point_v * point_dv_dx);
					dtg_dx = pow(ug, -2.) * (point_v * point_du_dx - point_u * point_dv_dx);
					dtg_dy = pow(ug, -2.) * (point_v * point_du_dy - point_u * point_dv_dy);

					//update
					res_vol->point_du_dx[i_res] = dug_dx * sin((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg)
													+ ug * cos((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg) * dtg_dx;
					res_vol->point_du_dy[i_res] = dug_dy * sin((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg)
													+ ug * cos((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg) * dtg_dy;

					res_vol->point_dv_dx[i_res] = dug_dx * cos((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg + M_PI)
													- ug * sin((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg + M_PI) * dtg_dx;
					res_vol->point_dv_dy[i_res] = dug_dy * cos((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg + M_PI)
													- ug * sin((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg + M_PI) * dtg_dy;
				}

				//update
				res_vol->point_u[i_res] = ug * sin((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg);
				res_vol->point_v[i_res] = ug * cos((res_vol->fcoriolis * res_vol->point_t[i_res]) - tg + M_PI);
			}
		}
		

		if (todo->calc_point_radial_vel_derivatives == 1) {
			res_vol->point_dradial_vel_ms_dx[i_res] = 
					(res_vol->point_du_dx[i_res] * res_vol->point_beamdir_enu_dx[i_res]) +
					(res_vol->point_dv_dx[i_res] * res_vol->point_beamdir_enu_dy[i_res]) +
					(res_vol->point_dw_dx[i_res] * res_vol->point_beamdir_enu_dz[i_res]);
			res_vol->point_dradial_vel_ms_dy[i_res] = 
					(res_vol->point_du_dy[i_res] * res_vol->point_beamdir_enu_dx[i_res]) +
					(res_vol->point_dv_dy[i_res] * res_vol->point_beamdir_enu_dy[i_res]) +
					(res_vol->point_dw_dy[i_res] * res_vol->point_beamdir_enu_dz[i_res]);
			res_vol->point_dradial_vel_ms_dz[i_res] = 
					(res_vol->point_du_dz[i_res] * res_vol->point_beamdir_enu_dx[i_res]) +
					(res_vol->point_dv_dz[i_res] * res_vol->point_beamdir_enu_dy[i_res]) +
					(res_vol->point_dw_dz[i_res] * res_vol->point_beamdir_enu_dz[i_res]);
			res_vol->point_dradial_vel_ms_dt[i_res] = 
					(res_vol->point_du_dt[i_res] * res_vol->point_beamdir_enu_dx[i_res]) +
					(res_vol->point_dv_dt[i_res] * res_vol->point_beamdir_enu_dy[i_res]) +
					(res_vol->point_dw_dt[i_res] * res_vol->point_beamdir_enu_dz[i_res]);				
		}
		
		if (todo->calc_point_shear == 1) {
			//calculate shear
			if (calcderivatives == 1) {
				res_vol->point_shear[i_res] =
					res_vol->point_beamdir_enu_dx[i_res] * 
						(	res_vol->point_du_dx[i_res] * res_vol->point_beamdir_enu_dx[i_res] +
							res_vol->point_dv_dx[i_res] * res_vol->point_beamdir_enu_dy[i_res] +
							res_vol->point_dw_dx[i_res] * res_vol->point_beamdir_enu_dz[i_res] 
						)
					+
					res_vol->point_beamdir_enu_dy[i_res] * 
						(	res_vol->point_du_dy[i_res] * res_vol->point_beamdir_enu_dx[i_res] +
							res_vol->point_dv_dy[i_res] * res_vol->point_beamdir_enu_dy[i_res] +
							res_vol->point_dw_dy[i_res] * res_vol->point_beamdir_enu_dz[i_res] 
						)
					+
					res_vol->point_beamdir_enu_dz[i_res] * 
						(	res_vol->point_du_dz[i_res] * res_vol->point_beamdir_enu_dx[i_res] +
							res_vol->point_dv_dz[i_res] * res_vol->point_beamdir_enu_dy[i_res] +
							res_vol->point_dw_dz[i_res] * res_vol->point_beamdir_enu_dz[i_res] 
						);
			}
		}
	}
	
	
	
		//res_vol->point_Z[i_res] = pow(10., res_vol->point_dBZ[i_res] / 10.);
	
	
	
	
	
	
	
	
	
	
	//apply noise if requested
	if (todo->apply_additive_noise == 1) {
		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			//apply additive noise
			tmpadj = 1. / (sqrt(res_vol->n));	//weighting factor, such that SNR applies to res vol weighted Z
			tmpZ = res_vol->point_Z[i_res];
			//tmpSNR = radarmeasurements->radar_constant * tmpZ / (radarmeasurements->noise_power * pow(res_vol->point_azel_range_m[i_res],2.));
			tmpSNR = res_vol->SNR;
			
			tmpR = ltqnorm(uniform_random());					//arbitrary number from Gaussian distribution
			tmpZ = fabs(tmpZ * (1. + (tmpadj * tmpR / tmpSNR )));		//apply

			//update
			res_vol->point_Z[i_res] 			= tmpZ;
			res_vol->point_dBZ[i_res]			= 10. * log10(tmpZ);

			tmpR = uniform_random();	//arbitrary number from Uniform distribution					
			res_vol->point_radial_vel_ms[i_res] = 
				((tmpSNR / (tmpSNR + tmpadj)) * res_vol->point_radial_vel_ms[i_res])
				+
				(tmpadj / (tmpSNR + tmpadj)) * res_vol->vmax * ((2. * tmpR) - 1.);
		}
	}

	if (todo->apply_multiplicative_noise == 1) {
		//apply multiplicative noise
		//radarmeasurements->Z[i_m] = radarmeasurements->Z[i_m] * (1. + (R / sqrt(n)))
		//not implemented yet
	}
	

	if (todo->calc_integrated_dBZ_derivatives == 1) {
		res_vol->integrated_ddBZ_dx	= 0.;
		res_vol->integrated_ddBZ_dy	= 0.;
		res_vol->integrated_ddBZ_dz	= 0.;
		res_vol->integrated_ddBZ_dt	= 0.;

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->integrated_ddBZ_dx += 
				res_vol->point_Z[i_res] * res_vol->point_ddBZ_dx[i_res];
			res_vol->integrated_ddBZ_dy += 
				res_vol->point_Z[i_res] * res_vol->point_ddBZ_dy[i_res];
			res_vol->integrated_ddBZ_dz += 
				res_vol->point_Z[i_res] * res_vol->point_ddBZ_dz[i_res];
			res_vol->integrated_ddBZ_dt += 
				res_vol->point_Z[i_res] * res_vol->point_ddBZ_dt[i_res];
		}
		
		res_vol->integrated_ddBZ_dx	/= res_vol->integrated_Z;
		res_vol->integrated_ddBZ_dy	/= res_vol->integrated_Z;
		res_vol->integrated_ddBZ_dz	/= res_vol->integrated_Z;
		res_vol->integrated_ddBZ_dt	/= res_vol->integrated_Z;
	}
	
	if (todo->calc_integrated_spectral_shape == 1) {
		//Calculate spectral shape
		mincdfP = 1.e-50;
		maxcdfP = 1. - 1.e-50;
		difcdfP = (maxcdfP - mincdfP) / res_vol->integrated_n_spectralshape;
		
		//walk through intervals
		for ( i_int = 0; i_int < res_vol->integrated_n_spectralshape; i_int++ ) {
			//unscaled
			lbound = ltqnorm(mincdfP + (i_int * difcdfP));
			ubound = ltqnorm(mincdfP + ((i_int + 1.)* difcdfP));
			
			//scaled
			lbound = (res_vol->integrated_radial_vel_width_ms * lbound) + res_vol->integrated_radial_vel_ms;
			ubound = (res_vol->integrated_radial_vel_width_ms * ubound) + res_vol->integrated_radial_vel_ms;
					
			this_int = 0.;
			
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				if ((lbound <= res_vol->point_radial_vel_ms[i_res]) &
					(res_vol->point_radial_vel_ms[i_res] < ubound)) {
						
					this_int += res_vol->point_Z[i_res];
				}
			}
			
			res_vol->integrated_spectral_shape[i_int] = this_int;
		}
	}


	if (todo->calc_integrated_u == 1) {
		res_vol->integrated_u 			= 0.;

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->integrated_u += 
				res_vol->point_Z[i_res] * res_vol->point_u[i_res];
		}

		res_vol->integrated_u /= res_vol->integrated_Z;		
	}
	
	if (todo->calc_integrated_u_derivatives == 1) {
		res_vol->integrated_du_dx 			= 0.;
		res_vol->integrated_du_dy 			= 0.;
		res_vol->integrated_du_dz 			= 0.;
		res_vol->integrated_du_dt 			= 0.;

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->integrated_du_dx += 
				res_vol->point_Z[i_res] * res_vol->point_du_dx[i_res];
			res_vol->integrated_du_dy += 
				res_vol->point_Z[i_res] * res_vol->point_du_dy[i_res];
			res_vol->integrated_du_dz += 
				res_vol->point_Z[i_res] * res_vol->point_du_dz[i_res];
			res_vol->integrated_du_dt += 
				res_vol->point_Z[i_res] * res_vol->point_du_dt[i_res];
		}

		res_vol->integrated_du_dx /= res_vol->integrated_Z;		
		res_vol->integrated_du_dy /= res_vol->integrated_Z;		
		res_vol->integrated_du_dz /= res_vol->integrated_Z;		
		res_vol->integrated_du_dt /= res_vol->integrated_Z;		
	}

	if (todo->calc_integrated_v == 1) {
		res_vol->integrated_v 			= 0.;

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->integrated_v += 
				res_vol->point_Z[i_res] * res_vol->point_v[i_res];
		}

		res_vol->integrated_v /= res_vol->integrated_Z;		
	}

	if (todo->calc_integrated_v_derivatives == 1) {
		res_vol->integrated_dv_dx 			= 0.;
		res_vol->integrated_dv_dy 			= 0.;
		res_vol->integrated_dv_dz 			= 0.;
		res_vol->integrated_dv_dt 			= 0.;

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->integrated_dv_dx += 
				res_vol->point_Z[i_res] * res_vol->point_dv_dx[i_res];
			res_vol->integrated_dv_dy += 
				res_vol->point_Z[i_res] * res_vol->point_dv_dy[i_res];
			res_vol->integrated_dv_dz += 
				res_vol->point_Z[i_res] * res_vol->point_dv_dz[i_res];
			res_vol->integrated_dv_dt += 
				res_vol->point_Z[i_res] * res_vol->point_dv_dt[i_res];
		}

		res_vol->integrated_dv_dx /= res_vol->integrated_Z;		
		res_vol->integrated_dv_dy /= res_vol->integrated_Z;		
		res_vol->integrated_dv_dz /= res_vol->integrated_Z;		
		res_vol->integrated_dv_dt /= res_vol->integrated_Z;		
	}
	
	if (todo->calc_integrated_w == 1) {
		res_vol->integrated_w 			= 0.;

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->integrated_w += 
				res_vol->point_Z[i_res] * res_vol->point_w[i_res];
		}

		res_vol->integrated_w /= res_vol->integrated_Z;		
	}

	if (todo->calc_integrated_w_derivatives == 1) {
		res_vol->integrated_dw_dx 			= 0.;
		res_vol->integrated_dw_dy 			= 0.;
		res_vol->integrated_dw_dz 			= 0.;
		res_vol->integrated_dw_dt 			= 0.;

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->integrated_dw_dx += 
				res_vol->point_Z[i_res] * res_vol->point_dw_dx[i_res];
			res_vol->integrated_dw_dy += 
				res_vol->point_Z[i_res] * res_vol->point_dw_dy[i_res];
			res_vol->integrated_dw_dz += 
				res_vol->point_Z[i_res] * res_vol->point_dw_dz[i_res];
			res_vol->integrated_dw_dt += 
				res_vol->point_Z[i_res] * res_vol->point_dw_dt[i_res];
		}

		res_vol->integrated_dw_dx /= res_vol->integrated_Z;		
		res_vol->integrated_dw_dy /= res_vol->integrated_Z;		
		res_vol->integrated_dw_dz /= res_vol->integrated_Z;		
		res_vol->integrated_dw_dt /= res_vol->integrated_Z;		
	}

	if (todo->calc_integrated_speed == 1) {
		res_vol->integrated_speed = sqrt( 
											  pow(res_vol->integrated_u,2.) 
											+ pow(res_vol->integrated_v,2.)
											+ pow(res_vol->integrated_w,2.)
										);
	}

	if (todo->calc_integrated_hspeed == 1) {
		tmp =  pow(res_vol->integrated_u,2.) + pow(res_vol->integrated_v,2.);
		if (tmp == 0.) {
			res_vol->integrated_hspeed = 0.;
		} else {
			res_vol->integrated_hspeed = sqrt(tmp);
		}
	}

	if (todo->calc_integrated_hdir == 1) {
		//meteorologic wind direction
		//u = -1. * H * Sin(Angle)
		//v = -1. * H * Cos(Angle)
		res_vol->integrated_hdir = M_PI + atan2(res_vol->integrated_u, res_vol->integrated_v);
	}

	//griddep calculation
	if (todo->calc_point_wind_griddep == 1) {
		//allocate memory
		if (res_vol->point_wind_griddep == NULL) {
			res_vol->point_wind_griddep = malloc(res_vol->n * windfield->field[0]->n * sizeof(double));}

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			outx[0] = res_vol->point_enu_x[i_res];
			outx[1] = res_vol->point_enu_y[i_res];
			outx[2] = res_vol->point_enu_z[i_res];
			outx[3] = res_vol->point_t[i_res];

			special_bilint(	
				&windfield->field[0]->n_dim,
				&windfield->field[0]->n,
				windfield->field[0]->int_ax_count,
				windfield->field[0]->int_ax_values,
				outx,
				&dummy_one,
				res_vol->point_wind_griddep + i_res * windfield->field[0]->n	//memory arithmetic
				);
		}		
	}









	int i_res;
	int i;
	
	double outx[4]; 		
	double derivatives[4];
	double u_edr, v_edr, w_edr;

	int bilint_special = 0;		//normal linear interpolation
	int calcderivatives;
	
	double point_enu_x_adj, point_enu_y_adj, point_enu_z_adj;
	double point_enu_speed;
	
	double tmpadj, tmpZ, tmpSNR, tmpR;
	double mincdfP, maxcdfP, difcdfP;
	double lbound;
	double ubound;
	double this_int;
	int i_int;
	double dummy_one = 1.0;
	double ug, tg;
	double u_adv, v_adv, w_adv, dt;
	double du_adv_dx, dv_adv_dy, dw_adv_dz;
	double point_u, point_v;
	double point_du_dx, point_du_dy, point_dv_dx, point_dv_dy;
	double dug_dx, dug_dy, dtg_dx, dtg_dy;
	double tmp;
	

	//griddep calculation
	if (todo->calc_point_edr_griddep == 1) {
		//allocate memory
		if (res_vol->point_edr_griddep == NULL)	
			{res_vol->point_edr_griddep = malloc(res_vol->n * edrfield->field->n * sizeof(double));}


		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			outx[0] = res_vol->point_enu_x[i_res];
			outx[1] = res_vol->point_enu_y[i_res];
			outx[2] = res_vol->point_enu_z[i_res];
			outx[3] = res_vol->point_t[i_res];

			special_bilint(	
				&edrfield->field->n_dim,
				&edrfield->field->n,
				edrfield->field->int_ax_count,
				edrfield->field->int_ax_values,
				outx,
				&dummy_one,
				res_vol->point_edr_griddep + i_res * edrfield->field->n	//memory arithmetic
				);
		}
	}
	
	//griddep calculation
	if (todo->calc_point_scatterer_griddep == 1) {
		//allocate memory
		if (res_vol->point_scatterer_griddep == NULL) 	
			{res_vol->point_scatterer_griddep = malloc(res_vol->n * scattererfield->field->n * sizeof(double));}


		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			outx[0] = res_vol->point_enu_x[i_res];
			outx[1] = res_vol->point_enu_y[i_res];
			outx[2] = res_vol->point_enu_z[i_res];
			outx[3] = res_vol->point_t[i_res];

			special_bilint(	
				&scattererfield->field->n_dim,
				&scattererfield->field->n,
				scattererfield->field->int_ax_count,
				scattererfield->field->int_ax_values,
				outx,
				&dummy_one,
				res_vol->point_scatterer_griddep + i_res * scattererfield->field->n	//memory arithmetic
				);
		}		
	}
	
	
	//griddep calculation
	if (todo->calc_integrated_wind_griddep == 1) {
		//allocate memory
		if (res_vol->integrated_wind_griddep == NULL) 		
			{res_vol->integrated_wind_griddep = malloc(windfield->field[0]->n * sizeof(double));}

		for ( i = 0; i < windfield->field[0]->n; i++ ) {
			res_vol->integrated_wind_griddep[i] = 0.;
		}

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			for ( i = 0; i < windfield->field[0]->n; i++ ) {
				res_vol->integrated_wind_griddep[i] +=
					res_vol->point_Z[i_res] * res_vol->point_wind_griddep[i_res * windfield->field[0]->n + i];
			}
		}
		
		for ( i = 0; i < windfield->field[0]->n; i++ ) {
			res_vol->integrated_wind_griddep[i] /= res_vol->integrated_Z;
		}
	}
	
	//griddep calculation
	if (todo->calc_integrated_edr_griddep == 1) {
		//allocate memory
		if (res_vol->integrated_edr_griddep == NULL) 		
			{res_vol->integrated_edr_griddep = malloc(edrfield->field->n * sizeof(double));}

		for ( i = 0; i < edrfield->field->n; i++ ) {
			res_vol->integrated_edr_griddep[i] = 0.;
		}

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			for ( i = 0; i < edrfield->field->n; i++ ) {
				res_vol->integrated_edr_griddep[i] +=
					res_vol->point_Z[i_res] * res_vol->point_edr_griddep[i_res * edrfield->field->n + i];
			}
		}
		
		for ( i = 0; i < edrfield->field->n; i++ ) {
			res_vol->integrated_edr_griddep[i] /= res_vol->integrated_Z;
		}
	}
	
	//griddep calculation
	if (todo->calc_integrated_scatterer_griddep == 1) {
		//allocate memory
		if (res_vol->integrated_scatterer_griddep == NULL) 	
			{res_vol->integrated_scatterer_griddep = malloc(scattererfield->field->n * sizeof(double));}

		for ( i = 0; i < scattererfield->field->n; i++ ) {
			res_vol->integrated_scatterer_griddep[i] = 0.;
		}

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			for ( i = 0; i < scattererfield->field->n; i++ ) {
				res_vol->integrated_scatterer_griddep[i] +=
					res_vol->point_Z[i_res] * res_vol->point_scatterer_griddep[i_res * scattererfield->field->n + i];
			}
		}
		
		for ( i = 0; i < scattererfield->field->n; i++ ) {
			res_vol->integrated_scatterer_griddep[i] /= res_vol->integrated_Z;
		}
	}


	//extra things
	if (todo->calc_integrated_u_derivatives == 1) {
		if (todo->calc_integrated_dBZ_derivatives == 1) {
			res_vol->integrated_ddBZ_du = res_vol->integrated_ddBZ_dx * -1. * outx[3];
		}
		if (todo->calc_integrated_radial_vel_ms_derivatives == 1) {
				res_vol->integrated_dradial_vel_ms_du = res_vol->integrated_dradial_vel_ms_dx * -1. * outx[3]; //TBD
		}
	}
	if (todo->calc_integrated_v_derivatives == 1) {
		if (todo->calc_integrated_dBZ_derivatives == 1) {
			res_vol->integrated_ddBZ_dv = res_vol->integrated_ddBZ_dy * -1. * outx[3];
		}
		if (todo->calc_integrated_radial_vel_ms_derivatives == 1) {
			res_vol->integrated_dradial_vel_ms_dv = res_vol->integrated_dradial_vel_ms_dy * -1. * outx[3];
		}
	}
	if (todo->calc_integrated_w_derivatives == 1) {
		if (todo->calc_integrated_dBZ_derivatives == 1) {
			res_vol->integrated_ddBZ_dw = res_vol->integrated_ddBZ_dz * -1. * outx[3];
		}
		if (todo->calc_integrated_radial_vel_ms_derivatives == 1) {
			res_vol->integrated_dradial_vel_ms_dw =  res_vol->integrated_dradial_vel_ms_dz * -1. * outx[3]; //TBD
		}
	}


	if (debug) diff = clock() - start;
	if (debug) printf("calc_res_vol_wind_edr: %.6e seconds\n", 1000. * diff / CLOCKS_PER_SEC);


	for ( i_int = 0; i_int < radarmeasurements->n_spectralshape; i_int++ ) {
		radarmeasurements->radial_vel_shape[radarmeasurements->n_spectralshape * i_m + i_int] = res_vol->integrated_spectral_shape[i_int];
	}		

*/




//radarfilter_initialize_resolution_volume
//- allocate memory
//i_mode,	//0 = simulation mode, 1 = retrieval mode

void radarfilter_initialize_resolution_volume(t_zephyros_config *cfg, int i_mode, t_radarfilter_res_vol** pres_vol, t_radarfilter_todolist *todo)
{
	t_radarfilter_res_vol *res_vol = calloc(1, sizeof(t_radarfilter_res_vol));
	int i_res, i_psd, i_par;
	int i_parmod;
	int i;
	
	t_zephyros_radarfilter	*rcfg;
	t_zephyros_windfield	*mywindfield;
	t_zephyros_scattererfield	*myscattererfield;
	if (i_mode == 0) {
		rcfg = cfg->simulation->radarfilter;
		mywindfield = cfg->simulation->windfield;
		myscattererfield = cfg->simulation->scattererfield;
		}
	if (i_mode == 1) {
		rcfg = cfg->retrieval->radarfilter;
		mywindfield = cfg->retrieval->post_windfield;
		myscattererfield = cfg->retrieval->post_scattererfield;		
		}

	res_vol->n_beam_range 			= rcfg->n_beam_range;
	res_vol->n_beam_theta 			= rcfg->n_beam_theta;
	res_vol->n_beam_phi 			= rcfg->n_beam_phi;
	res_vol->n_t					= rcfg->n_t;
	res_vol->n						= res_vol->n_beam_range * res_vol->n_beam_theta * res_vol->n_beam_phi * res_vol->n_t;

	//check configuration to see whether there is parametric turbulence
	res_vol->parametric_turbulence = 0;
	for ( i = 0; i < mywindfield->nturbulences; i++ ) {
		if ((mywindfield->turbulence[i] != NULL) & (mywindfield->turbulence[i]->type == 5)) {
			res_vol->parametric_turbulence = 1;
		}
	}
	
	if (res_vol->parametric_turbulence == 1) {
		res_vol->n_parmod_az			= rcfg->n_parmod_az;
		res_vol->n_parmod_el			= 2 * rcfg->n_parmod_el; //both positive and negative
		res_vol->n_parmod				= res_vol->n_parmod_el * res_vol->n_parmod_az;
		
		if (todo->der_edr13) {
			res_vol->n_parmod_plus = 2 * res_vol->n_parmod;
		} else {
			res_vol->n_parmod_plus = res_vol->n_parmod;			
		}
	} else {
		res_vol->n_parmod 				= 1;
		res_vol->n_parmod_plus			= 1;
	}
		
	res_vol->n_psd					= myscattererfield->npsd;
	res_vol->n_diameters			= calloc(myscattererfield->npsd, sizeof(int));
	for ( i_psd = 0; i_psd < myscattererfield->npsd; i_psd++ ) {
		if (myscattererfield->psd[i_psd] != NULL) {
			res_vol->n_diameters[i_psd] = myscattererfield->psd[i_psd]->n_diameters;
		} else {
			res_vol->n_diameters[i_psd] = 0;
		}
	}

	//coordinates
	res_vol->subvolume_coor		= calloc(res_vol->n, sizeof(t_zephyros_coordinates*));
	for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
		coordinates_initialize_coor(res_vol->subvolume_coor + i_res);
	}
	
	//subvolume scattering widget
	res_vol->subvolume_scat	= calloc(res_vol->n_psd, sizeof(t_zephyros_particles_widget**));
	for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
		res_vol->subvolume_scat[i_psd] = calloc(res_vol->n_diameters[i_psd], sizeof(t_zephyros_particles_widget*));
		for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++) {
			res_vol->subvolume_scat[i_psd][i_par] = calloc(1, sizeof(t_zephyros_particles_widget));
		}
	}

	//subvolume particle number density
	if (todo->calc_number_density) res_vol->subvolume_number_density_m3		= calloc(res_vol->n, sizeof(double**));
	if (todo->calc_number_density) res_vol->subvolume_ln_number_density_m3_der	= calloc(res_vol->n, sizeof(double***));
	for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
		if (todo->calc_number_density) res_vol->subvolume_number_density_m3[i_res]			= calloc(res_vol->n_psd, sizeof(double*));
		if (todo->calc_number_density) res_vol->subvolume_ln_number_density_m3_der[i_res]		= calloc(res_vol->n_psd, sizeof(double**));
		for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
			if (todo->calc_number_density) res_vol->subvolume_number_density_m3[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd], sizeof(double));
			if (todo->calc_number_density) res_vol->subvolume_ln_number_density_m3_der[i_res][i_psd]	= calloc(res_vol->n_diameters[i_psd], sizeof(double*));
			for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
				if (todo->calc_number_density) res_vol->subvolume_ln_number_density_m3_der[i_res][i_psd][i_par] = calloc(4, sizeof(double));
			}
		}
	}

	//subvolume eta_hh, eta_hv, eta_vv
	if (todo->calc_eta_hh) res_vol->subvolume_eta_hh	= calloc(res_vol->n, sizeof(double***));
	if (todo->calc_eta_hv) res_vol->subvolume_eta_hv	= calloc(res_vol->n, sizeof(double***));
	if (todo->calc_eta_vh) res_vol->subvolume_eta_vh	= calloc(res_vol->n, sizeof(double***));
	if (todo->calc_eta_vv) res_vol->subvolume_eta_vv	= calloc(res_vol->n, sizeof(double***));
	if (todo->calc_eta_ShhSvvc) res_vol->subvolume_eta_ShhSvvc	= calloc(res_vol->n, sizeof(double***));
	if (todo->calc_eta_ShhShvc) res_vol->subvolume_eta_ShhShvc	= calloc(res_vol->n, sizeof(double***));
	if (todo->calc_eta_SvvSvhc) res_vol->subvolume_eta_SvvSvhc	= calloc(res_vol->n, sizeof(double***));
	for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
		if (todo->calc_eta_hh) res_vol->subvolume_eta_hh[i_res]	= calloc(res_vol->n_psd, sizeof(double**));
		if (todo->calc_eta_hv) res_vol->subvolume_eta_hv[i_res]	= calloc(res_vol->n_psd, sizeof(double**));
		if (todo->calc_eta_vh) res_vol->subvolume_eta_vh[i_res]	= calloc(res_vol->n_psd, sizeof(double**));
		if (todo->calc_eta_vv) res_vol->subvolume_eta_vv[i_res]	= calloc(res_vol->n_psd, sizeof(double**));
		if (todo->calc_eta_ShhSvvc) res_vol->subvolume_eta_ShhSvvc[i_res]	= calloc(res_vol->n_psd, sizeof(double**));
		if (todo->calc_eta_ShhShvc) res_vol->subvolume_eta_ShhShvc[i_res]	= calloc(res_vol->n_psd, sizeof(double**));
		if (todo->calc_eta_SvvSvhc) res_vol->subvolume_eta_SvvSvhc[i_res]	= calloc(res_vol->n_psd, sizeof(double**));
		for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
			if (todo->calc_eta_hh) res_vol->subvolume_eta_hh[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd], sizeof(double*));
			if (todo->calc_eta_hv) res_vol->subvolume_eta_hv[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd], sizeof(double*));
			if (todo->calc_eta_vh) res_vol->subvolume_eta_vh[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd], sizeof(double*));
			if (todo->calc_eta_vv) res_vol->subvolume_eta_vv[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd], sizeof(double*));
			if (todo->calc_eta_ShhSvvc) res_vol->subvolume_eta_ShhSvvc[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd], sizeof(double*));
			if (todo->calc_eta_ShhShvc) res_vol->subvolume_eta_ShhShvc[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd], sizeof(double*));
			if (todo->calc_eta_SvvSvhc) res_vol->subvolume_eta_SvvSvhc[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd], sizeof(double*));
			for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
				if (todo->calc_eta_hh) res_vol->subvolume_eta_hh[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus, sizeof(double));
				if (todo->calc_eta_hv) res_vol->subvolume_eta_hv[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus, sizeof(double));
				if (todo->calc_eta_vh) res_vol->subvolume_eta_vh[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus, sizeof(double));
				if (todo->calc_eta_vv) res_vol->subvolume_eta_vv[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus, sizeof(double));
				if (todo->calc_eta_ShhSvvc) res_vol->subvolume_eta_ShhSvvc[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus, sizeof(double));
				if (todo->calc_eta_ShhShvc) res_vol->subvolume_eta_ShhShvc[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus, sizeof(double));
				if (todo->calc_eta_SvvSvhc) res_vol->subvolume_eta_SvvSvhc[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus, sizeof(double));
			}	
		}
	}
	
	if (todo->calc_air_velocity) {
		res_vol->subvolume_air_u		= calloc(res_vol->n, sizeof(double));
		res_vol->subvolume_air_v		= calloc(res_vol->n, sizeof(double));
		res_vol->subvolume_air_w		= calloc(res_vol->n, sizeof(double));

		res_vol->subvolume_air_u_der 	= calloc(res_vol->n, sizeof(double*));
		res_vol->subvolume_air_v_der 	= calloc(res_vol->n, sizeof(double*));
		res_vol->subvolume_air_w_der 	= calloc(res_vol->n, sizeof(double*));
		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {	
			res_vol->subvolume_air_u_der[i_res]	= calloc(4, sizeof(double));
			res_vol->subvolume_air_v_der[i_res]	= calloc(4, sizeof(double));
			res_vol->subvolume_air_w_der[i_res]	= calloc(4, sizeof(double));
		}	
	}
	
	if (todo->calc_particle_velocity) {
		res_vol->subvolume_particle_u	= calloc(res_vol->n , sizeof(double***));
		res_vol->subvolume_particle_v	= calloc(res_vol->n , sizeof(double***));
		res_vol->subvolume_particle_w	= calloc(res_vol->n , sizeof(double***));
		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->subvolume_particle_u[i_res]	= calloc(res_vol->n_psd , sizeof(double**));
			res_vol->subvolume_particle_v[i_res]	= calloc(res_vol->n_psd , sizeof(double**));
			res_vol->subvolume_particle_w[i_res]	= calloc(res_vol->n_psd , sizeof(double**));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				res_vol->subvolume_particle_u[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd] , sizeof(double*));
				res_vol->subvolume_particle_v[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd] , sizeof(double*));
				res_vol->subvolume_particle_w[i_res][i_psd]		= calloc(res_vol->n_diameters[i_psd] , sizeof(double*));			
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					res_vol->subvolume_particle_u[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus , sizeof(double));
					res_vol->subvolume_particle_v[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus , sizeof(double));
					res_vol->subvolume_particle_w[i_res][i_psd][i_par]	= calloc(res_vol->n_parmod_plus , sizeof(double));
				}
			}
		}
	}
	
	//TBD: calc_particle_velocity_der
	
	if (todo->calc_particle_direction) {
		res_vol->subvolume_particle_dir	= calloc(res_vol->n, sizeof(double****));
		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->subvolume_particle_dir[i_res]	= calloc(res_vol->n_psd, sizeof(double***));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				res_vol->subvolume_particle_dir[i_res][i_psd]	= calloc(res_vol->n_diameters[i_psd], sizeof(double**));
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					res_vol->subvolume_particle_dir[i_res][i_psd][i_par] = calloc(res_vol->n_parmod_plus, sizeof(double*));
					for ( i_parmod = 0; i_parmod < res_vol->n_parmod_plus; i_parmod++ ) {
						res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod] = calloc(4, sizeof(double));
					}
				}
			}
		}
	}
		
	//subvolume particle Doppler velocity
	if (todo->calc_Doppler_mean_velocity) {
		res_vol->subvolume_particle_Doppler_velocity	= calloc(res_vol->n, sizeof(double***));
		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			res_vol->subvolume_particle_Doppler_velocity[i_res]	= calloc(res_vol->n_psd, sizeof(double**));
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd] = calloc(res_vol->n_diameters[i_psd], sizeof(double*));
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par] = calloc(res_vol->n_parmod_plus, sizeof(double));
				}
			}
		}	
	}
	
/*
	if (todo->calc_edr) {
		res_vol->subvolume_edr				= malloc(res_vol->n * sizeof(double));
		res_vol->subvolume_edr13			= malloc(res_vol->n * sizeof(double));
	}
	
	if (todo->calc_maxL) {
		res_vol->subvolume_maxL				= malloc(res_vol->n * sizeof(double));
	}
*/
	//~ res_vol->point_Z				= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dBZ				= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_attenuation		= malloc(res_vol->n * sizeof(double));
//~ 
	//~ res_vol->point_df_adv_x_dx		= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_df_adv_y_dy		= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_df_adv_z_dz		= malloc(res_vol->n * sizeof(double));
//~ 
	//~ res_vol->point_ddBZ_dx			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_ddBZ_dy			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_ddBZ_dz			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_ddBZ_dt			= malloc(res_vol->n * sizeof(double));
//~ 
	//~ res_vol->point_dattenuation_dx	= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dattenuation_dy	= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dattenuation_dz	= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dattenuation_dt	= malloc(res_vol->n * sizeof(double));
//~ 
	//~ res_vol->point_du_dx			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_du_dy			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_du_dz			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_du_dt			= malloc(res_vol->n * sizeof(double));
	//~ 
	//~ res_vol->point_dv_dx			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dv_dy			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dv_dz			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dv_dt			= malloc(res_vol->n * sizeof(double));
	//~ 
	//~ res_vol->point_dw_dx			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dw_dy			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dw_dz			= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dw_dt			= malloc(res_vol->n * sizeof(double));
//~ 
	//~ res_vol->point_shear			= malloc(res_vol->n * sizeof(double));
//~ 
	//~ res_vol->point_dradial_vel_ms_dx= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dradial_vel_ms_dy= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dradial_vel_ms_dz= malloc(res_vol->n * sizeof(double));
	//~ res_vol->point_dradial_vel_ms_dt= malloc(res_vol->n * sizeof(double));
//~ 

	*pres_vol = res_vol;
}

void radarfilter_free_resolution_volume(t_zephyros_config *cfg, int i_mode, t_radarfilter_res_vol **pres_vol, t_radarfilter_todolist *todo)
{
	t_radarfilter_res_vol *res_vol = *pres_vol;
	int i_res, i_psd, i_par;
	int i_parmod;
	
	if (res_vol != NULL) {
		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			coordinates_free_coor(&(res_vol->subvolume_coor[i_res]));
		}
		
		//subvolume scattering widget
		for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
			for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++) {
				util_safe_free(&(res_vol->subvolume_scat[i_psd][i_par]));
			}
			util_safe_free(&(res_vol->subvolume_scat[i_psd]));
		}
		util_safe_free(&(res_vol->subvolume_scat));
		

		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					if (todo->calc_number_density) {free(res_vol->subvolume_ln_number_density_m3_der[i_res][i_psd][i_par]);}				
				}
				if (todo->calc_number_density) {free(res_vol->subvolume_number_density_m3[i_res][i_psd]);}
				if (todo->calc_number_density) {free(res_vol->subvolume_ln_number_density_m3_der[i_res][i_psd]);}
			}
			if (todo->calc_number_density) {free(res_vol->subvolume_number_density_m3[i_res]);}
			if (todo->calc_number_density) {free(res_vol->subvolume_ln_number_density_m3_der[i_res]);}
		}
		if (todo->calc_number_density) {free(res_vol->subvolume_number_density_m3);}
		if (todo->calc_number_density) {free(res_vol->subvolume_ln_number_density_m3_der);}


		if (todo->calc_air_velocity) {
			util_safe_free(&res_vol->subvolume_air_u);
			util_safe_free(&res_vol->subvolume_air_v);
			util_safe_free(&res_vol->subvolume_air_w);
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {			
				util_safe_free(&(res_vol->subvolume_air_u_der[i_res]));
				util_safe_free(&(res_vol->subvolume_air_v_der[i_res]));
				util_safe_free(&(res_vol->subvolume_air_w_der[i_res]));
			}
			util_safe_free(&res_vol->subvolume_air_u_der);
			util_safe_free(&res_vol->subvolume_air_v_der);
			util_safe_free(&res_vol->subvolume_air_w_der);
		}
		
		if (todo->calc_particle_velocity) {
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
					for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
						util_safe_free(&(res_vol->subvolume_particle_u[i_res][i_psd][i_par]));
						util_safe_free(&(res_vol->subvolume_particle_v[i_res][i_psd][i_par]));
						util_safe_free(&(res_vol->subvolume_particle_w[i_res][i_psd][i_par]));
					}
					util_safe_free(&(res_vol->subvolume_particle_u[i_res][i_psd]));
					util_safe_free(&(res_vol->subvolume_particle_v[i_res][i_psd]));
					util_safe_free(&(res_vol->subvolume_particle_w[i_res][i_psd]));
				}
				util_safe_free(&(res_vol->subvolume_particle_u[i_res]));
				util_safe_free(&(res_vol->subvolume_particle_v[i_res]));
				util_safe_free(&(res_vol->subvolume_particle_w[i_res]));
			}
			util_safe_free(&(res_vol->subvolume_particle_u));
			util_safe_free(&(res_vol->subvolume_particle_v));
			util_safe_free(&(res_vol->subvolume_particle_w));
		}

		if (todo->calc_particle_direction) {
			//subvolume particle uvw
			//subvolume particle dir
			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
					for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
						for ( i_parmod = 0; i_parmod < res_vol->n_parmod_plus; i_parmod++ ) {
							util_safe_free(&(res_vol->subvolume_particle_dir[i_res][i_psd][i_par][i_parmod]));			
						}
						util_safe_free(&(res_vol->subvolume_particle_dir[i_res][i_psd][i_par]));			
					}
					util_safe_free(&(res_vol->subvolume_particle_dir[i_res][i_psd]));
				}
				util_safe_free(&(res_vol->subvolume_particle_dir[i_res]));
			}
			util_safe_free(&(res_vol->subvolume_particle_dir));
		}

		//subvolume eta hh, hv, vv
		for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
			for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
					if (todo->calc_eta_hh) {util_safe_free(&(res_vol->subvolume_eta_hh[i_res][i_psd][i_par]));}
					if (todo->calc_eta_hv) {util_safe_free(&(res_vol->subvolume_eta_hv[i_res][i_psd][i_par]));}
					if (todo->calc_eta_vh) {util_safe_free(&(res_vol->subvolume_eta_vh[i_res][i_psd][i_par]));}
					if (todo->calc_eta_vv) {util_safe_free(&(res_vol->subvolume_eta_vv[i_res][i_psd][i_par]));}
					if (todo->calc_eta_ShhSvvc) {util_safe_free(&(res_vol->subvolume_eta_ShhSvvc[i_res][i_psd][i_par]));}
					if (todo->calc_eta_ShhShvc) {util_safe_free(&(res_vol->subvolume_eta_ShhShvc[i_res][i_psd][i_par]));}
					if (todo->calc_eta_SvvSvhc) {util_safe_free(&(res_vol->subvolume_eta_SvvSvhc[i_res][i_psd][i_par]));}
				}
				if (todo->calc_eta_hh) {util_safe_free(&(res_vol->subvolume_eta_hh[i_res][i_psd]));}
				if (todo->calc_eta_hv) {util_safe_free(&(res_vol->subvolume_eta_hv[i_res][i_psd]));}
				if (todo->calc_eta_vh) {util_safe_free(&(res_vol->subvolume_eta_vh[i_res][i_psd]));}
				if (todo->calc_eta_vv) {util_safe_free(&(res_vol->subvolume_eta_vv[i_res][i_psd]));}
				if (todo->calc_eta_ShhSvvc) {util_safe_free(&(res_vol->subvolume_eta_ShhSvvc[i_res][i_psd]));}
				if (todo->calc_eta_ShhShvc) {util_safe_free(&(res_vol->subvolume_eta_ShhShvc[i_res][i_psd]));}
				if (todo->calc_eta_SvvSvhc) {util_safe_free(&(res_vol->subvolume_eta_SvvSvhc[i_res][i_psd]));}
			}
			if (todo->calc_eta_hh) {util_safe_free(&(res_vol->subvolume_eta_hh[i_res]));}
			if (todo->calc_eta_hv) {util_safe_free(&(res_vol->subvolume_eta_hv[i_res]));}
			if (todo->calc_eta_vh) {util_safe_free(&(res_vol->subvolume_eta_vh[i_res]));}
			if (todo->calc_eta_vv) {util_safe_free(&(res_vol->subvolume_eta_vv[i_res]));}
			if (todo->calc_eta_ShhSvvc) {util_safe_free(&(res_vol->subvolume_eta_ShhSvvc[i_res]));}
			if (todo->calc_eta_ShhShvc) {util_safe_free(&(res_vol->subvolume_eta_ShhShvc[i_res]));}
			if (todo->calc_eta_SvvSvhc) {util_safe_free(&(res_vol->subvolume_eta_SvvSvhc[i_res]));}
		}
		if (todo->calc_eta_hh) {util_safe_free(&(res_vol->subvolume_eta_hh));}
		if (todo->calc_eta_hv) {util_safe_free(&(res_vol->subvolume_eta_hv));}
		if (todo->calc_eta_vh) {util_safe_free(&(res_vol->subvolume_eta_vh));}
		if (todo->calc_eta_vv) {util_safe_free(&(res_vol->subvolume_eta_vv));}
		if (todo->calc_eta_ShhSvvc) {util_safe_free(&(res_vol->subvolume_eta_ShhSvvc));}
		if (todo->calc_eta_ShhShvc) {util_safe_free(&(res_vol->subvolume_eta_ShhShvc));}
		if (todo->calc_eta_SvvSvhc) {util_safe_free(&(res_vol->subvolume_eta_SvvSvhc));}
		
		//subvolume particle Doppler velocity
		if (todo->calc_Doppler_mean_velocity) {

			for ( i_res = 0; i_res < res_vol->n; i_res++ ) {
				for ( i_psd = 0; i_psd < res_vol->n_psd; i_psd++ ) {
					for ( i_par = 0; i_par < res_vol->n_diameters[i_psd]; i_par++ ) {
						util_safe_free(&(res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd][i_par]));
					}
					util_safe_free(&(res_vol->subvolume_particle_Doppler_velocity[i_res][i_psd]));
				}
				util_safe_free(&(res_vol->subvolume_particle_Doppler_velocity[i_res]));
			}
			util_safe_free(&(res_vol->subvolume_particle_Doppler_velocity));
		}

		//if (res_vol->subvolume_particle_u != NULL) {free(res_vol->subvolume_particle_u);}
		//if (res_vol->subvolume_particle_v != NULL) {free(res_vol->subvolume_particle_v);}
		//if (res_vol->subvolume_particle_v != NULL) {free(res_vol->subvolume_particle_v);}
		
		//if (res_vol->subvolume_particle_radial_vel != NULL) {free(res_vol->subvolume_particle_radial_vel);}
		
		//if (res_vol->subvolume_edr != NULL) {free(res_vol->subvolume_edr);}
		//if (res_vol->subvolume_edr13 != NULL) {free(res_vol->subvolume_edr13);}
		//if (res_vol->subvolume_maxL != NULL) {free(res_vol->subvolume_maxL);}
	}
	res_vol = NULL;
}











void radarfilter_initialize_todolist(t_zephyros_config *cfg, int i_mode, t_radarfilter_todolist **ptodolist)
{
	t_radarfilter_todolist *todolist = malloc(sizeof(t_radarfilter_todolist));

	t_zephyros_radarfilter	*rcfg;
	if (i_mode == 0) rcfg = cfg->simulation->radarfilter; 
	if (i_mode == 1) rcfg = cfg->retrieval->radarfilter; 
	
	//reset
	todolist->calc_center_coordinates = 1;
	todolist->calc_number_density = 1;
	todolist->calc_number_density_der = 0;
	todolist->calc_cross_sections = 1;
	todolist->calc_eta_hh = 0;
	todolist->calc_eta_hv = 0;
	todolist->calc_eta_vh = 0;
	todolist->calc_eta_vv = 0;
	todolist->calc_eta_ShhSvvc = 0;
	todolist->calc_eta_ShhShvc = 0;
	todolist->calc_eta_SvvSvhc = 0;
	todolist->calc_air_velocity = 1;
	todolist->calc_air_velocity_der = 0;
	todolist->calc_terminal_fall_speeds = 1;
	todolist->calc_particle_velocity = 1;
	todolist->calc_particle_velocity_der = 0;
	todolist->calc_particle_direction = 1;

	todolist->calc_Doppler_mean_velocity = 0;
	todolist->calc_Doppler_spectrum_widths = 0;
	todolist->calc_Doppler_spectrum = 0;
	todolist->calc_Doppler_spectrum_ShhSvvc = 0;
	todolist->calc_Doppler_spectrum_ShhShvc = 0;
	todolist->calc_Doppler_spectrum_SvvSvhc = 0;
	
	todolist->calc_eta_i_hh = 0;
	todolist->calc_eta_i_hv = 0;
	todolist->calc_eta_i_vh = 0;
	todolist->calc_eta_i_vv = 0;
	
	todolist->calc_spectrum_eta_i_hh = 0;
	todolist->calc_spectrum_eta_i_hv = 0;
	todolist->calc_spectrum_eta_i_vh = 0;
	todolist->calc_spectrum_eta_i_vv = 0;
	
	todolist->der_edr13 = 0;
		
	radarfilter_prepare_todolist(cfg, i_mode, todolist);

    //res_vol->advection_integration_steps 							= 5;
    *ptodolist = todolist;
}

void radarfilter_prepare_todolist(t_zephyros_config *cfg, int i_mode, t_radarfilter_todolist *todolist)
{
	t_zephyros_radarfilter	*rcfg;
	if (i_mode == 0) rcfg = cfg->simulation->radarfilter; 
	if (i_mode == 1) rcfg = cfg->retrieval->radarfilter; 
	
	//follow some logic so that everything is correctly calculated and initialized
	if (rcfg->filter_dBZ_hh) {
		todolist->calc_eta_hh = 1;
	}
	if (rcfg->filter_dBZ_hv) {
		todolist->calc_eta_hv = 1;
	}
	if (rcfg->filter_dBZ_vh) {
		todolist->calc_eta_vh = 1;
	}
	if (rcfg->filter_dBZ_vv) {
		todolist->calc_eta_vv = 1;
	}
	if (rcfg->filter_dBZdr) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_vv = 1;
	}
	if (rcfg->filter_dBLdr) {
		todolist->calc_eta_hv = 1;
		todolist->calc_eta_vv = 1;
	}
	if (rcfg->filter_rho_co) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_vv = 1;
		todolist->calc_eta_ShhSvvc = 1;
	}
	if (rcfg->filter_rho_cxh) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_hv = 1;
		todolist->calc_eta_ShhShvc = 1;
	}
	if (rcfg->filter_rho_cxv) {
		todolist->calc_eta_vh = 1;
		todolist->calc_eta_vv = 1;
		todolist->calc_eta_SvvSvhc = 1;

	}
	if (rcfg->filter_Doppler_velocity_hh_ms) {
		todolist->calc_eta_hh = 1;
		todolist->calc_Doppler_mean_velocity = 1;
	}
	if (rcfg->filter_Doppler_velocity_hv_ms) {
		todolist->calc_eta_hv = 1;
		todolist->calc_Doppler_mean_velocity = 1;
	}
	if (rcfg->filter_Doppler_velocity_vh_ms) {
		todolist->calc_eta_vh = 1;
		todolist->calc_Doppler_mean_velocity = 1;
	}
	if (rcfg->filter_Doppler_velocity_vv_ms) {
		todolist->calc_eta_vv = 1;
		todolist->calc_Doppler_mean_velocity = 1;
	}
	if (rcfg->filter_Doppler_spectralwidth_hh_ms) {
		todolist->calc_eta_hh = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
	}
	if (rcfg->filter_Doppler_spectralwidth_hv_ms) {
		todolist->calc_eta_hv = 1;
		rcfg->filter_Doppler_velocity_hv_ms = 1;
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
	}
	if (rcfg->filter_Doppler_spectralwidth_vh_ms) {
		todolist->calc_eta_vh = 1;
		rcfg->filter_Doppler_velocity_vh_ms = 1;
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
	}
	if (rcfg->filter_Doppler_spectralwidth_vv_ms) {
		todolist->calc_eta_vv = 1;
		rcfg->filter_Doppler_velocity_vv_ms = 1;
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
	}
	if (rcfg->filter_Doppler_spectrum_dBZ_hh) {
		todolist->calc_eta_hh = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		rcfg->filter_Doppler_spectralwidth_hh_ms = 1;
		
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
		todolist->calc_Doppler_spectrum = 1;
	}
	if (rcfg->filter_Doppler_spectrum_dBZ_hv) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_hv = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		rcfg->filter_Doppler_spectralwidth_hh_ms = 1;

		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
		todolist->calc_Doppler_spectrum = 1;
	}
	if (rcfg->filter_Doppler_spectrum_dBZ_vh) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_vh = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		rcfg->filter_Doppler_spectralwidth_hh_ms = 1;

		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
		todolist->calc_Doppler_spectrum = 1;
	}
	if (rcfg->filter_Doppler_spectrum_dBZ_vv) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_vv = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		rcfg->filter_Doppler_spectralwidth_hh_ms = 1;

		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
		todolist->calc_Doppler_spectrum = 1;
	}
	if (rcfg->filter_specific_dBZdr) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_vv = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		rcfg->filter_Doppler_spectralwidth_hh_ms = 1;
		
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
		todolist->calc_Doppler_spectrum = 1;
		
		rcfg->filter_Doppler_spectrum_dBZ_hh = 1;
		rcfg->filter_Doppler_spectrum_dBZ_vv = 1;
	}
	if (rcfg->filter_specific_dBLdr) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_hv = 1;
		todolist->calc_eta_vv = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		rcfg->filter_Doppler_spectralwidth_hh_ms = 1;
		
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
		todolist->calc_Doppler_spectrum = 1;

		rcfg->filter_Doppler_spectrum_dBZ_hv = 1;
		rcfg->filter_Doppler_spectrum_dBZ_vv = 1;
	}
	if (rcfg->filter_specific_rho_co) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_vv = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		rcfg->filter_Doppler_spectralwidth_hh_ms = 1;
		
		todolist->calc_Doppler_spectrum_ShhSvvc = 1;
		rcfg->filter_Doppler_spectrum_dBZ_hh = 1;
		rcfg->filter_Doppler_spectrum_dBZ_vv = 1;
		
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
		todolist->calc_Doppler_spectrum = 1;
	}
	if (rcfg->filter_specific_rho_cxh) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_hv = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		rcfg->filter_Doppler_spectralwidth_hh_ms = 1;
		
		todolist->calc_Doppler_spectrum_ShhShvc = 1;
		rcfg->filter_Doppler_spectrum_dBZ_hh = 1;
		rcfg->filter_Doppler_spectrum_dBZ_hv = 1;
		
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
		todolist->calc_Doppler_spectrum = 1;
	}
	if (rcfg->filter_specific_rho_cxv) {
		todolist->calc_eta_hh = 1;
		todolist->calc_eta_vh = 1;
		todolist->calc_eta_vv = 1;
		rcfg->filter_Doppler_velocity_hh_ms = 1;
		rcfg->filter_Doppler_spectralwidth_hh_ms = 1;
		
		todolist->calc_Doppler_spectrum_SvvSvhc = 1;
		rcfg->filter_Doppler_spectrum_dBZ_vh = 1;
		rcfg->filter_Doppler_spectrum_dBZ_vv = 1;
		
		todolist->calc_Doppler_mean_velocity = 1;
		todolist->calc_Doppler_spectrum_widths = 1;
		todolist->calc_Doppler_spectrum = 1;
	}

	if (todolist->calc_eta_i_hh) todolist->calc_eta_hh = 1;
	if (todolist->calc_eta_i_hv) todolist->calc_eta_hv = 1;
	if (todolist->calc_eta_i_vh) todolist->calc_eta_vh = 1;
	if (todolist->calc_eta_i_vv) todolist->calc_eta_vv = 1;
	if (todolist->calc_spectrum_eta_i_hh) todolist->calc_eta_hh = 1;
	if (todolist->calc_spectrum_eta_i_hv) todolist->calc_eta_hv = 1;
	if (todolist->calc_spectrum_eta_i_vh) todolist->calc_eta_vh = 1;
	if (todolist->calc_spectrum_eta_i_vv) todolist->calc_eta_vv = 1;	
}

void radarfilter_free_todolist(t_radarfilter_todolist **ptodolist)
{
	t_radarfilter_todolist *todolist = *ptodolist;
	if (todolist != NULL) {
		free(todolist);
		todolist = NULL;
	}
}


void radarfilter_initialize_radarmeasurement(int n_measurements, t_radarmeasurement ***pradarmeasurement)
{
	t_radarmeasurement **radarmeasurement = calloc(n_measurements, sizeof(t_radarmeasurement*));
	int i;	
	for (i=0; i < n_measurements; i++) {
		radarmeasurement[i] = calloc(1, sizeof(t_radarmeasurement));
		
		coordinates_initialize_coor(&(radarmeasurement[i]->center_coor));

		radarmeasurement[i]->spectrum_lbound = NULL;
		radarmeasurement[i]->spectrum_ubound = NULL;
		radarmeasurement[i]->spectrum_center = NULL;

		radarmeasurement[i]->Doppler_spectrum_dBZ_hh = NULL;
		radarmeasurement[i]->Doppler_spectrum_dBZ_hh_err = NULL;
		radarmeasurement[i]->Doppler_spectrum_dBZ_hv = NULL;
		radarmeasurement[i]->Doppler_spectrum_dBZ_hv_err = NULL;
		radarmeasurement[i]->Doppler_spectrum_dBZ_vh = NULL;
		radarmeasurement[i]->Doppler_spectrum_dBZ_vh_err = NULL;
		radarmeasurement[i]->Doppler_spectrum_dBZ_vv = NULL;
		radarmeasurement[i]->Doppler_spectrum_dBZ_vv_err = NULL;
		radarmeasurement[i]->specific_dBZdr = NULL;
		radarmeasurement[i]->specific_dBZdr_err = NULL;
		radarmeasurement[i]->specific_dBLdr = NULL;
		radarmeasurement[i]->specific_dBLdr_err = NULL;
		
		radarmeasurement[i]->specific_rho_co = NULL;
		radarmeasurement[i]->specific_rho_co_err = NULL;
		radarmeasurement[i]->specific_rho_cxh = NULL;
		radarmeasurement[i]->specific_rho_cxh_err = NULL;
		radarmeasurement[i]->specific_rho_cxv = NULL;
		radarmeasurement[i]->specific_rho_cxv_err = NULL;
		
		radarmeasurement[i]->n_psd			= 0;
		radarmeasurement[i]->n_diameters	= NULL;
		radarmeasurement[i]->eta_i_hh = NULL;
		radarmeasurement[i]->eta_i_hv = NULL;
		radarmeasurement[i]->eta_i_vh = NULL;
		radarmeasurement[i]->eta_i_vv = NULL;
		radarmeasurement[i]->spectrum_eta_i_hh = NULL;
		radarmeasurement[i]->spectrum_eta_i_hv = NULL;
		radarmeasurement[i]->spectrum_eta_i_vh = NULL;
		radarmeasurement[i]->spectrum_eta_i_vv = NULL;
		
		radarmeasurement[i]->der_edr13_Doppler_spectrum_dBZ_hh = NULL;
		radarmeasurement[i]->der_edr13_Doppler_spectrum_dBZ_hv = NULL;
		radarmeasurement[i]->der_edr13_Doppler_spectrum_dBZ_vh = NULL;
		radarmeasurement[i]->der_edr13_Doppler_spectrum_dBZ_vv = NULL;
		
		//analysis
		radarmeasurement[i]->analysis = calloc(1, sizeof(t_radarmeasurement_analysis));
		radarmeasurement[i]->analysis->discrete_D_equiv_mm 		= NULL;
		radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_mean 		= NULL;
		radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_variance 	= NULL;
	}
	
	*pradarmeasurement = radarmeasurement;
}	 

void radarfilter_prepare_model_radarmeasurement(int n_measurements, t_radarmeasurement ***pdst, t_radarmeasurement **src)
{
	t_radarmeasurement **dst;
	int i,j;	
		
	radarfilter_initialize_radarmeasurement(n_measurements, pdst);
	dst = *pdst;
	
	for (i=0; i < n_measurements; i++) {
		dst[i]->azel_r1_m = src[i]->azel_r1_m;
		dst[i]->azel_r2_m = src[i]->azel_r2_m;
		dst[i]->azel_alpha_rad = src[i]->azel_alpha_rad;
		dst[i]->azel_gamma_rad = src[i]->azel_gamma_rad;
		dst[i]->beam_FWHM0_rad = src[i]->beam_FWHM0_rad;
		dst[i]->beam_FWHM1_rad = src[i]->beam_FWHM1_rad;
		dst[i]->dt = src[i]->dt;
		dst[i]->center_coor->enu_radar_location_xyzt[0] = src[i]->center_coor->enu_radar_location_xyzt[0];
		dst[i]->center_coor->enu_radar_location_xyzt[1] = src[i]->center_coor->enu_radar_location_xyzt[1];
		dst[i]->center_coor->enu_radar_location_xyzt[2] = src[i]->center_coor->enu_radar_location_xyzt[2];
		dst[i]->center_coor->enu_radar_location_xyzt[3] = src[i]->center_coor->enu_radar_location_xyzt[3];
		
		dst[i]->n_spectrum = src[i]->n_spectrum;
		dst[i]->spectrum_lbound = malloc(dst[i]->n_spectrum * sizeof(double));
		for (j=0; j < src[i]->n_spectrum; j++) 
			dst[i]->spectrum_lbound[j] = src[i]->spectrum_lbound[j];
		dst[i]->spectrum_ubound = malloc(dst[i]->n_spectrum * sizeof(double));
		for (j=0; j < src[i]->n_spectrum; j++) 
			dst[i]->spectrum_ubound[j] = src[i]->spectrum_ubound[j];
		dst[i]->spectrum_center = malloc(dst[i]->n_spectrum * sizeof(double));
		for (j=0; j < src[i]->n_spectrum; j++) 
			dst[i]->spectrum_center[j] = src[i]->spectrum_center[j];
	}
}

void radarfilter_free_radarmeasurement(int n_measurements, t_radarmeasurement ***pradarmeasurement)
{
	t_radarmeasurement **radarmeasurement = *pradarmeasurement;
	int i;
	int i_psd;
	int i_int;
	int i_par;
	
	for (i=0; i < n_measurements; i++) {
		if (radarmeasurement[i] != NULL) {
			coordinates_free_coor(&radarmeasurement[i]->center_coor);
			util_safe_free(&(radarmeasurement[i]->spectrum_lbound));
			util_safe_free(&(radarmeasurement[i]->spectrum_ubound));
			util_safe_free(&(radarmeasurement[i]->spectrum_center));
			util_safe_free(&(radarmeasurement[i]->Doppler_spectrum_dBZ_hh));
			util_safe_free(&(radarmeasurement[i]->Doppler_spectrum_dBZ_hh_err));
			util_safe_free(&(radarmeasurement[i]->Doppler_spectrum_dBZ_hv));
			util_safe_free(&(radarmeasurement[i]->Doppler_spectrum_dBZ_hv_err));
			util_safe_free(&(radarmeasurement[i]->Doppler_spectrum_dBZ_vh));
			util_safe_free(&(radarmeasurement[i]->Doppler_spectrum_dBZ_vh_err));
			util_safe_free(&(radarmeasurement[i]->Doppler_spectrum_dBZ_vv));
			util_safe_free(&(radarmeasurement[i]->Doppler_spectrum_dBZ_vv_err));
			util_safe_free(&(radarmeasurement[i]->specific_dBZdr));
			util_safe_free(&(radarmeasurement[i]->specific_dBZdr_err));
			util_safe_free(&(radarmeasurement[i]->specific_dBLdr));
			util_safe_free(&(radarmeasurement[i]->specific_dBLdr_err));
			util_safe_free(&(radarmeasurement[i]->specific_rho_co));
			util_safe_free(&(radarmeasurement[i]->specific_rho_co_err));
			util_safe_free(&(radarmeasurement[i]->specific_rho_cxh));
			util_safe_free(&(radarmeasurement[i]->specific_rho_cxh_err));
			util_safe_free(&(radarmeasurement[i]->specific_rho_cxv));
			util_safe_free(&(radarmeasurement[i]->specific_rho_cxv_err));
		
			for ( i_psd = 0; i_psd < radarmeasurement[i]->n_psd; i_psd++ ) {
				for ( i_par = 0; i_par < radarmeasurement[i]->n_diameters[i_psd]; i_par++ ) {
					if (radarmeasurement[i]->spectrum_eta_i_hh != NULL) util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_hh[i_psd][i_par]));
					if (radarmeasurement[i]->spectrum_eta_i_hv != NULL) util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_hv[i_psd][i_par]));
					if (radarmeasurement[i]->spectrum_eta_i_vh != NULL) util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_vh[i_psd][i_par]));
					if (radarmeasurement[i]->spectrum_eta_i_vv != NULL) util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_vv[i_psd][i_par]));
				}
				if (radarmeasurement[i]->spectrum_eta_i_hh != NULL) util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_hh[i_psd]));
				if (radarmeasurement[i]->spectrum_eta_i_hv != NULL) util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_hv[i_psd]));
				if (radarmeasurement[i]->spectrum_eta_i_vh != NULL) util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_vh[i_psd]));
				if (radarmeasurement[i]->spectrum_eta_i_vv != NULL) util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_vv[i_psd]));
				if (radarmeasurement[i]->eta_i_hh != NULL) util_safe_free(&(radarmeasurement[i]->eta_i_hh[i_psd]));
				if (radarmeasurement[i]->eta_i_hv != NULL) util_safe_free(&(radarmeasurement[i]->eta_i_hv[i_psd]));
				if (radarmeasurement[i]->eta_i_vh != NULL) util_safe_free(&(radarmeasurement[i]->eta_i_vh[i_psd]));
				if (radarmeasurement[i]->eta_i_vv != NULL) util_safe_free(&(radarmeasurement[i]->eta_i_vv[i_psd]));
			}
			util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_hh));
			util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_hv));
			util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_vh));
			util_safe_free(&(radarmeasurement[i]->spectrum_eta_i_vv));
			util_safe_free(&(radarmeasurement[i]->eta_i_hh));
			util_safe_free(&(radarmeasurement[i]->eta_i_hv));
			util_safe_free(&(radarmeasurement[i]->eta_i_vh));
			util_safe_free(&(radarmeasurement[i]->eta_i_vv));

			//analysis
			if (radarmeasurement[i]->analysis != NULL) {
				for ( i_psd = 0; i_psd < radarmeasurement[i]->n_psd; i_psd++ ) {
					if (radarmeasurement[i]->analysis->discrete_D_equiv_mm != NULL) util_safe_free(&(radarmeasurement[i]->analysis->discrete_D_equiv_mm[i_psd]));
					if (radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_mean != NULL) util_safe_free(&(radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_mean[i_psd]));
					if (radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_variance != NULL) util_safe_free(&(radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_variance[i_psd]));
				}
				util_safe_free(&(radarmeasurement[i]->analysis->discrete_D_equiv_mm));
				util_safe_free(&(radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_mean));
				util_safe_free(&(radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_variance));
				util_safe_free(&(radarmeasurement[i]->analysis));
			} 
			util_safe_free(&(radarmeasurement[i]->n_diameters));
			util_safe_free(&(radarmeasurement[i]));
		}
	}
	util_safe_free(&radarmeasurement);
}


void radarfilter_read_measurements(int *pn_measurements, t_radarmeasurement ***pradarmeasurement, char measurements_filename[8192])
{
	t_radarfilter_readout_widget *rwg = malloc(sizeof(t_radarfilter_readout_widget));

	char dummy[8192];
	int i;
	
	*pradarmeasurement = NULL;
	
	rwg->fp = fopen(measurements_filename,"r"); // read mode
	printf("Opening file: %s\n",measurements_filename); fflush(stdout);
	if( rwg->fp == NULL )
	{
		perror("Error while opening the file.\n");
		exit(EXIT_FAILURE);
	}

	while (1) {
		fgetpos(rwg->fp, &rwg->pos_thisline);								//store position of where this line starts
		if (fgets (rwg->line, sizeof(rwg->line), rwg->fp) == NULL) {break;} 	//read full line
		fgetpos(rwg->fp, &rwg->pos_nextline);						//store position of where next line starts

		//check if there is data here
		if ((rwg->line[0] == '!') & (rwg->line[1] == '!')) {			
			fsetpos(rwg->fp, &rwg->pos_thisline);
			fscanf(rwg->fp, "%s %s %i", dummy, rwg->identifier, &rwg->ndim); //puts !! in identifier


			//obtain dimensions
			for ( i = 0; i < rwg->ndim; i++ ) {
				fscanf(rwg->fp, "%i", rwg->dim + i);
			}

			//debug
			//printf("line = %s\n", rwg->line); fflush(stdout);
			//printf("rwg->ndim = %i\n\n\n", rwg->ndim); fflush(stdout);
			//printf("rwg->dim[0] = %i\n\n\n", rwg->dim[0]); fflush(stdout);

			//make sure measurments are initialized
			if (*pradarmeasurement == NULL) {
				*pn_measurements = rwg->dim[0];				
				radarfilter_initialize_radarmeasurement(*pn_measurements, pradarmeasurement);
			}

			if (*pradarmeasurement != NULL) {
				//read out measurements
				for (i=0; i < *pn_measurements; i++) {
					radarfilter_readout(rwg, (*pradarmeasurement)[i]);
				}
			}
		}
	}
}




void radarfilter_readout(t_radarfilter_readout_widget *rwg, t_radarmeasurement *radarmeasurement)
{
	int i;

	if (strcmp(rwg->identifier,"enu_radar_location_x") == 0 ) fscanf(rwg->fp, "%lf", radarmeasurement->center_coor->enu_radar_location_xyzt);
	if (strcmp(rwg->identifier,"enu_radar_location_y") == 0 ) fscanf(rwg->fp, "%lf", radarmeasurement->center_coor->enu_radar_location_xyzt + 1);
	if (strcmp(rwg->identifier,"enu_radar_location_z") == 0 ) fscanf(rwg->fp, "%lf", radarmeasurement->center_coor->enu_radar_location_xyzt + 2);
	if (strcmp(rwg->identifier,"enu_radar_time_t") == 0 ) fscanf(rwg->fp, "%lf", radarmeasurement->center_coor->enu_radar_location_xyzt + 3);
	if (strcmp(rwg->identifier,"azel_r1_m") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->azel_r1_m);
	if (strcmp(rwg->identifier,"azel_r2_m") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->azel_r2_m);
	if (strcmp(rwg->identifier,"azel_alpha_rad") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->azel_alpha_rad);
	if (strcmp(rwg->identifier,"azel_gamma_rad") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->azel_gamma_rad);
	if (strcmp(rwg->identifier,"beam_FWHM0_rad") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->beam_FWHM0_rad);
	if (strcmp(rwg->identifier,"beam_FWHM1_rad") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->beam_FWHM1_rad);
	if (strcmp(rwg->identifier,"dt") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dt);
	if (strcmp(rwg->identifier,"dBZ_hh") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZ_hh);
	if (strcmp(rwg->identifier,"dBZ_hh_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZ_hh_err);
	if (strcmp(rwg->identifier,"dBZ_hv") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZ_hv);
	if (strcmp(rwg->identifier,"dBZ_hv_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZ_hv_err);
	if (strcmp(rwg->identifier,"dBZ_vh") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZ_vh);
	if (strcmp(rwg->identifier,"dBZ_vh_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZ_vh_err);
	if (strcmp(rwg->identifier,"dBZ_vv") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZ_vv);
	if (strcmp(rwg->identifier,"dBZ_vv_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZ_vv_err);
	if (strcmp(rwg->identifier,"dBZdr") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZdr);
	if (strcmp(rwg->identifier,"dBZdr_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBZdr_err);
	if (strcmp(rwg->identifier,"dBLdr") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBLdr);
	if (strcmp(rwg->identifier,"dBLdr_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->dBLdr_err);
	if (strcmp(rwg->identifier,"rho_co") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->rho_co);
	if (strcmp(rwg->identifier,"rho_co_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->rho_co_err);
	if (strcmp(rwg->identifier,"rho_cxh") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->rho_cxh);
	if (strcmp(rwg->identifier,"rho_cxh_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->rho_cxh_err);
	if (strcmp(rwg->identifier,"rho_cxv") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->rho_cxv);
	if (strcmp(rwg->identifier,"rho_cxv_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->rho_cxv_err);
	if (strcmp(rwg->identifier,"KDP") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->KDP);
	if (strcmp(rwg->identifier,"KDP_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->KDP_err);
	if (strcmp(rwg->identifier,"Doppler_velocity_hh_ms") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_velocity_hh_ms);
	if (strcmp(rwg->identifier,"Doppler_velocity_hh_ms_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_velocity_hh_ms_err);
	if (strcmp(rwg->identifier,"Doppler_velocity_hv_ms") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_velocity_hv_ms);
	if (strcmp(rwg->identifier,"Doppler_velocity_hv_ms_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_velocity_hv_ms_err);
	if (strcmp(rwg->identifier,"Doppler_velocity_vh_ms") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_velocity_vh_ms);
	if (strcmp(rwg->identifier,"Doppler_velocity_vh_ms_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_velocity_vh_ms_err);
	if (strcmp(rwg->identifier,"Doppler_velocity_vv_ms") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_velocity_vv_ms);
	if (strcmp(rwg->identifier,"Doppler_velocity_vv_ms_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_velocity_vv_ms_err);
	if (strcmp(rwg->identifier,"Doppler_spectral_width_hh_ms") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_spectral_width_hh_ms);
	if (strcmp(rwg->identifier,"Doppler_spectral_width_hh_ms_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_spectral_width_hh_ms_err);
	if (strcmp(rwg->identifier,"Doppler_spectral_width_hv_ms") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_spectral_width_hv_ms);
	if (strcmp(rwg->identifier,"Doppler_spectral_width_hv_ms_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_spectral_width_hv_ms_err);
	if (strcmp(rwg->identifier,"Doppler_spectral_width_vh_ms") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_spectral_width_vh_ms);
	if (strcmp(rwg->identifier,"Doppler_spectral_width_vh_ms_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_spectral_width_vh_ms_err);
	if (strcmp(rwg->identifier,"Doppler_spectral_width_vv_ms") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_spectral_width_vv_ms);
	if (strcmp(rwg->identifier,"Doppler_spectral_width_vv_ms_err") == 0 ) fscanf(rwg->fp, "%lf", &radarmeasurement->Doppler_spectral_width_vv_ms_err);

	if (strcmp(rwg->identifier,"spectrum_velocity_lbound") == 0 ) {
		//allocate and read out
		radarmeasurement->n_spectrum = rwg->dim[1];
		radarmeasurement->spectrum_lbound = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->spectrum_lbound + i);
		}
	}
	if (strcmp(rwg->identifier,"spectrum_velocity_center") == 0 ) {
		//allocate and read out
		radarmeasurement->n_spectrum = rwg->dim[1];
		radarmeasurement->spectrum_center = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->spectrum_center + i);
		}
	}
	if (strcmp(rwg->identifier,"spectrum_velocity_ubound") == 0 ) {
		//allocate and read out
		radarmeasurement->n_spectrum = rwg->dim[1];
		radarmeasurement->spectrum_ubound = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->spectrum_ubound + i);
		}
	}
	
	if (strcmp(rwg->identifier,"Doppler_spectrum_dBZ_hh") == 0 ) {
		//allocate and read out
		radarmeasurement->Doppler_spectrum_dBZ_hh = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->Doppler_spectrum_dBZ_hh + i);
		}
	}
	if (strcmp(rwg->identifier,"Doppler_spectrum_dBZ_hh_err") == 0 ) {
		//allocate and read out
		radarmeasurement->Doppler_spectrum_dBZ_hh_err = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->Doppler_spectrum_dBZ_hh_err + i);
		}
	}
	if (strcmp(rwg->identifier,"Doppler_spectrum_dBZ_hv") == 0 ) {
		//allocate and read out
		radarmeasurement->Doppler_spectrum_dBZ_hv = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->Doppler_spectrum_dBZ_hv + i);
		}
	}
	if (strcmp(rwg->identifier,"Doppler_spectrum_dBZ_hv_err") == 0 ) {
		//allocate and read out
		radarmeasurement->Doppler_spectrum_dBZ_hv_err = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->Doppler_spectrum_dBZ_hv_err + i);
		}
	}
	if (strcmp(rwg->identifier,"Doppler_spectrum_dBZ_vh") == 0 ) {
		//allocate and read out
		radarmeasurement->Doppler_spectrum_dBZ_vh = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->Doppler_spectrum_dBZ_vh + i);
		}
	}
	if (strcmp(rwg->identifier,"Doppler_spectrum_dBZ_vh_err") == 0 ) {
		//allocate and read out
		radarmeasurement->Doppler_spectrum_dBZ_vh_err = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->Doppler_spectrum_dBZ_vh_err + i);
		}
	}
	if (strcmp(rwg->identifier,"Doppler_spectrum_dBZ_vv") == 0 ) {
		//allocate and read out
		radarmeasurement->Doppler_spectrum_dBZ_vv = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->Doppler_spectrum_dBZ_vv + i);
		}
	}
	if (strcmp(rwg->identifier,"Doppler_spectrum_dBZ_vv_err") == 0 ) {
		//allocate and read out
		radarmeasurement->Doppler_spectrum_dBZ_vv_err = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->Doppler_spectrum_dBZ_vv_err + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_dBZdr") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_dBZdr = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_dBZdr + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_dBZdr_err") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_dBZdr_err = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_dBZdr_err + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_dBLdr") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_dBLdr = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_dBLdr + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_dBLdr_err") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_dBLdr_err = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_dBLdr_err + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_rho_co") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_rho_co = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_rho_co + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_rho_co_err") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_rho_co_err = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_rho_co_err + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_rho_cxh") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_rho_cxh = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_rho_cxh + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_rho_cxh_err") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_rho_cxh_err = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_rho_cxh_err + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_rho_cxv") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_rho_cxv = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_rho_cxv + i);
		}
	}
	if (strcmp(rwg->identifier,"specific_rho_cxv_err") == 0 ) {
		//allocate and read out
		radarmeasurement->specific_rho_cxv_err = malloc(rwg->dim[1] * sizeof(double));
		for ( i = 0; i < rwg->dim[1]; i++ ) {
			fscanf(rwg->fp, "%lf", radarmeasurement->specific_rho_cxv_err + i);
		}
	}
}



void radarfilter_write_measurements(t_zephyros_config *cfg, int i_mode, int n_measurements, t_radarmeasurement **radarmeasurement, FILE *fp)
{
	int i, j;
	int i_psd, i_par;
	char myname[100];

	t_zephyros_radarfilter	*rcfg;

	//point radar filter configuration in the right way
	if (i_mode == 0) rcfg = cfg->simulation->radarfilter;
	if (i_mode == 1) rcfg = cfg->retrieval->radarfilter;
	
	fprintf(fp, "\n\n");
	
	//enu_radar_location_x
	fprintf(fp, "!! %-30s %-15i %-15i", "enu_radar_location_x", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->center_coor->enu_radar_location_xyzt[0]);
	fprintf(fp, " \n");
	//enu_radar_location_y
	fprintf(fp, "!! %-30s %-15i %-15i", "enu_radar_location_y", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->center_coor->enu_radar_location_xyzt[1]);
	fprintf(fp, " \n");
	//enu_radar_location_z
	fprintf(fp, "!! %-30s %-15i %-15i", "enu_radar_location_z", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->center_coor->enu_radar_location_xyzt[2]);
	fprintf(fp, " \n");
	//enu_radar_time_t
	fprintf(fp, "!! %-30s %-15i %-15i", "enu_radar_time_t", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->center_coor->enu_radar_location_xyzt[3]);
	fprintf(fp, " \n");
	
	//azel_r1_m
	fprintf(fp, "!! %-30s %-15i %-15i", "azel_r1_m", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->azel_r1_m);
	fprintf(fp, " \n");
	//azel_r2_m
	fprintf(fp, "!! %-30s %-15i %-15i", "azel_r2_m", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->azel_r2_m);
	fprintf(fp, " \n");
	//azel_alpha_rad
	fprintf(fp, "!! %-30s %-15i %-15i", "azel_alpha_rad", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->azel_alpha_rad);
	fprintf(fp, " \n");
	//azel_gamma_rad
	fprintf(fp, "!! %-30s %-15i %-15i", "azel_gamma_rad", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->azel_gamma_rad);
	fprintf(fp, " \n");
	//beam_FWHM0_rad
	fprintf(fp, "!! %-30s %-15i %-15i", "beam_FWHM0_rad", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->beam_FWHM0_rad);
	fprintf(fp, " \n");
	//beam_FWHM1_rad
	fprintf(fp, "!! %-30s %-15i %-15i", "beam_FWHM1_rad", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->beam_FWHM1_rad);
	fprintf(fp, " \n");
	//dt
	fprintf(fp, "!! %-30s %-15i %-15i", "dt", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dt);
	fprintf(fp, " \n");
	
	//center x
	fprintf(fp, "!! %-30s %-15i %-15i", "center_x", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->center_coor->enu_xyzt[0]);
	fprintf(fp, " \n");

	//center y
	fprintf(fp, "!! %-30s %-15i %-15i", "center_y", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->center_coor->enu_xyzt[1]);
	fprintf(fp, " \n");

	//center z
	fprintf(fp, "!! %-30s %-15i %-15i", "center_z", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->center_coor->enu_xyzt[2]);
	fprintf(fp, " \n");

	//center t
	fprintf(fp, "!! %-30s %-15i %-15i", "center_t", 1, n_measurements);
	for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->center_coor->enu_xyzt[3]);
	fprintf(fp, " \n");

	if (rcfg->filter_dBZ_hh) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZ_hh", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZ_hh);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_dBZ_hh) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZ_hh_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZ_hh_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_dBZ_hv) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZ_hv", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZ_hv);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_dBZ_hv) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZ_hv_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZ_hv_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_dBZ_vh) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZ_vh", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZ_vh);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_dBZ_vh) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZ_vh_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZ_vh_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_dBZ_vv) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZ_vv", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZ_vv);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_dBZ_vv) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZ_vv_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZ_vv_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_dBZdr) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZdr", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZdr);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_dBZdr) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBZdr_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBZdr_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_dBLdr) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBLdr", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBLdr);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_dBLdr) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "dBLdr_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->dBLdr_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_rho_co) {
		fprintf(fp, "!! %-30s %-15i %-15i", "rho_co", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->rho_co);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_rho_co) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "rho_co_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->rho_co_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_rho_cxh) {
		fprintf(fp, "!! %-30s %-15i %-15i", "rho_cxh", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->rho_cxh);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_rho_cxh) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "rho_cxh_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->rho_cxh_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_rho_cxv) {
		fprintf(fp, "!! %-30s %-15i %-15i", "rho_cxv", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->rho_cxv);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_rho_cxv) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "rho_cxv_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->rho_cxv_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_velocity_hh_ms) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_velocity_hh_ms", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_velocity_hh_ms);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_velocity_hh_ms) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_velocity_hh_ms_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_velocity_hh_ms_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_velocity_hv_ms) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_velocity_hv_ms", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_velocity_hv_ms);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_velocity_hv_ms) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_velocity_hv_ms_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_velocity_hv_ms_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_velocity_vh_ms) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_velocity_vh_ms", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_velocity_hv_ms);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_velocity_vh_ms) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_velocity_vh_ms_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_velocity_hv_ms_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_velocity_vv_ms) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_velocity_vv_ms", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_velocity_vv_ms);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_velocity_vv_ms) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_velocity_vv_ms_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_velocity_vv_ms_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_spectralwidth_hh_ms) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_spectral_width_hh_ms", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectral_width_hh_ms);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_spectralwidth_hh_ms) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_spectral_width_hh_ms_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectral_width_hh_ms_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_spectralwidth_hv_ms) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_spectral_width_hv_ms", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectral_width_hv_ms);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_spectralwidth_hv_ms) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_spectral_width_hv_ms_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectral_width_hv_ms_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_spectralwidth_vh_ms) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_spectral_width_vh_ms", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectral_width_hv_ms);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_spectralwidth_vh_ms) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_spectral_width_vh_ms_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectral_width_hv_ms_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_spectralwidth_vv_ms) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_spectral_width_vv_ms", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectral_width_vv_ms);
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_spectralwidth_vv_ms) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "Doppler_spectral_width_vv_ms_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectral_width_vv_ms_err);
		fprintf(fp, " \n");
	}
	if (rcfg->filter_KDP) {
		fprintf(fp, "!! %-30s %-15i %-15i", "KDP", 1, n_measurements);
		for (i=0; i < n_measurements; i++) {
			fprintf(fp, " %-15.3e", radarmeasurement[i]->KDP);
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_KDP) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i", "KDP_err", 1, n_measurements);
		for (i=0; i < n_measurements; i++) {
			fprintf(fp, " %-15.3e", radarmeasurement[i]->KDP_err);
		}
		fprintf(fp, " \n");
	}
	
	//*****
	//spectra
	if ( 	(rcfg->filter_Doppler_spectrum_dBZ_hh) |
			(rcfg->filter_Doppler_spectrum_dBZ_hv) |
			(rcfg->filter_Doppler_spectrum_dBZ_vh) |
			(rcfg->filter_Doppler_spectrum_dBZ_vv) |
			(rcfg->filter_specific_dBZdr) |
			(rcfg->filter_specific_dBLdr) |
			(rcfg->filter_specific_rho_co) |
			(rcfg->filter_specific_rho_cxh) |
			(rcfg->filter_specific_rho_cxv)	) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "spectrum_velocity_lbound", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->spectrum_lbound[j]);
			}
		}
		fprintf(fp, " \n");
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "spectrum_velocity_center", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->spectrum_center[j]);
			}
		}
		fprintf(fp, " \n");
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "spectrum_velocity_ubound", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->spectrum_ubound[j]);
			}
		}
		fprintf(fp, " \n");
	}

	if (rcfg->filter_Doppler_spectrum_dBZ_hh) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "Doppler_spectrum_dBZ_hh", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectrum_dBZ_hh[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_spectrum_dBZ_hh) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "Doppler_spectrum_dBZ_hh_err", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectrum_dBZ_hh_err[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_spectrum_dBZ_hv) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "Doppler_spectrum_dBZ_hv", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectrum_dBZ_hv[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_spectrum_dBZ_hv) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "Doppler_spectrum_dBZ_hv_err", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectrum_dBZ_hv_err[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_spectrum_dBZ_vh) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "Doppler_spectrum_dBZ_vh", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectrum_dBZ_vh[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_spectrum_dBZ_vh) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "Doppler_spectrum_dBZ_vh_err", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectrum_dBZ_vh_err[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if (rcfg->filter_Doppler_spectrum_dBZ_vv) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "Doppler_spectrum_dBZ_vv", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectrum_dBZ_vv[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_Doppler_spectrum_dBZ_vv) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "Doppler_spectrum_dBZ_vv_err", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->Doppler_spectrum_dBZ_vv_err[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if (rcfg->filter_specific_dBZdr) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_dBZdr", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_dBZdr[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_specific_dBZdr) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_dBZdr_err", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_dBZdr_err[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if (rcfg->filter_specific_dBLdr) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_dBLdr", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_dBLdr[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_specific_dBLdr) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_dBLdr_err", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_dBLdr_err[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if (rcfg->filter_specific_rho_co) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_rho_co", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_rho_co[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_specific_rho_co) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_rho_co_err", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_rho_co_err[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if (rcfg->filter_specific_rho_cxh) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_rho_cxh", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_rho_cxh[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_specific_rho_cxh) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_rho_cxh_err", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_rho_cxh_err[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if (rcfg->filter_specific_rho_cxv) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_rho_cxv", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_rho_cxv[j]);
			}
		}
		fprintf(fp, " \n");
	}
	if ((rcfg->filter_specific_rho_cxv) & (rcfg->filter_errors)) {
		fprintf(fp, "!! %-30s %-15i %-15i %-15i", "specific_rho_cxv_err", 2, n_measurements, radarmeasurement[0]->n_spectrum);
		for (i=0; i < n_measurements; i++) {
			for (j=0; j < radarmeasurement[i]->n_spectrum; j++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->specific_rho_cxv_err[j]);
			}
		}
		fprintf(fp, " \n");
	}
	
	fflush(stdout);
}

void radarfilter_write_measurements_detailed_analysis(t_zephyros_config *cfg, int i_mode, int n_measurements, t_radarmeasurement **radarmeasurement, FILE *fp)
{
	int i, j;
	int i_psd, i_par;
	char myname[100];

	for (i=0; i < n_measurements; i++) {
		fprintf(fp, "\n");
		fprintf(fp, "\n");
		fprintf(fp, "Detailed analysis for measurement %i \n", i);	fflush(stdout);
			
		for ( i_psd = 0; i_psd < radarmeasurement[i]->n_psd; i_psd++ ) {
			fprintf(fp, "\n");
			
			sprintf(myname, "analysis_psd%i_meas%i_discrete_D_equiv_mm", i_psd, i);
			fprintf(fp, "!! %-80s %-15i %-15i", myname, 1, radarmeasurement[i]->n_diameters[i_psd]);
			for (i_par=0; i_par < radarmeasurement[i]->n_diameters[i_psd]; i_par++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->analysis->discrete_D_equiv_mm[i_psd][i_par]);
			}
			fprintf(fp, " \n"); fflush(stdout);
		
			sprintf(myname, "analysis_psd%i_meas%i_unweighted_canting_angle_wrt_vertical_mean", i_psd, i);
			fprintf(fp, "!! %-80s %-15i %-15i", myname, 1, radarmeasurement[i]->n_diameters[i_psd]);
			for (i_par=0; i_par < radarmeasurement[i]->n_diameters[i_psd]; i_par++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_mean[i_psd][i_par]);
			}
			fprintf(fp, " \n"); fflush(stdout);
			
			sprintf(myname, "analysis_psd%i_meas%i_unweighted_canting_angle_wrt_vertical_variance", i_psd, i);
			fprintf(fp, "!! %-80s %-15i %-15i", myname, 1, radarmeasurement[i]->n_diameters[i_psd]);
			for (i_par=0; i_par < radarmeasurement[i]->n_diameters[i_psd]; i_par++) {
				fprintf(fp, " %-15.3e", radarmeasurement[i]->analysis->unweighted_canting_angle_wrt_vertical_variance[i_psd][i_par]);
			}
			fprintf(fp, " \n");	fflush(stdout);
		}
	}
	
	fflush(stdout);
}
