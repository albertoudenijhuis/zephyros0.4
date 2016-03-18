/*
Description: 
	Coordinate transformations

Revision History:
	2014

Functions:
	Coordinate transformations.

Author:
	Albert Oude Nijhuis <albertoudenijhuis@gmail.com>

Institute:
	Delft University of Technology
	
Zephyros version:
	0.2

Project:
	EU FP7 program, the UFO project

Dissemination:
	Confidential, only for members of the UFO project. Potentially public in the future.

Acknowledgement and citation:
	Whenever this code is used for publication of scientific results,
	the code writer should be informed, acknowledged and referenced.

Note:
	If you have any suggestions for improvements or amendments, please inform the author of this code.

*/


/* System include files */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "coordinates.h"

#define wgs84_finv 		298.257223563
#define wgs84_b 		6356752.31425
#define earth_a 		6371009.0
#define wgs84_a_mul_esq 285.834445797
#define earth_a_e 		8494678.66667
#define wgs84_a 		6378137.0
#define wgs84_f 		0.00335281066475
#define wgs84_e 		0.00669437999014
#define wgs84_esq 		4.48147234524e-05
#define wgs84_esqsq 	2.00835943812e-09
#define wgs84_asq 		4.06806315908e+13

void coordinates_initialize_coor(t_zephyros_coordinates **pcoor)
{
	t_zephyros_coordinates *coor = malloc(sizeof(t_zephyros_coordinates));

	//arrays
	coor->ecef_xyzt 					= calloc(4, sizeof(double));
	coor->ecef_radar_location_xyzt 		= calloc(4, sizeof(double));
	coor->enu_xyzt 						= calloc(4, sizeof(double));
	coor->enu_radar_location_xyzt 		= calloc(4, sizeof(double));
	coor->radar_enu_dir 				= calloc(4, sizeof(double));
	coor->radar_pol_hor_dir				= calloc(4, sizeof(double));
	coor->radar_pol_ver_dir 			= calloc(4, sizeof(double));

	coor->wgs84_phirad 						= 0.;
	coor->wgs84_lambdarad 					= 0.;
	coor->wgs84_h			 				= 0.;
	coor->wgs84_radar_location_phirad 		= 0.;
	coor->wgs84_radar_location_lambdarad 	= 0.;
	coor->wgs84_radar_location_h 			= 0.;
	
	coor->ecef_xyzt[0] = 0.;
	coor->ecef_xyzt[1] = 0.;
	coor->ecef_xyzt[2] = 0.;
	coor->ecef_xyzt[3] = 0.;
	coor->ecef_radar_location_xyzt[0] = 0.;
	coor->ecef_radar_location_xyzt[1] = 0.;
	coor->ecef_radar_location_xyzt[2] = 0.;
	coor->ecef_radar_location_xyzt[3] = 0.;
	
	coor->ref_wgs84_phirad = 0.;
	coor->ref_wgs84_lambdarad = 0.;
	coor->ref_wgs84_h = 0.;
	coor->enu_xyzt[0] = 0.;
	coor->enu_xyzt[1] = 0.;
	coor->enu_xyzt[2] = 0.;
	coor->enu_xyzt[3] = 0.;
	coor->enu_radar_location_xyzt[0] = 0.;
	coor->enu_radar_location_xyzt[1] = 0.;
	coor->enu_radar_location_xyzt[2] = 0.;
	coor->enu_radar_location_xyzt[3] = 0.;

	coor->radar_effearthcorrection 		= 1;
	coor->radar_azel_alpha	 			= 0.;
	coor->radar_azel_gamma	 			= 0.;
	coor->radar_azel_gamma_cor	 		= 0.;
	coor->radar_range	 				= 0.;
	
	coor->radar_enu_dir[0] = 0.;
	coor->radar_enu_dir[1] = 0.;
	coor->radar_enu_dir[2] = 0.;
	coor->radar_enu_dir[3] = 0.;
	
	coor->radar_beam_theta	 		= 0.;
	coor->radar_beam_phi	 		= 0.;

	coor->radar_pol_hor_dir[0] = 0.;
	coor->radar_pol_hor_dir[1] = 0.;
	coor->radar_pol_hor_dir[2] = 0.;
	coor->radar_pol_hor_dir[3] = 0.;
	coor->radar_pol_ver_dir[0] = 0.;
	coor->radar_pol_ver_dir[1] = 0.;
	coor->radar_pol_ver_dir[2] = 0.;
	coor->radar_pol_ver_dir[3] = 0.;

	strcpy(coor->name, "Untitled");
	coor->initialized = 1;
	
	*pcoor = coor;
}

void coordinates_assert_initialized(t_zephyros_coordinates *coor)
{
	if (coor == NULL) {
		printf("Coordinates was NULL pointers. Exiting\n");
		exit(0);
	}
	if (coor->initialized != 1) {
		printf("Coordinates were not initialized. Exiting.\n");
		exit(0);
	}
}

void coordinates_free_coor(t_zephyros_coordinates **pcoor)
{
	t_zephyros_coordinates *coor = *pcoor;
	
	if (coor != NULL) {
		coordinates_assert_initialized(coor);
	
		free(coor->ecef_xyzt);
		free(coor->ecef_radar_location_xyzt);
		free(coor->enu_xyzt);
		free(coor->enu_radar_location_xyzt);
		free(coor->radar_enu_dir);
		free(coor->radar_pol_hor_dir);
		free(coor->radar_pol_ver_dir);

		free(coor);
		*pcoor = NULL;
	}
}


void coordinates_geodetic2ecef(t_zephyros_coordinates *coor)
{              
	double N;
	
	coordinates_assert_initialized(coor);

	N = wgs84_a / sqrt(1. - wgs84_esq * pow(sin(coor->wgs84_phirad),2));
	coor->ecef_xyzt[0] = (N + coor->wgs84_h) * cos(coor->wgs84_phirad) * cos(coor->wgs84_lambdarad);
	coor->ecef_xyzt[1] = (N + coor->wgs84_h) * cos(coor->wgs84_phirad) * sin(coor->wgs84_lambdarad);
	coor->ecef_xyzt[2] = (N * (1. - wgs84_esq) + coor->wgs84_h) * sin(coor->wgs84_phirad);
}

void coordinates_ecef2geodetic(t_zephyros_coordinates *coor, int method)
{
	double 	psq;
	double 	p;
	int		nsteps = 3;
	double	kappa[nsteps+1];
	double	c[nsteps];
	int		i;
	
	double  zeta, rho, s, t, u, v, w;

	coordinates_assert_initialized(coor);
		
	coor->wgs84_lambdarad 	= atan2(coor->ecef_xyzt[1], coor->ecef_xyzt[0]);

	psq			= pow(coor->ecef_xyzt[0],2.) + pow(coor->ecef_xyzt[1],2.);
	p 			= sqrt(psq);
	kappa[0] 	= 1. / (1. - wgs84_esq);

	if (method == 1) {
		/*method 1: Bowring's method  */
		for ( i = 0; i < nsteps; i++ ) {
			c[i] 		= pow(psq + (1. - wgs84_esq) * pow(coor->ecef_xyzt[2],2.) * pow(kappa[i],2.) , 3./2.) / wgs84_a_mul_esq;
			kappa[i+1]	= (c[i] + (1. - wgs84_esq) * pow(coor->ecef_xyzt[2],2.) * pow(kappa[i],3.)) / (c[i] - psq);
		}
		
		coor->wgs84_h 		= (pow(kappa[nsteps], -1.) - pow(kappa[0], -1.)) * sqrt(psq + pow(coor->ecef_xyzt[2] * kappa[nsteps],2.)) / wgs84_esq;
		coor->wgs84_phirad 	= atan2(kappa[nsteps] * coor->ecef_xyzt[2] , p);
	} 

	if (method == 2) {
		/* Ferrari's solution */
		zeta 			= (1. - wgs84_esq) * pow(coor->ecef_xyzt[2],2.) / wgs84_asq;
		rho				= ((psq / wgs84_asq) + zeta - wgs84_esqsq) / 6.;
		s 				= wgs84_esqsq * zeta * psq / (4. * wgs84_asq);
		t 				= pow(pow(rho,3.) + s + sqrt(s * (s + 2. * pow(rho,3.))),1./3.);
		u 				= rho + t + (pow(rho,2.) / t);
		v				= sqrt(pow(u,2.) + wgs84_esqsq * zeta);
		w				= wgs84_esq * (u + v - zeta) / (2. * v);
		kappa[1]		= 1. + wgs84_esq * (sqrt(u + v + pow(w,2.)) + w) / (u + v);
		
		coor->wgs84_h 		= (pow(kappa[1],-1.) - pow(kappa[0],-1.)) * sqrt(psq + pow(coor->ecef_xyzt[2] * kappa[1], 2.)) / wgs84_esq;
		coor->wgs84_phirad 	= atan2(kappa[1] * coor->ecef_xyzt[2] , p);
	}
}

void coordinates_ecef2enu(t_zephyros_coordinates *coor)
{
	t_zephyros_coordinates *refcoor	= malloc(sizeof(t_zephyros_coordinates));

	double ecef_dx;
	double ecef_dy;
	double ecef_dz;

	coordinates_assert_initialized(coor);
		
	refcoor->wgs84_phirad 		= coor->ref_wgs84_phirad;
	refcoor->wgs84_lambdarad 	= coor->ref_wgs84_lambdarad;
	refcoor->wgs84_h 			= coor->ref_wgs84_h;
	coordinates_geodetic2ecef(refcoor);
	
	ecef_dx = coor->ecef_xyzt[0] - refcoor->ecef_xyzt[0];
	ecef_dy = coor->ecef_xyzt[1] - refcoor->ecef_xyzt[1];
	ecef_dz = coor->ecef_xyzt[2] - refcoor->ecef_xyzt[2];
	
	coor->enu_xyzt[0] = - 	sin(coor->ref_wgs84_lambdarad) * ecef_dx 
			+ 	cos(coor->ref_wgs84_lambdarad) * ecef_dy;
	coor->enu_xyzt[1] = - 	sin(coor->ref_wgs84_phirad) * cos(coor->ref_wgs84_lambdarad) * ecef_dx
			- 	sin(coor->ref_wgs84_phirad) * sin(coor->ref_wgs84_lambdarad) * ecef_dy
			+ 	cos(coor->ref_wgs84_phirad) * ecef_dz;
	coor->enu_xyzt[2] = 	cos(coor->ref_wgs84_phirad) * cos(coor->ref_wgs84_lambdarad) * ecef_dx
			+ 	cos(coor->ref_wgs84_phirad) * sin(coor->ref_wgs84_lambdarad) * ecef_dy
			+ 	sin(coor->ref_wgs84_phirad) * ecef_dz;

	free(refcoor);
}

void coordinates_enu2ecef(t_zephyros_coordinates *coor)
{
	t_zephyros_coordinates *refcoor	= malloc(sizeof(t_zephyros_coordinates));

	coordinates_assert_initialized(coor);
	
	refcoor->wgs84_phirad 			= coor->ref_wgs84_phirad;
	refcoor->ref_wgs84_lambdarad 	= coor->ref_wgs84_lambdarad;
	refcoor->ref_wgs84_h 			= coor->ref_wgs84_h;
	coordinates_geodetic2ecef(refcoor);
		
	coor->ecef_xyzt[0] = refcoor->ecef_xyzt[0] 
			- 	sin(coor->ref_wgs84_lambdarad) * coor->enu_xyzt[0] 
			-	sin(coor->ref_wgs84_phirad) * cos(coor->ref_wgs84_lambdarad) * coor->enu_xyzt[1]
			+	cos(coor->ref_wgs84_phirad) * cos(coor->ref_wgs84_lambdarad) * coor->enu_xyzt[2];
	coor->ecef_xyzt[1] = refcoor->ecef_xyzt[1]
			+	cos(coor->ref_wgs84_lambdarad) * coor->enu_xyzt[0]
			- 	sin(coor->ref_wgs84_phirad) * sin(coor->ref_wgs84_lambdarad) * coor->enu_xyzt[1]
			+ 	cos(coor->ref_wgs84_phirad) * sin(coor->ref_wgs84_lambdarad) * coor->enu_xyzt[2];
	coor->ecef_xyzt[2] = refcoor->ecef_xyzt[2]
			+	cos(coor->ref_wgs84_phirad) * coor->enu_xyzt[1]
			+ 	sin(coor->ref_wgs84_phirad) * coor->enu_xyzt[2];

	free(refcoor); 
}

void coordinates_radar_azel2enu(t_zephyros_coordinates *coor)
{
	coordinates_assert_initialized(coor);
	
	if (coor->radar_effearthcorrection == 1) {
		/* effective earth correction */
		coor->radar_azel_gamma_cor = coor->radar_azel_gamma + atan((coor->radar_range * cos(coor->radar_azel_gamma)) / (earth_a_e + coor->radar_range * sin(coor->radar_azel_gamma)));
	} else {
		coor->radar_azel_gamma_cor = coor->radar_azel_gamma;
	}
	
	coor->enu_xyzt[0] = coor->enu_radar_location_xyzt[0] + (coor->radar_range * cos(coor->radar_azel_gamma_cor) * sin(coor->radar_azel_alpha));
	coor->enu_xyzt[1] = coor->enu_radar_location_xyzt[1] + (coor->radar_range * cos(coor->radar_azel_gamma_cor) * cos(coor->radar_azel_alpha));
	coor->enu_xyzt[2] = coor->enu_radar_location_xyzt[2] + (coor->radar_range * sin(coor->radar_azel_gamma_cor));
	coor->enu_xyzt[3] = coor->enu_radar_location_xyzt[3];
}

void coordinates_radar_azelrangedir2enu(t_zephyros_coordinates *coor)
{
	double dgamma_dr;

	coordinates_assert_initialized(coor);
		
	if (coor->radar_effearthcorrection == 1) {
		/* effective earth correction */
		coor->radar_azel_gamma_cor = coor->radar_azel_gamma + atan((coor->radar_range * cos(coor->radar_azel_gamma)) /
												( earth_a_e + (coor->radar_range * sin(coor->radar_azel_gamma))  ));
		dgamma_dr = (earth_a_e * cos(coor->radar_azel_gamma)) / 
			(pow(earth_a_e + coor->radar_range * sin(coor->radar_azel_gamma), 2.) + pow(coor->radar_range * cos(coor->radar_azel_gamma), 2.));
		coor->radar_enu_dir[0] = cos(coor->radar_azel_gamma_cor) * sin(coor->radar_azel_alpha)	+ 
					(-1. * coor->radar_range * dgamma_dr * sin(coor->radar_azel_gamma_cor) * sin(coor->radar_azel_alpha));
		coor->radar_enu_dir[1] = cos(coor->radar_azel_gamma_cor) * cos(coor->radar_azel_alpha)	+ 
					(-1. * coor->radar_range * dgamma_dr * sin(coor->radar_azel_gamma_cor) * cos(coor->radar_azel_alpha));
		coor->radar_enu_dir[2] = sin(coor->radar_azel_gamma_cor)						+
					(1. * coor->radar_range * dgamma_dr * cos(coor->radar_azel_gamma_cor));

	} else {
		coor->radar_enu_dir[0] = cos(coor->radar_azel_gamma) * sin(coor->radar_azel_alpha);
		coor->radar_enu_dir[1] = cos(coor->radar_azel_gamma) * cos(coor->radar_azel_alpha);
		coor->radar_enu_dir[2] = sin(coor->radar_azel_gamma);
	}
	coor->radar_enu_dir[3] = 0.;
}


void coordinates_radar_enu2azel(t_zephyros_coordinates *coor)
{
	coordinates_assert_initialized(coor);
	
	coor->radar_range = sqrt(		pow(coor->enu_xyzt[0] - coor->enu_radar_location_xyzt[0],2.)
									+ pow(coor->enu_xyzt[1] - coor->enu_radar_location_xyzt[1],2.)
									+ pow(coor->enu_xyzt[2] - coor->enu_radar_location_xyzt[2],2.));

	if (coor->radar_range == 0.) {
		coor->radar_azel_alpha = 0.;
		coor->radar_azel_gamma = 0.;
	} else {
		coor->radar_azel_alpha		= atan2(coor->enu_xyzt[0] - coor->enu_radar_location_xyzt[0] , coor->enu_xyzt[1] - coor->enu_radar_location_xyzt[1]);
		if (coor->radar_effearthcorrection == 1) {
			/* not implemented yet */
			coor->radar_azel_gamma = asin((coor->enu_xyzt[2] - coor->enu_radar_location_xyzt[2])/ coor->radar_range);
		} else {
			coor->radar_azel_gamma = asin((coor->enu_xyzt[2] - coor->enu_radar_location_xyzt[2]) / coor->radar_range);
		}
	}
}

void coordinates_radar_beam2azel(t_zephyros_coordinates *coor)
{
	coordinates_assert_initialized(coor);
	
	coor->radar_azel_alpha = atan(tan(coor->radar_beam_theta) * cos(coor->radar_beam_phi));
	coor->radar_azel_gamma = asin(sin(coor->radar_beam_theta) * sin(coor->radar_beam_phi));
}

void coordinates_radar_azel2beam(t_zephyros_coordinates *coor)
{
	coordinates_assert_initialized(coor);
	
	coor->radar_beam_theta		= acos(cos(coor->radar_azel_gamma) * cos(coor->radar_azel_alpha));
	coor->radar_beam_phi		= atan(tan(coor->radar_azel_gamma) * sin(coor->radar_azel_alpha));
}

void coordinates_radar_pol_dir(t_zephyros_coordinates *coor)
{
	coordinates_assert_initialized(coor);
	
	//taken over from DeWolf (1990)
	//direction of H polarization
	//coor->radar_pol_hor_dir[0] = -1. * cos((M_PI/2.) - coor->radar_azel_gamma_cor) * sin(coor->radar_azel_alpha);
	//coor->radar_pol_hor_dir[1] = cos(coor->radar_azel_alpha);
	//coor->radar_pol_hor_dir[2] = -1. * sin((M_PI/2.) - coor->radar_azel_gamma_cor) * sin(coor->radar_azel_alpha);
	//coor->radar_pol_hor_dir[3] = 0.;
	
	//direction of V polarization
	//coor->radar_pol_ver_dir[0] = cos((M_PI/2.) - coor->radar_azel_gamma_cor) * cos(coor->radar_azel_alpha);
	//coor->radar_pol_ver_dir[1] = sin(coor->radar_azel_alpha);
	//coor->radar_pol_ver_dir[2] = sin((M_PI/2.) - coor->radar_azel_gamma_cor) * cos(coor->radar_azel_alpha);
	//coor->radar_pol_ver_dir[3] = 0.;


	//From Bringi
	//translation of radar angles to spherical angles
	//theta_i = pi/2 - theta_e
	//phi_i = pi/2 - phi_r
	
	//direction of H polarization
	coor->radar_pol_hor_dir[0] = -1. * sin((M_PI / 2.) - coor->radar_azel_alpha);
	coor->radar_pol_hor_dir[1] = cos((M_PI / 2.) - coor->radar_azel_alpha);
	coor->radar_pol_hor_dir[2] = 0.;
	coor->radar_pol_hor_dir[3] = 0.;
	
	//direction of V polarization
	coor->radar_pol_ver_dir[0] = cos((M_PI / 2.) - coor->radar_azel_gamma_cor) * cos((M_PI / 2.) - coor->radar_azel_alpha);
	coor->radar_pol_ver_dir[1] = cos((M_PI / 2.) - coor->radar_azel_gamma_cor) * sin((M_PI / 2.) - coor->radar_azel_alpha);
	coor->radar_pol_ver_dir[2] = -1. * sin((M_PI / 2.) - coor->radar_azel_gamma_cor);
	coor->radar_pol_ver_dir[3] = 0.;
}



void coordinates_print(t_zephyros_coordinates *coor)
{
	coordinates_assert_initialized(coor);
	
	//incomplete, for debugging purposes
	printf("\n\n\n");
	printf("%-30s = %.2e\n", "wgs84_phirad", coor->wgs84_phirad);	
	printf("%-30s = %.2e\n", "wgs84_lambdarad", coor->wgs84_lambdarad);	
	printf("%-30s = %.2e\n", "wgs84_h", coor->wgs84_h);	
	printf("%-30s = %.2e\n", "wgs84_radar_location_phirad", coor->wgs84_radar_location_phirad);	
	printf("%-30s = %.2e\n", "wgs84_radar_location_lambdarad", coor->wgs84_radar_location_lambdarad);	
	printf("%-30s = %.2e\n", "wgs84_radar_location_h", coor->wgs84_radar_location_h);	
	printf("\n");
	
	printf("%-30s = (%.2e, %.2e, %.2e, %.2e)\n", "ecef_xyzt", coor->ecef_xyzt[0], coor->ecef_xyzt[1], coor->ecef_xyzt[2], coor->ecef_xyzt[3]);
	printf("%-30s = (%.2e, %.2e, %.2e, %.2e)\n", "ecef_radar_location_xyzt", coor->ecef_radar_location_xyzt[0], coor->ecef_radar_location_xyzt[1], coor->ecef_radar_location_xyzt[2], coor->ecef_radar_location_xyzt[3]);
	printf("\n");

	printf("%-30s = %.2e\n", "ref_wgs84_phirad", coor->ref_wgs84_phirad);	
	printf("%-30s = %.2e\n", "ref_wgs84_lambdarad", coor->ref_wgs84_lambdarad);	
	printf("%-30s = %.2e\n", "ref_wgs84_h", coor->ref_wgs84_h);	
	printf("%-30s = (%.2e, %.2e, %.2e, %.2e)\n", "enu_xyzt", coor->enu_xyzt[0], coor->enu_xyzt[1], coor->enu_xyzt[2], coor->enu_xyzt[3]);
	printf("%-30s = (%.2e, %.2e, %.2e, %.2e)\n", "enu_radar_location_xyzt", coor->enu_radar_location_xyzt[0], coor->enu_radar_location_xyzt[1], coor->enu_radar_location_xyzt[2], coor->enu_radar_location_xyzt[3]);
	printf("\n");

	printf("%-30s = %i\n", "radar_effearthcorrection", coor->radar_effearthcorrection);	
	printf("%-30s = %.2e\n", "radar_azel_alpha", coor->radar_azel_alpha);	
	printf("%-30s = %.2e\n", "radar_azel_gamma", coor->radar_azel_gamma);	
	printf("%-30s = %.2e\n", "radar_azel_gamma_cor", coor->radar_azel_gamma_cor);	
	printf("%-30s = %.2e\n", "radar_range", coor->radar_range);	
	printf("\n");

	printf("%-30s = (%.2e, %.2e, %.2e, %.2e)\n", "radar_enu_dir", coor->radar_enu_dir[0], coor->radar_enu_dir[1], coor->radar_enu_dir[2], coor->radar_enu_dir[3]);
	printf("\n");

	printf("%-30s = %.2e\n", "radar_beam_theta", coor->radar_beam_theta);	
	printf("%-30s = %.2e\n", "radar_beam_phi", coor->radar_beam_phi);	

	printf("%-30s = (%.2e, %.2e, %.2e, %.2e)\n", "radar_pol_hor_dir", coor->radar_pol_hor_dir[0], coor->radar_pol_hor_dir[1], coor->radar_pol_hor_dir[2], coor->radar_pol_hor_dir[3]);
	printf("%-30s = (%.2e, %.2e, %.2e, %.2e)\n", "radar_pol_ver_dir", coor->radar_pol_ver_dir[0], coor->radar_pol_ver_dir[1], coor->radar_pol_ver_dir[2], coor->radar_pol_ver_dir[3]);

}
