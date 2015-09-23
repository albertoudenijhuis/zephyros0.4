/*
Description: 
	TBD

Revision History:
	2014

Functions:
	TBD
	
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

#include "turbulence.h"
#include "ltqnorm.h"
#include "edr.h"
#include "util.h"
#include "func.h"


#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif	
	
void turbulence_initialize_widget(t_zephyros_turbulence_widget **pcfg)
{
	t_zephyros_turbulence_widget *cfg = malloc(sizeof(t_zephyros_turbulence_widget));
	
	cfg->field = NULL;
	cfg->grid_edr = NULL;
	cfg->grid_edr13 = NULL;
	cfg->lut_edr13 = NULL;
	cfg->freqx = NULL;
	cfg->freqy = NULL;
	cfg->freqz = NULL;
	cfg->grid_karman_L = NULL;
	cfg->grid_kolmogorov_constant = NULL;
	cfg->lut_karman_L = NULL;
	cfg->lut_kolmogorov_constant = NULL;
	cfg->random_fourier_cx = NULL;
	cfg->random_fourier_cy = NULL;
	cfg->random_fourier_cz = NULL;
	cfg->random_shift_x = NULL;
	cfg->random_shift_y = NULL;
	cfg->random_shift_z = NULL;
	cfg->random_fourier_c0 = NULL;
	cfg->random_fourier_c1 = NULL;
	cfg->random_fourier_c2 = NULL;
	cfg->u_fourier = NULL;
	cfg->v_fourier = NULL;
	cfg->w_fourier = NULL;
	cfg->random_fourier_beta = NULL;
	cfg->random_a = NULL;
	cfg->random_b = NULL;
	cfg->lambdak = NULL;
	cfg->alphak = NULL;
	cfg->turb_field = NULL;
	cfg->turb_u = NULL;
	cfg->turb_v = NULL;
	cfg->turb_w = NULL;
	cfg->turb_lut_u = NULL;
	cfg->turb_lut_v = NULL;
	cfg->turb_lut_w = NULL;

	cfg->fit_edr13 = 0;
	cfg->grid_edr13_err = NULL;
	cfg->edr13_ecm.mat 		= NULL;
	cfg->edr13_ecm.mat_S 	= NULL;
	cfg->edr13_ecm.mat_N 	= NULL;
	*pcfg = cfg;
}



void turbulence_prepare_widget(t_zephyros_turbulence_widget *cfg)
{
	double RN1, RN2;
	double coef;
	int ix, iy, iz;
	int ik, im;
	int i;
	int tmpi;
	FILE *fp;
	int *tau;
	double deltax;
	double deltay;
	double deltaz;
	double *dummy;
	double *tmpxyzt;
	
	t_turbulence_karmanspec *karmanspec = malloc(sizeof(t_turbulence_karmanspec));
	
	cfg->calibration_factor = 1.;	
		
	//calculate edr13 for interpolation
	cfg->grid_edr13 = malloc(cfg->field->n * sizeof(double));
	for (i=0; i<cfg->field->n; i++) {	
		if (cfg->grid_edr[i] == 0.) {
			cfg->grid_edr13[i] = 0.;
		} else {
			cfg->grid_edr13[i] = pow(cfg->grid_edr[i], 1./3.);
		}
	}
	
	//prepare luts
	util_field2lut(cfg->field, cfg->grid_edr13, 0, &cfg->lut_edr13);

	if ((cfg->type == 1) | (cfg->type == 2)) {	//mann1998, ctm
		util_field2lut(cfg->field, cfg->grid_karman_L, 0, &cfg->lut_karman_L);
		util_field2lut(cfg->field, cfg->grid_kolmogorov_constant, 0, &cfg->lut_kolmogorov_constant);
	}
	
	//set frequencies
	if ((cfg->type == 1) |
		(cfg->type == 2) |
		(cfg->type == 3))
	 {
		cfg->freqx = malloc(cfg->Nx * sizeof(double));
		deltax = cfg->Lx / cfg->Nx;
		f_fftfreq( &cfg->Nx, &deltax, cfg->freqx);
		for (ix=0; ix<cfg->Nx; ix++) {cfg->freqx[ix] *= 2. * M_PI;}
					
		cfg->freqy = malloc(cfg->Ny * sizeof(double));
		deltay = cfg->Ly / cfg->Ny;
		f_fftfreq( &cfg->Ny, &deltay, cfg->freqy);
		for (iy=0; iy<cfg->Ny; iy++) {cfg->freqy[iy] *= 2. * M_PI;}
	}
	if ((cfg->type == 1) | (cfg->type == 2)) {
		cfg->freqz = malloc(cfg->Nz * sizeof(double));
		deltaz = cfg->Lz / cfg->Nz;
		f_fftfreq( &cfg->Nz, &deltaz, cfg->freqz);
		for (iz=0; iz<cfg->Nz; iz++) {cfg->freqz[iz] *= 2. * M_PI;}
	}
	
	if (cfg->type == 1) {
		//mann1998 ...
		//memory allocation
		cfg->random_fourier_c0 = malloc(cfg->Nx * sizeof(double complex**));
		cfg->random_fourier_c1 = malloc(cfg->Nx * sizeof(double complex**));
		cfg->random_fourier_c2 = malloc(cfg->Nx * sizeof(double complex**));
		for (ix=0; ix<cfg->Nx; ix++) {
			cfg->random_fourier_c0[ix] = malloc(cfg->Ny * sizeof(double complex*));
			cfg->random_fourier_c1[ix] = malloc(cfg->Ny * sizeof(double complex*));
			cfg->random_fourier_c2[ix] = malloc(cfg->Ny * sizeof(double complex*));
			for (iy=0; iy<cfg->Ny; iy++) {
				cfg->random_fourier_c0[ix][iy] = malloc(cfg->Nz * sizeof(double complex));
				cfg->random_fourier_c1[ix][iy] = malloc(cfg->Nz * sizeof(double complex));
				cfg->random_fourier_c2[ix][iy] = malloc(cfg->Nz * sizeof(double complex));			
			}
		}

		if (cfg->random_numbers == 0) {
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					for (iz=0; iz<cfg->Nz; iz++) {
						RN1 = ltqnorm(uniform_random()); RN2 = uniform_random();
						cfg->random_fourier_c0[ix][iy][iz] = RN1 * cexp( I * (2.  * M_PI * RN2) );
						RN1 = ltqnorm(uniform_random()); RN2 = uniform_random();
						cfg->random_fourier_c1[ix][iy][iz] = RN1 * cexp( I * (2.  * M_PI * RN2) );
						RN1 = ltqnorm(uniform_random()); RN2 = uniform_random();
						cfg->random_fourier_c2[ix][iy][iz] = RN1 * cexp( I * (2.  * M_PI * RN2) );
					}
				}
			}
		} 
		
		if (cfg->random_numbers == 1) {
			fp = fopen(cfg->rn_file,"r"); // read mode
			printf("Reading random numbers file: %s\n",cfg->rn_file);
			if( fp == NULL )
			{
				perror("Error while opening the configuration file.\n");
				exit(EXIT_FAILURE);
			}
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					for (iz=0; iz<cfg->Nz; iz++) {
						RN1 = ltqnorm(read_uniform_random(fp)); RN2 = read_uniform_random(fp);
						cfg->random_fourier_c0[ix][iy][iz] = RN1 * cexp( I * (2.  * M_PI * RN2) );
						RN1 = ltqnorm(read_uniform_random(fp)); RN2 = read_uniform_random(fp);
						cfg->random_fourier_c1[ix][iy][iz] = RN1 * cexp( I * (2.  * M_PI * RN2) );
						RN1 = ltqnorm(read_uniform_random(fp)); RN2 = read_uniform_random(fp);
						cfg->random_fourier_c2[ix][iy][iz] = RN1 * cexp( I * (2.  * M_PI * RN2) );
					}
				}
			}
			fclose(fp);
		}

		cfg->u_fourier = malloc(cfg->Nx * sizeof(double complex**));
		cfg->v_fourier = malloc(cfg->Nx * sizeof(double complex**));
		cfg->w_fourier = malloc(cfg->Nx * sizeof(double complex**));
		for (ix=0; ix<cfg->Nx; ix++) {
			cfg->u_fourier[ix] = malloc(cfg->Ny * sizeof(double complex*));
			cfg->v_fourier[ix] = malloc(cfg->Ny * sizeof(double complex*));
			cfg->w_fourier[ix] = malloc(cfg->Ny * sizeof(double complex*));
			for (iy=0; iy<cfg->Ny; iy++) {
				cfg->u_fourier[ix][iy] = malloc(cfg->Nz * sizeof(double complex));
				cfg->v_fourier[ix][iy] = malloc(cfg->Nz * sizeof(double complex));
				cfg->w_fourier[ix][iy] = malloc(cfg->Nz * sizeof(double complex));			
			}
		}

		karmanspec->a = 1.;
		karmanspec->L = 1.e10;
		karmanspec->edr = 1.;
		
		for (ix=0; ix<cfg->Nx; ix++) {
			for (iy=0; iy<cfg->Ny; iy++) {
				for (iz=0; iz<cfg->Nz; iz++) {
					//u
					cfg->u_fourier[ix][iy][iz]  = turbulence_mann1998_C(karmanspec, cfg->freqx[ix], cfg->freqy[iy], cfg->freqz[iz], 0, 0, cfg->Lx / cfg->Nx, cfg->Ly / cfg->Ny, cfg->Lz / cfg->Nz) * cfg->random_fourier_c0[ix][iy][iz];
					cfg->u_fourier[ix][iy][iz] += turbulence_mann1998_C(karmanspec, cfg->freqx[ix], cfg->freqy[iy], cfg->freqz[iz], 0, 1, cfg->Lx / cfg->Nx, cfg->Ly / cfg->Ny, cfg->Lz / cfg->Nz) * cfg->random_fourier_c1[ix][iy][iz];
					cfg->u_fourier[ix][iy][iz] += turbulence_mann1998_C(karmanspec, cfg->freqx[ix], cfg->freqy[iy], cfg->freqz[iz], 0, 2, cfg->Lx / cfg->Nx, cfg->Ly / cfg->Ny, cfg->Lz / cfg->Nz) * cfg->random_fourier_c2[ix][iy][iz];

					//v
					cfg->v_fourier[ix][iy][iz]  = turbulence_mann1998_C(karmanspec, cfg->freqx[ix], cfg->freqy[iy], cfg->freqz[iz], 1, 0, cfg->Lx / cfg->Nx, cfg->Ly / cfg->Ny, cfg->Lz / cfg->Nz) * cfg->random_fourier_c0[ix][iy][iz];
					cfg->v_fourier[ix][iy][iz] += turbulence_mann1998_C(karmanspec, cfg->freqx[ix], cfg->freqy[iy], cfg->freqz[iz], 1, 1, cfg->Lx / cfg->Nx, cfg->Ly / cfg->Ny, cfg->Lz / cfg->Nz) * cfg->random_fourier_c1[ix][iy][iz];
					cfg->v_fourier[ix][iy][iz] += turbulence_mann1998_C(karmanspec, cfg->freqx[ix], cfg->freqy[iy], cfg->freqz[iz], 1, 2, cfg->Lx / cfg->Nx, cfg->Ly / cfg->Ny, cfg->Lz / cfg->Nz) * cfg->random_fourier_c2[ix][iy][iz];

					//w
					cfg->w_fourier[ix][iy][iz]  = turbulence_mann1998_C(karmanspec, cfg->freqx[ix], cfg->freqy[iy], cfg->freqz[iz], 2, 0, cfg->Lx / cfg->Nx, cfg->Ly / cfg->Ny, cfg->Lz / cfg->Nz) * cfg->random_fourier_c0[ix][iy][iz];
					cfg->w_fourier[ix][iy][iz] += turbulence_mann1998_C(karmanspec, cfg->freqx[ix], cfg->freqy[iy], cfg->freqz[iz], 2, 1, cfg->Lx / cfg->Nx, cfg->Ly / cfg->Ny, cfg->Lz / cfg->Nz) * cfg->random_fourier_c1[ix][iy][iz];
					cfg->w_fourier[ix][iy][iz] += turbulence_mann1998_C(karmanspec, cfg->freqx[ix], cfg->freqy[iy], cfg->freqz[iz], 2, 2, cfg->Lx / cfg->Nx, cfg->Ly / cfg->Ny, cfg->Lz / cfg->Nz) * cfg->random_fourier_c2[ix][iy][iz];
				}
			}
		}

		//transfer to field
		//reduces calculation time
		if (1) {			
			fields_initialize(&(cfg->turb_field));
			strcpy(cfg->turb_field->name, "turb_field");		
			
			cfg->turb_field->n_x = cfg->Nx;
			cfg->turb_field->n_y = cfg->Ny;
			cfg->turb_field->n_z = cfg->Nz;
			cfg->turb_field->n_t = 1;
			cfg->turb_field->vec_x = malloc(cfg->turb_field->n_x*sizeof(double));
			cfg->turb_field->vec_y = malloc(cfg->turb_field->n_y*sizeof(double));
			cfg->turb_field->vec_z = malloc(cfg->turb_field->n_z*sizeof(double));
			cfg->turb_field->vec_t = malloc(cfg->turb_field->n_t*sizeof(double));
			
			for (ix=0; ix<cfg->Nx; ix++) 
				cfg->turb_field->vec_x[ix] = 1. * ix * cfg->Lx / cfg->Nx;
			for (iy=0; iy<cfg->Ny; iy++) 
				cfg->turb_field->vec_y[iy] = 1. * iy * cfg->Ly / cfg->Ny;
			for (iz=0; iz<cfg->Nz; iz++) 
				cfg->turb_field->vec_z[iz] = 1. * iz * cfg->Lz / cfg->Nz;
			cfg->turb_field->vec_t[0] = 0.;

			cfg->turb_field->n = cfg->turb_field->n_x * cfg->turb_field->n_y * cfg->turb_field->n_z * cfg->turb_field->n_t;
			cfg->turb_u = malloc(cfg->turb_field->n * sizeof(double));
			cfg->turb_v = malloc(cfg->turb_field->n * sizeof(double));
			cfg->turb_w = malloc(cfg->turb_field->n * sizeof(double));
			
			tmpxyzt = malloc(4 * sizeof(double));
			i= -1;

			for (iz=0; iz<cfg->Nz; iz++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					for (ix=0; ix<cfg->Nx; ix++) {						
						i++;
						
						//print status
						tmpi = (cfg->turb_field->n / 10); if (tmpi == 0) tmpi = 1;
						if (	(i < 10) |
								( (i % tmpi) == 0)) {
							printf("prepare turbulence field %5i/%5i\n", i + 1, cfg->turb_field->n); fflush(stdout);
						}
						
						tmpxyzt[0] = cfg->turb_field->vec_x[ix];
						tmpxyzt[1] = cfg->turb_field->vec_y[iy];
						tmpxyzt[2] = cfg->turb_field->vec_z[iz];
						tmpxyzt[3] = 0.;
						turbulence_mann1998_uvw(cfg, tmpxyzt, 0, cfg->turb_u + i, 0, dummy);
						turbulence_mann1998_uvw(cfg, tmpxyzt, 1, cfg->turb_v + i, 0, dummy);
						turbulence_mann1998_uvw(cfg, tmpxyzt, 2, cfg->turb_w + i, 0, dummy);
					}
				}
			}
			free(tmpxyzt);

			//in principle memory can be safed by freeing mann variables that are not necessary anymore.
						
			util_field2lut(cfg->turb_field, cfg->turb_u, 0, &cfg->turb_lut_u);
			util_field2lut(cfg->turb_field, cfg->turb_v, 0, &cfg->turb_lut_v);
			util_field2lut(cfg->turb_field, cfg->turb_w, 0, &cfg->turb_lut_w);
			
			//make periodic
			cfg->turb_lut_u->periodic_L = malloc(4 * sizeof(double));
			cfg->turb_lut_u->periodic_L[0] = cfg->Lx;
			cfg->turb_lut_u->periodic_L[1] = cfg->Ly;
			cfg->turb_lut_u->periodic_L[2] = cfg->Lz;
			cfg->turb_lut_u->periodic_L[3] = 1.;
			cfg->turb_lut_u->periodic = 1;

			cfg->turb_lut_v->periodic_L = malloc(4 * sizeof(double));
			cfg->turb_lut_v->periodic_L[0] = cfg->Lx;
			cfg->turb_lut_v->periodic_L[1] = cfg->Ly;
			cfg->turb_lut_v->periodic_L[2] = cfg->Lz;
			cfg->turb_lut_v->periodic_L[3] = 1.;
			cfg->turb_lut_v->periodic = 1;
			
			cfg->turb_lut_w->periodic_L = malloc(4 * sizeof(double));
			cfg->turb_lut_w->periodic_L[0] = cfg->Lx;
			cfg->turb_lut_w->periodic_L[1] = cfg->Ly;
			cfg->turb_lut_w->periodic_L[2] = cfg->Lz;
			cfg->turb_lut_w->periodic_L[3] = 1.;
			cfg->turb_lut_w->periodic = 1;
			
		}
		
		//make the coefficients for a real signal!
		/*
		for (ix=0; ix<cfg->Nx; ix++) {
			for (iy=0; iy<cfg->Ny; iy++) {
				for (iz=0; iz<cfg->Nz; iz++) {
					if ((ix == ((-ix -1) % cfg->Nx)) && (iy == ((-iy -1) % cfg->Ny)) && (iz == ((-iz-1) % cfg->Nz))) {
						cfg->random_fourier_c0[ix][iy][iz] = 
							cabs(cfg->random_fourier_c0[cfg->Nx-ix-1][cfg->Ny-iy-1][cfg->Nz-iz-1]);
						cfg->random_fourier_c1[ix][iy][iz] = 
							cabs(cfg->random_fourier_c1[cfg->Nx-ix-1][cfg->Ny-iy-1][cfg->Nz-iz-1]);
						cfg->random_fourier_c2[ix][iy][iz] = 
							cabs(cfg->random_fourier_c2[cfg->Nx-ix-1][cfg->Ny-iy-1][cfg->Nz-iz-1]);
					} else {
						cfg->random_fourier_c0[ix][iy][iz] =
							conj(cfg->random_fourier_c0[cfg->Nx-ix-1][cfg->Ny-iy-1][cfg->Nz-iz-1]);
						cfg->random_fourier_c1[ix][iy][iz] =
							conj(cfg->random_fourier_c1[cfg->Nx-ix-1][cfg->Ny-iy-1][cfg->Nz-iz-1]);
						cfg->random_fourier_c2[ix][iy][iz] =
							conj(cfg->random_fourier_c2[cfg->Nx-ix-1][cfg->Ny-iy-1][cfg->Nz-iz-1]);
					}
				}
			}
		}
		*/
	}
		
	//type 2, ctm
	if (cfg->type == 2) {
		cfg->random_fourier_cx = malloc(cfg->Nx * sizeof(double complex));
		cfg->random_fourier_cy = malloc(cfg->Ny * sizeof(double complex));
		cfg->random_fourier_cz = malloc(cfg->Nz * sizeof(double complex));
		
		cfg->random_shift_x	= malloc(100 * sizeof(double));
		cfg->random_shift_y	= malloc(100 * sizeof(double));
		cfg->random_shift_z	= malloc(100 * sizeof(double));
		
		if (cfg->random_numbers == 1) {
			fp = fopen(cfg->rn_file,"r"); // read mode
			printf("Reading random numbers file: %s\n",cfg->rn_file);
			if( fp == NULL )
			{
				perror("Error while opening the configuration file.\n");
				exit(EXIT_FAILURE);
			}
			for (ix=0; ix<cfg->Nx; ix++) {
				RN1 = ltqnorm(read_uniform_random(fp)); RN2 = read_uniform_random(fp);
				cfg->random_fourier_cx[ix] = RN1 * cexp( I * (2.  * M_PI * RN2) );
			}
			for (iy=0; iy<cfg->Ny; iy++) {
				RN1 = ltqnorm(read_uniform_random(fp)); RN2 = read_uniform_random(fp);
				cfg->random_fourier_cy[iy] = RN1 * cexp( I * (2.  * M_PI * RN2) );
			}
			for (iz=0; iz<cfg->Nz; iz++) {
				RN1 = ltqnorm(read_uniform_random(fp)); RN2 = read_uniform_random(fp);
				cfg->random_fourier_cz[iz] = RN1 * cexp( I * (2.  * M_PI * RN2) );
			}
			for (i=0; i<100; i++) {
				RN1 = read_uniform_random(fp);
				cfg->random_shift_x[i] = RN1;
				RN1 = read_uniform_random(fp);
				cfg->random_shift_y[i] = RN1;
				RN1 = read_uniform_random(fp);
				cfg->random_shift_z[i] = RN1;
			}
			fclose(fp);
		}

		if (cfg->random_numbers == 0) {
			for (ix=0; ix<cfg->Nx; ix++) {
				RN1 = ltqnorm(uniform_random()); RN2 = uniform_random();
				cfg->random_fourier_cx[ix] = RN1 * cexp( I * (2.  * M_PI * RN2) );
			}
			for (iy=0; iy<cfg->Ny; iy++) {
				RN1 = ltqnorm(uniform_random()); RN2 = uniform_random();
				cfg->random_fourier_cy[iy] = RN1 * cexp( I * (2.  * M_PI * RN2) );
			}
			for (iz=0; iz<cfg->Nz; iz++) {
				RN1 = ltqnorm(uniform_random()); RN2 = uniform_random();
				cfg->random_fourier_cz[iz] = RN1 * cexp( I * (2.  * M_PI * RN2) );
			}
			for (i=0; i<100; i++) {
				RN1 = uniform_random();
				cfg->random_shift_x[i] = RN1;
				RN1 = uniform_random();
				cfg->random_shift_y[i] = RN1;
				RN1 = uniform_random();
				cfg->random_shift_z[i] = RN1;
			}
		}

		//Delete all even frequencies
		tau = malloc(cfg->Nx * sizeof(int));
		f_tau(&cfg->Nx, tau);
		for (ix=0; ix<cfg->Nx; ix++) {
			if (tau[ix] % 2 == 0) {
				cfg->random_fourier_cx[ix] = 0.;
			}
		}
		free(tau);

		tau = malloc(cfg->Ny * sizeof(int));
		f_tau(&cfg->Ny, tau);
		for (iy=0; iy<cfg->Ny; iy++) {
			if (tau[iy] % 2 == 0) {
				cfg->random_fourier_cy[iy] = 0.;
			}
		}
		free(tau);
		
		tau = malloc(cfg->Nz * sizeof(int));
		f_tau(&cfg->Nz, tau);
		for (iz=0; iz<cfg->Nz; iz++) {
			if (tau[iz] % 2 == 0) {
				cfg->random_fourier_cz[iz] = 0.;
			}
		}
		free(tau);
	}

	if (cfg->type == 3) {
		//careta1993
		//memory allocation
		cfg->random_fourier_beta = malloc(cfg->Nx * sizeof(double complex*));

		for (ix=0; ix<cfg->Nx; ix++) {
			cfg->random_fourier_beta[ix] = malloc(cfg->Ny * sizeof(double complex));
		}
		
		if (cfg->random_numbers == 0) {
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					RN1 = ltqnorm(uniform_random()); RN2 = uniform_random();
					cfg->random_fourier_beta[ix][iy] = RN1 * cexp( I * (2.  * M_PI * RN2) );
				}
			}
		} 
		
		if (cfg->random_numbers == 1) {
			fp = fopen(cfg->rn_file,"r"); // read mode
			printf("Reading random numbers file: %s\n",cfg->rn_file);
			if( fp == NULL )
			{
				perror("Error while opening the configuration file.\n");
				exit(EXIT_FAILURE);
			}
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					RN1 = ltqnorm(read_uniform_random(fp)); RN2 = read_uniform_random(fp);
					cfg->random_fourier_beta[ix][iy] = RN1 * cexp( I * (2.  * M_PI * RN2) );
				}
			}
			fclose(fp);
		}
	}

	if (cfg->type == 4) {
		//pinsky2006
		//memory allocation
		cfg->random_a = malloc(cfg->K * sizeof(double*));
		cfg->random_b = malloc(cfg->K * sizeof(double*));
		for (ik=0; ik<cfg->K; ik++) {
			cfg->random_a[ik] = malloc(cfg->M * sizeof(double));
			cfg->random_b[ik] = malloc(cfg->M * sizeof(double));
		}

		if (cfg->random_numbers == 0) {
			for (ik=0; ik<cfg->K; ik++) {
				for (im=0; im<cfg->M; im++) {
					RN1 = ltqnorm(uniform_random());
					cfg->random_a[ik][im] = RN1;
					RN1 = ltqnorm(uniform_random());
					cfg->random_b[ik][im] = RN1;
				}
			}
		} 
		
		if (cfg->random_numbers == 1) {
			fp = fopen(cfg->rn_file,"r"); // read mode
			printf("Reading random numbers file: %s\n",cfg->rn_file);
			if( fp == NULL )
			{
				perror("Error while opening the configuration file.\n");
				exit(EXIT_FAILURE);
			}
			for (ik=0; ik<cfg->K; ik++) {
				for (im=0; im<cfg->M; im++) {
					RN1 = ltqnorm(read_uniform_random(fp));
					cfg->random_a[ik][im] = RN1;
					RN1 = ltqnorm(read_uniform_random(fp));
					cfg->random_b[ik][im] = RN1;
				}
			}
			fclose(fp);
		}
		
		//pinsky2006 file
		cfg->lambdak 	= malloc(cfg->K * sizeof(double));
		cfg->alphak 	= malloc(cfg->K * sizeof(double));
		fp = fopen(cfg->pinsky2006_file,"r"); // read mode
		printf("Reading pinsky2006 file: %s\n",cfg->pinsky2006_file);
		if( fp == NULL )
		{
			perror("Error while opening the configuration file.\n");
			exit(EXIT_FAILURE);
		}
		for (ik=0; ik<cfg->K; ik++) {
			cfg->lambdak[ik] = read_uniform_random(fp); //function name misused
		}
		for (ik=0; ik<cfg->K; ik++) {
			cfg->alphak[ik] = read_uniform_random(fp); //function name misused
		}
		fclose(fp);
	}
	
	if (cfg->type == 5) {	//parametric turbulence
		util_field2lut(cfg->field, cfg->grid_kolmogorov_constant, 0, &cfg->lut_kolmogorov_constant);
		
		if (cfg->fit_edr13) {
			//prepare error coviance matrix
			util_initialize_field_err_cov_matrix(cfg->field, cfg->grid_edr13_err, &cfg->edr13_ecm);
		}
	}
	
	//calibrate turbulence
	if ((cfg->calibration_method == 1)
			| (cfg->calibration_method == 2)
			| (cfg->calibration_method == 3)
			| (cfg->calibration_method == 4)) {
		turbulence_calibrate(cfg);
	}

	free(karmanspec);
}

void turbulence_free_widget(t_zephyros_turbulence_widget **pcfg)
{
	int ix, iy;
	int i;
	t_zephyros_turbulence_widget *cfg = *pcfg;
	
	if (cfg != NULL) {
		fields_free(&cfg->field);
		util_safe_free(&cfg->grid_edr);
		util_safe_free(&cfg->grid_edr13);
		interpolation_free_lut(&cfg->lut_edr13);
		util_safe_free(&cfg->freqx);
		util_safe_free(&cfg->freqy);
		util_safe_free(&cfg->freqz);
		util_safe_free(&cfg->grid_karman_L);
		util_safe_free(&cfg->grid_kolmogorov_constant);
		interpolation_free_lut(&cfg->lut_karman_L);
		interpolation_free_lut(&cfg->lut_kolmogorov_constant);
		util_safe_free(&cfg->random_fourier_cx);
		util_safe_free(&cfg->random_fourier_cy);
		util_safe_free(&cfg->random_fourier_cz);

		util_safe_free(&cfg->random_shift_x);
		util_safe_free(&cfg->random_shift_y);
		util_safe_free(&cfg->random_shift_z);
		
		if (cfg->random_fourier_c0 != NULL) {
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					free(cfg->random_fourier_c0[ix][iy]);
				}
				free(cfg->random_fourier_c0[ix]);
			}
			free(cfg->random_fourier_c0);
		}
		
		if (cfg->random_fourier_c1 != NULL) {
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					free(cfg->random_fourier_c1[ix][iy]);
				}
				free(cfg->random_fourier_c1[ix]);
			}
			free(cfg->random_fourier_c1);
		}
		
		if (cfg->random_fourier_c2 != NULL) {	
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					free(cfg->random_fourier_c2[ix][iy]);
				}
				free(cfg->random_fourier_c2[ix]);
			}
			free(cfg->random_fourier_c2);
		}


		if (cfg->u_fourier != NULL) {
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					free(cfg->u_fourier[ix][iy]);
				}
				free(cfg->u_fourier[ix]);
			}
			free(cfg->u_fourier);
		}

		if (cfg->v_fourier != NULL) {
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					free(cfg->v_fourier[ix][iy]);
				}
				free(cfg->v_fourier[ix]);
			}
			free(cfg->v_fourier);
		}
		if (cfg->w_fourier != NULL) {
			for (ix=0; ix<cfg->Nx; ix++) {
				for (iy=0; iy<cfg->Ny; iy++) {
					free(cfg->w_fourier[ix][iy]);
				}
				free(cfg->w_fourier[ix]);
			}
			free(cfg->w_fourier);
		}		
		if (cfg->random_fourier_beta != NULL) {
			for (i=0; i<cfg->Nx; i++) 
				free(cfg->random_fourier_beta[i]);
			free(cfg->random_fourier_beta);
		}
		if (cfg->random_a != NULL) {
			for (i=0; i<cfg->K; i++) 
				free(cfg->random_a[i]);
			free(cfg->random_a);
		}
		if (cfg->random_b != NULL) {
			for (i=0; i<cfg->K; i++) 
				free(cfg->random_b[i]);
			free(cfg->random_b);
		}
						
		util_safe_free(&cfg->lambdak);
		util_safe_free(&cfg->alphak);
		util_safe_free(&cfg->turb_u);
		util_safe_free(&cfg->turb_v);
		util_safe_free(&cfg->turb_w);		
		interpolation_free_lut(&cfg->turb_lut_u);
		interpolation_free_lut(&cfg->turb_lut_v);
		interpolation_free_lut(&cfg->turb_lut_w);

		if (cfg->edr13_ecm.mat != NULL) cs_spfree(cfg->edr13_ecm.mat);
		if (cfg->edr13_ecm.mat_S != NULL) cs_sfree(cfg->edr13_ecm.mat_S);
		if (cfg->edr13_ecm.mat_N != NULL) cs_nfree(cfg->edr13_ecm.mat_N);

		free(cfg);
		cfg = NULL;
	}
}
	
void turbulence_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives)
{
	double edr13;
	double *dummy;
	int calcderivatives = 0;
	int bilint_special = 0;
	
	interpolation_bilint(cfg->lut_edr13,
			xyzt,
			&edr13,
			0, //no derivatives
			dummy);
	
	if (cfg->type == 1) {
		//mann1998
		turbulence_mann1998_uvw(cfg, xyzt, i, ans, uvw_calcderivatives, uvw_derivatives);
	}
	if (cfg->type == 2) {
		//ctm
		turbulence_ctm_uvw(cfg, xyzt, i, ans, uvw_calcderivatives, uvw_derivatives);
	}
	if (cfg->type == 3) {
		//careta93
		turbulence_careta93_uvw(cfg, xyzt, i, ans, uvw_calcderivatives, uvw_derivatives);
	}
	if (cfg->type == 4) {
		//pinsky2006
		turbulence_pinsky2006_uvw(cfg, xyzt, i, ans, uvw_calcderivatives, uvw_derivatives);
	}

	//calibration factor
	*ans *= edr13 * cfg->calibration_factor;
	if (uvw_calcderivatives) {
		uvw_derivatives[0] *= edr13 * cfg->calibration_factor;
		uvw_derivatives[1] *= edr13 * cfg->calibration_factor;
		uvw_derivatives[2] *= edr13 * cfg->calibration_factor;
		uvw_derivatives[3] *= edr13 * cfg->calibration_factor;
	}
}

/*
void turbulence_mann1998_uvw_particle(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives,
	double 	*D_maj_mm)
{
	//test modus
	//.....
	//in development ...
	//.....

	int ix, iy, iz;
	double complex phase;
	double complex *amp;
	t_turbulence_karmanspec *karmanspec = malloc(sizeof(t_turbulence_karmanspec));
	int calcderivatives, bilint_special;
	double x = xyzt[0];
	double y = xyzt[1];
	double z = xyzt[2];
	double t = xyzt[3];
	
	double u;
	double v;
	double w;

	double 	T_K				= 288.15;
	double 	pair_hPa		= 1013.25;
	double 	pvapor_hPa		= 0.;
	int		particle_type	= 2;
	double  terminal_fall_speed;
	double  inertial_distance_maj;
	double  inertial_distance_min;

	double inertial_fraction;
	double inertial_phaselag;
	double inertial_omega;
	double *dummy;

	//obtain some fall speed and intertial parameters
	terminal_fall_speed_khvorostyanov2005(
		&T_K,
		&pair_hPa,
		&pvapor_hPa,
		D_maj_mm,
		&particle_type,
		&terminal_fall_speed,
		&inertial_distance_maj,
		&inertial_distance_min
	);
	
	//Calculate u, v ,w
	mann1998_uvw(cfg, xyzt, 0, &u, 0, tmpderivatives);
	mann1998_uvw(cfg, xyzt, 1, &v, 1, tmpderivatives);
	mann1998_uvw(cfg, xyzt, 2, &w, 2, tmpderivatives);
		
	calcderivatives = 0;
	bilint_special = 0;

	interpolation_bilint(cfg->lut_kolmogorov_constant,
			xyzt,
			&karmanspec->a,
			0, //no derivatives
			dummy);
	interpolation_bilint(cfg->lut_karman_L,
			xyzt,
			&karmanspec->L,
			0, //no derivatives
			dummy);			
	
	karmanspec->edr = 1.;

	//TBD account for spectrum derivatives

	*ans = 0.;
	if (uvw_calcderivatives) {
		uvw_derivatives[0] = 0.;
		uvw_derivatives[1] = 0.;
		uvw_derivatives[2] = 0.;
		uvw_derivatives[3] = 0.;
	}
	
	for (ix=0; ix<cfg->Nx; ix++) {
		for (iy=0; iy<cfg->Ny; iy++) {
			for (iz=0; iz<cfg->Nz; iz++) {
				phase = I * ((cfg->freqx[ix] * x) + (cfg->freqy[iy] * y) + (cfg->freqz[iz] * z));

				inertial_omega  = sqrt(pow(u * cfg->freqx[ix], 2.) + pow(v * cfg->freqy[iy], 2.) + pow((w - terminal_fall_speed) * cfg->freqz[iz], 2.));
				
				if ((i == 0) | (i == 1)) {
					inertiamodel_fraction(&inertial_omega, &inertial_distance_min, &inertial_fraction);
					inertiamodel_phaselag(&inertial_omega, &inertial_distance_min, &inertial_phaselag);
				}
				if (i == 2) {
					inertiamodel_fraction(&inertial_omega, &inertial_distance_maj, &inertial_fraction);
					inertiamodel_phaselag(&inertial_omega, &inertial_distance_maj, &inertial_phaselag);
				}

				if (i == 0) {
					//u
					amp = &cfg->u_fourier[ix][iy][iz];
				}
				if (i == 1) {
					//v
					amp = &cfg->v_fourier[ix][iy][iz];
				}
				if (i == 2) {
					//w
					amp = &cfg->w_fourier[ix][iy][iz];
				}
				*ans += creal(inertia_fraction * *amp * cexp(phase - (I * inertia_phaselag)));
				
				if (uvw_calcderivatives) {
					uvw_derivatives[0] = creal((I * cfg->freqx[ix]) * inertia_fraction * *amp * cexp(phase - (I * inertia_phaselag)));
					uvw_derivatives[1] = creal((I * cfg->freqy[iy]) * inertia_fraction * *amp * cexp(phase - (I * inertia_phaselag)));
					uvw_derivatives[2] = creal((I * cfg->freqz[iz]) * inertia_fraction * *amp * cexp(phase - (I * inertia_phaselag)));
				}
			}
		}
	}
	
	free(karmanspec);
}
*/

void turbulence_mann1998_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives)
{
	int ix, iy, iz;
	double complex phase;
	double complex *amp;
	t_turbulence_karmanspec *karmanspec; 
	int calcderivatives, bilint_special;
	double x = xyzt[0];
	double y = xyzt[1];
	double z = xyzt[2];
	double t = xyzt[3];
	double *dummy;
	int done;
	
	done = 0;
	
	//fast interpolation
	if ((i == 0) & (cfg->turb_lut_u != NULL)) {
		interpolation_bilint(cfg->turb_lut_u,
				xyzt,
				ans,
				uvw_calcderivatives, 
				uvw_derivatives);
		done = 1;
	}
	if ((i == 1) & (cfg->turb_lut_v != NULL)) {
		interpolation_bilint(cfg->turb_lut_v,
				xyzt,
				ans,
				uvw_calcderivatives, 
				uvw_derivatives);
		done = 1;
	}
	if ((i == 2) & (cfg->turb_lut_w != NULL)) {
		interpolation_bilint(cfg->turb_lut_w,
				xyzt,
				ans,
				uvw_calcderivatives, 
				uvw_derivatives);
		done = 1;
	}
	
	if (done == 0) {
		karmanspec= malloc(sizeof(t_turbulence_karmanspec));	
		interpolation_bilint(cfg->lut_kolmogorov_constant,
				xyzt,
				&karmanspec->a,
				0, //no derivatives
				dummy);
		interpolation_bilint(cfg->lut_karman_L,
				xyzt,
				&karmanspec->L,
				0, //no derivatives
				dummy);			
		karmanspec->edr = 1.;


		*ans = 0.;
		if (uvw_calcderivatives) {
			uvw_derivatives[0] = 0.;
			uvw_derivatives[1] = 0.;
			uvw_derivatives[2] = 0.;
			uvw_derivatives[3] = 0.;
		}
		for (ix=0; ix<cfg->Nx; ix++) {
			for (iy=0; iy<cfg->Ny; iy++) {
				for (iz=0; iz<cfg->Nz; iz++) {
					phase = I * ((cfg->freqx[ix] * x) + (cfg->freqy[iy] * y) + (cfg->freqz[iz] * z));
					if (i == 0) {
						//u
						amp = &(cfg->u_fourier[ix][iy][iz]);
					}
					if (i == 1) {
						//v
						amp = &(cfg->v_fourier[ix][iy][iz]);
					}
					if (i == 2) {
						//w
						amp = &(cfg->w_fourier[ix][iy][iz]);
					}

					*ans += creal(*amp * cexp(phase));
					
					if (uvw_calcderivatives) {
						uvw_derivatives[0] = creal((I * cfg->freqx[ix]) * *amp * cexp(phase));
						uvw_derivatives[1] = creal((I * cfg->freqy[iy]) * *amp * cexp(phase));
						uvw_derivatives[2] = creal((I * cfg->freqz[iz]) * *amp * cexp(phase));
					}
				}
			}
		}
		free(karmanspec);
	}
}

void turbulence_ctm_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives)
{
	t_turbulence_karmanspec *karmanspec = malloc(sizeof(t_turbulence_karmanspec));
	int calcderivatives, bilint_special;
	
	double x = xyzt[0];
	double y = xyzt[1];
	double z = xyzt[2];
	double t = xyzt[3];

	double x2;
	double y2;
	double z2;
	
	int		j;
	int 	jmax = 100;

	double 	L0, L1;
	int ix, iy, iz;
	double complex phase;
	double amp;
	double fct;
	double *dummy;
	
	calcderivatives = 0;
	bilint_special = 0;
	
	interpolation_bilint(cfg->lut_kolmogorov_constant,
			xyzt,
			&karmanspec->a,
			0, //no derivatives
			dummy);
	interpolation_bilint(cfg->lut_karman_L,
			xyzt,
			&karmanspec->L,
			0, //no derivatives
			dummy);			
	karmanspec->edr = 1.;

	*ans = 0.;
	
	for (j = 0; j < jmax; j++) {
		fct = pow(2. , -j);

		if ((i == 0) & (((cfg->Lx / cfg->Nx) * fct)) < (cfg->minL_div_maxL * cfg->Lx)) break;
		if ((i == 1) & (((cfg->Ly / cfg->Ny) * fct)) < (cfg->minL_div_maxL * cfg->Ly)) break;
		if ((i == 2) & (((cfg->Lz / cfg->Nz) * fct)) < (cfg->minL_div_maxL * cfg->Lz)) break;
		
		if (i == 0) {
			x2 = x + cfg->random_shift_x[j % 100] * cfg->Lx * fct;
			for (ix=0; ix<cfg->Nx; ix++) {		 
				phase = I * (cfg->freqx[ix] / fct) * x2;
				amp = turbulence_ctm_C(karmanspec, fabs(cfg->freqx[ix]) / fct, cfg->Nx, cfg->Lx) * cfg->random_fourier_cx[ix];
				*ans += creal(amp * cexp(phase));
			}
		}
		if (i == 1) {
			y2 = y + cfg->random_shift_y[j % 100] * cfg->Ly * fct;
			for (iy=0; iy<cfg->Ny; iy++) {		 
				phase = I * (cfg->freqx[iy] / fct) * y2;
				amp = turbulence_ctm_C(karmanspec, fabs(cfg->freqy[iy]) / fct, cfg->Ny, cfg->Ly) * cfg->random_fourier_cy[iy];
				*ans += creal(amp * cexp(phase));
			}
		}
		if (i == 2) {
			z2 = z + cfg->random_shift_z[j % 100] * cfg->Lz * fct;
			for (iz=0; iz<cfg->Nz; iz++) {		 
				phase = I * (cfg->freqx[iz] / fct) * z2;
				amp = turbulence_ctm_C(karmanspec, fabs(cfg->freqz[iz]) / fct, cfg->Nz, cfg->Lz) * cfg->random_fourier_cz[iz];
				*ans += creal(amp * cexp(phase));
			}
		}
	}
}




void turbulence_careta93_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives)
{
	int ix, iy, iz;
	double complex phase;
	double amp;
	double x = xyzt[0];
	double y = xyzt[1];
	double z = xyzt[2];
	double t = xyzt[3];
	double eps, coef;
		

	eps = 10.; //irrelevant, when calibration turned on.
	coef = sqrt(8. * pow(M_PI, 2.) * 1. * eps);

	*ans = 0.;
	if (uvw_calcderivatives) {
		uvw_derivatives[0] = 0.;
		uvw_derivatives[1] = 0.;
		uvw_derivatives[2] = 0.;
		uvw_derivatives[3] = 0.;
	}
	


	//velocity field is in the xy-plane
	for (ix=0; ix<cfg->Nx; ix++) {
		for (iy=0; iy<cfg->Ny; iy++) {
			amp = coef * turbulence_careta1993_Q(cfg, cfg->freqx[ix], cfg->freqy[iy]) * cfg->random_fourier_beta[ix][iy];
			//phase = 2.0I * M_PI * ((x * ix / cfg->Lx) + (y * iy / cfg->Ly));
			phase = I * ((cfg->freqx[ix] * x) + (cfg->freqy[iy] * y)); 
			
			//eta = amp * phase. eta is describing the stream function
			//u = - deta/dy
			//v = deta/dx
			
			//u
			if (i == 0) {
				*ans += -1. * creal((I * cfg->freqy[iy]) * amp * cexp(phase));
			}
			
			//v
			if (i == 1) {
				*ans += creal((I * cfg->freqx[ix]) * amp * cexp(phase));
			}
			
			if (uvw_calcderivatives) {
				if (i == 0) {
					uvw_derivatives[0] += -1. * creal((I * cfg->freqx[ix]) * (I * cfg->freqy[iy]) * amp * cexp(phase));
					uvw_derivatives[1] += -1. * creal((I * cfg->freqy[iy]) * (I * cfg->freqy[iy]) * amp * cexp(phase));
				}
				if (i == 1) {
					uvw_derivatives[0] += creal((I * cfg->freqx[ix]) * (I * cfg->freqx[ix]) * amp * cexp(phase));
					uvw_derivatives[1] += creal((I * cfg->freqy[iy]) * (I * cfg->freqx[ix]) * amp * cexp(phase));
				}
			}
		}
	}
}

void turbulence_pinsky2006_uvw(
	t_zephyros_turbulence_widget *cfg,
	double *xyzt,
	int i,				//u (i=0), v (i=1), w (i=2)
	double *ans,
	int uvw_calcderivatives,
	double *uvw_derivatives)
{
	double x = xyzt[0];
	double y = xyzt[1];
	double z = xyzt[2];
	double t = xyzt[3];
	int im, ik;
	double phim;
	double tmp;
	
	*ans = 0.;
	if (uvw_calcderivatives) {
		uvw_derivatives[0] = 0.;
		uvw_derivatives[1] = 0.;
		uvw_derivatives[2] = 0.;
		uvw_derivatives[3] = 0.;
	}

	for (im = 1; im <= cfg->M; im++) {
		phim = 2. * M_PI * im / cfg->M;
		for (ik = 0; ik < cfg->K; ik++) {	
			tmp = (  (cfg->random_a[ik][im-1] * cos(cfg->lambdak[ik] * ((x * cos(phim)) + (y * sin(phim))))) +
						(cfg->random_b[ik][im-1] * sin(cfg->lambdak[ik] * ((x * cos(phim)) + (y * sin(phim)))))
						);
						
			//u
			if (i == 0) {
				*ans += cfg->alphak[ik] * cfg->lambdak[ik] * sin(phim) * tmp;
			}
			
			//v
			if (i == 1) {
				*ans += cfg->alphak[ik] * cfg->lambdak[ik] * cos(phim) * tmp;
			}
			
			//derivatives
			//TBD		

		}
	}

}
    
double turbulence_mann1998_C(
	t_turbulence_karmanspec *sp,
	double kx,
	double ky,
	double kz,
	int i,
	int j,
	double deltax,
	double deltay,
	double deltaz)
{
    return (pow(2. * M_PI,3) * (1. / deltax) * (1. / deltay) * (1. / deltaz)) * turbulence_mann1998_A(sp, kx, ky, kz, i, j);
}

double turbulence_mann1998_A(
	t_turbulence_karmanspec *sp,
	double kx,
	double ky,
	double kz,
	int i,
	int j)
{
	double kabs;
	double Ek;
    double fct;
    //i is matrix row, j = matrix column, note that indices start at 0 (i.p.v. 1)
    
    if (i == j) return 0;
    
    kabs = pow(kx,2) + pow(ky,2) + pow(kz,2);
    if (kabs != 0.) {kabs = sqrt(kabs);}

    Ek = turbulence_mann1998_Ek(sp, kabs);
    
    if (Ek == 0) {
        fct = 0.;
    } else {
        fct = sqrt(Ek) / ( sqrt(4. * M_PI) * pow(kabs,2.)  );
    }
    
    if ((i == 0) && (j == 1)) return fct * kz;
    if ((i == 0) && (j == 2)) return -fct * ky;
    if ((i == 1) && (j == 0)) return -fct * kz;
    if ((i == 1) && (j == 2)) return fct * kx;
    if ((i == 2) && (j == 0)) return fct * ky;
    if ((i == 2) && (j == 1)) return -fct * kx;
}

//Eq. (9)
double turbulence_mann1998_Ek(
	t_turbulence_karmanspec *sp,
	double kabs)
{
    return  (   sp->a * pow(sp->edr, 2./3.) * pow(sp->L, 5./3.) *
                (pow(sp->L, 4.) * pow(kabs, 4.) /
                pow((1. + (pow(sp->L, 2.) * pow(kabs, 2.))), 17. / 6.))
            );
}

double turbulence_ctm_C(
	t_turbulence_karmanspec *sp,
	double kabs,
	int 	Nx,
	double Lmax)
{
	double Ek;
	double deltak = (2. * M_PI / Lmax);
	double kabs2 = fabs(kabs);

    if (kabs2 == 0.) {
		return 0.;
	} else {
		//Kolmogorov spectrum
		return sqrt((3./2.) * sp->a * pow(sp->edr, 2./3.) * (pow(kabs2, -5./3.)));
	}
}

double turbulence_careta1993_Q(
	t_zephyros_turbulence_widget *cfg,
	double kx,
	double ky)
{
	double deltax = cfg->Lx / cfg->Nx;
	double deltay = cfg->Ly / cfg->Ny;
	
	//careta1993, spectrum 2
	return  pow(   1. + 
				( 
					(2. * (cfg->lambdax * cfg->lambday) / (deltax * deltay))
					* 
					(2. - cos(deltax * kx) - cos(deltay * ky) )
				)
			,-7. / 6.);	

}



void turbulence_calibrate(t_zephyros_turbulence_widget *cfg)
{
	int 	n = cfg->calibration_n;
	double 	retr_C = cfg->calibration_C;
	int 	periodic = cfg->calibration_periodic;
	int		method = cfg->calibration_method;
	int 	retr_nintervals = cfg->calibration_nint;
	
	int retr_domain = 0; //space

	double *u 				= malloc(n * sizeof(double));
	double *t				= malloc(n * sizeof(double));
	double *l				= malloc(n * sizeof(double));
	double *v2 				= malloc(n * sizeof(double));
	
	double *retr_edr 		= malloc(n * sizeof(double));
	double *retr_edr13_err 	= malloc(n * sizeof(double));
	double *retr_av		 	= malloc(n * sizeof(double));
	double *retr_var	 	= malloc(n * sizeof(double));
	int i;
	
	double *tmp_edr13;
	t_zephyros_interpolation_bilint_lut *original_lut_edr13;
	
	double *xyzt = malloc(4 * sizeof(double));
	double *dummy;
	
	int 	calcderivatives;
	
	//set calibration factor 1
	cfg->calibration_factor = 1.;	
	
	//store memory allocation of original lut.
	original_lut_edr13 = cfg->lut_edr13;
		
	//make temporary lut with edr = 1.
	tmp_edr13 = malloc(cfg->field->n * sizeof(double));
	for (i = 0; i < cfg->field->n; i++) {
		tmp_edr13[i] = 1.;
	}
	util_field2lut(cfg->field, tmp_edr13, 0, &cfg->lut_edr13);

	//generate l, and also u
	for (i = 0; i < n; i++) {
		xyzt[0] = 0.;
		xyzt[1] = 0.;
		xyzt[2] = 0.;
		xyzt[3] = 0.;

		if (cfg->calibration_dir == 0) {
			l[i] = (cfg->calibration_L / n) * i;
			xyzt[0] = l[i];
		}
		if (cfg->calibration_dir == 1) {
			l[i] = (cfg->calibration_L / n) * i;
			xyzt[1] = l[i];
		}
		if (cfg->calibration_dir == 2) {
			l[i] = (cfg->calibration_L / n) * i;
			xyzt[2] = l[i];
		}
	
		turbulence_uvw(cfg, xyzt, cfg->calibration_dir, u + i, 0, dummy);
	}
	
	//apply retrieval and update edr
	if (method == 1) {
		retr_edr_variance(&n, u, t, l, v2, &retr_domain, &retr_C, &n, retr_av, retr_var, retr_edr, retr_edr13_err);
	}

	if (method == 2) {	
		retr_edr_powerspectrum(&n, u, t, l, v2, &retr_domain, &retr_C, &n, &periodic, &retr_nintervals, retr_edr, retr_edr13_err);
	}
	
	if (method == 3) {
		retr_edr_2nd_order_structure_function(&n, u, t, l, v2, &retr_domain, &retr_C, &n, &periodic, retr_edr,	retr_edr13_err);
	}

	cfg->calibration_factor = pow(*retr_edr, -1./3.);
	
	free(u);
	free(t);
	free(l);
	free(v2);
	
	free(retr_edr);
	free(retr_edr13_err);
	free(retr_av);
	free(retr_var);

	//place original edr and edr13 back
	interpolation_free_lut(&cfg->lut_edr13);
	cfg->lut_edr13 = original_lut_edr13;
	free(xyzt);
}
