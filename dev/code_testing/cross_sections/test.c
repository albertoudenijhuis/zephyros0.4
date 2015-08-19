#include "particles.h"
#include "zephyros_config.h"


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

void bhmie_(double *x, double *refrel_r, double *refrel_i, int *nang_in, double *qext0, double *qsca, double *qback, double *gsca);

int main () {
	t_zephyros_config   		*cfg = NULL;
	char cfg_filename[8192];
	char additional_output_filename[8192];
	double tmp;
	double *dummy;

	t_zephyros_particles_widget *scat = malloc(sizeof(t_zephyros_particles_widget));
	t_zephyros_coordinates *coor;
	int i,j,n;
	int i_abc;
	double res[100];
	FILE *fp;
	
	double kmsq;
	
	//needed for fortran Mie code
	double x;	
	double refrel_r, refrel_i;
	int nang_in;
	double qext0, qsca, qback, gsca;

	double eh_u3_sq, ev_u3_sq;


	strcpy(cfg_filename, "../../../input_files/general/standard_output.cfg;../../../input_files/general/water_refractive_index_segelstein.cfg;../../../input_files/general/atmosphere/US1976.cfg;../../../input_files/general/instruments/tara.cfg;../../../input_files/simulation/scatterers/droplet_5mm.cfg;../../../input_files/simulation/wind/vector.cfg;");
	strcpy(additional_output_filename,"./additional_output.zout");

	//load configuration
	zephyros_config_read(cfg_filename, additional_output_filename, &cfg);
	
	//initilize coordinates
	coordinates_initialize_coor(&coor);

	//make up some coordinates
	coor->radar_range 		= 1.e3;
	coor->radar_azel_alpha 	= 0.; //rad
	coor->radar_azel_gamma 	= M_PI / 4.; //rad
	coor->enu_radar_location_xyzt[3] = 0.;

	coordinates_radar_azel2enu(coor);
	coordinates_radar_azelrangedir2enu(coor);
	coordinates_radar_pol_dir(coor);

	//make up some scattering stuff
	scat->air_temperature_K = 285.;
	scat->air_drypressure_hPa = 1013.;
	scat->air_vaporpressure_hPa = 0.;
		
	particles_air_parameters(scat, coor);

	scat->particle_type			= 2;
	scat->particle_D_eqvol_mm	= 5.;
	particles_spheroid_geometry_beard1987(scat);
	particles_terminal_fall_speed_khvorostyanov2005(scat);
	

	
	//*****
	//print some things.
	if (0) {
		scat->radar_wavelength_m = 9.090e-02;	//TARA wavelength
		scat->particle_refractive_index = 0.;

		interpolation_bilint(
			cfg->general->water_refractive_index->lut_realindex,
			&scat->radar_wavelength_m,
			&tmp,
			0,
			dummy);
		scat->particle_refractive_index += tmp;
			
		interpolation_bilint(
			cfg->general->water_refractive_index->lut_imagindex,
			&scat->radar_wavelength_m,
			&tmp,
			0,
			dummy);
		scat->particle_refractive_index -= (tmp * 1.i);

		//set coordinates    
		coor->radar_range 		= 1.e3;
		coor->radar_azel_alpha 	= M_PI / 4.; //rad
		coor->radar_azel_gamma 	= M_PI / 4.; //rad
		coor->enu_radar_location_xyzt[3] = 0.;
		coordinates_radar_azel2enu(coor);
		coordinates_radar_azelrangedir2enu(coor);
		coordinates_radar_pol_dir(coor);

		//particle direction (has to be unit vector)
		scat->particle_dir[0] = 0.;
		scat->particle_dir[1] = 0.;
		scat->particle_dir[2] = 1.;
		
		n = 2;
		for ( i = 0; i < n; i++ ) {
			if (i == 0) {
				scat->particle_D_eqvol_mm = 8.;
			}
			
			if (i == 1) {
				scat->particle_D_eqvol_mm = 9.86;
			}
			//do calculations
			scat->particle_type			= 2;	//spheroid
			particles_spheroid_geometry_beard1987(scat);
			particles_terminal_fall_speed_khvorostyanov2005(scat);
			//particles_cross_sections_mischenko2000(scat, coor, 0);			
			particles_cross_sections_dewolf1990(scat, coor);			


			
			
			printf("\n\n\n\n\n\n\n\nDe Wolf (1990) cross sections for specific configuration and equivolumetric diameter of %.2e mm.\n\n", scat->particle_D_eqvol_mm );

			//inproducts
			eh_u3_sq = pow(		(coor->radar_pol_hor_dir[0] * scat->particle_dir[0])
								+ (coor->radar_pol_hor_dir[1] * scat->particle_dir[1])
								+ (coor->radar_pol_hor_dir[2] * scat->particle_dir[2]), 2.);
			ev_u3_sq = pow(		(coor->radar_pol_ver_dir[0] * scat->particle_dir[0])
								+ (coor->radar_pol_ver_dir[1] * scat->particle_dir[1])
								 + (coor->radar_pol_ver_dir[2] * scat->particle_dir[2]), 2.);
			
					


			printf("eh_u3_sq = %.5e\n", eh_u3_sq);
			printf("ev_u3_sq = %.5e\n", ev_u3_sq);
			printf("pow(scat->dewolf1990_shapefactor_biglambda12, 2.) + (pow(scat->dewolf1990_shapefactor_biglambda12 - scat->dewolf1990_shapefactor_biglambda3, 2.) * pow(ev_u3_sq , 2.)) = %.10e\n", pow(scat->dewolf1990_shapefactor_biglambda12, 2.) + (pow(scat->dewolf1990_shapefactor_biglambda12 - scat->dewolf1990_shapefactor_biglambda3, 2.) * pow(ev_u3_sq , 2.)) );
			printf("(2. * scat->dewolf1990_shapefactor_biglambda12 * (scat->dewolf1990_shapefactor_biglambda12 - scat->dewolf1990_shapefactor_biglambda3) * ev_u3_sq) = %.10e\n", (2. * scat->dewolf1990_shapefactor_biglambda12 * (scat->dewolf1990_shapefactor_biglambda12 - scat->dewolf1990_shapefactor_biglambda3) * ev_u3_sq));
					
			//print coordinates
			coordinates_print(coor);
			
			//print particle widget
			particle_print_widget(scat);
		}
	}


	
	
	
	
	
	//*****
	//output.txt
	if (1) {
	   fp = fopen("output.txt", "w");

		scat->radar_wavelength_m = 9.090e-02;	
		scat->particle_refractive_index = 0.;

		interpolation_bilint(
			cfg->general->water_refractive_index->lut_realindex,
			&scat->radar_wavelength_m,
			&tmp,
			0,
			dummy);
		scat->particle_refractive_index += tmp;
			
		interpolation_bilint(
			cfg->general->water_refractive_index->lut_imagindex,
			&scat->radar_wavelength_m,
			&tmp,
			0,
			dummy);
		scat->particle_refractive_index -= (tmp * 1.i);


		//particle direction (has to be unit vector)
		scat->particle_dir[0] = 0.;
		scat->particle_dir[1] = .2;
		scat->particle_dir[2] = .8;
		
		tmp = sqrt(pow(scat->particle_dir[0], 2.) + pow(scat->particle_dir[1], 2.) + pow(scat->particle_dir[2], 2.));
		scat->particle_dir[0] /= tmp;
		scat->particle_dir[1] /= tmp;
		scat->particle_dir[2] /= tmp;
		
		//calculate cross sections
		particles_cross_sections_dewolf1990(scat, coor);			
		//printf("\n\n\n\nDe Wolf (1990) cross sections \n");
		//particle_print_widget(scat);
		
		particles_cross_sections_mischenko2000(scat, coor, 0);			
		//printf("\n\n\n\nMishchenko (2000) cross sections \n");
		//particle_print_widget(scat);
			
		
		n = 100;
		for ( i = 0; i < n; i++ ) {
			for ( j = 0; j < n; j++ ) {
				//update coordinates
				coor->radar_azel_alpha 	= (2. * M_PI * i) / n; //rad
				res[0] = coor->radar_azel_alpha;
				coor->radar_azel_gamma = ((M_PI/2.) * j) / n; //rad
				res[1] = coor->radar_azel_gamma;
				
				coordinates_radar_azel2enu(coor);
				coordinates_radar_azelrangedir2enu(coor);
				coordinates_radar_pol_dir(coor);

				//De Wolf
				particles_cross_sections_dewolf1990(scat, coor);			
				res[2] = scat->particle_sigma_hh;
				res[3] = scat->particle_sigma_hv;
				res[4] = scat->particle_sigma_vv;
				
				//Mishchenko
				particles_cross_sections_mischenko2000_update_coor(scat, coor, 0);
				res[5] = scat->particle_sigma_hh;
				res[6] = scat->particle_sigma_hv;
				res[7] = scat->particle_sigma_vv;
				
				fprintf(fp, "%.2e %.2e     %.2e %.2e %.2e     %.2e %.2e %.2e\n", res[0], res[1], res[2], res[3], res[4], res[5], res[6], res[7]);
			}
		}

		fclose(fp);
	}


	//*****
	//output2.txt
	if (1) {
		//cross sections as a function of wavelength
	   fp = fopen("output2.txt", "w");

		//particle direction (has to be unit vector)
		scat->particle_dir[0] = 0.;
		scat->particle_dir[1] = .2;
		scat->particle_dir[2] = .8;
		
		tmp = sqrt(pow(scat->particle_dir[0], 2.) + pow(scat->particle_dir[1], 2.) + pow(scat->particle_dir[2], 2.));
		scat->particle_dir[0] /= tmp;
		scat->particle_dir[1] /= tmp;
		scat->particle_dir[2] /= tmp;

		//set coordinates
		coor->radar_azel_alpha 	= (M_PI/8.);
		coor->radar_azel_gamma = (M_PI/4.); //pi / 4 rad = 45 deg
		coordinates_radar_azel2enu(coor);
		coordinates_radar_azelrangedir2enu(coor);
		coordinates_radar_pol_dir(coor);
		
		n = 500;
		for ( i = 0; i < n; i++ ) {
			//update wavelength
			scat->radar_wavelength_m = 
				pow(10., -3. + ((1. - -3.) * (i / (n-1.))));	
			res[0] = scat->radar_wavelength_m;
					
			scat->particle_refractive_index = 0.;
			interpolation_bilint(
				cfg->general->water_refractive_index->lut_realindex,
				&scat->radar_wavelength_m,
				&tmp,
				0,
				dummy);
			scat->particle_refractive_index += tmp;
				
			interpolation_bilint(
				cfg->general->water_refractive_index->lut_imagindex,
				&scat->radar_wavelength_m,
				&tmp,
				0,
				dummy);
			scat->particle_refractive_index -= (tmp * 1.i);
			
			
			//De Wolf
			particles_cross_sections_dewolf1990(scat, coor);			
			res[1] = 10. * log10(scat->particle_sigma_hh);
			res[2] = 10. * log10(scat->particle_sigma_hv);
			res[3] = 10. * log10(scat->particle_sigma_vv);
			res[7] = scat->particle_sigma_hh;
			res[8] = scat->particle_sigma_hv;
			res[9] = scat->particle_sigma_vv;
			
			//Mishchenko
			particles_cross_sections_mischenko2000(scat, coor, 0);			
			res[4] = 10. * log10(scat->particle_sigma_hh);
			res[5] = 10. * log10(scat->particle_sigma_hv);
			res[6] = 10. * log10(scat->particle_sigma_vv);
			res[10] = scat->particle_sigma_hh;
			res[11] = scat->particle_sigma_hv;
			res[12] = scat->particle_sigma_vv;
				
			//classical expression, (Eq. 3.6 from Doviak).
			kmsq = pow(cabs((cpow(scat->particle_refractive_index, 2.)- 1.) / (cpow(scat->particle_refractive_index, 2.) + 2.)), 2.);
			res[13] = (pow(M_PI, 5.) / pow(scat->radar_wavelength_m, 4.)) * kmsq * pow(scat->particle_D_eqvol_mm * 1.e-3, 6.);
			res[14] = 10. * log10(res[13]);

			//Mie Code for spherical droplets (borhen hufman)
			x = 2. * M_PI * (scat->particle_D_eqvol_mm * 1.e-3 / 2.) / scat->radar_wavelength_m;
			refrel_r = creal(scat->particle_refractive_index);
			refrel_i = cimag(scat->particle_refractive_index);
			nang_in = 2;	
			bhmie_(&x, &refrel_r, &refrel_i, &nang_in, &qext0, &qsca, &qback, &gsca);
			res[15] = (M_PI * 4.) * (M_PI / 4.) *  pow(scat->particle_D_eqvol_mm * 1.e-3, 2.) * qback;
			res[16] =  10. * log10(res[15]);
				
			for ( j = 0; j < 17; j++ ) {
				fprintf(fp, "%.2e   ", res[j]);
			}
			fprintf(fp, "\n");
		}

		fclose(fp);
	}



	//*****
	//output3_a.txt, output3_b.txt, output3_c.txt
	
	if (1) {

		//comparison of spherical vs spheroid cross sections
		//spheroid cross sections, for simple case
		for ( i_abc = 0; i_abc < 3; i_abc++ ) {
			if (i_abc == 0) fp = fopen("output3_a.txt", "w");
			if (i_abc == 1) fp = fopen("output3_b.txt", "w");
			if (i_abc == 2) fp = fopen("output3_c.txt", "w");

			//wavelength
			if (i_abc == 0) scat->radar_wavelength_m = 9.090e-02;	//TARA wavelength
			if (i_abc == 1) scat->radar_wavelength_m = 0.0316;		//IDRA wavelength
			if (i_abc == 2) scat->radar_wavelength_m = 0.00856;		//35 GHz, cloud radar

			//set coordinates    
			coor->radar_range 		= 1.e3;
			coor->radar_azel_alpha 	= M_PI / 4.; //rad
			coor->radar_azel_gamma 	= M_PI / 4.; //rad
			coor->enu_radar_location_xyzt[3] = 0.;
			coordinates_radar_azel2enu(coor);
			coordinates_radar_azelrangedir2enu(coor);
			coordinates_radar_pol_dir(coor);

			//particle direction (has to be unit vector)
			scat->particle_dir[0] = 0.;
			scat->particle_dir[1] = 0.;
			scat->particle_dir[2] = 1.;


			//interpolate refractive index
			scat->particle_refractive_index = 0.;
			interpolation_bilint(
				cfg->general->water_refractive_index->lut_realindex,
				&scat->radar_wavelength_m,
				&tmp,
				0,
				dummy);
			scat->particle_refractive_index += tmp;
				
			interpolation_bilint(
				cfg->general->water_refractive_index->lut_imagindex,
				&scat->radar_wavelength_m,
				&tmp,
				0,
				dummy);
			scat->particle_refractive_index -= (tmp * 1.i);
				
					
			particle_print_widget(scat);
			
			n = 500;
			for ( i = 0; i < n; i++ ) {
				//update diameter
				scat->particle_D_eqvol_mm = 
					pow(10., 0. + ((1. - 0.) * (i / (n-1.))));	
				res[0] = scat->particle_D_eqvol_mm;

				//************************
				//spherical cross sections
				scat->particle_type			= 1;	//spherical
				particles_spheroid_geometry_beard1987(scat);
				particles_terminal_fall_speed_khvorostyanov2005(scat);


				//classical expression, (Eq. 3.6 from Doviak).
				kmsq = pow(cabs((cpow(scat->particle_refractive_index, 2.)- 1.) / (pow(scat->particle_refractive_index, 2.) + 2.)), 2.);
				res[1] = (pow(M_PI, 5.) / pow(scat->radar_wavelength_m, 4.)) * kmsq * pow(scat->particle_D_eqvol_mm * 1.e-3, 6.);
				res[2] = 10. * log10(res[1]);

				//De Wolf
				particles_cross_sections_dewolf1990(scat, coor);			
				res[3] = scat->particle_sigma_hh;
				res[4] = scat->particle_sigma_hv;
				res[5] = scat->particle_sigma_vv;
				res[6] = 10. * log10(res[3]);
				res[7] = 10. * log10(res[4]);
				res[8] = 10. * log10(res[5]);
				
				//Mishchenko
				particles_cross_sections_mischenko2000(scat, coor, 0);			
				res[9] = scat->particle_sigma_hh;
				res[10] = scat->particle_sigma_hv;
				res[11] = scat->particle_sigma_vv;
				res[12] = 10. * log10(res[9]);
				res[13] = 10. * log10(res[10]);
				res[14] = 10. * log10(res[11]);
					

				//************************
				//spheroid cross sections
				scat->particle_type			= 2;	//spheroid
				particles_spheroid_geometry_beard1987(scat);
				particles_terminal_fall_speed_khvorostyanov2005(scat);

				//De Wolf
				particles_cross_sections_dewolf1990(scat, coor);			
				res[15] = scat->particle_sigma_hh;
				res[16] = scat->particle_sigma_hv;
				res[17] = scat->particle_sigma_vv;
				res[18] = 10. * log10(res[15]);
				res[19] = 10. * log10(res[16]);
				res[20] = 10. * log10(res[17]);
				
				//Mishchenko
				particles_cross_sections_mischenko2000(scat, coor, 0);			
				res[21] = scat->particle_sigma_hh;
				res[22] = scat->particle_sigma_hv;
				res[23] = scat->particle_sigma_vv;
				res[24] = 10. * log10(res[21]);
				res[25] = 10. * log10(res[22]);
				res[26] = 10. * log10(res[23]);
				
				res[27] = scat->particle_axis_ratio;

				

				//Mie Code for spherical droplets (borhen hufman)
				x = 2. * M_PI * (scat->particle_D_eqvol_mm * 1.e-3 / 2.) / scat->radar_wavelength_m;
				refrel_r = creal(scat->particle_refractive_index);
				refrel_i = cimag(scat->particle_refractive_index);
				nang_in = 2;	
				bhmie_(&x, &refrel_r, &refrel_i, &nang_in, &qext0, &qsca, &qback, &gsca);
				res[28] = (M_PI * 4.) * (M_PI / 4.) *  pow(scat->particle_D_eqvol_mm * 1.e-3, 2.) * qback;
				res[29] =  10. * log10(res[28]);
				
				for ( j = 0; j < 30; j++ ) {
					fprintf(fp, "%.2e   ", res[j]);
				}
				fprintf(fp, "\n");
			}

			fclose(fp);
		}
	}







	//clean up.
	coordinates_free_coor(&coor);
	free(scat);
	
	zephyros_config_free(&cfg);	

    return 0;
}
