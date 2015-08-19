#include "particles.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

int main () {
	t_zephyros_particles_widget *scat = malloc(sizeof(t_zephyros_particles_widget));
	t_zephyros_coordinates *coor;
	int i,j,n;
	double res[100];
	FILE *fp;




	//initilize coordinates
	coordinates_initialize_coor(&coor);

	//make up some coordinates
	coor->radar_range 		= 1.e3;
	coor->radar_azel_alpha 	= 0.; //rad
	coor->radar_azel_gamma 	= M_PI / 4.; //rad

	coordinates_radar_azel2enu(coor);
	coordinates_radar_azelrangedir2enu(coor);
	coordinates_radar_pol_dir(coor);

	//make up some scattering stuff
	scat->air_temperature_K =  288.15;
	scat->air_drypressure_hPa = 1013.25;
	scat->air_vaporpressure_hPa = 0.;
		
	particles_air_parameters(scat, coor);

	

	

	
	
	
	
	//*****
	//output.txt
	if (1) {
	   fp = fopen("output.txt", "w");

		n = 100;
		for ( i = 0; i < n; i++ ) {
			scat->particle_type			= 2; //spheroid
			scat->particle_D_eqvol_mm	= pow(10., -2. + 3. * i / (n-1) );
			particles_spheroid_geometry_beard1987(scat);
			particles_terminal_fall_speed_khvorostyanov2005(scat);

			res[0] = scat->particle_D_eqvol_mm;
			res[1] = scat->particle_axis_ratio;
			res[2] = scat->particle_terminal_fall_speed;
			res[3] = scat->particle_inertial_distance_z;
			res[4] = scat->particle_inertial_distance_z / scat->particle_terminal_fall_speed;
			res[5] = scat->particle_inertial_distance_xy;
			res[6] = scat->particle_inertial_distance_xy / 1.;


			scat->particle_type			= 1; //spherical
			scat->particle_D_eqvol_mm	= pow(10., -2. + 3. * i / (n-1) );
			particles_spheroid_geometry_beard1987(scat);
			particles_terminal_fall_speed_khvorostyanov2005(scat);

			res[7] = scat->particle_terminal_fall_speed;
			res[8] = scat->particle_inertial_distance_z;
			res[9] = scat->particle_inertial_distance_z /  scat->particle_terminal_fall_speed;
			res[10] = scat->particle_inertial_distance_xy;
			res[11] = scat->particle_inertial_distance_xy / 1.;


			for ( j = 0; j < 12; j++ ) {
				fprintf(fp, "%.2e   ", res[j]);
			}
			fprintf(fp, "\n");	
		}
		fclose(fp);
	}

	//clean up.
	coordinates_free_coor(&coor);
	free(scat);
	

    return 0;
}
