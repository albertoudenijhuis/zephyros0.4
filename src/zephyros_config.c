/*
Description: 
	Zephyros configuration read code

Revision History:
	2014

Functions:
	Zephyros configuration read code
	
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
#include <string.h> 
#include <math.h>

#include "zephyros_config.h"
#include "specialfunctions.h"

//uncomment next statement for debug mode
//#define _ZEPHYROS_CONFIG_DEBUG

void zephyros_config_read(char file_name[8192], char additional_output_filename[8192], t_zephyros_config **pcfg)
{
	t_zephyros_config_read_widget *rwg = malloc(sizeof(t_zephyros_config_read_widget));
	t_zephyros_config *cfg;
	
	double dummy1,dummy2,dummy3; //TBD
	int dummyint;
	int tmpn;
	double tmp;
	int i, j;
	
	int filename_i0;
	int filename_i1;
	char this_file_name[8192];
	char this_dir_name[8192];
	char this_tmp_char[8192];
	
	t_zephyros_radarfilter *myradarfilter;
	t_zephyros_windfield *mywindfield;
	t_zephyros_scattererfield *myscattererfield;

	t_zephyros_config_read_widget_indices *myindices;

	zephyros_config_initialize(pcfg);
	cfg = *pcfg;
		
	//initialize read widget, rwg
	rwg->commentSign[0] 									= '#';
	strcpy(rwg->section, "");		
	strcpy(rwg->subsection, "");		
	strcpy(rwg->identifier, "");		
	strcpy(rwg->line, "");		
	rwg->simulation_indices.windfield_grid 					= 0;
	rwg->simulation_indices.windfield_wave 					= 0;
	rwg->simulation_indices.windfield_vortex 				= 0;
	rwg->simulation_indices.windfield_turbulence 			= 0;
	rwg->simulation_indices.scattererfield_psd 				= 0;
	rwg->simulation_indices.scattererfield_psd_diameter_i 	= 0;
	rwg->simulation_indices.algorithm_run 					= 0;
	rwg->retrieval_indices.windfield_grid 					= 0;
	rwg->retrieval_indices.windfield_wave 					= 0;
	rwg->retrieval_indices.windfield_vortex 				= 0;
	rwg->retrieval_indices.windfield_turbulence 			= 0;
	rwg->retrieval_indices.scattererfield_psd 				= 0;
	rwg->retrieval_indices.scattererfield_psd_diameter_i 	= 0;
	rwg->retrieval_indices.algorithm_run 					= 0;
	
	//update file_name in configuration
	strcpy(cfg->general->overall->file_name,file_name);

	filename_i0 = 0;
			
	while (1) {
		if (filename_i0 >= strlen(file_name)) break;
		
		//obtain another configuration filename
		for ( filename_i1 = filename_i0; filename_i1 < strlen(file_name); filename_i1++ ) {
			if (*(file_name + filename_i1) == ';') {
				break;
			}
		}

		memset(&this_file_name[0], 0, sizeof(this_file_name));
		strncpy(this_file_name, file_name + filename_i0, filename_i1 - filename_i0);
		this_file_name[filename_i1 - filename_i0] = '\0';

		//initialize position for next configuration file
		filename_i0 = filename_i1 +1;
								
		if (strlen(this_file_name) > 3) {
			rwg->fp = fopen(this_file_name,"r"); // read mode

			printf("Opening configuration file: %s\n",this_file_name); fflush(stdout);
			if( rwg->fp == NULL )
			{
				perror("Error while opening the configuration file.\n");
				exit(EXIT_FAILURE);
			}

			//get this_dir_name .... 
			zephyros_config_get_dirname(this_file_name, this_dir_name);
		 
			while (1) {
				fgetpos(rwg->fp, &rwg->pos_thisline);		//store position of where this line starts
				
				//read full line
				if (fgets (rwg->line, sizeof(rwg->line), rwg->fp) == NULL) {break;}

				//debug, print every line
				//printf("%s\n", rwg->line); fflush(stdout);
						
				fgetpos(rwg->fp, &rwg->pos_nextline);	//store position of where next line starts
				rwg->firstchar[0] = rwg->line[0];		//get first char
								
				if (strncmp(rwg->firstchar,rwg->commentSign, 1) == 0) {
					//comment sign, nothing to be done
				} else {
					//read identifier
					fsetpos(rwg->fp, &rwg->pos_thisline);
					fscanf(rwg->fp, "%s", rwg->identifier);			

					if ( strcmp(rwg->identifier,"section") == 0 ) {
						//read section identifier
						fsetpos(rwg->fp, &rwg->pos_thisline);
						fscanf(rwg->fp, "%s %s", rwg->identifier, rwg->section);			
					}

					if ( strcmp(rwg->identifier,"subsection") == 0 ) {
						//read subsection identifier
						fsetpos(rwg->fp, &rwg->pos_thisline);
						fscanf(rwg->fp, "%s %s", rwg->identifier, rwg->subsection);
					}
				
					//general
					if ( strcmp(rwg->section,"general") == 0 ) {
						zephyros_config_read_________s(rwg, "general", "overall", "version_number", cfg->general->overall->versionnumber);

						//additional output	
						zephyros_config_read_________i(rwg, "general", "additional_output", "print_configuration", &cfg->general->additional_output->print_configuration);
						zephyros_config_read_________i(rwg, "general", "additional_output", "print_detailed_analysis", &cfg->general->additional_output->print_detailed_analysis);

						//atmosphere
						zephyros_config_read__lf_array(rwg, "general", "atmosphere", "vec_x", &cfg->general->atmosphere->field->n_x, &cfg->general->atmosphere->field->vec_x);
						zephyros_config_read__lf_array(rwg, "general", "atmosphere", "vec_y", &cfg->general->atmosphere->field->n_y, &cfg->general->atmosphere->field->vec_y);
						zephyros_config_read__lf_array(rwg, "general", "atmosphere", "vec_z", &cfg->general->atmosphere->field->n_z, &cfg->general->atmosphere->field->vec_z);
						zephyros_config_read__lf_array(rwg, "general", "atmosphere", "vec_t", &cfg->general->atmosphere->field->n_t, &cfg->general->atmosphere->field->vec_t);
						zephyros_config_read__lf_array(rwg, "general", "atmosphere", "grid_T_K", &cfg->general->atmosphere->field->n, &cfg->general->atmosphere->grid_T_K);
						zephyros_config_read__lf_array(rwg, "general", "atmosphere", "grid_pair_hPa", &cfg->general->atmosphere->field->n, &cfg->general->atmosphere->grid_pair_hPa);
						zephyros_config_read__lf_array(rwg, "general", "atmosphere", "grid_pvapor_hPa", &cfg->general->atmosphere->field->n, &cfg->general->atmosphere->grid_pvapor_hPa);

						//instrument
						zephyros_config_read_________i(rwg, "general", "instrument", "type", &cfg->general->instrument->type);
						zephyros_config_read________lf(rwg, "general", "instrument", "transmit_power_watt", &cfg->general->instrument->transmit_power_watt);
						zephyros_config_read________lf(rwg, "general", "instrument", "k_nonuniformity", &cfg->general->instrument->k_nonuniformity);
						zephyros_config_read________lf(rwg, "general", "instrument", "power_gain", &cfg->general->instrument->power_gain);
						zephyros_config_read________lf(rwg, "general", "instrument", "lr_receiver_loss_factor", &cfg->general->instrument->lr_receiver_loss_factor);
						zephyros_config_read________lf(rwg, "general", "instrument", "central_frequency_hz", &cfg->general->instrument->central_frequency_hz);
						zephyros_config_read________lf(rwg, "general", "instrument", "prf_hz", &cfg->general->instrument->prf_hz);
						zephyros_config_read________lf(rwg, "general", "instrument", "temperature_K", &cfg->general->instrument->temperature_K);

						//water_refractive_index
						zephyros_config_read__lf_array(rwg, "general", "water_refractive_index", "wavelength_m", &cfg->general->water_refractive_index->n, &cfg->general->water_refractive_index->wavelength_m);
						zephyros_config_read__lf_array(rwg, "general", "water_refractive_index", "realindex", &cfg->general->water_refractive_index->n, &cfg->general->water_refractive_index->realindex);
						zephyros_config_read__lf_array(rwg, "general", "water_refractive_index", "imagindex", &cfg->general->water_refractive_index->n, &cfg->general->water_refractive_index->imagindex);

						//white1999_integral
						zephyros_config_read__lf_array(rwg, "general", "white1999_integral", "vec_ln_a", &cfg->general->white1999_integral->na, &cfg->general->white1999_integral->vec_ln_a);
						zephyros_config_read__lf_array(rwg, "general", "white1999_integral", "vec_ln_b", &cfg->general->white1999_integral->nb, &cfg->general->white1999_integral->vec_ln_b);
						zephyros_config_read__lf_array(rwg, "general", "white1999_integral", "vec_ln_L", &cfg->general->white1999_integral->nL, &cfg->general->white1999_integral->vec_ln_L);					
						zephyros_config_read__lf_array(rwg, "general", "white1999_integral", "integral_sqrt", &cfg->general->white1999_integral->n, &cfg->general->white1999_integral->integral_sqrt);					
					}

					//simulation or retrieval
					if (( strcmp(rwg->section,"simulation") == 0 ) | ( strcmp(rwg->section,"retrieval") == 0 )) {
						//initialize
						if ( strcmp(rwg->section,"simulation") == 0 ) zephyros_config_initialize_simulation(cfg);
						if ( strcmp(rwg->section,"retrieval") == 0 ) zephyros_config_initialize_retrieval(cfg);
	
						//myindices
						if ( strcmp(rwg->section,"simulation") == 0 ) myindices = &(rwg->simulation_indices);
						if ( strcmp(rwg->section,"retrieval") == 0 ) myindices = &(rwg->retrieval_indices);
							
						//radar filter
						if ( strcmp(rwg->subsection,"radarfilter") == 0 ) {
							if ( strcmp(rwg->section,"simulation") == 0 ) myradarfilter = cfg->simulation->radarfilter;
							if ( strcmp(rwg->section,"retrieval") == 0 ) myradarfilter = cfg->retrieval->radarfilter;
														
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "n_beam_range", &myradarfilter->n_beam_range);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "n_beam_theta", &myradarfilter->n_beam_theta);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "n_beam_phi", &myradarfilter->n_beam_phi);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "n_parmod_az", &myradarfilter->n_parmod_az);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "n_parmod_el", &myradarfilter->n_parmod_el);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "n_t", &myradarfilter->n_t);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "n_spectrum", &myradarfilter->n_spectrum);
							
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_dBZ_hh", &myradarfilter->filter_dBZ_hh);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_dBZ_hv", &myradarfilter->filter_dBZ_hv);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_dBZ_vh", &myradarfilter->filter_dBZ_vh);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_dBZ_vv", &myradarfilter->filter_dBZ_vv);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_dBZdr", &myradarfilter->filter_dBZdr);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_dBLdr", &myradarfilter->filter_dBLdr);
							
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_rho_co", &myradarfilter->filter_rho_co);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_rho_cxh", &myradarfilter->filter_rho_cxh);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_rho_cxv", &myradarfilter->filter_rho_cxv);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_KDP", &myradarfilter->filter_KDP);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_velocity_hh_ms", &myradarfilter->filter_Doppler_velocity_hh_ms);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_velocity_hv_ms", &myradarfilter->filter_Doppler_velocity_hv_ms);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_velocity_vh_ms", &myradarfilter->filter_Doppler_velocity_vh_ms);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_velocity_vv_ms", &myradarfilter->filter_Doppler_velocity_vv_ms);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_spectralwidth_hh_ms", &myradarfilter->filter_Doppler_spectralwidth_hh_ms);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_spectralwidth_hv_ms", &myradarfilter->filter_Doppler_spectralwidth_hv_ms);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_spectralwidth_vh_ms", &myradarfilter->filter_Doppler_spectralwidth_vh_ms);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_spectralwidth_vv_ms", &myradarfilter->filter_Doppler_spectralwidth_vv_ms);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_spectrum_dBZ_hh", &myradarfilter->filter_Doppler_spectrum_dBZ_hh);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_spectrum_dBZ_hv", &myradarfilter->filter_Doppler_spectrum_dBZ_hv);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_spectrum_dBZ_vh", &myradarfilter->filter_Doppler_spectrum_dBZ_vh);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_Doppler_spectrum_dBZ_vv", &myradarfilter->filter_Doppler_spectrum_dBZ_vv);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_specific_dBZdr", &myradarfilter->filter_specific_dBZdr);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_specific_dBLdr", &myradarfilter->filter_specific_dBLdr);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_specific_rho_co", &myradarfilter->filter_specific_rho_co);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_specific_rho_cxh", &myradarfilter->filter_specific_rho_cxh);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_specific_rho_cxv", &myradarfilter->filter_specific_rho_cxv);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "filter_errors", &myradarfilter->filter_errors);
							
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "inertia_effect", &myradarfilter->inertia_effect);
							
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "effective_earth_correction", &myradarfilter->effective_earth_correction);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "additive_noise", &myradarfilter->additive_noise);
							zephyros_config_read_________i(rwg, rwg->section, "radarfilter", "multiplicative_noise", &myradarfilter->multiplicative_noise);
						}
						
						//wind field
						if (	( strcmp(rwg->subsection,"windfield") == 0 ) |
								( strcmp(rwg->subsection,"prior_windfield") == 0 ) ) 
							{
							if ( strcmp(rwg->section,"simulation") == 0 ) mywindfield = cfg->simulation->windfield;
							if ( strcmp(rwg->section,"retrieval") == 0 ) mywindfield = cfg->retrieval->prior_windfield;
							
							if (strcmp(rwg->subsection,"windfield") == 0 )			mywindfield->type = 0;
							if (strcmp(rwg->subsection,"prior_windfield") == 0 )	mywindfield->type = 1;
							
							//initialize wind field grid
							if (zephyros_config_read_position(rwg, rwg->section, rwg->subsection, "grid")) {
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "grid", &myindices->windfield_grid);
								if ((myindices->windfield_grid + 1) > mywindfield->nfields) {
									mywindfield->nfields = myindices->windfield_grid + 1;
								}
								if (mywindfield->field[myindices->windfield_grid] == NULL) {
									fields_initialize(&(mywindfield->field[myindices->windfield_grid]));
									strcpy(mywindfield->field[myindices->windfield_grid]->name, "windfield");		

								}
							}

							if (mywindfield->field[myindices->windfield_grid] != NULL) {
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "vec_x", &mywindfield->field[myindices->windfield_grid]->n_x, &mywindfield->field[myindices->windfield_grid]->vec_x);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "vec_y", &mywindfield->field[myindices->windfield_grid]->n_y, &mywindfield->field[myindices->windfield_grid]->vec_y);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "vec_z", &mywindfield->field[myindices->windfield_grid]->n_z, &mywindfield->field[myindices->windfield_grid]->vec_z);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "vec_t", &mywindfield->field[myindices->windfield_grid]->n_t, &mywindfield->field[myindices->windfield_grid]->vec_t);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "grid_u", &mywindfield->field[myindices->windfield_grid]->n, &mywindfield->grid_u[myindices->windfield_grid]);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "grid_v", &mywindfield->field[myindices->windfield_grid]->n, &mywindfield->grid_v[myindices->windfield_grid]);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "grid_w", &mywindfield->field[myindices->windfield_grid]->n, &mywindfield->grid_w[myindices->windfield_grid]);
							}
							
							//initialize the wind field wave
							if (zephyros_config_read_position(rwg, rwg->section, rwg->subsection, "wave")) {
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "wave", &myindices->windfield_wave);
								if ((myindices->windfield_wave + 1) > mywindfield->nwaves) {
									mywindfield->nwaves = myindices->windfield_wave + 1;
								}
								if (mywindfield->wave[myindices->windfield_wave] == NULL) {
									mywindfield->wave[myindices->windfield_wave] = malloc(sizeof(t_zephyros_wave));
								}
							}

							if (mywindfield->wave[myindices->windfield_wave] != NULL) {
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "wave_amplitude_msinv", &dummyint, &mywindfield->wave[myindices->windfield_wave]->amplitude_msinv);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "wave_phi0_rad",  &mywindfield->wave[myindices->windfield_wave]->phi0_rad);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "wave_k_minv", &dummyint, &mywindfield->wave[myindices->windfield_wave]->k_minv);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "wave_f_sinv",  &mywindfield->wave[myindices->windfield_wave]->f_sinv);
							}
						
							//initialize the wind field vortex
							if (zephyros_config_read_position(rwg, rwg->section, rwg->subsection, "vortex")) {
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "vortex", &myindices->windfield_vortex);
								if ((myindices->windfield_vortex + 1) > mywindfield->nvortices) {
									mywindfield->nvortices = myindices->windfield_vortex + 1;
								}
								if (mywindfield->vortex[myindices->windfield_vortex] == NULL) {
									mywindfield->vortex[myindices->windfield_vortex] = malloc(sizeof(t_zephyros_vortex));
								}
							}

							if (mywindfield->vortex[myindices->windfield_vortex] != NULL) {
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "vortex_type", &mywindfield->vortex[myindices->windfield_vortex]->type);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "vortex_xyz0_m", &dummyint, &mywindfield->vortex[myindices->windfield_vortex]->xyz0_m);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "vortex_xyz_scaling_m", &dummyint, &mywindfield->vortex[myindices->windfield_vortex]->xyz_scaling_m);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "vortex_rotation_direction", &dummyint, &mywindfield->vortex[myindices->windfield_vortex]->rotation_direction);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "vortex_vmax_msinv",  &mywindfield->vortex[myindices->windfield_vortex]->vmax_msinv);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "vortex_maxr_m",  &mywindfield->vortex[myindices->windfield_vortex]->maxr_m);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "vortex_r_m",  &mywindfield->vortex[myindices->windfield_vortex]->r_m);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "vortex_nu",  &mywindfield->vortex[myindices->windfield_vortex]->nu);
							}

							//initialize the wind field turbulence
							if (zephyros_config_read_position(rwg, rwg->section, rwg->subsection, "turbulence")) {
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence", &myindices->windfield_turbulence);
								if ((myindices->windfield_turbulence + 1) > mywindfield->nturbulences) {
									mywindfield->nturbulences = myindices->windfield_turbulence + 1;
								}
								if (mywindfield->turbulence[myindices->windfield_turbulence] == NULL) {
									mywindfield->turbulence[myindices->windfield_turbulence] = malloc(sizeof(t_zephyros_turbulence_widget));
									fields_initialize(&(mywindfield->turbulence[myindices->windfield_turbulence]->field));
									strcpy(mywindfield->turbulence[myindices->windfield_turbulence]->field->name, "turbulence");		

								}
							}

							if (mywindfield->turbulence[myindices->windfield_turbulence] != NULL) {
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_type", &mywindfield->turbulence[myindices->windfield_turbulence]->type);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "turbulence_vec_x", &mywindfield->turbulence[myindices->windfield_turbulence]->field->n_x, &mywindfield->turbulence[myindices->windfield_turbulence]->field->vec_x);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "turbulence_vec_y", &mywindfield->turbulence[myindices->windfield_turbulence]->field->n_y, &mywindfield->turbulence[myindices->windfield_turbulence]->field->vec_y);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "turbulence_vec_z", &mywindfield->turbulence[myindices->windfield_turbulence]->field->n_z, &mywindfield->turbulence[myindices->windfield_turbulence]->field->vec_z);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "turbulence_vec_t", &mywindfield->turbulence[myindices->windfield_turbulence]->field->n_t, &mywindfield->turbulence[myindices->windfield_turbulence]->field->vec_t);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "turbulence_grid_edr", &mywindfield->turbulence[myindices->windfield_turbulence]->field->n, &mywindfield->turbulence[myindices->windfield_turbulence]->grid_edr);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "turbulence_grid_karman_L", &mywindfield->turbulence[myindices->windfield_turbulence]->field->n, &mywindfield->turbulence[myindices->windfield_turbulence]->grid_karman_L);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "turbulence_grid_kolmogorov_constant", &mywindfield->turbulence[myindices->windfield_turbulence]->field->n, &mywindfield->turbulence[myindices->windfield_turbulence]->grid_kolmogorov_constant);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "turbulence_Lx", &mywindfield->turbulence[myindices->windfield_turbulence]->Lx);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "turbulence_Ly", &mywindfield->turbulence[myindices->windfield_turbulence]->Ly);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "turbulence_Lz", &mywindfield->turbulence[myindices->windfield_turbulence]->Lz);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_Nx", &mywindfield->turbulence[myindices->windfield_turbulence]->Nx);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_Ny", &mywindfield->turbulence[myindices->windfield_turbulence]->Ny);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_Nz", &mywindfield->turbulence[myindices->windfield_turbulence]->Nz);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_K", &mywindfield->turbulence[myindices->windfield_turbulence]->K);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_M", &mywindfield->turbulence[myindices->windfield_turbulence]->M);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "turbulence_lambdax", &mywindfield->turbulence[myindices->windfield_turbulence]->lambdax);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "turbulence_lambday", &mywindfield->turbulence[myindices->windfield_turbulence]->lambday);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "turbulence_minL_div_maxL", &mywindfield->turbulence[myindices->windfield_turbulence]->minL_div_maxL);

								if (zephyros_config_read_position(rwg, rwg->section, rwg->subsection, "turbulence_pinsky2006_file")) {
									strcpy(mywindfield->turbulence[myindices->windfield_turbulence]->pinsky2006_file, this_dir_name);
									zephyros_config_read_________s(rwg, rwg->section, rwg->subsection, "turbulence_pinsky2006_file", this_tmp_char);
									//add filename
									strcat(mywindfield->turbulence[myindices->windfield_turbulence]->pinsky2006_file, this_tmp_char);
								}
										
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_random_numbers", &mywindfield->turbulence[myindices->windfield_turbulence]->random_numbers);

								if (zephyros_config_read_position(rwg, rwg->section, rwg->subsection, "turbulence_rn_file")) {
									strcpy(mywindfield->turbulence[myindices->windfield_turbulence]->rn_file, this_dir_name);
									zephyros_config_read_________s(rwg, rwg->section, rwg->subsection, "turbulence_rn_file", this_tmp_char);
									//add filename
									strcat(mywindfield->turbulence[myindices->windfield_turbulence]->rn_file, this_tmp_char);
								}
										
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_calibration_method", &mywindfield->turbulence[myindices->windfield_turbulence]->calibration_method);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_calibration_n", &mywindfield->turbulence[myindices->windfield_turbulence]->calibration_n);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "turbulence_calibration_C", &mywindfield->turbulence[myindices->windfield_turbulence]->calibration_C);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_calibration_periodic", &mywindfield->turbulence[myindices->windfield_turbulence]->calibration_periodic);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_calibration_nint", &mywindfield->turbulence[myindices->windfield_turbulence]->calibration_nint);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "turbulence_calibration_dir", &mywindfield->turbulence[myindices->windfield_turbulence]->calibration_dir);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "turbulence_calibration_L", &mywindfield->turbulence[myindices->windfield_turbulence]->calibration_L);
							}
						}
						
						
						//scattererfield
						if (	( strcmp(rwg->subsection,"scattererfield") == 0 ) |
								( strcmp(rwg->subsection,"prior_scattererfield") == 0 ) )
								{
							if ( strcmp(rwg->section,"simulation") == 0 ) myscattererfield = cfg->simulation->scattererfield;
							if ( strcmp(rwg->section,"retrieval") == 0 ) myscattererfield = cfg->retrieval->prior_scattererfield;
								
							if (strcmp(rwg->subsection,"scattererfield") == 0 )			myscattererfield->type = 0;
							if (strcmp(rwg->subsection,"prior_scattererfield") == 0 )	myscattererfield->type = 1;
								
							//initialize psd
							if (zephyros_config_read_position(rwg, rwg->section, rwg->subsection, "psd")) {
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "psd", &myindices->scattererfield_psd);
								if ((myindices->scattererfield_psd + 1) > myscattererfield->npsd) {
									myscattererfield->npsd = myindices->scattererfield_psd + 1;
								}
								if (myscattererfield->psd[myindices->scattererfield_psd] == NULL) {
									util_initialize_psd(&(myscattererfield->psd[myindices->scattererfield_psd]));									
								}
							}

							if (myscattererfield->psd[myindices->scattererfield_psd] != NULL) {
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "psd_distribution_type", &myscattererfield->psd[myindices->scattererfield_psd]->distribution_type);
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "psd_particle_type", &myscattererfield->psd[myindices->scattererfield_psd]->particle_type);

								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_vec_x", &myscattererfield->psd[myindices->scattererfield_psd]->field->n_x, &myscattererfield->psd[myindices->scattererfield_psd]->field->vec_x);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_vec_y", &myscattererfield->psd[myindices->scattererfield_psd]->field->n_y, &myscattererfield->psd[myindices->scattererfield_psd]->field->vec_y);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_vec_z", &myscattererfield->psd[myindices->scattererfield_psd]->field->n_z, &myscattererfield->psd[myindices->scattererfield_psd]->field->vec_z);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_vec_t", &myscattererfield->psd[myindices->scattererfield_psd]->field->n_t, &myscattererfield->psd[myindices->scattererfield_psd]->field->vec_t);
								
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_grid_lwc_gm3", &myscattererfield->psd[myindices->scattererfield_psd]->field->n, &myscattererfield->psd[myindices->scattererfield_psd]->grid_lwc_gm3);
								
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_grid_gammadistribution_mu", &myscattererfield->psd[myindices->scattererfield_psd]->field->n, &myscattererfield->psd[myindices->scattererfield_psd]->grid_gammadistribution_mu);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_grid_gammadistribution_D0_mm", &myscattererfield->psd[myindices->scattererfield_psd]->field->n, &myscattererfield->psd[myindices->scattererfield_psd]->grid_gammadistribution_D0_mm);

								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_discrete_D_equiv_mm", &myscattererfield->psd[myindices->scattererfield_psd]->n_diameters, &myscattererfield->psd[myindices->scattererfield_psd]->discrete_D_equiv_mm);
								if (zephyros_config_read_position(rwg, rwg->section, rwg->subsection, "psd_discrete_D_equiv_mm")) {
									myscattererfield->psd[myindices->scattererfield_psd]->grid_number_density_m3 = 
										malloc(myscattererfield->psd[myindices->scattererfield_psd]->n_diameters*sizeof(double*));
									myscattererfield->psd[myindices->scattererfield_psd]->grid_dBnumber_density_err_m3 = 
										malloc(myscattererfield->psd[myindices->scattererfield_psd]->n_diameters*sizeof(double*));
									for (i = 0; i < myscattererfield->psd[myindices->scattererfield_psd]->n_diameters; i++ ) {
										myscattererfield->psd[myindices->scattererfield_psd]->grid_number_density_m3[i] = NULL;
										myscattererfield->psd[myindices->scattererfield_psd]->grid_dBnumber_density_err_m3[i] = NULL;
									}
									myindices->scattererfield_psd_diameter_i = -1;
								}
								
								if (zephyros_config_read_position(rwg, rwg->section, rwg->subsection, "psd_grid_number_density_m3")) {
									myindices->scattererfield_psd_diameter_i += 1;
									zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_grid_number_density_m3", &myscattererfield->psd[myindices->scattererfield_psd]->field->n, &myscattererfield->psd[myindices->scattererfield_psd]->grid_number_density_m3[myindices->scattererfield_psd_diameter_i]);
								}
									zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_grid_dBnumber_density_err_m3", &myscattererfield->psd[myindices->scattererfield_psd]->field->n, &myscattererfield->psd[myindices->scattererfield_psd]->grid_dBnumber_density_err_m3[myindices->scattererfield_psd_diameter_i]);

								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_grid_gammadistribution_N0", &myscattererfield->psd[myindices->scattererfield_psd]->field->n, &myscattererfield->psd[myindices->scattererfield_psd]->grid_gammadistribution_N0);
								zephyros_config_read__lf_array(rwg, rwg->section, rwg->subsection, "psd_grid_gammadistribution_N0_err", &myscattererfield->psd[myindices->scattererfield_psd]->field->n, &myscattererfield->psd[myindices->scattererfield_psd]->grid_gammadistribution_N0_err);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "psd_gammadistribution_mu", &myscattererfield->psd[myindices->scattererfield_psd]->gammadistribution_mu);							
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "psd_gammadistribution_mu_err", &myscattererfield->psd[myindices->scattererfield_psd]->gammadistribution_mu_err);							
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "psd_gammadistribution_D0_mm", &myscattererfield->psd[myindices->scattererfield_psd]->gammadistribution_D0_mm);							
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "psd_gammadistribution_D0_err_mm", &myscattererfield->psd[myindices->scattererfield_psd]->gammadistribution_D0_err_mm);							
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "psd_gammadistribution_dmin_mm", &myscattererfield->psd[myindices->scattererfield_psd]->gammadistribution_dmin_mm);							
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "psd_gammadistribution_dmax_mm", &myscattererfield->psd[myindices->scattererfield_psd]->gammadistribution_dmax_mm);	
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "psd_n_diameters", &myscattererfield->psd[myindices->scattererfield_psd]->n_diameters);	
								
								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "fit_dBlwc", &myscattererfield->psd[myindices->scattererfield_psd]->fit_dBlwc);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "cl_x_m_dBlwc", &myscattererfield->psd[myindices->scattererfield_psd]->dBlwc_ecm.cl_x_m);	
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "cl_y_m_dBlwc", &myscattererfield->psd[myindices->scattererfield_psd]->dBlwc_ecm.cl_y_m);	
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "cl_z_m_dBlwc", &myscattererfield->psd[myindices->scattererfield_psd]->dBlwc_ecm.cl_z_m);	
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "cl_t_s_dBlwc", &myscattererfield->psd[myindices->scattererfield_psd]->dBlwc_ecm.cl_t_s);	
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "c_threshold_dBlwc", &myscattererfield->psd[myindices->scattererfield_psd]->dBlwc_ecm.c_threshold);	

								zephyros_config_read_________i(rwg, rwg->section, rwg->subsection, "fit_dBN", &myscattererfield->psd[myindices->scattererfield_psd]->fit_dBN);
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "cl_x_m_dBN", &myscattererfield->psd[myindices->scattererfield_psd]->dBN_ecm.cl_x_m);	
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "cl_y_m_dBN", &myscattererfield->psd[myindices->scattererfield_psd]->dBN_ecm.cl_y_m);	
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "cl_z_m_dBN", &myscattererfield->psd[myindices->scattererfield_psd]->dBN_ecm.cl_z_m);	
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "cl_t_s_dBN", &myscattererfield->psd[myindices->scattererfield_psd]->dBN_ecm.cl_t_s);	
								zephyros_config_read________lf(rwg, rwg->section, rwg->subsection, "c_threshold_dBN", &myscattererfield->psd[myindices->scattererfield_psd]->dBN_ecm.c_threshold);	

							}
						}
					}


					//simulation
					//if ( strcmp(rwg->section,"simulation") == 0 ) {
					//	zephyros_config_initialize_simulation(cfg);
					//
					//}

					//retrieval
					if ( strcmp(rwg->section,"retrieval") == 0 ) {
						//ensure that retrieval is initialized
						zephyros_config_initialize_retrieval(cfg);
						myindices = &rwg->retrieval_indices;
						
						zephyros_config_read_________i(rwg, "retrieval", "algorithm", "run", &myindices->algorithm_run);
						
						//initialize
						if (zephyros_config_read_position(rwg, "retrieval", "algorithm", "type")) {
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "type", cfg->retrieval->algorithm->type + myindices->algorithm_run);							
							if (cfg->retrieval->algorithm->type[myindices->algorithm_run] == 1) {
								if (cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run] == NULL) {
									cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run] = 
										malloc(sizeof(t_zephyros_config_retrieval_lwm_cfg));
									cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->vd_zcfg = (void*) cfg;	
								}
							}
							if (cfg->retrieval->algorithm->type[myindices->algorithm_run] == 2) {
								if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run] == NULL) {
									cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run] = 
										malloc(sizeof(t_zephyros_config_retrieval_fdvar_cfg));
								}
							}
						}

						//lwm cfg
						if (cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run] != NULL) {
							zephyros_config_read__lf_array(rwg, "retrieval", "algorithm", "xvec_m", &dummyint, &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->xvec_m);
							zephyros_config_read__lf_array(rwg, "retrieval", "algorithm", "yvec_m", &dummyint, &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->yvec_m);
							zephyros_config_read__lf_array(rwg, "retrieval", "algorithm", "zvec_m", &dummyint, &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->zvec_m);
							zephyros_config_read__lf_array(rwg, "retrieval", "algorithm", "tvec_s", &dummyint, &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->tvec_s);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_u0", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_u0);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_u_x", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_u_x);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_u_z", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_u_z);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_v0", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_v0);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_v_y", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_v_y);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_v_z", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_v_z);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_u_y_plus_v_x", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_u_y_plus_v_x);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_w0", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_w0);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_w_x", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_w_x);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_w_y", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_w_y);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_w_z", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_w_z);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "fit_u_t_plus_v_t_plus_w_t", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->fit_u_t_plus_v_t_plus_w_t);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "apply_weights", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->apply_weights);
							zephyros_config_read________lf(rwg, "retrieval", "algorithm", "maximum_time_s", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->maximum_time_s);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "extra_points_n", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->extra_points_n);
							zephyros_config_read________lf(rwg, "retrieval", "algorithm", "extra_points_dx", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->extra_points_dx);
							zephyros_config_read________lf(rwg, "retrieval", "algorithm", "extra_points_dy", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->extra_points_dy);
							zephyros_config_read________lf(rwg, "retrieval", "algorithm", "extra_points_dz", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->extra_points_dz);
							zephyros_config_read________lf(rwg, "retrieval", "algorithm", "extra_points_dt", &cfg->retrieval->algorithm->lwm_cfg[myindices->algorithm_run]->extra_points_dt);
						}

						//fdvar cfg
						if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run] != NULL) {
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "costfunction_dBZ_hh", &cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_dBZ_hh);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "costfunction_dBZdr", &cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_dBZdr);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "costfunction_dBLdr", &cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_dBLdr);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "costfunction_Doppler_velocity_hh_ms", &cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_Doppler_velocity_hh_ms);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "costfunction_Doppler_spectral_width_hh_ms", &cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_Doppler_spectral_width_hh_ms);
							
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "costfunction_Doppler_spectrum_dBZ_hh", &cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_Doppler_spectrum_dBZ_hh);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "costfunction_specific_dBZdr", &cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_specific_dBZdr);
							zephyros_config_read_________i(rwg, "retrieval", "algorithm", "costfunction_specific_dBLdr", &cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_specific_dBLdr);
							
							zephyros_config_read________lf(rwg, "retrieval", "algorithm", "maximum_time_s", &cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->maximum_time_s);
						}
					}
					
					
				}

				fsetpos(rwg->fp, &rwg->pos_nextline);		//set position to next line	
			}
		
			fclose(rwg->fp);
		}
   }
     	
   //open file pointer for additional output
   printf("Opening additional output file: %s\n",additional_output_filename); fflush(stdout);
   cfg->fp_ao = fopen(additional_output_filename, "w");
	if( cfg->fp_ao == NULL )
	{
		perror("Error while opening the additional output file.\n");
		exit(EXIT_FAILURE);
	}


	//make lut for atmosphere
	if (cfg->general->atmosphere->grid_T_K != NULL) {
		fields_prepare(cfg->general->atmosphere->field);
		
		//add small number, so that interpolation of log values goes ok.
		for ( i = 0; i < cfg->general->atmosphere->field->n; i++ ) {		
			cfg->general->atmosphere->grid_pair_hPa[i] += 1.e-5;
			cfg->general->atmosphere->grid_pvapor_hPa[i] += 1.e-5;
		}
		
		util_field2lut(cfg->general->atmosphere->field, cfg->general->atmosphere->grid_T_K, 0, &cfg->general->atmosphere->lut_T_K);
		util_field2lut(cfg->general->atmosphere->field, cfg->general->atmosphere->grid_pair_hPa, 3, &cfg->general->atmosphere->lut_ln_grid_pair_hPa);
		util_field2lut(cfg->general->atmosphere->field, cfg->general->atmosphere->grid_pvapor_hPa, 3, &cfg->general->atmosphere->lut_ln_grid_pvapor_hPa);
	}

	//Make lut for refractive index!
	if (cfg->general->water_refractive_index->wavelength_m != NULL) {
		interpolation_linearinterpolation2lut(
			cfg->general->water_refractive_index->n,
			cfg->general->water_refractive_index->wavelength_m,
			cfg->general->water_refractive_index->realindex,
			0,
			&(cfg->general->water_refractive_index->lut_realindex));
		interpolation_linearinterpolation2lut(
			cfg->general->water_refractive_index->n,
			cfg->general->water_refractive_index->wavelength_m,
			cfg->general->water_refractive_index->imagindex,
			0,
			&(cfg->general->water_refractive_index->lut_imagindex));
	}
	
	//Make lut for white1999_integral!
	if (cfg->general->white1999_integral->vec_ln_a != NULL) {
		interpolation_bilint3D2lut(
			cfg->general->white1999_integral->na,
			cfg->general->white1999_integral->vec_ln_a,
			cfg->general->white1999_integral->nb,
			cfg->general->white1999_integral->vec_ln_b,
			cfg->general->white1999_integral->nL,
			cfg->general->white1999_integral->vec_ln_L,
			cfg->general->white1999_integral->integral_sqrt,
			0,
			&(cfg->general->white1999_integral->lut_integral_sqrt));
	}

#ifdef _ZEPHYROS_CONFIG_DEBUG
	printf("initialize simulation\n"); fflush(stdout);
#endif 

	//initialize simulation
	if (cfg->simulation != NULL) {	
		util_prepare_windfield(cfg->simulation->windfield);
		util_prepare_scattererfield(cfg->simulation->scattererfield);
	}

#ifdef _ZEPHYROS_CONFIG_DEBUG
	printf("initialize retrieval\n"); fflush(stdout);
#endif 

	//initialize retrieval
	if (cfg->retrieval != NULL) {		
		util_prepare_windfield(cfg->retrieval->prior_windfield);
		util_prepare_scattererfield(cfg->retrieval->prior_scattererfield);		

		//logic: if something is fitted, it has to be in the filter.
		for ( i = 0; i <= 100; i++ ) {		
			if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run] != NULL) {
				if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_dBZ_hh) 
					cfg->retrieval->radarfilter->filter_dBZ_hh = 1;
				if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_dBZdr) 
					cfg->retrieval->radarfilter->filter_dBZdr = 1;
				if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_dBLdr) 
					cfg->retrieval->radarfilter->filter_dBLdr = 1;
				if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_Doppler_velocity_hh_ms) 
					cfg->retrieval->radarfilter->filter_Doppler_velocity_hh_ms = 1;
				if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_Doppler_spectral_width_hh_ms) 
					cfg->retrieval->radarfilter->filter_Doppler_spectralwidth_hh_ms = 1;
				if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_Doppler_spectrum_dBZ_hh) 
					cfg->retrieval->radarfilter->filter_Doppler_spectrum_dBZ_hh = 1;
				if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_specific_dBZdr) 
					cfg->retrieval->radarfilter->filter_specific_dBZdr = 1;
				if (cfg->retrieval->algorithm->fdvar_cfg[myindices->algorithm_run]->costfunction_specific_dBLdr) 
					cfg->retrieval->radarfilter->filter_specific_dBLdr = 1;
			}
		}
	}

#ifdef _ZEPHYROS_CONFIG_DEBUG
	printf("derive quantities\n"); fflush(stdout);
#endif 

	//derive quantities
   zephyros_config_derive_quantities_zephyros_config(cfg);

	//print configuration to additional output if aksed for
	if (cfg->general->additional_output->print_configuration) {
		zephyros_config_print(cfg, cfg->fp_ao);
	}

#ifdef _ZEPHYROS_CONFIG_DEBUG
	printf("prepare retrieval\n"); fflush(stdout);
#endif 

	//prepare retrieval
	if (cfg->retrieval != NULL) {		
		util_prepare_post_windfield(&cfg->retrieval->post_windfield, cfg->retrieval->prior_windfield);
		util_prepare_post_scattererfield(&cfg->retrieval->post_scattererfield, cfg->retrieval->prior_scattererfield);
	}
	
#ifdef _ZEPHYROS_CONFIG_DEBUG
	printf("free read widget\n"); fflush(stdout);
#endif 
   //free read widget
   	free(rwg);      
}

int zephyros_config_read_position(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192]) 
{
	int result;
	result = ((strcmp(rwg->section,section) == 0 ) & (strcmp(rwg->subsection,subsection) == 0) & (strcmp(rwg->identifier,identifier) == 0));
	return result;
}

void zephyros_config_read_________i(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192],
	int *destination)
{
	char dummy[8192];
	if (zephyros_config_read_position(rwg, section, subsection, identifier)) {
		fsetpos(rwg->fp, &rwg->pos_thisline);
		fscanf(rwg->fp, "%s %i", dummy, destination);
	}
}

void zephyros_config_read________lf(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192],
	double *destination)
{
	char dummy[8192];
	if (zephyros_config_read_position(rwg, section, subsection, identifier)) {
		fsetpos(rwg->fp, &rwg->pos_thisline);
		fscanf(rwg->fp, "%s %lf", dummy, destination);
	}
}

void zephyros_config_read_________s(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192],
	char *destination)
{
	char dummy[8192];
	if (zephyros_config_read_position(rwg, section, subsection, identifier)) {
		fsetpos(rwg->fp, &rwg->pos_thisline);
		fscanf(rwg->fp, "%s %s", dummy, destination);
	}
}

void zephyros_config_read__lf_array(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192],
	int  *arraysize,
	double **destination)
{
	char dummy[8192];
	int i;
	double *newarray;
	
	if (zephyros_config_read_position(rwg, section, subsection, identifier)) {		
		fsetpos(rwg->fp, &rwg->pos_thisline);
		fscanf(rwg->fp, "%s %i", dummy, arraysize );
		//allocate
		newarray = malloc(*arraysize * sizeof(double));
		for ( i = 0; i < *arraysize; i++ ) {
			fscanf(rwg->fp, "%lf", newarray + i);
		}
		*destination = newarray;
	}
}





void zephyros_config_derive_quantities_zephyros_config(t_zephyros_config *cfg)
{
	double c0 = 299792458.;
	double *dummy;
	double tmp;

	/*
	double Ksq_water = 0.93;
	double Ksq_ice	 = 0.197;
	double SB_k = 1.38e-23;	 //Joules / K
	double radar_resolution_m;
	double radar_maximum_range_m;
	*/
	
	cfg->derived_quantities->central_wavelength_m 	= c0 / cfg->general->instrument->central_frequency_hz;
	cfg->derived_quantities->radar_vmax_ms			= cfg->derived_quantities->central_wavelength_m * cfg->general->instrument->prf_hz / 4.;

	//interpolate water refractive index to radar wavelength
	cfg->derived_quantities->radar_water_refractive_index = 0.;
	interpolation_bilint(
		cfg->general->water_refractive_index->lut_realindex,
		&cfg->derived_quantities->central_wavelength_m,
		&tmp,
		0,
		dummy);

	cfg->derived_quantities->radar_water_refractive_index += tmp;
	
	interpolation_bilint(
		cfg->general->water_refractive_index->lut_imagindex,
		&cfg->derived_quantities->central_wavelength_m,
		&tmp,
		0,
		dummy);
	cfg->derived_quantities->radar_water_refractive_index -= (tmp * 1.j);
}

void zephyros_config_print(t_zephyros_config *cfg, FILE *fp)
{
	double x;
	int i, j;
	
	fprintf(fp, "# configuration file for Zephyros\n");
	fprintf(fp, "# '# ' #-sign + blank space at the beginning of a line means a comment and these lines are skipped\n");
	fprintf(fp, "# comments can also be put in parentheses (comment)\n");
	fprintf(fp, "#\n");
	fprintf(fp, "#configuration read from file: %s\n", cfg->general->overall->file_name);
	fprintf(fp, "#\n");
	fprintf(fp, "section general\n");

	fflush(fp);

	fprintf(fp, "subsection additional_output\n");
	fprintf(fp, "%-30s %-15i\n",
				"print_configuration",
				cfg->general->additional_output->print_configuration);
	fprintf(fp, "%-30s %-15i\n",
				"print_detailed_analysis",
				cfg->general->additional_output->print_detailed_analysis);
	fprintf(fp, "\n");

	fflush(fp);

	//atmosphere
	fprintf(fp, "subsection atmosphere\n");
	fprintf(fp, "%-30s %-15i", 
				"vec_x",
				cfg->general->atmosphere->field->n_x);
	fprintf_array(fp, cfg->general->atmosphere->field->n_x, cfg->general->atmosphere->field->vec_x);
	fprintf(fp, "\n");
	fprintf(fp, "%-30s %-15i", 
				"vec_y",
				cfg->general->atmosphere->field->n_y);
	fprintf_array(fp, cfg->general->atmosphere->field->n_y, cfg->general->atmosphere->field->vec_y);
	fprintf(fp, "\n");
	fprintf(fp, "%-30s %-15i", 
				"vec_z",
				cfg->general->atmosphere->field->n_z);
	fprintf_array(fp, cfg->general->atmosphere->field->n_z, cfg->general->atmosphere->field->vec_z);
	fprintf(fp, "\n");
	fprintf(fp, "%-30s %-15i", 
				"vec_t",
				cfg->general->atmosphere->field->n_t);
	fprintf_array(fp, cfg->general->atmosphere->field->n_t, cfg->general->atmosphere->field->vec_t);
	fprintf(fp, "\n");

	fprintf(fp, "%-30s %-15i", 
				"grid_T_K",
				cfg->general->atmosphere->field->n);
	fprintf_array(fp, cfg->general->atmosphere->field->n, cfg->general->atmosphere->grid_T_K);
	fprintf(fp, "\n");
	fprintf(fp, "%-30s %-15i", 
				"grid_pair_hPa",
				cfg->general->atmosphere->field->n);
	fprintf_array(fp, cfg->general->atmosphere->field->n, cfg->general->atmosphere->grid_pair_hPa);
	fprintf(fp, "\n");
	fprintf(fp, "%-30s %-15i", 
				"grid_pvapor_hPa",
				cfg->general->atmosphere->field->n);
	fprintf_array(fp, cfg->general->atmosphere->field->n, cfg->general->atmosphere->grid_pvapor_hPa);
	fprintf(fp, "\n");			
		
	fflush(fp);

	fprintf(fp, "subsection instrument\n");
	fprintf(fp, "%-30s %-15i\n",
				"type",
				cfg->general->instrument->type);
	fprintf(fp, "%-30s %-15.3e\n",
			"transmit_power_watt",
			cfg->general->instrument->transmit_power_watt);
	fprintf(fp, "%-30s %-15.3e\n",
			"k_nonuniformity",
			cfg->general->instrument->k_nonuniformity);
	fprintf(fp, "%-30s %-15.3e\n",
			"power_gain",
			cfg->general->instrument->power_gain);
	fprintf(fp, "%-30s %-15.3e\n",
			"lr_receiver_loss_factor",
			cfg->general->instrument->lr_receiver_loss_factor);
	fprintf(fp, "%-30s %-15.3e\n",
			"central_frequency_hz",
			cfg->general->instrument->central_frequency_hz);
	fprintf(fp, "%-30s %-15.3e\n",
			"prf_hz",
			cfg->general->instrument->prf_hz);
	fprintf(fp, "%-30s %-15.3e\n",
			"temperature_K",
			cfg->general->instrument->temperature_K);
		
	fflush(fp);
				
	fprintf(fp, "subsection overall\n");
	fprintf(fp, "%-30s%-10s\n",
				"version_number",
				zephyros_version);
	fprintf(fp, "\n");
					
	fprintf(fp, "subsection water_refractive_index\n");
	if (cfg->general->water_refractive_index->wavelength_m != NULL) {
		fprintf(fp, "%-30s %-15i", 
					"wavelength_m",
					cfg->general->water_refractive_index->n);
		fprintf_array(fp, cfg->general->water_refractive_index->n, cfg->general->water_refractive_index->wavelength_m);
		fprintf(fp, "\n");
		}
	if (cfg->general->water_refractive_index->realindex != NULL) {
		fprintf(fp, "%-30s %-15i", 
					"realindex",
					cfg->general->water_refractive_index->n);
		fprintf_array(fp, cfg->general->water_refractive_index->n, cfg->general->water_refractive_index->realindex);
		fprintf(fp, "\n");
	}
	if (cfg->general->water_refractive_index->imagindex != NULL) {
		fprintf(fp, "%-30s %-15i", 
					"imagindex",
					cfg->general->water_refractive_index->n);
		fprintf_array(fp, cfg->general->water_refractive_index->n, cfg->general->water_refractive_index->imagindex);
		fprintf(fp, "\n");
	}			
	fprintf(fp, "\n");

	fflush(fp);

	fprintf(fp, "subsection white1999_integral\n");
	if (cfg->general->white1999_integral->vec_ln_a != NULL) {
		fprintf(fp, "%-30s %-15i", 
					"vec_ln_a",
					cfg->general->white1999_integral->na);
		fprintf_array(fp, cfg->general->white1999_integral->na, cfg->general->white1999_integral->vec_ln_a);
		fprintf(fp, "\n");
	}
	if (cfg->general->white1999_integral->vec_ln_b != NULL) {
		fprintf(fp, "%-30s %-15i", 
					"vec_ln_b",
					cfg->general->white1999_integral->nb);
		fprintf_array(fp, cfg->general->white1999_integral->nb, cfg->general->white1999_integral->vec_ln_b);
		fprintf(fp, "\n");
	}
	if (cfg->general->white1999_integral->vec_ln_L != NULL) {
		fprintf(fp, "%-30s %-15i", 
					"vec_ln_L",
					cfg->general->white1999_integral->nL);
		fprintf_array(fp, cfg->general->white1999_integral->nL, cfg->general->white1999_integral->vec_ln_L);
		fprintf(fp, "\n");
	}
	if (cfg->general->white1999_integral->integral_sqrt != NULL) {
		fprintf(fp, "%-30s %-15i", 
					"integral_sqrt",
					cfg->general->white1999_integral->n);
		fprintf_array(fp, cfg->general->white1999_integral->n, cfg->general->white1999_integral->integral_sqrt);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	fflush(fp);

	if (cfg->simulation != NULL) {
		fprintf(fp, "section simulation\n");
		fprintf(fp, "subsection scattererfield\n");
		zephyros_config_print_scattererfield(cfg->simulation->scattererfield, fp);
		
		fprintf(fp, "subsection windfield\n");
		zephyros_config_print_windfield(cfg->simulation->windfield, fp);
		
		zephyros_config_print_radarfilter(cfg->simulation->radarfilter, fp);
	}



	if (cfg->retrieval != NULL) {
		fprintf(fp, "section retrieval\n");
		
		fprintf(fp, "subsection prior_scattererfield\n");
		zephyros_config_print_scattererfield(cfg->retrieval->prior_scattererfield, fp);

		fprintf(fp, "subsection prior_windfield\n");
		zephyros_config_print_windfield(cfg->retrieval->prior_windfield, fp);
		
		
		zephyros_config_print_radarfilter(cfg->retrieval->radarfilter, fp);

		fprintf(fp, "subsection algorithm\n");
		for ( i = 0; i <= 100; i++ ) {
			if (cfg->retrieval->algorithm->type[i] == 1) {
				fprintf(fp, "%-30s %-15i\n",
					"run",
					i);
				fprintf(fp, "%-30s %-15i\n",
					"type",
					cfg->retrieval->algorithm->type[i]);
				fprintf(fp, "%-30s %-15i %-15.3e %-15.3e %-15.3e\n", 
							"xvec_m", 3,
							cfg->retrieval->algorithm->lwm_cfg[i]->xvec_m[0],
							cfg->retrieval->algorithm->lwm_cfg[i]->xvec_m[1],
							cfg->retrieval->algorithm->lwm_cfg[i]->xvec_m[2]);
				fprintf(fp, "%-30s %-15i %-15.3e %-15.3e %-15.3e\n",
							"yvec_m", 3,
							cfg->retrieval->algorithm->lwm_cfg[i]->yvec_m[0],
							cfg->retrieval->algorithm->lwm_cfg[i]->yvec_m[1],
							cfg->retrieval->algorithm->lwm_cfg[i]->yvec_m[2]);
				fprintf(fp, "%-30s %-15i %-15.3e %-15.3e %-15.3e\n",
							"zvec_m", 3,
							cfg->retrieval->algorithm->lwm_cfg[i]->zvec_m[0],
							cfg->retrieval->algorithm->lwm_cfg[i]->zvec_m[1],
							cfg->retrieval->algorithm->lwm_cfg[i]->zvec_m[2]);
				fprintf(fp, "%-30s %-15i %-15.3e %-15.3e %-15.3e\n",
							"tvec_s", 3,
							cfg->retrieval->algorithm->lwm_cfg[i]->tvec_s[0],
							cfg->retrieval->algorithm->lwm_cfg[i]->tvec_s[1],
							cfg->retrieval->algorithm->lwm_cfg[i]->tvec_s[2]);
				fprintf(fp, "%-30s %-15i\n",
					"fit_u0",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_u0);
				fprintf(fp, "%-30s %-15i\n",
					"fit_u_x",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_u_x);
				fprintf(fp, "%-30s %-15i\n",
					"fit_u_z",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_u_z);
				fprintf(fp, "%-30s %-15i\n",
					"fit_v0",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_v0);
				fprintf(fp, "%-30s %-15i\n",
					"fit_v_y",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_v_y);
				fprintf(fp, "%-30s %-15i\n",
					"fit_v_z",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_v_z);
				fprintf(fp, "%-30s %-15i\n",
					"fit_u_y_plus_v_x",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_u_y_plus_v_x);
				fprintf(fp, "%-30s %-15i\n",
					"fit_w0",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_w0);
				fprintf(fp, "%-30s %-15i\n",
					"fit_w_x",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_w_x);
				fprintf(fp, "%-30s %-15i\n",
					"fit_w_y",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_w_y);
				fprintf(fp, "%-30s %-15i\n",
					"fit_w_z",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_w_z);
				fprintf(fp, "%-30s %-15i\n",
					"fit_u_t_plus_v_t_plus_w_t",
					cfg->retrieval->algorithm->lwm_cfg[i]->fit_u_t_plus_v_t_plus_w_t);
				fprintf(fp, "%-30s %-15i\n",
					"apply_weights",
					cfg->retrieval->algorithm->lwm_cfg[i]->apply_weights);
				fprintf(fp, "%-30s %-15.3e\n",
					"maximum_time_s",
					cfg->retrieval->algorithm->lwm_cfg[i]->maximum_time_s);
				fprintf(fp, "%-30s %-15i\n",
					"extra_points_n",
					cfg->retrieval->algorithm->lwm_cfg[i]->extra_points_n);					
				fprintf(fp, "%-30s %-15.3e\n",
					"extra_points_dx",
					cfg->retrieval->algorithm->lwm_cfg[i]->extra_points_dx);
				fprintf(fp, "%-30s %-15.3e\n",
					"extra_points_dy",
					cfg->retrieval->algorithm->lwm_cfg[i]->extra_points_dy);
				fprintf(fp, "%-30s %-15.3e\n",
					"extra_points_dz",
					cfg->retrieval->algorithm->lwm_cfg[i]->extra_points_dz);
				fprintf(fp, "%-30s %-15.3e\n",
					"extra_points_dt",
					cfg->retrieval->algorithm->lwm_cfg[i]->extra_points_dt);					
			}

			if (cfg->retrieval->algorithm->type[i] == 2) {
				fprintf(fp, "%-30s %-15i\n",
					"run",
					i);
				fprintf(fp, "%-30s %-15i\n",
					"type",
					cfg->retrieval->algorithm->type[i]);
				fprintf(fp, "%-30s %-15i\n",
							"costfunction_dBZ_hh",
							cfg->retrieval->algorithm->fdvar_cfg[i]->costfunction_dBZ_hh);
				fprintf(fp, "%-30s %-15i\n",
							"costfunction_dBZdr",
							cfg->retrieval->algorithm->fdvar_cfg[i]->costfunction_dBZdr);
				fprintf(fp, "%-30s %-15i\n",
							"costfunction_dBLdr",
							cfg->retrieval->algorithm->fdvar_cfg[i]->costfunction_dBLdr);
				fprintf(fp, "%-30s %-15i\n",
							"costfunction_Doppler_velocity_hh_ms",
							cfg->retrieval->algorithm->fdvar_cfg[i]->costfunction_Doppler_velocity_hh_ms);
				fprintf(fp, "%-30s %-15i\n",
							"costfunction_Doppler_spectral_width_hh_ms",
							cfg->retrieval->algorithm->fdvar_cfg[i]->costfunction_Doppler_spectral_width_hh_ms);
				fprintf(fp, "%-30s %-15i\n",
							"costfunction_Doppler_spectrum_dBZ_hh",
							cfg->retrieval->algorithm->fdvar_cfg[i]->costfunction_Doppler_spectrum_dBZ_hh);
				fprintf(fp, "%-30s %-15i\n",
							"costfunction_specific_dBZdr",
							cfg->retrieval->algorithm->fdvar_cfg[i]->costfunction_specific_dBZdr);
				fprintf(fp, "%-30s %-15i\n",
							"costfunction_specific_dBLdr",
							cfg->retrieval->algorithm->fdvar_cfg[i]->costfunction_specific_dBLdr);
							
							
				fprintf(fp, "%-30s %-15.3e\n",
							"maximum_time_s",
							cfg->retrieval->algorithm->fdvar_cfg[i]->maximum_time_s);
			}
		}
		fflush(fp);
	}



	fprintf(fp, "derived quantities\n");
	fprintf(fp, "central_wavelength_m            %-15.3e\n", cfg->derived_quantities->central_wavelength_m);
	fprintf(fp, "radar_vmax_ms                   %-15.3e\n", cfg->derived_quantities->radar_vmax_ms);
	fprintf(fp, "radar_water_refractive_index_real    %-15.3e\n", creal(cfg->derived_quantities->radar_water_refractive_index));
	fprintf(fp, "radar_water_refractive_index_imag    %-15.3e\n", cimag(cfg->derived_quantities->radar_water_refractive_index));
	//fprintf(fp, "noise_power                     %-15.3e\n", cfg->derived_quantities->noise_power);
	fprintf(fp, "\n\n\n");

	
	fprintf(fp, "\n\n\n");

}



void zephyros_config_print_radarfilter(t_zephyros_radarfilter *myradarfilter, FILE *fp)
{
	fprintf(fp, "subsection radarfilter\n");

	fprintf(fp, "%-30s %-15i\n",
			"n_beam_range",
			myradarfilter->n_beam_range);
	fprintf(fp, "%-30s %-15i\n",
			"n_beam_theta",
			myradarfilter->n_beam_theta);
	fprintf(fp, "%-30s %-15i\n",
			"n_beam_phi",
			myradarfilter->n_beam_phi);
	fprintf(fp, "%-30s %-15i\n",
			"n_t",
			myradarfilter->n_t );
	fprintf(fp, "%-30s %-15i\n",
			"n_parmod_az",
			myradarfilter->n_parmod_az );
	fprintf(fp, "%-30s %-15i\n",
			"n_parmod_el",
			myradarfilter->n_parmod_el );
	fprintf(fp, "%-30s %-15i\n",
			"n_spectrum",
			myradarfilter->n_spectrum );
	fprintf(fp, "\n");
		
	fprintf(fp, "%-30s %-15i\n", "filter_dBZ_hh", myradarfilter->filter_dBZ_hh);
	fprintf(fp, "%-30s %-15i\n", "filter_dBZ_hv", myradarfilter->filter_dBZ_hv);
	fprintf(fp, "%-30s %-15i\n", "filter_dBZ_vh", myradarfilter->filter_dBZ_vh);
	fprintf(fp, "%-30s %-15i\n", "filter_dBZ_vv", myradarfilter->filter_dBZ_vv);
	fprintf(fp, "%-30s %-15i\n", "filter_dBZdr", myradarfilter->filter_dBZdr);
	fprintf(fp, "%-30s %-15i\n", "filter_dBLdr", myradarfilter->filter_dBLdr);
	fprintf(fp, "%-30s %-15i\n", "filter_rho_co", myradarfilter->filter_rho_co);
	fprintf(fp, "%-30s %-15i\n", "filter_rho_cxh", myradarfilter->filter_rho_cxh);
	fprintf(fp, "%-30s %-15i\n", "filter_rho_cxv", myradarfilter->filter_rho_cxv);
	fprintf(fp, "%-30s %-15i\n", "filter_KDP", myradarfilter->filter_KDP);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_velocity_hh_ms", myradarfilter->filter_Doppler_velocity_hh_ms);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_velocity_hv_ms", myradarfilter->filter_Doppler_velocity_hv_ms);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_velocity_vh_ms", myradarfilter->filter_Doppler_velocity_vh_ms);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_velocity_vv_ms", myradarfilter->filter_Doppler_velocity_vv_ms);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_spectralwidth_hh_ms", myradarfilter->filter_Doppler_spectralwidth_hh_ms);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_spectralwidth_hv_ms", myradarfilter->filter_Doppler_spectralwidth_hv_ms);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_spectralwidth_vh_ms", myradarfilter->filter_Doppler_spectralwidth_vh_ms);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_spectralwidth_vv_ms", myradarfilter->filter_Doppler_spectralwidth_vv_ms);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_spectrum_dBZ_hh", myradarfilter->filter_Doppler_spectrum_dBZ_hh);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_spectrum_dBZ_hv", myradarfilter->filter_Doppler_spectrum_dBZ_hv);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_spectrum_dBZ_vh", myradarfilter->filter_Doppler_spectrum_dBZ_vh);
	fprintf(fp, "%-30s %-15i\n", "filter_Doppler_spectrum_dBZ_vv", myradarfilter->filter_Doppler_spectrum_dBZ_vv);
	fprintf(fp, "%-30s %-15i\n", "filter_specific_dBZdr", myradarfilter->filter_specific_dBZdr);
	fprintf(fp, "%-30s %-15i\n", "filter_specific_dBLdr", myradarfilter->filter_specific_dBLdr);
	fprintf(fp, "%-30s %-15i\n", "filter_specific_rho_co", myradarfilter->filter_specific_rho_co);
	fprintf(fp, "%-30s %-15i\n", "filter_specific_rho_cxh", myradarfilter->filter_specific_rho_cxh);
	fprintf(fp, "%-30s %-15i\n", "filter_specific_rho_cxv", myradarfilter->filter_specific_rho_cxv);
	fprintf(fp, "%-30s %-15i\n", "filter_errors", myradarfilter->filter_errors);
	
	fprintf(fp, "%-30s %-15i\n",
			"additive_noise",
			myradarfilter->additive_noise);
	fprintf(fp, "%-30s %-15i\n",
			"multiplicative_noise",
			myradarfilter->multiplicative_noise);
	fprintf(fp, "%-30s %-15i\n",
			"inertia_effect",
			myradarfilter->inertia_effect);
	fprintf(fp, "%-30s %-15i\n",
			"effective_earth_correction",
			myradarfilter->effective_earth_correction);
	fprintf(fp, "\n");
	
	
	fprintf(fp, "\n\n\n");

	fflush(fp);
}

void zephyros_config_print_windfield(t_zephyros_windfield *mywindfield, FILE *fp)
{
	int i, j;
	
	for ( i = 0; i <= 100; i++ ) {
		if (mywindfield->field[i] != NULL) {
			fprintf(fp, "%-30s %-15i\n", 
					"grid",
					i);
		
			fprintf(fp, "%-30s %-15i", 
						"vec_x",
						mywindfield->field[i]->n_x);
			fprintf_array(fp, mywindfield->field[i]->n_x, mywindfield->field[i]->vec_x);
			fprintf(fp, "\n");

			fprintf(fp, "%-30s %-15i", 
						"vec_y",
						mywindfield->field[i]->n_y);
			fprintf_array(fp, mywindfield->field[i]->n_y, mywindfield->field[i]->vec_y);
			fprintf(fp, "\n");
			
			fprintf(fp, "%-30s %-15i", 
						"vec_z",
						mywindfield->field[i]->n_z);
			fprintf_array(fp, mywindfield->field[i]->n_z, mywindfield->field[i]->vec_z);
			fprintf(fp, "\n");
			
			fprintf(fp, "%-30s %-15i", 
						"vec_t",
						mywindfield->field[i]->n_t);
			fprintf_array(fp, mywindfield->field[i]->n_t, mywindfield->field[i]->vec_t);
			fprintf(fp, "\n");
			
			fprintf(fp, "%-30s %-15i", 
						"grid_u",
						mywindfield->field[i]->n);
			fprintf_array(fp, mywindfield->field[i]->n, mywindfield->grid_u[i]);
			fprintf(fp, "\n");
			
			fprintf(fp, "%-30s %-15i", 
						"grid_v",
						mywindfield->field[i]->n);
			fprintf_array(fp, mywindfield->field[i]->n, mywindfield->grid_v[i]);
			fprintf(fp, "\n");
			
			fprintf(fp, "%-30s %-15i", 
						"grid_w",
						mywindfield->field[i]->n);
			fprintf_array(fp, mywindfield->field[i]->n, mywindfield->grid_w[i]);
			fprintf(fp, "\n");
		}
	}

	for ( i = 0; i <= 100; i++ ) {
		if (mywindfield->wave[i] != NULL) {	
			fprintf(fp, "%-30s %-15i\n", 
					"wave",
					i);

			fprintf(fp, "%-30s %-15i", 
						"wave_amplitude_msinv",
						3);
			fprintf_array(fp, 3, mywindfield->wave[i]->amplitude_msinv);
			fprintf(fp, "\n");

			fprintf(fp, "%-30s %-15.3e\n",
					"wave_phi0_rad",
					mywindfield->wave[i]->phi0_rad);

			fprintf(fp, "%-30s %-15i", 
						"wave_k_minv",
						3);
			fprintf_array(fp, 3, mywindfield->wave[i]->k_minv);
			fprintf(fp, "\n");

			fprintf(fp, "%-30s %-15.3e\n",
					"wave_f_sinv",
					mywindfield->wave[i]->f_sinv);
		}
	}

	for ( i = 0; i <= 100; i++ ) {
		if (mywindfield->vortex[i] != NULL) {
			fprintf(fp, "%-30s %-15i\n", 
					"vortex",
					i);
			fprintf(fp, "%-30s %-15i\n", 
					"vortex_type",
					mywindfield->vortex[i]->type);

			fprintf(fp, "%-30s %-15i", 
						"vortex_xyz0_m",
						3);
			fprintf_array(fp, 3, mywindfield->vortex[i]->xyz0_m);
			fprintf(fp, "\n");

			fprintf(fp, "%-30s %-15i", 
						"vortex_xyz_scaling_m",
						3);
			fprintf_array(fp, 3, mywindfield->vortex[i]->xyz_scaling_m);
			fprintf(fp, "\n");

			fprintf(fp, "%-30s %-15i", 
						"vortex_rotation_direction",
						3);
			fprintf_array(fp, 3, mywindfield->vortex[i]->rotation_direction);
			fprintf(fp, "\n");

			fprintf(fp, "%-30s %-15.3e\n",
					"vortex_vmax_msinv",
					mywindfield->vortex[i]->vmax_msinv);
			fprintf(fp, "%-30s %-15.3e\n",
					"vortex_maxr_m",
					mywindfield->vortex[i]->maxr_m);
			if (mywindfield->vortex[i]->type == 1) {
				fprintf(fp, "%-30s %-15.3e\n",
						"vortex_r_m",
						mywindfield->vortex[i]->r_m);
			}
	
			if (mywindfield->vortex[i]->type == 2) {
				fprintf(fp, "%-30s %-15.3e\n",
						"vortex_nu",
						mywindfield->vortex[i]->nu);
				fprintf(fp, "%-30s %-15.3e\n",
						"vortex_r_m",
						mywindfield->vortex[i]->r_m);
			}
		}
	}


	for ( i = 0; i <= 100; i++ ) {
		if (mywindfield->turbulence[i] != NULL) {				
			fprintf(fp, "%-30s %-15i\n", 
					"turbulence",
					i);
					
			fprintf(fp, "%-30s %-15i\n",
					"turbulence_type",
					mywindfield->turbulence[i]->type);
					
			fprintf(fp, "%-30s %-15i", 
						"turbulence_vec_x",
						mywindfield->turbulence[i]->field->n_x);
			for (j=0; j<mywindfield->turbulence[i]->field->n_x; j++) {
				fprintf(fp, " %-15.3e", mywindfield->turbulence[i]->field->vec_x[j]);					
			}
			fprintf(fp, "\n");

			fprintf(fp, "%-30s %-15i", 
						"turbulence_vec_y",
						mywindfield->turbulence[i]->field->n_y);
			for (j=0; j<mywindfield->turbulence[i]->field->n_y; j++) {
				fprintf(fp, " %-15.3e", mywindfield->turbulence[i]->field->vec_y[j]);					
			}
			fprintf(fp, "\n");
			
			fprintf(fp, "%-30s %-15i", 
						"turbulence_vec_z",
						mywindfield->turbulence[i]->field->n_z);
			for (j=0; j<mywindfield->turbulence[i]->field->n_z; j++) {
				fprintf(fp, " %-15.3e", mywindfield->turbulence[i]->field->vec_z[j]);					
			}
			fprintf(fp, "\n");
			
			fprintf(fp, "%-30s %-15i", 
						"turbulence_vec_t",
						mywindfield->turbulence[i]->field->n_t);
			for (j=0; j<mywindfield->turbulence[i]->field->n_t; j++) {
				fprintf(fp, " %-15.3e", mywindfield->turbulence[i]->field->vec_t[j]);					
			}
			fprintf(fp, "\n");
			
			fprintf(fp, "%-30s %-15i", 
						"turbulence_grid_edr",
						mywindfield->turbulence[i]->field->n);
			for (j=0; j<mywindfield->turbulence[i]->field->n; j++) {
				fprintf(fp, " %-15.3e", mywindfield->turbulence[i]->grid_edr[j]);					
			}
			fprintf(fp, "\n");
			
			if ((mywindfield->turbulence[i]->type == 1) |  (mywindfield->turbulence[i]->type == 2)) {
				fprintf(fp, "%-30s %-15i", 
							"turbulence_grid_karman_L",
							mywindfield->turbulence[i]->field->n);
				for (j=0; j<mywindfield->turbulence[i]->field->n; j++) {
					fprintf(fp, " %-15.3e", mywindfield->turbulence[i]->grid_karman_L[j]);					
				}
				fprintf(fp, "\n");
			}
			
			fprintf(fp, "%-30s %-15i", 
						"turbulence_grid_kolmogorov_constant",
						mywindfield->turbulence[i]->field->n);
			for (j=0; j<mywindfield->turbulence[i]->field->n; j++) {
				fprintf(fp, " %-15.3e", mywindfield->turbulence[i]->grid_kolmogorov_constant[j]);					
			}
			fprintf(fp, "\n");

			if (
				(mywindfield->turbulence[i]->type == 1) |
				(mywindfield->turbulence[i]->type == 2) |
				(mywindfield->turbulence[i]->type == 3)) {
				fprintf(fp, "%-30s %-15.3e\n",
						"turbulence_Lx",
						mywindfield->turbulence[i]->Lx);
				fprintf(fp, "%-30s %-15.3e\n",
						"turbulence_Ly",
						mywindfield->turbulence[i]->Ly);
			}
			if ((mywindfield->turbulence[i]->type == 1) |  (mywindfield->turbulence[i]->type == 2)) {
				fprintf(fp, "%-30s %-15.3e\n",
						"turbulence_Lz",
						mywindfield->turbulence[i]->Lz);
			}
			
			if (
				(mywindfield->turbulence[i]->type == 1) |
				(mywindfield->turbulence[i]->type == 2) |
				(mywindfield->turbulence[i]->type == 3)) {
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_Nx",
						mywindfield->turbulence[i]->Nx);
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_Ny",
						mywindfield->turbulence[i]->Ny);
			}
			if ((mywindfield->turbulence[i]->type == 1) |  (mywindfield->turbulence[i]->type == 2)) {
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_Nz",
						mywindfield->turbulence[i]->Nz);
			}
			
			if (mywindfield->turbulence[i]->type == 4) {
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_K",
						mywindfield->turbulence[i]->K);
			}
			if (mywindfield->turbulence[i]->type == 4) {
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_M",
						mywindfield->turbulence[i]->M);
			}
			if (mywindfield->turbulence[i]->type == 3) {				
			fprintf(fp, "%-30s %-15.3e\n",
					"turbulence_lambdax",
					mywindfield->turbulence[i]->lambdax);
			fprintf(fp, "%-30s %-15.3e\n",
					"turbulence_lambday",
					mywindfield->turbulence[i]->lambday);
			}
			if (mywindfield->turbulence[i]->type == 4) {				
				fprintf(fp, "%-30s %-15s\n",
						"turbulence_pinsky2006_file",
						mywindfield->turbulence[i]->pinsky2006_file);						
			}
			
			if (mywindfield->turbulence[i]->type == 2) {	
				fprintf(fp, "%-30s %-15.3e\n",
						"turbulence_minL_div_maxL",
						mywindfield->turbulence[i]->minL_div_maxL);
			}
			
			if ((1 <= mywindfield->turbulence[i]->type) 
					& (mywindfield->turbulence[i]->type <= 4)) {	
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_random_numbers",
						mywindfield->turbulence[i]->random_numbers);
				fprintf(fp, "%-30s %-15s\n",
						"turbulence_rn_file",
						mywindfield->turbulence[i]->rn_file);

				fprintf(fp, "%-30s %-15i\n",
						"turbulence_calibration_method",
						mywindfield->turbulence[i]->calibration_method);
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_calibration_n",
						mywindfield->turbulence[i]->calibration_n);
				fprintf(fp, "%-30s %-15.3e\n",
						"turbulence_calibration_C",
						mywindfield->turbulence[i]->calibration_C);
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_calibration_periodic",
						mywindfield->turbulence[i]->calibration_periodic);
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_calibration_nint",
						mywindfield->turbulence[i]->calibration_nint);
				fprintf(fp, "%-30s %-15i\n",
						"turbulence_calibration_dir",
						mywindfield->turbulence[i]->calibration_dir);
				fprintf(fp, "%-30s %-15.3e\n",
						"turbulence_calibration_L",
						mywindfield->turbulence[i]->calibration_L);
			}
		}
	}
	
	fprintf(fp, "\n\n\n");

	fflush(fp);
}

void zephyros_config_print_scattererfield(t_zephyros_scattererfield *myscattererfield, FILE *fp)
{
	int i, j;
	
	
	for ( i = 0; i <= 100; i++ ) {
		if (myscattererfield->psd[i] != NULL) {		
			fprintf(fp, "%-30s %-15i\n", 
					"psd",
					i);
			fprintf(fp, "%-30s %-15i\n",
					"psd_distribution_type",
					myscattererfield->psd[i]->distribution_type);
			fprintf(fp, "%-30s %-15i\n",
					"psd_particle_type",
					myscattererfield->psd[i]->particle_type);

			fprintf(fp, "%-30s %-15i", 
						"vec_x",
						myscattererfield->psd[i]->field->n_x);
			fprintf_array(fp, myscattererfield->psd[i]->field->n_x, myscattererfield->psd[i]->field->vec_x);
			fprintf(fp, "\n");

			fprintf(fp, "%-30s %-15i", 
						"vec_y",
						myscattererfield->psd[i]->field->n_y);
			fprintf_array(fp, myscattererfield->psd[i]->field->n_y, myscattererfield->psd[i]->field->vec_y);
			fprintf(fp, "\n");
			
			fprintf(fp, "%-30s %-15i", 
						"vec_z",
						myscattererfield->psd[i]->field->n_z);
			fprintf_array(fp, myscattererfield->psd[i]->field->n_z, myscattererfield->psd[i]->field->vec_z);
			fprintf(fp, "\n");

			fprintf(fp, "%-30s %-15i", 
						"vec_t",
						myscattererfield->psd[i]->field->n_t);
			fprintf_array(fp, myscattererfield->psd[i]->field->n_t, myscattererfield->psd[i]->field->vec_t);
			fprintf(fp, "\n");
			
			if (myscattererfield->psd[i]->grid_lwc_gm3 != NULL) {
				fprintf(fp, "%-30s %-15i", 
							"psd_grid_lwc_gm3",
							myscattererfield->psd[i]->field->n);
				fprintf_array(fp, myscattererfield->psd[i]->field->n, myscattererfield->psd[i]->grid_lwc_gm3);
				fprintf(fp, "\n");
			}
			if (myscattererfield->psd[i]->grid_dBlwc_err_gm3 != NULL) {
				fprintf(fp, "%-30s %-15i", 
							"psd_grid_dBlwc_err_gm3",
							myscattererfield->psd[i]->field->n);
				fprintf_array(fp, myscattererfield->psd[i]->field->n, myscattererfield->psd[i]->grid_dBlwc_err_gm3);
				fprintf(fp, "\n");
			}
			
			if (myscattererfield->psd[i]->distribution_type == 1) {
				fprintf(fp, "%-30s %-15i", 
							"psd_grid_gammadistribution_N0",
							myscattererfield->psd[i]->field->n);
				fprintf_array(fp, myscattererfield->psd[i]->field->n, myscattererfield->psd[i]->grid_gammadistribution_N0);
				fprintf(fp, "\n");
				

				if (myscattererfield->type == 1) {
					if (myscattererfield->psd[i]->grid_gammadistribution_N0_err != NULL) {
						fprintf(fp, "%-30s %-15i", 
									"psd_grid_gammadistribution_N0_err",
									myscattererfield->psd[i]->field->n);
						fprintf_array(fp, myscattererfield->psd[i]->field->n, myscattererfield->psd[i]->grid_gammadistribution_N0_err);
						fprintf(fp, "\n");
					}
				}
				
				fprintf(fp, "%-30s %-15.3e\n",
						"psd_gammadistribution_mu",
						myscattererfield->psd[i]->gammadistribution_mu);
				if (myscattererfield->type == 1) {
					fprintf(fp, "%-30s %-15.3e\n",
							"psd_gammadistribution_mu_err",
							myscattererfield->psd[i]->gammadistribution_mu_err);
				}
				
				fprintf(fp, "%-30s %-15.3e\n",
						"psd_gammadistribution_D0_mm",
						myscattererfield->psd[i]->gammadistribution_D0_mm);
				if (myscattererfield->type == 1) {
					fprintf(fp, "%-30s %-15.3e\n",
							"psd_gammadistribution_D0_err_mm",
							myscattererfield->psd[i]->gammadistribution_D0_err_mm);
				}
				 
				fprintf(fp, "%-30s %-15.3e\n",
						"psd_gammadistribution_dmin_mm",
						myscattererfield->psd[i]->gammadistribution_dmin_mm);
				fprintf(fp, "%-30s %-15.3e\n",
						"psd_gammadistribution_dmax_mm",
						myscattererfield->psd[i]->gammadistribution_dmax_mm);
						
				fprintf(fp, "#\n");

				if (myscattererfield->psd[i]->grid_gammadistribution_mu != NULL) {
					fprintf(fp, "%-30s %-15i", 
								"grid_gammadistribution_mu",
								myscattererfield->psd[i]->field->n);
					fprintf_array(fp, myscattererfield->psd[i]->field->n, myscattererfield->psd[i]->grid_gammadistribution_mu);
					fprintf(fp, "\n");
				}				
				if (myscattererfield->psd[i]->grid_gammadistribution_D0_mm != NULL) {
					fprintf(fp, "%-30s %-15i", 
								"grid_gammadistribution_D0_mm",
								myscattererfield->psd[i]->field->n);
					fprintf_array(fp, myscattererfield->psd[i]->field->n, myscattererfield->psd[i]->grid_gammadistribution_D0_mm);
					fprintf(fp, "\n");
				}				
			}

			fprintf(fp, "%-30s %-15i\n",
					"psd_n_diameters",
					myscattererfield->psd[i]->n_diameters);

			if ((myscattererfield->type == 1) | (myscattererfield->type == 2)) {
				fprintf(fp, "%-30s %-15i\n",
						"fit_dBlwc",
						myscattererfield->psd[i]->fit_dBlwc);
			}
			
			if (myscattererfield->type == 1) {
				if (myscattererfield->psd[i]->fit_dBlwc) {
					fprintf(fp, "%-30s %-15.3e\n",
							"cl_x_m_dBlwc",
							myscattererfield->psd[i]->dBlwc_ecm.cl_x_m);
					fprintf(fp, "%-30s %-15.3e\n",
							"cl_y_m_dBlwc",
							myscattererfield->psd[i]->dBlwc_ecm.cl_y_m);
					fprintf(fp, "%-30s %-15.3e\n",
							"cl_z_m_dBlwc",
							myscattererfield->psd[i]->dBlwc_ecm.cl_z_m);
					fprintf(fp, "%-30s %-15.3e\n",
							"cl_t_s_dBlwc",
							myscattererfield->psd[i]->dBlwc_ecm.cl_t_s);
					fprintf(fp, "%-30s %-15.3e\n",
							"c_threshold_dBlwc",
							myscattererfield->psd[i]->dBlwc_ecm.c_threshold);
				}
			}

			if ((myscattererfield->type == 1) | (myscattererfield->type == 2)) {
				fprintf(fp, "%-30s %-15i\n",
						"fit_dBN",
						myscattererfield->psd[i]->fit_dBN);
			}
			
			if (myscattererfield->type == 1) {
				if (myscattererfield->psd[i]->fit_dBN) {
					fprintf(fp, "%-30s %-15.3e\n",
							"cl_x_m_dBN",
							myscattererfield->psd[i]->dBN_ecm.cl_x_m);
					fprintf(fp, "%-30s %-15.3e\n",
							"cl_y_m_dBN",
							myscattererfield->psd[i]->dBN_ecm.cl_y_m);
					fprintf(fp, "%-30s %-15.3e\n",
							"cl_z_m_dBN",
							myscattererfield->psd[i]->dBN_ecm.cl_z_m);
					fprintf(fp, "%-30s %-15.3e\n",
							"cl_t_s_dBN",
							myscattererfield->psd[i]->dBN_ecm.cl_t_s);
					fprintf(fp, "%-30s %-15.3e\n",
							"c_threshold_dBN",
							myscattererfield->psd[i]->dBN_ecm.c_threshold);
				}
			}
			
			if (myscattererfield->psd[i]->distribution_type == 1) {
				fprintf(fp, "#    Equivalent representation of the gamma distribution in discrete pdf\n");
				fprintf(fp, "#    \n");
				fprintf(fp, "#    %-30s %-15i\n",
						"psd_distribution_type",
						0);
			}

			if (myscattererfield->psd[i]->distribution_type == 1) fprintf(fp, "#    ");
			fprintf(fp, "%-30s %-15i", 
						"psd_discrete_D_equiv_mm",
						myscattererfield->psd[i]->n_diameters);
			fprintf_array(fp, myscattererfield->psd[i]->n_diameters, myscattererfield->psd[i]->discrete_D_equiv_mm);
			fprintf(fp, "\n");
			
			for (j=0; j<myscattererfield->psd[i]->n_diameters; j++) {
				if (myscattererfield->psd[i]->distribution_type == 1) fprintf(fp, "#    ");	
				fprintf(fp, "%-30s %-15i", 
							"psd_grid_number_density_m3",
							myscattererfield->psd[i]->field->n);
				fprintf_array(fp, myscattererfield->psd[i]->field->n, myscattererfield->psd[i]->grid_number_density_m3[j]);
				fprintf(fp, "\n");
				
				if (myscattererfield->psd[i]->distribution_type == 1) fprintf(fp, "#    ");	
				fprintf(fp, "%-30s %-15i", 
							"psd_grid_dBnumber_density_err_m3",
							myscattererfield->psd[i]->field->n);
				fprintf_array(fp, myscattererfield->psd[i]->field->n, myscattererfield->psd[i]->grid_dBnumber_density_err_m3[j]);
				fprintf(fp, "\n");
			}
			
		}
	}
	
	fprintf(fp, "\n\n\n");

	fflush(fp);
}





void fprintf_array(FILE *fp, int n, double *variable)
{
	int j;
	for (j=0; j<n; j++) {
		fprintf(fp, " %-15.3e", variable[j]);					
	}
}




void zephyros_config_initialize(t_zephyros_config **pcfg)
{
	t_zephyros_config *cfg 			= malloc(sizeof(t_zephyros_config));
	
	//general
	cfg->general 									= malloc(sizeof(t_zephyros_config_general));
	cfg->general->additional_output					= malloc(sizeof(t_zephyros_config_general_additional_output));
	cfg->general->atmosphere 						= malloc(sizeof(t_zephyros_atmosphere));
	fields_initialize(&cfg->general->atmosphere->field);
	strcpy(cfg->general->atmosphere->field->name, "general->atmosphere");		

	cfg->general->atmosphere->grid_T_K 				= NULL;
	cfg->general->atmosphere->grid_pair_hPa 		= NULL;
	cfg->general->atmosphere->grid_pvapor_hPa 		= NULL;
	cfg->general->instrument 						= malloc(sizeof(t_zephyros_instrument));

	cfg->general->overall									= malloc(sizeof(t_zephyros_config_general_overall));
	cfg->general->water_refractive_index					= malloc(sizeof(t_zephyros_config_general_water_refractive_index));
	cfg->general->white1999_integral						= malloc(sizeof(t_zephyros_config_general_white1999_integral));
	cfg->general->water_refractive_index->wavelength_m 		= NULL;
	cfg->general->water_refractive_index->realindex 		= NULL;
	cfg->general->water_refractive_index->imagindex 		= NULL;
	cfg->general->water_refractive_index->lut_realindex 	= NULL;
	cfg->general->water_refractive_index->lut_imagindex 	= NULL;
	cfg->general->white1999_integral->vec_ln_a 				= NULL;
	cfg->general->white1999_integral->vec_ln_b 				= NULL;
	cfg->general->white1999_integral->vec_ln_L 				= NULL;
	cfg->general->white1999_integral->integral_sqrt			= NULL;
	cfg->general->white1999_integral->lut_integral_sqrt		= NULL;

	//simulation, optional
	cfg->simulation					= NULL;
	
	//retrieval, optional
	cfg->retrieval					= NULL;
	
	cfg->derived_quantities			= malloc(sizeof(t_zephyros_config_derived_quantities));
	
	*pcfg = cfg;
}

void zephyros_config_initialize_simulation(t_zephyros_config *cfg)
{
	if (cfg->simulation == NULL) {
		cfg->simulation										= malloc(sizeof(t_zephyros_config_simulation));	
		cfg->simulation->radarfilter						= calloc(1, sizeof(t_zephyros_radarfilter));
		util_initialize_scattererfield(&(cfg->simulation->scattererfield));
		util_initialize_windfield(&(cfg->simulation->windfield));
	}
}

void zephyros_config_initialize_retrieval(t_zephyros_config *cfg)
{
	int i;
	if (cfg->retrieval == NULL) {
		cfg->retrieval					= malloc(sizeof(t_zephyros_config_retrieval));	
		cfg->retrieval->radarfilter		= calloc(1,sizeof(t_zephyros_radarfilter));
		util_initialize_scattererfield(&(cfg->retrieval->prior_scattererfield));
		cfg->retrieval->post_scattererfield = NULL;
		util_initialize_windfield(&(cfg->retrieval->prior_windfield));
		cfg->retrieval->post_windfield = NULL;
		
		cfg->retrieval->algorithm		= malloc(sizeof(t_zephyros_config_retrieval_algorithm));
		for ( i = 0; i <= 100; i++ ) {
			cfg->retrieval->algorithm->type[i] = 0;
			cfg->retrieval->algorithm->fdvar_cfg[i] = NULL;
			cfg->retrieval->algorithm->lwm_cfg[i] = NULL;
		}
	}
}

void zephyros_config_free(t_zephyros_config **pcfg)
{
	int i;
	t_zephyros_config *cfg = *pcfg;
	
	if (cfg != NULL) {
		free(cfg->general->additional_output);

		fields_free(&cfg->general->atmosphere->field);
		if (cfg->general->atmosphere->grid_T_K != NULL) {free(cfg->general->atmosphere->grid_T_K); cfg->general->atmosphere->grid_T_K = NULL;}
		if (cfg->general->atmosphere->grid_pair_hPa != NULL) {free(cfg->general->atmosphere->grid_pair_hPa); cfg->general->atmosphere->grid_pair_hPa = NULL;}
		if (cfg->general->atmosphere->grid_pvapor_hPa != NULL) {free(cfg->general->atmosphere->grid_pvapor_hPa); cfg->general->atmosphere->grid_pvapor_hPa = NULL;}
		free(cfg->general->atmosphere);
		
		free(cfg->general->instrument);
		
		free(cfg->general->overall);

		if (cfg->general->water_refractive_index->wavelength_m != NULL) {free(cfg->general->water_refractive_index->wavelength_m);}
		if (cfg->general->water_refractive_index->realindex != NULL) {free(cfg->general->water_refractive_index->realindex);}
		if (cfg->general->water_refractive_index->imagindex != NULL) {free(cfg->general->water_refractive_index->imagindex);}
		interpolation_free_lut(&cfg->general->water_refractive_index->lut_realindex);
		interpolation_free_lut(&cfg->general->water_refractive_index->lut_imagindex);
		free(cfg->general->water_refractive_index);
		
		if (cfg->general->white1999_integral->vec_ln_a != NULL) 		{free(cfg->general->white1999_integral->vec_ln_a);}
		if (cfg->general->white1999_integral->vec_ln_b != NULL) 		{free(cfg->general->white1999_integral->vec_ln_b);}
		if (cfg->general->white1999_integral->vec_ln_L != NULL) 		{free(cfg->general->white1999_integral->vec_ln_L);}
		if (cfg->general->white1999_integral->integral_sqrt != NULL) {free(cfg->general->white1999_integral->integral_sqrt);}
		interpolation_free_lut(&cfg->general->white1999_integral->lut_integral_sqrt);
		free(cfg->general->white1999_integral);
		
		free(cfg->general);

		if (cfg->simulation != NULL) {	
			free(cfg->simulation->radarfilter);
		
			util_free_scattererfield(&cfg->simulation->scattererfield);
			util_free_windfield(&cfg->simulation->windfield);
			
			free(cfg->simulation); cfg->simulation = NULL;
		}
		
		if (cfg->retrieval != NULL) {	
			free(cfg->retrieval->radarfilter);
			
			util_free_scattererfield(&cfg->retrieval->prior_scattererfield);
			util_free_scattererfield(&cfg->retrieval->post_scattererfield);
			util_free_windfield(&cfg->retrieval->prior_windfield);
			util_free_windfield(&cfg->retrieval->post_windfield);
			
			for ( i = 0; i <= 100; i++ ) {
				if (cfg->retrieval->algorithm->fdvar_cfg[i] != NULL) {
					free(cfg->retrieval->algorithm->fdvar_cfg[i]);
				}
				if (cfg->retrieval->algorithm->lwm_cfg[i] != NULL) {
					free(cfg->retrieval->algorithm->lwm_cfg[i]->xvec_m);
					free(cfg->retrieval->algorithm->lwm_cfg[i]->yvec_m);
					free(cfg->retrieval->algorithm->lwm_cfg[i]->zvec_m);
					free(cfg->retrieval->algorithm->lwm_cfg[i]->tvec_s);
					free(cfg->retrieval->algorithm->lwm_cfg[i]);
					cfg->retrieval->algorithm->lwm_cfg[i] = NULL;
				}
			}
			free(cfg->retrieval->algorithm);
			free(cfg->retrieval); cfg->retrieval = NULL;
		}

		free(cfg->derived_quantities);
		
		fclose(cfg->fp_ao);
		free(cfg);
		cfg = NULL;
	}
}


void zephyros_config_get_dirname(char filename[8192], char dirname[8192])
{
	int i;
	
	memset(&dirname[0], 0, 8192);
	for ( i = strlen(filename) - 1; i > -1; i--) {
		if ((*(filename + i) == '/') | (*(filename + i) == '\\')) {
			break;
		}
	}
	strncpy(dirname, filename, i+1);	
}


