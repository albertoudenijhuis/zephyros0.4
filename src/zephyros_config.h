#ifndef _ZEPHYROS_CONFIG
#define _ZEPHYROS_CONFIG

#define zephyros_version "0.4"

#include <stdio.h>
#include <complex.h>

#include "fields.h"
#include "interpolation.h"
#include "util.h"

typedef struct t_zephyros_config_general_additional_output
{
int print_configuration;
int print_detailed_analysis;
} t_zephyros_config_general_additional_output;

typedef struct t_zephyros_config_general_overall
{
	char versionnumber[8192];
	char file_name[8192];
} t_zephyros_config_general_overall;

typedef struct t_zephyros_config_general_water_refractive_index
{
int n;
double *wavelength_m;
double *realindex;
double *imagindex;
t_zephyros_interpolation_bilint_lut	*lut_realindex;
t_zephyros_interpolation_bilint_lut	*lut_imagindex;
} t_zephyros_config_general_water_refractive_index;

typedef struct t_zephyros_config_general_white1999_integral
{
int na;
double *vec_ln_a;
int nb;
double *vec_ln_b;
int nL;
double *vec_ln_L;
int n;
double *integral_sqrt;
t_zephyros_interpolation_bilint_lut	*lut_integral_sqrt;
} t_zephyros_config_general_white1999_integral;

typedef struct t_zephyros_config_general
{
	t_zephyros_config_general_additional_output *additional_output;
	t_zephyros_atmosphere							*atmosphere;
	t_zephyros_instrument							*instrument;
	t_zephyros_config_general_overall *overall;
	t_zephyros_config_general_water_refractive_index *water_refractive_index;
	t_zephyros_config_general_white1999_integral *white1999_integral;
} t_zephyros_config_general;


typedef struct t_zephyros_config_simulation
{
	t_zephyros_radarfilter 							*radarfilter;
	t_zephyros_scattererfield						*scattererfield;
	t_zephyros_windfield			 				*windfield;
} t_zephyros_config_simulation;

typedef struct t_zephyros_config_retrieval_fdvar_cfg
{
	//void		*vd_zcfg;
	//double		fcoriolis;
    int			costfunction_dBZ_hh;
    int			costfunction_dBZdr;
    int			costfunction_dBLdr;
    int			costfunction_Doppler_velocity_hh_ms;
    int			costfunction_Doppler_spectral_width_hh_ms;
    
    //spectra
    int			costfunction_Doppler_spectrum_dBZ_hh;
    int			costfunction_specific_dBZdr;
    int			costfunction_specific_dBLdr;

	//special stuff
    int			n_active_windfield_grid_nrs;
    int			*active_windfield_grid_nrs;
    int			n_active_windfield_turbulence_nrs;
    int			*active_windfield_turbulence_nrs;
    int			n_active_scattererfield_nrs;
    int			*active_scattererfield_nrs;
    
    int			n_cast_windfield_grid_nrs;
    int			*cast_windfield_grid_nrs;

    double		update_windfield_hspeed_err;		
    double		update_windfield_hdir_err;		
    double		update_windfield_u_err;		
    double		update_windfield_v_err;		
    double		update_windfield_w_err;		

    
    double		maximum_time_s;		
    
    
    
    //fit h_speed --> in windfield
    //fit h_dir --> in windifeld
    
    //int			fit_w;
//    int			apply_advection;
//    int			apply_geostrophic_correction;
//    double		fit_maximum_time_s;
//    int			advection_integration_steps;
} t_zephyros_config_retrieval_fdvar_cfg;
void zephyros_config_initialize_retrieval_fdvar_cfg(t_zephyros_config_retrieval_fdvar_cfg **pcfg);
void zephyros_config_free_retrieval_fdvar_cfg(t_zephyros_config_retrieval_fdvar_cfg **pcfg);

typedef struct t_zephyros_config_retrieval_lwm_cfg
{
	void		*vd_zcfg;
	double		*xvec_m;
	double		*yvec_m;
	double		*zvec_m;
	double		*tvec_s;	
	int			fit_u0;
	int			fit_u_x;
	int			fit_u_z;
	int			fit_v0;
	int			fit_v_y;
	int			fit_v_z;
	int			fit_u_y_plus_v_x;
	int			fit_w0;
	int			fit_w_x;
	int			fit_w_y;
	int			fit_w_z;
	int			fit_u_t_plus_v_t_plus_w_t;
	int			apply_weights;
	double 		maximum_time_s;
	int	 		extra_points_n;
	double 		extra_points_dx;
	double 		extra_points_dy;
	double 		extra_points_dz;
	double 		extra_points_dt;
} t_zephyros_config_retrieval_lwm_cfg;
void zephyros_config_initialize_retrieval_lwm_cfg(t_zephyros_config_retrieval_lwm_cfg **pcfg);

typedef struct t_zephyros_config_retrieval_algorithm
{
	int		type[101];
	t_zephyros_config_retrieval_fdvar_cfg			*fdvar_cfg[101];
	t_zephyros_config_retrieval_lwm_cfg				*lwm_cfg[101];
} t_zephyros_config_retrieval_algorithm;
typedef struct t_zephyros_config_retrieval
{
	t_zephyros_radarfilter 							*radarfilter;
	t_zephyros_config_retrieval_algorithm	 		*algorithm;
	t_zephyros_scattererfield						*prior_scattererfield;	
	t_zephyros_scattererfield						*post_scattererfield;	
	t_zephyros_windfield			 				*prior_windfield;		
	t_zephyros_windfield			 				*post_windfield;		
} t_zephyros_config_retrieval;



typedef struct t_zephyros_config_derived_quantities
{
	double central_wavelength_m;
	//double noise_power;
	double radar_vmax_ms;
	
	double complex radar_water_refractive_index;
} t_zephyros_config_derived_quantities;


typedef struct t_zephyros_config
{
	t_zephyros_config_general 				*general;
	t_zephyros_config_simulation 			*simulation;
	t_zephyros_config_retrieval 			*retrieval;
	t_zephyros_config_derived_quantities   	*derived_quantities;
		
	FILE *fp_ao; //file pointer for additional output
} t_zephyros_config;

typedef struct t_zephyros_config_read_widget_indices
{
	int windfield_grid;
	int windfield_wave;
	int windfield_vortex;
	int windfield_turbulence;
	int scattererfield_psd;
	int scattererfield_psd_diameter_i;

	int algorithm_run;
} t_zephyros_config_read_widget_indices;

typedef struct t_zephyros_config_read_widget
{
	FILE *fp;
	fpos_t pos_thisline;
	fpos_t pos_nextline;
	
	char identifier[8192];
	char section[8192];
	char subsection[8192];
	char firstchar[1];
	char line[8192];
	char commentSign[1];

	t_zephyros_config_read_widget_indices simulation_indices;
	t_zephyros_config_read_widget_indices retrieval_indices;
} t_zephyros_config_read_widget;

//functions
void zephyros_config_read(char file_name[8192], char additional_output_filename[8192], t_zephyros_config **pcfg);

int zephyros_config_read_position(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192]);

void zephyros_config_read_________i(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192],
	int *destination);	
	
void zephyros_config_read________lf(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192],
	double *destination);
		
void zephyros_config_read_________s(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192],
	char *destination);
	
void zephyros_config_read__lf_array(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192],
	int  *arraysize,
	double **destination);	
	
void zephyros_config_read___i_array(
	t_zephyros_config_read_widget *rwg,
	char section[8192],
	char subsection[8192],
	char identifier[8192],
	int  *arraysize,
	int **destination);	
	
void zephyros_config_derive_quantities_zephyros_config(t_zephyros_config *cfg);
	
void zephyros_config_print(t_zephyros_config *cfg, FILE *fp);
void zephyros_config_print_radarfilter(t_zephyros_radarfilter *myradarfilter, FILE *fp);
void zephyros_config_print_windfield(t_zephyros_windfield *mywindfield, FILE *fp);
void zephyros_config_print_scattererfield(t_zephyros_scattererfield *myscattererfield, FILE *fp);

void fprintf_array(FILE *fp, int n, double *variable);
void fprinti_array(FILE *fp, int n, int *variable);
void zephyros_config_initialize(t_zephyros_config **pcfg);

void zephyros_config_initialize_simulation(t_zephyros_config *cfg);
void zephyros_config_initialize_retrieval(t_zephyros_config *cfg);

void zephyros_config_free(t_zephyros_config **pcfg);

void zephyros_config_get_dirname(char filename[8192], char dirname[8192]);



#endif
