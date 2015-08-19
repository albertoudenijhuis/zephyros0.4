#ifndef _ZEPHYROS_RADARFILTER
#define _ZEPHYROS_RADARFILTER

#include "zephyros_config.h"
#include "fields.h"
#include "coordinates.h"
#include "particles.h"

typedef struct st_radarmeasurement_analysis
{
    double  	**discrete_D_equiv_mm;
    double		**unweighted_canting_angle_wrt_vertical_mean;
    double		**unweighted_canting_angle_wrt_vertical_variance;
} t_radarmeasurement_analysis;

typedef struct st_radarmeasurement
{
	double		azel_r1_m;			/* radial distance to begin of the resolution volume [m] */
	double		azel_r2_m;			/* radial distance to end of the resolution volume [m] */
	double		azel_alpha_rad;		/* azimuth angle, w.r.t north, positive towards east [rad] */
	double		azel_gamma_rad;		/* elevation angle [rad] */
	double		beam_FWHM0_rad;		/* beam width FWHM [rad], azimutal direction */
	double		beam_FWHM1_rad;		/* beam width FWHM [rad], elevation direction */
	double		dt;					/* dt */
	
	t_zephyros_coordinates		*center_coor;	
		
	//needed for calculations
	double		eta_hh;
	double		eta_hv;
	double		eta_vh;
	double		eta_vv;
	
	double		eta_ShhSvvc;
	double		eta_ShhShvc;
	double		eta_SvvSvhc;
	
	double		n_Re_Shh_min_Svv;

	double		*Doppler_spectrum_ShhSvvc;
	double		*Doppler_spectrum_ShhShvc;
	double		*Doppler_spectrum_SvvSvhc;
			
	//output
	double		dBZ_hh;
	double		dBZ_hh_err;
	double		dBZ_hv;
	double		dBZ_hv_err;
	double		dBZ_vh;
	double		dBZ_vh_err;
	double		dBZ_vv;
	double		dBZ_vv_err;
	
	double		dBZdr;
	double		dBZdr_err;
	double		dBLdr;
	double		dBLdr_err;
	
	//correlation coefficients
	double		rho_co;
	double		rho_co_err;
	double		rho_cxh;
	double		rho_cxh_err;
	double		rho_cxv;
	double		rho_cxv_err;
	
	double		Doppler_velocity_hh_ms; 	/* reflectivity weighted radial velocity [m/s] */
	double		Doppler_velocity_hh_ms_err; 	/* reflectivity weighted radial velocity [m/s] */
	double		Doppler_velocity_hv_ms;
	double		Doppler_velocity_hv_ms_err;
	double		Doppler_velocity_vh_ms;
	double		Doppler_velocity_vh_ms_err;
	double		Doppler_velocity_vv_ms;
	double		Doppler_velocity_vv_ms_err;
	
	double		Doppler_spectral_width_hh_ms; /* reflectivity weighted radial velocity width (1 STD) 	[m/s] */
	double		Doppler_spectral_width_hh_ms_err; 
	double		Doppler_spectral_width_hv_ms;
	double		Doppler_spectral_width_hv_ms_err;
	double		Doppler_spectral_width_vh_ms;
	double		Doppler_spectral_width_vh_ms_err;
	double		Doppler_spectral_width_vv_ms;
	double		Doppler_spectral_width_vv_ms_err;
	
	//spectra
	int			n_spectrum;

	double 		*spectrum_lbound, *spectrum_ubound, *spectrum_center;
	
	double		*Doppler_spectrum_dBZ_hh;
	double		*Doppler_spectrum_dBZ_hh_err;
	double		*Doppler_spectrum_dBZ_hv;
	double		*Doppler_spectrum_dBZ_hv_err;
	double		*Doppler_spectrum_dBZ_vh;
	double		*Doppler_spectrum_dBZ_vh_err;
	double		*Doppler_spectrum_dBZ_vv;
	double		*Doppler_spectrum_dBZ_vv_err;

	double 		*specific_dBZdr;
	double 		*specific_dBZdr_err;
	double 		*specific_dBLdr;
	double 		*specific_dBLdr_err;

	//correlation coefficients
	double		*specific_rho_co;
	double		*specific_rho_co_err;
	double		*specific_rho_cxh;
	double		*specific_rho_cxh_err;
	double		*specific_rho_cxv;
	double		*specific_rho_cxv_err;
	
	//analysis
	t_radarmeasurement_analysis *analysis;
    	
	//differential phase
	double 		KDP;
	double 		KDP_err;
	
	
	
	//additional variables, needed to calculate derivatives
	int 		n_psd;
	int			*n_diameters;
	double		**eta_i_hh; //(psd, particle)
	double		**eta_i_hv; //(psd, particle)
	double		**eta_i_vh; //(psd, particle)
	double		**eta_i_vv; //(psd, particle)
	
	double		***spectrum_eta_i_hh; //(psd, particle, spectrum interval)
	double		***spectrum_eta_i_hv; //(psd, particle, spectrum interval)
	double		***spectrum_eta_i_vh; //(psd, particle, spectrum interval)
	double		***spectrum_eta_i_vv; //(psd, particle, spectrum interval)
	
	//derivatives to eddy dissipation rate ^ (1/3)
	double		der_edr13_dBZ_hh;
	double		der_edr13_dBZdr;
	double		der_edr13_dBLdr;
	double		der_edr13_Doppler_velocity_hh_ms;
	double		der_edr13_Doppler_spectral_width_hh_ms;

	/*
	double		*Doppler_spectrum_dBZ_hh;
	double 		*specific_dBZdr;
	double 		*specific_dBLdr;
	*/
} t_radarmeasurement;

typedef struct st_radarfilter_res_vol			/* resolution volume */
{
	int			n;					/* =n_r*n_theta*n_phi*n_t */
	int			n_beam_range;
	int			n_beam_theta;
	int			n_beam_phi;
	int			n_t;

	//parametric model, number of division points in azimuth and elevation
	int			parametric_turbulence; //(0 = no, 1 = parametric turbulence model used)
	int 		n_parmod_az;
	int 		n_parmod_el;
	int 		n_parmod;

	int			n_psd;
	int			*n_diameters;

	//cross sections widget. dimensions n_psd x n_diameters
	t_zephyros_particles_widget ***subvolume_scat;

	//subvolume center coordinates
	t_zephyros_coordinates **subvolume_coor;

	//resolution volume air velocity
	double		air_u;
	double		air_v;
	double		air_w;
	double		air_u_der[4];
	double		air_v_der[4];
	double		air_w_der[4];

	//subvolume air velocity
	double		*subvolume_air_u;
	double		*subvolume_air_v;
	double		*subvolume_air_w;
	double		**subvolume_air_u_der;
	double		**subvolume_air_v_der;
	double		**subvolume_air_w_der;

	//subvolume particle velocity. dimensions n_res x n_psd x n_diameters x n_parmod
	double		****subvolume_particle_u;
	double		****subvolume_particle_v;
	double		****subvolume_particle_w;

	//subvolume particle Doppler velocity, dimensions n_res x n_psd x n_diameters x n_parmod
	double		****subvolume_particle_Doppler_velocity;

	//subvolume particle orientation,  dimensions n_res x n_psd x n_diameters x n_parmod x 4
	double		*****subvolume_particle_dir;

	//number densities
	double		***subvolume_number_density_m3;
	double		****subvolume_ln_number_density_m3_der;
	
	//subvolume particle velocity, dimensions n_res x n_psd x n_diameters x n_parmod
	double		****subvolume_eta_hh;
	double		****subvolume_eta_hv;
	double		****subvolume_eta_vh;
	double		****subvolume_eta_vv;
	
	double		****subvolume_eta_ShhSvvc;
	double		****subvolume_eta_ShhShvc;
	double		****subvolume_eta_SvvSvhc;
	
	
	//double		*subvolume_edr;
	//double		*subvolume_edr13;
	//double		*subvolume_maxL;



	//double		*subvolume_particle_u_der[101][4];
	//double		*subvolume_particle_v_der[101][4];
	//double		*subvolume_particle_w_der[101][4];


	/*


	double*		subvolume_Z;
	double*		subvolume_dBZ;
	double*		subvolume_attenuation;

	double		subvolume_f_adv_x_der[4];

	double		subvolume_dBZ_der[4];

	double		subvolume_attenuation_der[4];

	
	double		subvolume_edr13_der[4];
	double		subvolume_maxL_der[4];

	double		subvolume_radial_vel_ms;
	double		subvolume_shear;
	
	double*		subvolume_radial_vel_ms_der;


	//special
	double		*subvolume_wind_griddep;
	double		*subvolume_edr_griddep;
	double		*subvolume_scatterer_griddep;
	double		*integrated_wind_griddep;
	double		*integrated_edr_griddep;
	double		*integrated_scatterer_griddep;
	
	double		fcoriolis;
	int			advection_integration_steps;
	
	*/
} t_radarfilter_res_vol;

typedef struct st_radarfilter_todolist			/* resolution volume to do list */
{
	int			calc_center_coordinates;
	
	int			calc_number_density;
	int			calc_number_density_der;

	int			calc_cross_sections; //one of eta calculations
	int			calc_eta_hh;
	int			calc_eta_hv;
	int			calc_eta_vh;
	int			calc_eta_vv;
	int			calc_eta_ShhSvvc;
	int			calc_eta_ShhShvc;
	int			calc_eta_SvvSvhc;
	
    int			calc_air_velocity;
    int			calc_air_velocity_der;

	int			calc_terminal_fall_speeds;

    int			calc_particle_velocity;
    int			calc_particle_velocity_der;

    int			calc_particle_direction;
    
    int			calc_Doppler_mean_velocity;
    int			calc_Doppler_spectrum_widths;



	int 		calc_Doppler_spectrum;	//one of those others
	int 		calc_Doppler_spectrum_ShhSvvc;
	int 		calc_Doppler_spectrum_ShhShvc;
	int 		calc_Doppler_spectrum_SvvSvhc;
	


	int 		calc_eta_i_hh;
	int 		calc_eta_i_hv;
	int 		calc_eta_i_vh;
	int 		calc_eta_i_vv;
	
	int 		calc_spectrum_eta_i_hh;
	int 		calc_spectrum_eta_i_hv;
	int 		calc_spectrum_eta_i_vh;
	int 		calc_spectrum_eta_i_vv;
	
	//int			apply_advection;
	//int			apply_geostrophic_correction;
} t_radarfilter_todolist;

void radarfilter_exec(
	t_zephyros_config 				*cfg,
	int								i_mode,	//0 = simulation mode, 1 = retrieval mode
	t_radarfilter_todolist			*todo,	
	int								n_measurements,
	t_radarmeasurement				**radarmeasurement
);

void radarfilter_initialize_resolution_volume(t_zephyros_config *cfg, int i_mode, t_radarfilter_res_vol** pres_vol, t_radarfilter_todolist *todo);
void radarfilter_free_resolution_volume(t_zephyros_config *cfg, int i_mode, t_radarfilter_res_vol **pres_vol, t_radarfilter_todolist *todo);

void radarfilter_initialize_todolist(t_zephyros_config *cfg, int i_mode, t_radarfilter_todolist **ptodolist);
void radarfilter_prepare_todolist(t_zephyros_config *cfg, int i_mode, t_radarfilter_todolist *todolist);
void radarfilter_free_todolist(t_radarfilter_todolist **ptodolist);

void radarfilter_initialize_radarmeasurement(int n_measurements, t_radarmeasurement ***pradarmeasurement);
void radarfilter_prepare_model_radarmeasurement(int n_measurements, t_radarmeasurement ***pdst, t_radarmeasurement **src);
void radarfilter_free_radarmeasurement(int n_measurements, t_radarmeasurement ***pradarmeasurement);

typedef struct t_radarfilter_readout_widget
{
	FILE *fp;
	fpos_t pos_thisline;
	fpos_t pos_nextline;
	
	char identifier[8192];
	char line[8192];
	
	int ndim;
	int dim[10];
} t_radarfilter_readout_widget;

void radarfilter_read_measurements(int *pn_measurements, t_radarmeasurement ***pradarmeasurement, char measurements_filename[8192]);
void radarfilter_readout(t_radarfilter_readout_widget *rwg, t_radarmeasurement *radarmeasurement);
void radarfilter_write_measurements(t_zephyros_config *cfg, int i_mode, int n_measurements, t_radarmeasurement **radarmeasurement, FILE *fp);
void radarfilter_write_measurements_detailed_analysis(t_zephyros_config *cfg, int i_mode, int n_measurements, t_radarmeasurement **radarmeasurement, FILE *fp);


#endif
