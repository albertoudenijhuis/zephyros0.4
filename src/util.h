#ifndef _ZEPHYROS_UTIL
#define _ZEPHYROS_UTIL

#include <complex.h>
#include "interpolation.h"
#include "particles.h"
#include "fields.h"
#include "cs.h"

//Zephyros types
typedef struct st_zephyros_field_errcovmat
{
	double cl_x_m;
	double cl_y_m;
	double cl_z_m;
	double cl_t_s;
	double c_threshold;
	
	int n;
	cs *mat;
    css *mat_S;
    csn *mat_N;      
} t_zephyros_field_errcovmat;

//radar filter configuration
typedef struct t_zephyros_radarfilter
{
	int 	n_beam_range;
	int 	n_beam_theta;
	int 	n_beam_phi;
	int 	n_t;
	
	int 	n_parmod_az;
	int 	n_parmod_el;
	
	int		n_spectrum;
	
	int 	filter_dBZ_hh;
	int 	filter_dBZ_hv;
	int 	filter_dBZ_vh;
	int 	filter_dBZ_vv;
	int 	filter_dBZdr;
	int 	filter_dBLdr;
	int 	filter_rho_co;
	int 	filter_rho_cxh;
	int 	filter_rho_cxv;
	int 	filter_KDP;
	int 	filter_Doppler_velocity_hh_ms;
	int 	filter_Doppler_velocity_hv_ms;
	int 	filter_Doppler_velocity_vh_ms;
	int 	filter_Doppler_velocity_vv_ms;
	int 	filter_Doppler_spectralwidth_hh_ms;
	int 	filter_Doppler_spectralwidth_hv_ms; 
	int 	filter_Doppler_spectralwidth_vh_ms; 
	int 	filter_Doppler_spectralwidth_vv_ms;
	int 	filter_Doppler_spectrum_dBZ_hh;
	int 	filter_Doppler_spectrum_dBZ_hv;
	int 	filter_Doppler_spectrum_dBZ_vh;
	int 	filter_Doppler_spectrum_dBZ_vv;
	int 	filter_specific_dBZdr;
	int 	filter_specific_dBLdr;
	int 	filter_specific_rho_co;
	int 	filter_specific_rho_cxh;
	int 	filter_specific_rho_cxv;
	
	int		filter_errors;
	
	int 	effective_earth_correction;
	
	int 	inertia_effect;
	
	int 	additive_noise;
	int 	multiplicative_noise;
} t_zephyros_radarfilter;

typedef struct t_zephyros_atmosphere
{
	t_zephyros_field		*field;
	double*		grid_T_K;
	double*		grid_pair_hPa;		/* air pressure */
	double*		grid_pvapor_hPa; 	/* vapor pressure */	
	
	t_zephyros_interpolation_bilint_lut	*lut_T_K;
	t_zephyros_interpolation_bilint_lut	*lut_ln_grid_pair_hPa;
	t_zephyros_interpolation_bilint_lut	*lut_ln_grid_pvapor_hPa;	
} t_zephyros_atmosphere;

typedef struct t_zephyros_instrument
{
	int		type;						//0 = radar
	double 	transmit_power_watt;		//Pt
	double 	k_nonuniformity;			//
	double 	power_gain;					//-
	double 	lr_receiver_loss_factor;	//-
	double 	central_frequency_hz;		//hz		-->> determines wavelength
	double 	prf_hz;						//hz		-->> determines vmax
	double 	temperature_K;				//
} t_zephyros_instrument;

typedef struct st_zephyros_wave
{
	double		*amplitude_msinv;
	double 		phi0_rad;
	double		*k_minv;	//wavenumbers
	double		f_sinv;		//frequency
} t_zephyros_wave;

typedef struct st_zephyros_vortex
{
	int			type;				//(1 = Rankine vortex, 2 = Lamb-Oseen vortex)
	double		*xyz0_m;
	double		t0_s;
	double 		*xyz_scaling_m;
	double 		*rotation_direction;
	double 		vmax_msinv;
	double 		maxr_m;
	double 		r_m;			//Rankine parameter, radius in m
	double 		nu;
} t_zephyros_vortex;


typedef struct t_zephyros_turbulence_widget
{
	int			type;				//(1 = Mann1998, 2 = CTM, 3 = Careta1993, 4 = Pinsky2006, 5 = parametric)
	t_zephyros_field 	*field;
	double 		*grid_edr;
	double 		*grid_edr13;
	t_zephyros_interpolation_bilint_lut	*lut_edr13;

	double 		Lx;							
	double 		Ly;							
	double 		Lz;							
	int			Nx;							
	int			Ny;							
	int			Nz;							
	double		*freqx;						
	double		*freqy;						
	double		*freqz;						
		
	double 		*grid_karman_L;						//mann1998, ctm
	double 		*grid_kolmogorov_constant;			//mann1998, ctm
	t_zephyros_interpolation_bilint_lut	*lut_karman_L;
	t_zephyros_interpolation_bilint_lut	*lut_kolmogorov_constant;
	
	double complex		*random_fourier_cx;		//ctm
	double complex		*random_fourier_cy;		//ctm
	double complex		*random_fourier_cz;		//ctm
	double 				*random_shift_x;		//ctm
	double 				*random_shift_y;		//ctm
	double 				*random_shift_z;		//ctm
		
	double complex		***random_fourier_c0;	//mann1998
	double complex 		***random_fourier_c1;	//mann1998
	double complex 		***random_fourier_c2;	//mann1998
	double complex 		***u_fourier;	//mann1998
	double complex 		***v_fourier;	//mann1998
	double complex 		***w_fourier;	//mann1998
	
	double complex		**random_fourier_beta;	//careta1993
	
	char		pinsky2006_file[8192];					//pinsky2006
	
	double 		minL_div_maxL;						//ctm
	
	int			random_numbers;					//mann1998, ctm
	char		rn_file[8192];					//mann1998, ctm

	double 		lambdax; 						//careta1993
	double 		lambday; 						//careta1993
	
	int			M;								//pinsky2006
	int			K;								//pinsky2006
	double 		**random_a;						//pinsky2006
	double 		**random_b;						//pinsky2006
	double 		*lambdak;						//pinsky2006
	double 		*alphak;						//pinsky2006
	
	double		calibration_factor;
	int			calibration_method;				
	int			calibration_n;
	double		calibration_C;
	int			calibration_periodic;
	int			calibration_nint;
	int			calibration_dir;
	double		calibration_L;
	
	//field representation
	t_zephyros_field 	*turb_field;
	double *turb_u;
	double *turb_v;
	double *turb_w;
	t_zephyros_interpolation_bilint_lut	*turb_lut_u;
	t_zephyros_interpolation_bilint_lut	*turb_lut_v;
	t_zephyros_interpolation_bilint_lut	*turb_lut_w;
	
	//fitting
	double 		*grid_edr13_err;
	int fit_edr13;
	int fit_edr13_Knr;
	t_zephyros_field_errcovmat edr13_ecm;
} t_zephyros_turbulence_widget;

typedef struct t_zephyros_windfield
{
	int					type;		//0 = normal, 1 = prior, 2 = post
	
    int			nfields;
	t_zephyros_field 	*field[101];
    double 		*grid_hspeed[101];
    double 		*grid_hdir[101];
    double 		*grid_u[101];
    double	 	*grid_v[101];
    double	 	*grid_w[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_hspeed[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_hdir[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_u[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_v[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_w[101];

    int								nvortices;
    t_zephyros_vortex				*vortex[101];
    int								nwaves;
    t_zephyros_wave					*wave[101];
    int								nturbulences;
    t_zephyros_turbulence_widget 	*turbulence[101];
    
    //for retrieval mode, prior parameters
    double 		*grid_hspeed_err[101];
    double 		*grid_hdir_err[101];
    double 		*grid_u_err[101];
    double	 	*grid_v_err[101];
    double	 	*grid_w_err[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_hspeed_err[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_hdir_err[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_u_err[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_v_err[101];
    t_zephyros_interpolation_bilint_lut	 	*lut_w_err[101];
    
	int fit_hspeed[101];
	int fit_hdir[101];
	int fit_u[101];
	int fit_v[101];
	int fit_w[101];
	
	int fit_hspeed_Knr[101];
	int fit_hdir_Knr[101];
	int fit_u_Knr[101];
	int fit_v_Knr[101];
	int fit_w_Knr[101];
	
	//error covariance matrices
	t_zephyros_field_errcovmat hspeed_ecm[101];
	t_zephyros_field_errcovmat hdir_ecm[101];
	t_zephyros_field_errcovmat u_ecm[101];
	t_zephyros_field_errcovmat v_ecm[101];
	t_zephyros_field_errcovmat w_ecm[101];
} t_zephyros_windfield;

void util_initialize_windfield(t_zephyros_windfield **pwindfield);
void util_prepare_windfield(t_zephyros_windfield *windfield);
void util_prepare_post_windfield(t_zephyros_windfield **pdst, t_zephyros_windfield *src);
void util_free_windfield(t_zephyros_windfield **pwindfield);


void util_windfield_wave_fuvw(t_zephyros_wave *wave, double *xyzt, int i_uvw, double *output, int calcderivatives, double *derivatives);
void util_windfield_vortex_fuvw(t_zephyros_vortex *vortex, double *xyzt, int i_uvw, double *output, int calcderivatives, double *derivatives);
void util_windfield_fuvw(t_zephyros_windfield *windfield, double* xyzt, int i_uvw, double *output, int calcderivatives, double *derivatives);

void util_field2lut(t_zephyros_field *field, double *variable, int special, t_zephyros_interpolation_bilint_lut **plut);






//particle size distribution
typedef struct st_zephyros_psd
{
	int			distribution_type;				//0 = discrete distribution, 1 = gamma-distribution
    int 		particle_type;					//1 = spherical droplet, 2 = speroid droplet

	t_zephyros_field 	*field;
   
    double  	*grid_lwc_gm3;
    double  	*grid_dBlwc_err_gm3;
    double  	*discrete_D_equiv_mm;
    double 		**grid_number_density_m3;
    double 		**grid_dBnumber_density_err_m3;
    t_zephyros_interpolation_bilint_lut	**lut_ln_number_density_m3;
    t_zephyros_interpolation_bilint_lut	**lut_dBnumber_density_err_m3;

	double *grid_gammadistribution_N0; 	//(gridded number density parameter of the gamma distribution)
	double *grid_gammadistribution_N0_err; 	//(gridded number density parameter of the gamma distribution)
	double gammadistribution_mu;
	double gammadistribution_mu_err;
	double gammadistribution_D0_mm;
	double gammadistribution_D0_err_mm;
	double gammadistribution_dmin_mm;
	double gammadistribution_dmax_mm; 
    int	    n_diameters;					//number of diameters used in discrete psd

	double *grid_gammadistribution_mu; 	//(gridded number density parameter of the gamma distribution)
	double *grid_gammadistribution_D0_mm; 	//(gridded number density parameter of the gamma distribution)

	int fit_dBlwc;			//fit liquid water content
	int fit_dBlwc_Knr;
	t_zephyros_field_errcovmat dBlwc_ecm;
	
	int fit_dBN;			//fit inidividual number densities
	int fit_dBN_Knr;		
	t_zephyros_field_errcovmat dBN_ecm;
	
	//int fit_N; TBD

	
	//error covariance matrices
} t_zephyros_psd;
typedef struct st_zephyros_scattererfield
{
	int					type;		//0 = normal, 1 = prior, 2 = post
    int					npsd;
    t_zephyros_psd 		*psd[101];		//particle size distribution
} t_zephyros_scattererfield;



void util_initialize_scattererfield(t_zephyros_scattererfield **pscattererfield);
void util_prepare_scattererfield(t_zephyros_scattererfield *scattererfield);
void util_prepare_post_scattererfield(t_zephyros_scattererfield **pdst, t_zephyros_scattererfield *src);
void util_free_scattererfield(t_zephyros_scattererfield **pscattererfield);

void util_initialize_psd(t_zephyros_psd **pthepsd);
void util_prepare_psd(t_zephyros_psd *thepsd, t_zephyros_scattererfield *scattererfield);
void util_prepare_post_psd(t_zephyros_psd **pdst, t_zephyros_psd *src);




double util_inverse_gammaDalpha_cdf(double mu, double D0, double Dmin, double Dmax, double cdfP);
double util_inverse_gammaDalpha_cdf_costf(double D, void *vp);
double util_inverse_gammaDalpha_cdf_costf_der(double D, void *vp);
double util_gammaDalpha_cdf(double D, void *vp) ;
double util_gamma_integral(double N0, double mu, double D0, double D1, double D2);

void util_initialize_field_err_cov_matrix(t_zephyros_field *field, double *err, t_zephyros_field_errcovmat *ecm);

void util_safe_free(void **ptr);

void util_copy_dbl_array(int n, double **pdst, double *src);

#endif
