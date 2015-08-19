/*
Description: 
	Scatterer functions

Revision History:
	2014

Functions:
	Particle functions

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

#include "particles.h"
#include "func.h"

/*
http://en.wikipedia.org/wiki/Viscosity
Sutherland's formula can be used to derive the dynamic viscosity of an ideal gas as a function of the temperature:[22]

 {\mu} = {\mu}_0 \frac {T_0+C} {T + C} \left (\frac {T} {T_0} \right )^{3/2}.
This in turn is equal to

\lambda\,\frac{T^{3/2}}{T+C}\,,  where  \lambda = \frac{\mu_0(T_0+C)}{T_0^{3/2}}\,  is a constant for the gas.
in Sutherland's formula:

μ = dynamic viscosity (Pa·s or μPa·s) at input temperature T,
μ0 = reference viscosity (in the same units as μ) at reference temperature T0,
T = input temperature (kelvin),
T0 = reference temperature (kelvin),
C = Sutherland's constant for the gaseous material in question.
Valid for temperatures between 0 < T < 555 K with an error due to pressure less than 10% below 3.45 MPa.

According to Sutherland's formula, if the absolute temperature is less than C, the relative change in viscosity for a small change in temperature is greater than the relative change in the absolute temperature, but it is smaller when T is above C. The kinematic viscosity though always increases faster than the temperature (that is, d log(ν)/d log(T) is greater than 1).

Sutherland's constant, reference values and λ values for some gases:
*/


//function: particles_air_parameters
//input:
//- scat->T_K
//- scat->air_drypressure_hPa
//- scat->air_vaporpressure_hPa
//
//output:
//- scat->air_density
//- scat->air_totalpressure_hPa
//- scat->air_dynamic_viscosity
//- scat->air_kinematic_viscosity
//- scat->air_gravity

void particles_air_parameters(t_zephyros_particles_widget *scat, t_zephyros_coordinates *coor)
{
	double viscosity_C 		= 120;			//K
	double viscosity_T0 	= 291.15;		//K
	double viscosity_mu0 	= 18.27;		//muPa s
	double viscosity_mu;	

	double specific_gas_constant_Rd = 287.058; //Specific gas constant for dry air, 287.058 J/(kg·K)
	double specific_gas_constant_Rv = 461.495; //Specific gas constant for water vapor, 461.495 J/(kg·K)

	double gravity_g0 = 9.80665; //[m s-2]
	double gravity_re = 6371.e3; //m

	//air density
	scat->air_density = ((100. * scat->air_drypressure_hPa) / (specific_gas_constant_Rd * scat->air_temperature_K)) + ((100. * scat->air_vaporpressure_hPa) / (specific_gas_constant_Rv * scat->air_temperature_K));

	//air total pressure
	scat->air_totalpressure_hPa = scat->air_drypressure_hPa + scat->air_vaporpressure_hPa;
	
	//air viscosity
	//typically, dynamic viscosity is mu, eta
	//typically, kinematic viscosity is nu
	scat->air_dynamic_viscosity 	=  1.e-6 * viscosity_mu0 * ((viscosity_T0 + viscosity_C) / (scat->air_temperature_K + viscosity_C)) * pow(scat->air_temperature_K / viscosity_T0, 3./2.); //Pa s = kg m-1 s-3  (S.I.)
	scat->air_kinematic_viscosity 	= scat->air_dynamic_viscosity / scat->air_density;
	
	//gravity, old way
	//rough estimate of the height, based on the pressure
	//(scat->air_density / (100. * scat->air_totalpressure_hPa * gravity_g0)) * log( scat->air_totalpressure_hPa / 1013.25);
	
	//gravity, new way
	scat->air_gravity = gravity_g0 * pow(gravity_re / (gravity_re + coor->enu_xyzt[2]), 2.);
}










//function: particles_spheroid_geometry_beard1987
//input:
//- scat->particle_type
//- scat->particle_D_eqvol_mm
//- scat->air_density
//- scat->air_gravity
//
//output:
//- scat->particle_density
//- scat->particle_axis_ratio
//- scat->particle_D_maj_mm
//- scat->particle_D_min_mm
//- scat->particle_mass_kg
//- scat->particle_A_maj
//- scat->particle_A_min
//- scat->particle_volume
void particles_spheroid_geometry_beard1987(t_zephyros_particles_widget *scat)
{	
	if ((0 < scat->particle_type) & (scat->particle_type < 10)) {
		//water density
		scat->particle_density = 998.2071; //[kg m-3]
	} else {
		//ice crystal density
		//rhoB = ...
		printf("Ice crystals not implemented yet! exiting...");
		exit(0);
	}
	
	//diameters
	if ((scat->particle_type == 2) | (scat->particle_type == 4)) {
		//spheroid water droplet
		scat->particle_axis_ratio 	= 
			1.0048 +
			(5.7e-4	* scat->particle_D_eqvol_mm) +
			(-2.628e-2 * pow(scat->particle_D_eqvol_mm, 2.)) +
			(3.682e-3 * pow(scat->particle_D_eqvol_mm, 3.)) +
			(-1.677e-4 * pow(scat->particle_D_eqvol_mm, 4.));

		scat->particle_D_maj_mm 		= scat->particle_D_eqvol_mm * pow(scat->particle_axis_ratio, -1./3.);
		scat->particle_D_min_mm 		= scat->particle_D_maj_mm * scat->particle_axis_ratio;		
	} else {
		//spherical water droplet
		scat->particle_axis_ratio 	= 1.;
		scat->particle_D_maj_mm 		= scat->particle_D_eqvol_mm;
		scat->particle_D_min_mm 		= scat->particle_D_eqvol_mm;
	}


	scat->particle_mass_kg 		= scat->particle_density * (4./3.) * M_PI * pow((scat->particle_D_eqvol_mm/2.) * 1.e-3,3.); //kg
	scat->particle_A_maj 		= M_PI * pow((scat->particle_D_maj_mm/2.) * 1.e-3, 2.);  //m2
	scat->particle_A_min 		= M_PI * scat->particle_D_maj_mm * scat->particle_D_min_mm * pow((1./2.) * 1.e-3, 2.);  //m2
	scat->particle_volume 		= scat->particle_mass_kg / scat->particle_density;
}






//function: particles_terminal_fall_speed_khvorostyanov2005
//input:
//- scat->particle_type
//- scat->particle_D_maj_mm
//- scat->particle_D_min_mm
//- scat->particle_A_maj
//- scat->particle_A_min
//- scat->particle_density
//- scat->particle_mass_kg
//- scat->air_kinematic_viscosity
//- scat->air_density
//- scat->air_gravity
//
//output:
//- scat->particle_terminal_fall_speed
//- scat->particle_inertial_time_maj
void particles_terminal_fall_speed_khvorostyanov2005(t_zephyros_particles_widget *scat)
{
	//KC02 constants both for drops and crystals
	double delta0 	= 9.06;			//KC02 parameter
	double C0		= 0.292;		//KC02 parameter
	double C1		= 4. / (pow(delta0,2.) * pow(C0, 1./2.));
	
	double X_maj, X_min;	//Davies (or Best) number	[-]
	double aRe_l_maj, aRe_l_min, bRe_l_maj, bRe_l_min;	//[-], [-]
	double delta_bRe_maj, delta_bRe_min;		//[-]
	double aRe_t_maj, aRe_t_min, bRe_t_maj, bRe_t_min;	//[-], [-]
	double X0;				//[-]
	double Ct;				//[-]
	double k;				//[-]
	double z_maj, z_min;	//[-]
	double psi_maj, psi_min;//[-]
	double zeta_maj, zeta_min; //[-]
	
	double Re_maj, Re_min;	//[-]
	double CD_maj, CD_min;	//[-]

	//Khvorostyanov (2005) uses the following expression (which is slightly better):
	X_maj = (2. * scat->particle_volume * fabs(scat->particle_density - scat->air_density) * scat->air_gravity * pow(scat->particle_D_maj_mm * 1.e-3, 2.)) / (scat->particle_A_maj * scat->air_density * pow(scat->air_kinematic_viscosity, 2.)); //S.I.
	X_min = (2. * scat->particle_volume * fabs(scat->particle_density - scat->air_density) * scat->air_gravity * scat->particle_D_maj_mm * scat->particle_D_min_mm * pow(1.e-3, 2.)) / (scat->particle_A_min * scat->air_density * pow(scat->air_kinematic_viscosity, 2.)); //S.I.
		
	bRe_l_maj = C1 * pow(X_maj, 1./2.) /
			(
			2. 
				* (pow(1 + (C1 * pow(X_maj, 1./2.)),1./2.) -1.) 
				* pow( 1. + (C1 * pow(X_maj, 1./2.)), 1./2.)
			);
	bRe_l_min = C1 * pow(X_min, 1./2.) /
			(
			2. 
				* (pow(1 + (C1 * pow(X_min, 1./2.)),1./2.) -1.) 
				* pow( 1. + (C1 * pow(X_min, 1./2.)), 1./2.)
			);

	aRe_l_maj = (pow(delta0, 2.) / 4.)
				* (pow(pow(1. + (C1 * pow(X_maj, 1./2.)), 1./2.) - 1., 2) / pow(X_maj, bRe_l_maj));
	aRe_l_min = (pow(delta0, 2.) / 4.)
				* (pow(pow(1. + (C1 * pow(X_min, 1./2.)), 1./2.) - 1., 2) / pow(X_min, bRe_l_min));

	//turbulence correction
	Ct 	= 1.6;
	k 	= 2.;

	if ((0 < scat->particle_type) & (scat->particle_type < 10)) {		
		//droplets
		X0 = 6.7e6;
	}
	//ice crystals
	//X0 = 2.8e6;	
	
	z_maj = X_maj / X0;
	z_min = X_min / X0;
	
	delta_bRe_maj = -1. *
					k * (Ct - 1) * pow(z_maj, k)
					/
					(
					2. * (1. + pow(z_maj, k)) * (1. + (Ct * pow(z_maj, k)))
					);
	delta_bRe_min = -1. *
					k * (Ct - 1) * pow(z_min, k)
					/
					(
					2. * (1. + pow(z_min, k)) * (1. + (Ct * pow(z_min, k)))
					);
					
	psi_maj = 	(1. + pow(z_maj, k)) 
			/ 
			(1. + (Ct * pow(z_maj, k)));
	psi_min = 	(1. + pow(z_min, k)) 
			/ 
			(1. + (Ct * pow(z_min, k)));

	zeta_maj = pow(psi_maj, 1./2.) / pow(X_maj, delta_bRe_maj);
	zeta_min = pow(psi_min, 1./2.) / pow(X_min, delta_bRe_min);
	
	//apply turbulence correction
	if ((scat->particle_type == 1) | (scat->particle_type == 2)) {
		aRe_t_maj = aRe_l_maj * zeta_maj;
		aRe_t_min = aRe_l_min * zeta_min;
		
		bRe_t_maj = bRe_l_maj + delta_bRe_maj;
		bRe_t_min = bRe_l_min + delta_bRe_min;
	}

	//do not apply turbulence correction
	if ((scat->particle_type == 3) | (scat->particle_type == 4)) {
		aRe_t_maj = aRe_l_maj;
		aRe_t_min = aRe_l_min;
		
		bRe_t_maj = bRe_l_maj;
		bRe_t_min = bRe_l_min;
	}
		
	scat->particle_terminal_fall_speed 
			= aRe_t_maj * pow(scat->air_kinematic_viscosity, 1. - (2. * bRe_t_maj))
			* pow( (2. * scat->particle_mass_kg * scat->air_gravity * (1. - scat->air_density / scat->particle_density)) / (scat->air_density * scat->particle_A_maj) , bRe_t_maj)
			* pow(scat->particle_D_maj_mm * 1.e-3, (2. * bRe_t_maj) - 1.);

	//calculate inertial times
	Re_maj	= aRe_t_maj * pow(X_maj , bRe_t_maj);
	Re_min	= aRe_t_min * pow(X_min , bRe_t_min);
	
	CD_maj	= C0 * pow(1. + (delta0 / pow(Re_maj , 1./2.)), 2.);
	CD_min	= C0 * pow(1. + (delta0 / pow(Re_min , 1./2.)), 2.);
	
	
	scat->particle_inertial_eta_z = (scat->air_density * scat->particle_A_maj * CD_maj) / (2. * scat->particle_mass_kg) ;
	scat->particle_inertial_eta_xy = (scat->air_density * scat->particle_A_min * CD_min) / (2. * scat->particle_mass_kg);
	
	scat->particle_inertial_distance_z = 2. / scat->particle_inertial_eta_z;
	scat->particle_inertial_distance_xy = (1. - exp(-1.)) / (exp(-1.) * scat->particle_inertial_eta_xy);
}

//function: particles_terminal_fall_speed_mitchell1996
//input:
//- scat->particle_D_maj_mm
//- scat->particle_D_min_mm
//- scat->particle_A_maj
//- scat->particle_A_min
//- scat->particle_mass_kg
//- scat->air_dynamic_viscosity
//- scat->air_density
//- scat->air_gravity
//
//output:
//- scat->particle_terminal_fall_speed
//- scat->particle_inertial_time_maj
void particles_terminal_fall_speed_mitchell1996(t_zephyros_particles_widget *scat)
{
	double delta0 	= 9.06;	//[-]
	double C0 		= 0.292;//[-]
		
	double aRe, bRe;		//[-],[-]
	double eta, nu;			//[kg m-1 s-1], [m2 s-1]
	double Re_maj, Re_min, CD_maj, CD_min;			//[-], [-]
	double X_maj, X_min;	//Davies (or Best) number
	
	double ptotal_hPa;		//[hPa]
	double rhoF, rhoB;		//[kg m-3], [kg m-3]
	double re, g0, g, h;	//[m], [m s-2], [m s-2], [m]
	
	double Rd, Rv;			//[J/(kg·K)], [J/(kg·K)]

	double D_eqvol_mm;
	double D_min_mm;

	//Mitchel (1996) uses the following expression:
	X_maj = (2. * scat->particle_mass_kg * scat->air_gravity * scat->air_density * pow(scat->particle_D_maj_mm * 1.e-3, 2.)) / (scat->particle_A_maj * pow(scat->air_dynamic_viscosity, 2.));
	X_min = (2. * scat->particle_mass_kg * scat->air_gravity * scat->air_density * pow(scat->particle_D_min_mm * 1.e-3, 2.)) / (scat->particle_A_min * pow(scat->air_dynamic_viscosity, 2.));
		
	Re_maj 	= 		(pow(delta0,2.) / 4.)
			*	pow(
					pow(1. + 
						(	
							(4. * pow(X_maj, 1./2.))
							/
							(pow(delta0, 2.) * pow(C0, 1./2.))
						)
					      , 1./2.)
					- 1.
				, 2.);

	Re_min 	= 		(pow(delta0,2.) / 4.)
			*	pow(
					pow(1. + 
						(	
							(4. * pow(X_min, 1./2.))
							/
							(pow(delta0, 2.) * pow(C0, 1./2.))
						)
					      , 1./2.)
					- 1.
				, 2.);
	
	CD_maj = C0 *
		pow(1. + (delta0 / pow(Re_maj, 1./2.)),2.);

	CD_min = C0 *
		pow(1. + (delta0 / pow(Re_min, 1./2.)),2.);
	
	scat->particle_terminal_fall_speed 
		= pow((2. * scat->particle_mass_kg *  scat->air_gravity) / (scat->air_density * scat->particle_A_maj * CD_maj)      , 1./2.);
	
	
	scat->particle_inertial_eta_z =  (scat->air_density * scat->particle_A_maj * CD_maj ) / (2. * scat->particle_mass_kg);
	scat->particle_inertial_eta_xy = (scat->air_density * scat->particle_A_min * CD_min ) / (2. * scat->particle_mass_kg);
	
	scat->particle_inertial_distance_z = 2. / scat->particle_inertial_eta_z;
	scat->particle_inertial_distance_xy = (1. - exp(-1.)) / (exp(-1.) * scat->particle_inertial_eta_xy);
}

//function: particles_fall_speed_atlas1973
//input:
//- scat->particle_D_eqvol_mm
//- scat->air_drypressure_hPa
//- scat->air_vaporpressure_hPa
//- scat->air_density

//output:
//- scat->particle_terminal_fall_speed
void particles_fall_speed_atlas1973(t_zephyros_particles_widget *scat)
{
	double alpha = 9.65;
	double beta = 10.3;	
	double specific_gas_constant_Rd = 287.058; //Specific gas constant for dry air, 287.058 J/(kg·K)
	double specific_gas_constant_Rv = 461.495; //Specific gas constant for water vapor, 461.495 J/(kg·K)
	double rho0;
	
	rho0 = ((100. * scat->air_drypressure_hPa) / (specific_gas_constant_Rd * 293.15)) + ((100. * scat->air_vaporpressure_hPa) / (specific_gas_constant_Rv * 293.15));
	scat->particle_terminal_fall_speed = alpha - (beta * exp(-0.6 * scat->particle_D_eqvol_mm));
	scat->particle_terminal_fall_speed *= pow(rho0/scat->air_density, 0.4);
}




//function: particles_cross_sections_dewolf1990
//input:
//- scat->particle_axis_ratio
//- scat->particle_refractive_index
//-	coor->radar_pol_hor_dir
//-	coor->radar_pol_ver_dir
//-	scat->particle_dir

//output:
//- scat->particle_terminal_fall_speed
void particles_cross_sections_dewolf1990(t_zephyros_particles_widget *scat, t_zephyros_coordinates *coor)
{
	double e, esq, f, fsq;
	double shapefactor_lambda12, shapefactor_lambda3;
	
	//epsilon_r defined as the real part of relative permitivity
	//relative permitivity = n^2, where n is the complex refractive index of water
	//+/- 80 for water droplets, lambda = 9 cm.
	double epsilon_r = creal(scat->particle_refractive_index * scat->particle_refractive_index);
	
	//solution taken over from Yanovsky, GRSPaper hr.pdf
	if (scat->particle_axis_ratio < 1.) {
		esq = 1. - pow(scat->particle_axis_ratio, 2.);
		e = sqrt(esq);
		shapefactor_lambda3 =
			((1. - esq) / esq) * (-1. + ((1. / (2. * e)) * log((1. + e)/ (1. - e))));
	} else if (scat->particle_axis_ratio == 1.) {
		shapefactor_lambda3 = 1./3.;
	} else {
		fsq = pow(scat->particle_axis_ratio, 2.) - 1.;
		f = sqrt(fsq);
		shapefactor_lambda3 =
			((1. + fsq) / fsq) * (1. - ((1. / f) * atan(f)));
	}

	shapefactor_lambda12	= (1. - shapefactor_lambda3) / 2.;

	scat->dewolf1990_shapefactor_biglambda12	= 1. / (1. + (shapefactor_lambda12 * (epsilon_r - 1.)));
	scat->dewolf1990_shapefactor_biglambda3		= 1. / (1. + (shapefactor_lambda3 * (epsilon_r - 1.)));
	scat->dewolf1990_fct 						= (pow(M_PI, 5.) * pow(scat->particle_D_eqvol_mm * 1.e-3, 6.) / (9. * pow(scat->radar_wavelength_m, 4.))) * pow(epsilon_r - 1., 2.);

	particles_cross_sections_dewolf1990_update_coor(scat, coor);
}

void particles_cross_sections_dewolf1990_update_coor(t_zephyros_particles_widget *scat, t_zephyros_coordinates *coor)
{
	double qhh, qvv, qhv;
	double eh_u3_sq, ev_u3_sq;
	
	//inproducts
	eh_u3_sq = pow(		(coor->radar_pol_hor_dir[0] * scat->particle_dir[0])
						+ (coor->radar_pol_hor_dir[1] * scat->particle_dir[1])
						+ (coor->radar_pol_hor_dir[2] * scat->particle_dir[2]), 2.);
	ev_u3_sq = pow(		(coor->radar_pol_ver_dir[0] * scat->particle_dir[0])
						+ (coor->radar_pol_ver_dir[1] * scat->particle_dir[1])
						 + (coor->radar_pol_ver_dir[2] * scat->particle_dir[2]), 2.);

	qhh = pow( (scat->dewolf1990_shapefactor_biglambda3 - scat->dewolf1990_shapefactor_biglambda12) * eh_u3_sq + scat->dewolf1990_shapefactor_biglambda12 , 2.);
	qvv = pow( (scat->dewolf1990_shapefactor_biglambda3 - scat->dewolf1990_shapefactor_biglambda12) * ev_u3_sq + scat->dewolf1990_shapefactor_biglambda12 , 2.);
	qhv =  eh_u3_sq * ev_u3_sq * pow(scat->dewolf1990_shapefactor_biglambda3 - scat->dewolf1990_shapefactor_biglambda12, 2.);

	scat->particle_sigma_hh = scat->dewolf1990_fct * qhh;
	scat->particle_sigma_vv = scat->dewolf1990_fct * qvv;
	scat->particle_sigma_hv = scat->dewolf1990_fct * qhv;
}














//function: particles_cross_sections_mischenko2000
//input:
//- scat->particle_axis_ratio
//- scat->particle_refractive_index
//-	coor->radar_pol_hor_dir
//-	coor->radar_pol_ver_dir
//-	scat->particle_dir

//output:
//- scat->particle_terminal_fall_speed
void particles_cross_sections_mischenko2000(t_zephyros_particles_widget *scat, t_zephyros_coordinates *coor, int i_backward_forward)
{
	//TBD: not sure about units and so on ...
	//here everything given in m.

	//calculate t-matrix
	scat->mischenko2000_wid.axi 	= 1.e-3 * scat->particle_D_eqvol_mm / 2.;		//equivalent-sphere radius   
	scat->mischenko2000_wid.rat 	= 1;
	scat->mischenko2000_wid.lam 	= scat->radar_wavelength_m;
	scat->mischenko2000_wid.mrr 	= creal(scat->particle_refractive_index);
	scat->mischenko2000_wid.mri 	= -1. * cimag(scat->particle_refractive_index);
	scat->mischenko2000_wid.eps 	= 1. / scat->particle_axis_ratio;
	scat->mischenko2000_wid.np 		= -1;	
	//scat->mischenko2000_wid.ddelt 	= .001;	//standard value
	scat->mischenko2000_wid.ddelt 	= .005;	
	scat->mischenko2000_wid.ndgs 	= 2;

	//calculate t-matrix
	mischenko2000_calc_tmatrix(&(scat->mischenko2000_wid));

	particles_cross_sections_mischenko2000_update_coor(scat, coor, i_backward_forward);
}

void particles_cross_sections_mischenko2000_update_coor(t_zephyros_particles_widget *scat, t_zephyros_coordinates *coor, int i_backward_forward)
{	
	//set angles
	//beta is the angle w.r.t. the z-axis
	//alpha is the angle w.r.t. the x-axis, positive in the direction of the y-axis
	scat->mischenko2000_wid.beta 	= fmod((180. / M_PI) * fabs(acos(scat->particle_dir[2])), 180.);	
	scat->mischenko2000_wid.alpha 	= fmod((180. / M_PI) * atan2(scat->particle_dir[1], scat->particle_dir[0]) + 360.,360.);	
	
	//theta is the angle w.r.t. the vertical, the z-axis
	//phi is the angle w.r.t. x-axis, positive in direction of the y-axis.
	scat->mischenko2000_wid.thet0 	= (180. / M_PI) * fabs(acos(coor->radar_enu_dir[2]));
	scat->mischenko2000_wid.phi0 	= fmod((180. / M_PI) * atan2(coor->radar_enu_dir[1], coor->radar_enu_dir[0]) + 360.,360.);  
	
	if (i_backward_forward == 0) {
		//for back scattering	
		scat->mischenko2000_wid.thet 	= (180. / M_PI) * fabs(acos(-coor->radar_enu_dir[2]));
		scat->mischenko2000_wid.phi 	= fmod((180. / M_PI) * atan2(-coor->radar_enu_dir[1], -coor->radar_enu_dir[0]) + 360.,360.);
	}
	
	if (i_backward_forward == 1) {
		//for forward scattering	
		scat->mischenko2000_wid.thet 	= (180. / M_PI) * fabs(acos(coor->radar_enu_dir[2]));
		scat->mischenko2000_wid.phi 	= fmod((180. / M_PI) * atan2(coor->radar_enu_dir[1], coor->radar_enu_dir[0]) + 360.,360.);
	}
	
	/*
	//extra check for angle
	if (scat->mischenko2000_wid.beta > 180. - delta) 	scat->mischenko2000_wid.beta -= delta;
	if (scat->mischenko2000_wid.beta < 0. + delta) 		scat->mischenko2000_wid.beta += delta;
	if (scat->mischenko2000_wid.alpha > 360. - delta) 	scat->mischenko2000_wid.alpha -= delta;
	if (scat->mischenko2000_wid.alpha < 0. + delta) 	scat->mischenko2000_wid.alpha += delta;
	
	if (scat->mischenko2000_wid.thet0 > 180. - delta) 	scat->mischenko2000_wid.thet0 -= delta;
	if (scat->mischenko2000_wid.thet0 < 0. + delta) 	scat->mischenko2000_wid.thet0 += delta;
	if (scat->mischenko2000_wid.phi0 > 360. - delta) 	scat->mischenko2000_wid.phi0 -= delta;
	if (scat->mischenko2000_wid.phi0 < 0. + delta) 		scat->mischenko2000_wid.phi0 += delta;
	
	if (scat->mischenko2000_wid.thet > 180. - delta) 	scat->mischenko2000_wid.thet -= delta;
	if (scat->mischenko2000_wid.thet < 0. + delta) 		scat->mischenko2000_wid.thet += delta;
	if (scat->mischenko2000_wid.phi > 360. - delta) 	scat->mischenko2000_wid.phi -= delta;
	if (scat->mischenko2000_wid.phi < 0. + delta) 		scat->mischenko2000_wid.phi += delta;
	*/
	
	//calculate amplitude and phase matrices
	mischenko2000_amp_phase_matrices(&(scat->mischenko2000_wid));

	//radar cross section 3.10.1 from Bringi
	//or 3.134
	//eta_hh = n * 4. * M_PI * pow(cabs(wid->s11), 2.);
	//eta_vv = n * 4. * M_PI * pow(cabs(wid->s22), 2.);
	//eta_hv = n * 4 * M_PI * pow(cabs(wid->s12), 2.);
	
	if (i_backward_forward == 0) {
		scat->particle_sigma_hh = 4. * M_PI * pow(cabs(scat->mischenko2000_wid.s22), 2.);
		scat->particle_sigma_vv = 4. * M_PI * pow(cabs(scat->mischenko2000_wid.s11), 2.);
		scat->particle_sigma_hv = 4. * M_PI * pow(cabs(scat->mischenko2000_wid.s21), 2.);
		scat->particle_sigma_vh = 4. * M_PI * pow(cabs(scat->mischenko2000_wid.s12), 2.);
		
		//extra to calculate correlation coefficients
		scat->particle_sigma_ShhSvvc = 4. * M_PI * cabs(scat->mischenko2000_wid.s22 * conj(scat->mischenko2000_wid.s11));
		scat->particle_sigma_ShhShvc = 4. * M_PI * cabs(scat->mischenko2000_wid.s22 * conj(scat->mischenko2000_wid.s21));
		scat->particle_sigma_SvvSvhc = 4. * M_PI * cabs(scat->mischenko2000_wid.s11 * conj(scat->mischenko2000_wid.s12));	
	}
	
	if (i_backward_forward == 1) {	
		//extra to calculate differential phase
		scat->particle_Re_Shh_min_Svv = creal(scat->mischenko2000_wid.s22 - scat->mischenko2000_wid.s11);
	}
}





void particles_inertiamodel_fraction(
	double *omega,
	double *etaI,
	double *fraction)
{
	 double A0 = 0.4346;
	 double A1 = 1.711;
	 double x = *omega * *etaI;
	 
	if (x == 0.) {
		*fraction = 0.;
	} else {
		*fraction = ((1. - (0.5 * tanh(1.))) * tanh(pow(A0 * x, -1.)))
                + (0.5 * tanh(exp(-1. * A1 * x)));
	}
}

void particles_inertiamodel_phaselag(
	double *omega,
	double *etaI,
	double *phaselag)
{
	double B0 = 1.252;
	double B1 = 0.3662;
	double B2 = 0.5390;
	double B3 = 0.05852;
	double B4 = 0.3518;
	double x = *omega * *etaI;
	
    if (x == 0.) {
        *phaselag = 0.;
    } else {
        *phaselag = (
                    (B0 * tanh(pow(B1 * x, B2)))
                    +
                    (((M_PI / 2.) - B0) * tanh(pow(B3 * x, B4)))
                    );
	}
}
    


void particle_print_widget(t_zephyros_particles_widget *scat)
{
	printf("\n\n\n");
	printf("%-30s = %.2e\n", "air_totalpressure_hPa", scat->air_totalpressure_hPa);
	printf("%-30s = %.2e\n", "air_drypressure_hPa", scat->air_drypressure_hPa);
	printf("%-30s = %.2e\n", "air_vaporpressure_hPa", scat->air_vaporpressure_hPa);
	printf("%-30s = %.2e\n", "air_temperature_K", scat->air_temperature_K);
	printf("%-30s = %.2e\n", "air_dynamic_viscosity", scat->air_dynamic_viscosity);
	printf("%-30s = %.2e\n", "air_kinematic_viscosity", scat->air_kinematic_viscosity);
	printf("%-30s = %.2e\n", "air_density", scat->air_density);
	printf("%-30s = %.2e\n", "air_gravity", scat->air_gravity);
	
	printf("%-30s = %i\n", "particle_type", scat->particle_type);
	printf("%-30s = %.2e\n", "particle_density", scat->particle_density);
	printf("%-30s = %.2e\n", "particle_axis_ratio", scat->particle_axis_ratio);
	printf("%-30s = %.2e\n", "particle_D_maj_mm", scat->particle_D_maj_mm);
	printf("%-30s = %.2e\n", "particle_D_eqvol_mm", scat->particle_D_eqvol_mm);
	printf("%-30s = %.2e\n", "particle_D_min_mm", scat->particle_D_min_mm);
	printf("%-30s = %.2e\n", "particle_mass_kg", scat->particle_mass_kg);
	printf("%-30s = %.2e\n", "particle_A_maj", scat->particle_A_maj);
	printf("%-30s = %.2e\n", "particle_A_min", scat->particle_A_min);
	printf("%-30s = %.2e\n", "particle_volume", scat->particle_volume);
	printf("%-30s = (%.2e, %.2e, %.2e)\n", "particle_dir", scat->particle_dir[0], scat->particle_dir[1], scat->particle_dir[2]);
	printf("%-30s = %.2e\n", "particle_terminal_fall_speed", scat->particle_terminal_fall_speed);
	printf("%-30s = %.2e\n", "particle_inertial_eta_z", scat->particle_inertial_eta_z);
	printf("%-30s = %.2e\n", "particle_inertial_eta_xy", scat->particle_inertial_eta_xy);
	printf("%-30s = %.2e\n", "particle_inertial_distance_z", scat->particle_inertial_distance_z);
	printf("%-30s = %.2e\n", "particle_inertial_distance_xy", scat->particle_inertial_distance_xy);
	
	
	printf("%-30s = %.2e + %.2ei\n", "particle_refractive_index", creal(scat->particle_refractive_index), cimag(scat->particle_refractive_index));
	printf("%-30s = %.2e\n", "particle_sigma_hh", scat->particle_sigma_hh);
	printf("%-30s = %.2e\n", "particle_sigma_hv", scat->particle_sigma_hv);
	printf("%-30s = %.2e\n", "particle_sigma_vv", scat->particle_sigma_vv);
	printf("%-30s = %.2e\n", "particle_attenuation", scat->particle_attenuation);
	printf("%-30s = %.2e\n", "dewolf1990_shapefactor_biglambda12", scat->dewolf1990_shapefactor_biglambda12);
	printf("%-30s = %.2e\n", "dewolf1990_shapefactor_biglambda3", scat->dewolf1990_shapefactor_biglambda3);
	printf("%-30s = %.2e\n", "dewolf1990_fct", scat->dewolf1990_fct);
	printf("%-30s = %.2e\n", "radar_wavelength_m", scat->radar_wavelength_m);
	
}

