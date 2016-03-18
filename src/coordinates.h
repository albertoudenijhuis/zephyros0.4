#ifndef _ZEPHYROS_COORDINATES
#define _ZEPHYROS_COORDINATES

/*
WGS84 (World Geodetic System of 1984)
wgs84_phirad  		latitude in rad
wgs84_lambdarad  	longitude in rad
wgs84_h 			height in m

ECEF system (Earth-centred earth-fixed)
ecef_x in m
ecef_y in m
ecef_z in m

ENU system (EAST NORHT UP)
enu_x in m
enu_y in m
enu_z in m
*/

typedef struct st_zephyros_coordinates			
{
	//WGS84
	double	wgs84_phirad;
	double	wgs84_lambdarad;
	double	wgs84_h;
	double 	wgs84_radar_location_phirad;
	double 	wgs84_radar_location_lambdarad;
	double 	wgs84_radar_location_h;

	//ECEF
	double	*ecef_xyzt;
	double 	*ecef_radar_location_xyzt;

	//ENU
	double	ref_wgs84_phirad;
	double	ref_wgs84_lambdarad;
	double	ref_wgs84_h;
	double	*enu_xyzt;
	double  *enu_radar_location_xyzt;
	
	//Radar coordinates, AZEL+range
	int		radar_effearthcorrection;
	double	radar_azel_alpha;		/* azimuth angle, w.r.t north, positive towards east [rad] */
	double	radar_azel_gamma;		/* elevation angle [rad] */
	double  radar_azel_gamma_cor; 	//required for effearthcorrection
	double	radar_range;

	//Radar coordinates, ENU direction
	double	*radar_enu_dir; //radar_enu_dx, radar_enu_dy, radar_enu_dz
	
	//Radar coordinates, BEAM coordinates
	double	radar_beam_theta;
	double	radar_beam_phi;
	
	//Radar, Polarization direction
	double	*radar_pol_hor_dir;
	double	*radar_pol_ver_dir;
	
	char		name[8192];	//for error reporting	
	int initialized;
} t_zephyros_coordinates;

void coordinates_initialize_coor(t_zephyros_coordinates **pcoor);
void coordinates_assert_initialized(t_zephyros_coordinates *coor);
void coordinates_free_coor(t_zephyros_coordinates **pcoor);

void coordinates_geodetic2ecef(t_zephyros_coordinates *coor);
void coordinates_ecef2geodetic(t_zephyros_coordinates *coor, int method);
void coordinates_ecef2enu(t_zephyros_coordinates *coor);
void coordinates_enu2ecef(t_zephyros_coordinates *coor);
void coordinates_radar_azel2enu(t_zephyros_coordinates *coor);
void coordinates_radar_azelrangedir2enu(t_zephyros_coordinates *coor);
void coordinates_radar_enu2azel(t_zephyros_coordinates *coor);
void coordinates_radar_beam2azel(t_zephyros_coordinates *coor);
void coordinates_radar_azel2beam(t_zephyros_coordinates *coor);
void coordinates_radar_pol_dir(t_zephyros_coordinates *coor);

void coordinates_print(t_zephyros_coordinates *coor);

#endif
