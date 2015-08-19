%module coordinates

%{
    #define SWIG_FILE_WITH_INIT
	#include "coordinates.h"
%}

%include "numpy.i"

%init %{
    import_array();
%}


%apply (int DIM1, double* IN_ARRAY1) {
						(int len100, double* azel_alpha100),
						(int len101, double* azel_gamma101),
						(int len102, double* azel_r102),
						(int len103, double* effearthcorrection103)
							};
%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len200, double* enu_x200),
						(int len201, double* enu_y201),
						(int len202, double* enu_z202)
						};

%rename (azel2enu) my_azel2enu;
%inline %{
    double my_azel2enu(
						int len100, double* azel_alpha100,
						int len101, double* azel_gamma101,
						int len102, double* azel_r102,
						int len103, double* effearthcorrection103,
						int len200, double* enu_x200,
						int len201, double* enu_y201,
						int len202, double* enu_z202
						) {
		int i;

		t_zephyros_coordinates *coor = malloc(sizeof(t_zephyros_coordinates));
		
		for ( i = 0; i < len100; i++ ) {
			coor->radar_azel_alpha 	= azel_alpha100[i];
			coor->radar_azel_gamma 	= azel_gamma101[i];
			coor->radar_range 		= azel_r102[i];
			coor->radar_effearthcorrection = (int) effearthcorrection103[i];
			coordinates_radar_azel2enu(coor);
			enu_x200[i] = coor->enu_xyzt[0];
			enu_y201[i] = coor->enu_xyzt[1];
			enu_z202[i] = coor->enu_xyzt[2];
		}
		
		free(coor);
	}
%}




%apply (int DIM1, double* IN_ARRAY1) {
						(int len300, double* azel_alpha300),
						(int len301, double* azel_gamma301),
						(int len302, double* azel_r302),
						(int len303, double* effearthcorrection303)
							};
%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len400, double* enu_dx400),
						(int len401, double* enu_dy401),
						(int len402, double* enu_dz402)
						};

%rename (azelrangedir2enu) my_azelrangedir2enu;
%inline %{
    double my_azelrangedir2enu(
						int len300, double* azel_alpha300,
						int len301, double* azel_gamma301,
						int len302, double* azel_r302,
						int len303, double* effearthcorrection303,
						int len400, double* enu_dx400,
						int len401, double* enu_dy401,
						int len402, double* enu_dz402
						) {
		int i;

		t_zephyros_coordinates *coor = malloc(sizeof(t_zephyros_coordinates));
		
		for ( i = 0; i < len300; i++ ) {
			coor->radar_azel_alpha 	= azel_alpha300[i];
			coor->radar_azel_gamma 	= azel_gamma301[i];
			coor->radar_range 		= azel_r302[i];
			coor->radar_effearthcorrection = (int) effearthcorrection303[i];
			coordinates_radar_azelrangedir2enu(coor);
			enu_dx400[i] = coor->radar_enu_dir[0];
			enu_dy401[i] = coor->radar_enu_dir[1];
			enu_dz402[i] = coor->radar_enu_dir[2];
		}
		
		free(coor);
	}
%}






%apply (int DIM1, double* IN_ARRAY1) {
						(int len500, double* enu_x500),
						(int len501, double* enu_y501),
						(int len502, double* enu_z502),
						(int len503, double* effearthcorrection503)
							};
%apply (int DIM1, double* ARGOUT_ARRAY1) {
						(int len600, double* azel_alpha600),
						(int len601, double* azel_gamma601),
						(int len602, double* azel_r602)
						};

%rename (enu2azel) my_enu2azel;
%inline %{
    double my_enu2azel(
						int len500, double* enu_x500,
						int len501, double* enu_y501,
						int len502, double* enu_z502,
						int len503, double* effearthcorrection503,
						int len600, double* azel_alpha600,
						int len601, double* azel_gamma601,
						int len602, double* azel_r602
						) {
		int i;

		t_zephyros_coordinates *coor = malloc(sizeof(t_zephyros_coordinates));
		
		for ( i = 0; i < len500; i++ ) {
			coor->enu_xyzt[0] = enu_x500[i];
			coor->enu_xyzt[1] = enu_y501[i];
			coor->enu_xyzt[2] = enu_z502[i];
			coor->radar_effearthcorrection = (int) effearthcorrection503[i];
			coordinates_radar_enu2azel(coor);
			azel_alpha600[i] = coor->radar_azel_alpha;
			azel_gamma601[i] = coor->radar_azel_gamma;
			azel_r602[i] = coor->radar_range;
		}
		
		free(coor);
	}
%}



