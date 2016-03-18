#include "zephyros_config.h"
#include "radarfilter.h"
#include "fields.h"

typedef struct st_lwm_o
{
    int 					n;					//number of observations
	t_radarmeasurement 		**radarmeasurement;	//radar observations
 
	//output
	double		*windvector_u;
	double		*windvector_v;
	double		*windvector_w;
	double		*windvector_u_err;
	double		*windvector_v_err;
	double		*windvector_w_err;

	double		*lwm_vr;
	
	int			*used_for_analysis;
	int			*in_analysis_volume;

	int			in_analysis_volume_n;
	int			used_for_analysis_n;
} t_lwm_o ;
void retrieval_lwm_initialize_o(t_lwm_o **p_lwm_o, char measurements_filename[8192]);
void retrieval_lwm_free_o(t_lwm_o **p_lwm_o);

typedef struct st_lwm_p
{
	int			i;			//current index of lwm retrieval
	
	t_zephyros_field 	*field;
	
	int 		Kn; 		//number of fitted parameters
	
	double		*u0;
	double		*u_x;
	double		*u_z;
	double		*v0;
	double		*v_y;
	double		*v_z;
	double		*u_y_plus_v_x;
	double		*w0;
	double		*w_x;
	double		*w_y;
	double		*w_z;
	double		*u_t_plus_v_t_plus_w_t;
	
	//observation space
	double		**coef_u0;
	double		**coef_u_x;
	double		**coef_u_z;
	double		**coef_v0;
	double		**coef_v_y;
	double		**coef_v_z;
	double		**coef_u_y_plus_v_x;
	double		**coef_w0;
	double		**coef_w_x;
	double		**coef_w_y;
	double		**coef_w_z;
	double		**coef_u_t_plus_v_t_plus_w_t;

	double		*center_coef_u0;
	double		*center_coef_u_x;
	double		*center_coef_u_z;
	double		*center_coef_v0;
	double		*center_coef_v_y;
	double		*center_coef_v_z;
	double		*center_coef_u_y_plus_v_x;
	double		*center_coef_w0;
	double		*center_coef_w_x;
	double		*center_coef_w_y;
	double		*center_coef_w_z;
	double		*center_coef_u_t_plus_v_t_plus_w_t;

	int			n_sv; //number of subvolume points
	t_zephyros_coordinates ***subvolume_coor;
	
	int			fit_Knr_u0;
	int			fit_Knr_u_x;
	int			fit_Knr_u_z;
	int			fit_Knr_v0;
	int			fit_Knr_v_y;
	int			fit_Knr_v_z;
	int			fit_Knr_u_y_plus_v_x;
	int			fit_Knr_w0;
	int			fit_Knr_w_x;
	int			fit_Knr_w_y;
	int			fit_Knr_w_z;
	int			fit_Knr_u_t_plus_v_t_plus_w_t;
} t_lwm_p ;

typedef struct st_lwm_opc
{
	t_lwm_o *o;
	t_lwm_p *p;
	t_zephyros_config_retrieval_lwm_cfg *c;
	t_zephyros_config *zc;
} t_lwm_opc;

void retrieval_lwm_initialize_p(t_lwm_opc *opc);
void retrieval_lwm_free_p(t_lwm_opc *opc);

void retrieval_lwm_apply(t_lwm_opc *opc);

int retrieval_lwm_minimize_chisquared(t_lwm_opc *opc);

double retrieval_lwm_chisquared_nlopt(
	unsigned n,
	const double *x,
	double *grad,
	void *params);
	
void retrieval_lwm_chisquared(
	void 	*vd_opc,		//optional arguments
	int		*n,				//number of parameters
	double 	*x,				//parameters
	double 	*f,				//function value
	double 	*fd,			//function derivatives
	int 	*iflag			
	);

void retrieval_lwm_calculate_xyz(t_lwm_opc *opc);
void retrieval_lwm_calculate_coefficients(t_lwm_opc *opc);
void retrieval_lwm_store_results(t_lwm_opc *opc);
void retrieval_lwm_additional_output(t_lwm_opc *opc, FILE *fp);

void retrieval_lwm_unpack_K(t_lwm_opc *opc, double *K);






