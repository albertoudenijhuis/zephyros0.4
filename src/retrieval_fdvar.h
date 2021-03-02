#include "zephyros_config.h"
#include "fields.h"
#include "radarfilter.h"


double prev_costfunction;

typedef struct st_fdvar_o
{
    int 					n;								//number of radar observations
	t_radarmeasurement 		**radarmeasurement;				//radar observations
	t_radarmeasurement 		**model_radarmeasurement;		//model radar observations
} t_fdvar_o ;

typedef struct st_fdvar_p
{
	int			Kn;			//number of parameters, typically 3 + 4 * nxyz;
	double*		Kubound;	//K upper bound
	double*		Klbound;	//K lower bound
	int			cost_no;	//total number of observations that are fitted, typically 2 x o->n (reflection and vr).

} t_fdvar_p ;

typedef struct st_fdvar_opc
{
	t_fdvar_o 								*o;
	t_fdvar_p 								*p;
	t_zephyros_config_retrieval_fdvar_cfg 	*c;
	t_zephyros_config 						*zc;
	
	t_radarfilter_res_vol 					*res_vol;
	t_radarfilter_todolist					*todo_no_derivatives;	
	t_radarfilter_todolist					*todo_with_derivatives;	
} t_fdvar_opc;




void retrieval_fdvar_apply(t_fdvar_opc *opc);

int retrieval_fdvar_minimize_cost_function(t_fdvar_opc *opc);

double retrieval_fdvar_cost_function_nlopt(
	unsigned n,
	const double *x,
	double *grad,
	void *params);
	
	
void retrieval_fdvar_cost_function(
	void 	*vd_opc,		//optional arguments
	int		*n,				//number of parameters
	double 	*x,				//parameters
	double 	*f,				//function value
	double 	*fd,			//function derivatives
	int 	*iflag			
	);
	



	
double retrieval_fdvar_K2x(t_fdvar_p *p, double *x, int Knr);
double retrieval_fdvar_x2K(t_fdvar_p *p, double *K, int Knr);

void retrieval_fdvar_init_x(void *vd_opc, double **x);
void retrieval_fdvar_unpack_x(t_fdvar_opc *opc, double *x);

void retrieval_fdvar_estimate_posterior_errors(t_fdvar_opc *opc);
void retrieval_fdvar_estimate_posterior_errors_solver1(cs *mat1_triplet, double *prior_errors, double *result_errors);

void retrieval_fdvar_cast_solutions(t_fdvar_opc *opc);

void retrieval_fdvar_additional_output(t_fdvar_opc *opc, FILE *fp);

void retrieval_fdvar_initialize_o(t_fdvar_o **p_fdvar_o, char measurements_filename[8192]);
void retrieval_fdvar_free_o(t_fdvar_o **p_fdvar_o);

void retrieval_fdvar_initialize_p(t_fdvar_opc *opc);
void retrieval_fdvar_free_p(t_fdvar_opc *opc);

void retrieval_fdvar_cs_qrsol(t_zephyros_field_errcovmat *ecm, double *b, double *x);



void cs_matrix_vector (
	cs *A,
	double* x,
	double* y,
	int n
	);


