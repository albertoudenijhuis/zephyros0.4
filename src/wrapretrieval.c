#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "retrieval_lwm.h"
#include "retrieval_fdvar.h"
#include "wrapretrieval.h"
#include "zephyros_config.h"

void wrapretrieval(char cfg_filename[8192], char additional_output_filename[8192], char measurements_filename[8192])

{
	int iretr;

	t_lwm_opc 	*lwm_opc		= malloc(sizeof(t_lwm_opc));
	t_fdvar_opc *fdvar_opc		= malloc(sizeof(t_fdvar_opc));
	t_lwm_opc 	*lwm_opc_prev	= malloc(sizeof(t_lwm_opc));
	t_fdvar_opc *fdvar_opc_prev	= malloc(sizeof(t_fdvar_opc));
	
	t_fdvar_o 	*fdvar_o;
	t_lwm_o 	*lwm_o;
	
	t_zephyros_config *cfg;

	
	//load configuration
	zephyros_config_read(cfg_filename, additional_output_filename, &cfg);
	lwm_opc->zc = cfg;
	lwm_opc_prev->zc = cfg;
	fdvar_opc->zc = cfg;
	fdvar_opc_prev->zc = cfg;
	
	//initialize & read observations
	retrieval_lwm_initialize_o(&lwm_o, measurements_filename);
	lwm_opc->o = lwm_o;
	lwm_opc_prev->o = lwm_o;

	retrieval_fdvar_initialize_o(&fdvar_o, measurements_filename);
	fdvar_opc->o = fdvar_o;
	fdvar_opc_prev->o = fdvar_o;

	//walk through retrievals
	lwm_opc_prev->p = NULL;
	fdvar_opc_prev->p = NULL;
	for ( iretr = 0; iretr <= 100; iretr++ ) {
		if (cfg->retrieval->algorithm->type[iretr] == 1) {
			//set configuration
			lwm_opc->c = cfg->retrieval->algorithm->lwm_cfg[iretr];
					
			//Cast previous solution (post) to current prior
			if (lwm_opc_prev->p != NULL) {printf("not implemented yet..."); exit(0);}
			if (fdvar_opc_prev->p != NULL) {printf("not implemented yet..."); exit(0);}

			//do the retrieval 
			retrieval_lwm_apply(lwm_opc);
			
			//cleanup
			if (lwm_opc_prev->p != NULL)  {retrieval_lwm_free_p(lwm_opc_prev);}
			if (fdvar_opc_prev->p != NULL)  {retrieval_fdvar_free_p(fdvar_opc_prev);}
			lwm_opc_prev->p = lwm_opc->p;
		}
	
		if (cfg->retrieval->algorithm->type[iretr] == 2) {
			//set configuration
			fdvar_opc->c = cfg->retrieval->algorithm->fdvar_cfg[iretr];
			
			//Cast previous solution (post) to current prior
			if (lwm_opc_prev->p != NULL) {printf("not implemented yet..."); exit(0);}
			if (fdvar_opc_prev->p != NULL) {printf("not implemented yet..."); exit(0);}
			//if (fdvar_opc_prev->p != NULL) {retrieval_fdvar_cast_p(fdvar_opc_prev, fdvar_opc);}

			//do the retrieval 
			retrieval_fdvar_apply(fdvar_opc);

			//cleanup
			if (lwm_opc_prev->p != NULL)  {retrieval_lwm_free_p(lwm_opc_prev);}
			if (fdvar_opc_prev->p != NULL)  {retrieval_fdvar_free_p(fdvar_opc_prev);}
			fdvar_opc_prev->p = fdvar_opc->p;
		}
	}

	//print additional output
	if (lwm_opc_prev->p != NULL)  {retrieval_lwm_additional_output(lwm_opc, cfg->fp_ao);}
	if (fdvar_opc_prev->p != NULL)  {retrieval_fdvar_additional_output(fdvar_opc, cfg->fp_ao);}
		
	//cleanup
	if (lwm_opc_prev->p != NULL)  {retrieval_lwm_free_p(lwm_opc_prev);}
	if (fdvar_opc_prev->p != NULL) {retrieval_fdvar_free_p(fdvar_opc_prev);}
	
	retrieval_lwm_free_o(lwm_o);
	//retrieval_fdvar_free_o(fdvar_o);
	
	free(lwm_opc);
	free(fdvar_opc);
		
	zephyros_config_free(&cfg);	
}

