#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "retrieval_lwm.h"
#include "retrieval_fdvar.h"
#include "wrapretrieval.h"
#include "zephyros_config.h"

//uncomment next statement for debug mode
#define _ZEPHYROS_WRAPRETRIEVAL_DEBUG
	
void wrapretrieval(char cfg_filename[8192], char additional_output_filename[8192], char measurements_filename[8192])

{
	int iretr;

	t_lwm_opc 	*lwm_opc		= malloc(sizeof(t_lwm_opc));
	t_fdvar_opc *fdvar_opc		= malloc(sizeof(t_fdvar_opc));
	
	t_fdvar_o 	*fdvar_o;
	t_lwm_o 	*lwm_o;
	
	t_zephyros_config *cfg;
	
	int lwm_retrieval_executed = 0;
	int fdvar_retrieval_executed = 0;
	
	//load configuration
	zephyros_config_read(cfg_filename, additional_output_filename, &cfg);
	lwm_opc->zc = cfg;
	fdvar_opc->zc = cfg;
	
	//initialize & read observations
	retrieval_lwm_initialize_o(&lwm_o, measurements_filename);
	lwm_opc->o = lwm_o;

	retrieval_fdvar_initialize_o(&fdvar_o, measurements_filename);
	fdvar_opc->o = fdvar_o;

	for ( iretr = 0; iretr <= 100; iretr++ ) {
		if (cfg->retrieval->algorithm->type[iretr] == 1) {
			//set configuration
			lwm_opc->c = cfg->retrieval->algorithm->lwm_cfg[iretr];
					
			//do the retrieval 
			retrieval_lwm_apply(lwm_opc);	
			
			lwm_retrieval_executed = 1;		
		}
	
		if (cfg->retrieval->algorithm->type[iretr] == 2) {
			//set configuration
			fdvar_opc->c = cfg->retrieval->algorithm->fdvar_cfg[iretr];
			
			//do the retrieval 
			retrieval_fdvar_apply(fdvar_opc);
			
			fdvar_retrieval_executed = 1;
		}
	}

	#ifdef _ZEPHYROS_WRAPRETRIEVAL_DEBUG
		printf("print additional output\n"); fflush(stdout);
	#endif 
	if (cfg->general->additional_output->print_configuration) {
		zephyros_config_print(cfg, cfg->fp_ao);
	}
	//print additional output
	if (lwm_retrieval_executed)  {retrieval_lwm_additional_output(lwm_opc, cfg->fp_ao);}
	if (fdvar_retrieval_executed)  {retrieval_fdvar_additional_output(fdvar_opc, cfg->fp_ao);}


	#ifdef _ZEPHYROS_WRAPRETRIEVAL_DEBUG
		printf("clean up\n"); fflush(stdout);
	#endif 
				
	retrieval_lwm_free_o(&lwm_o);
	retrieval_fdvar_free_o(&fdvar_o);
	
	free(lwm_opc);
	free(fdvar_opc);

	#ifdef _ZEPHYROS_WRAPRETRIEVAL_DEBUG
		printf("free zephyros config\n"); fflush(stdout);
	#endif 
		
	zephyros_config_free(&cfg);	
	
	#ifdef _ZEPHYROS_WRAPRETRIEVAL_DEBUG
		printf("Done.\n"); fflush(stdout);
	#endif 	
}

