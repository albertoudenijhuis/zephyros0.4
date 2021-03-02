#include <math.h>
#include <string.h> 
#include <stdlib.h>
#include <stdio.h>

#include "wrapwindfield.h"
#include "zephyros_config.h"

void wrapwindfield(char cfg_filename[8192], char additional_output_filename[8192], int n, double *x, double *y, double *z, double *t, double *u, double *v, double *w)
	{
	int i;
	double *dummy;
	double *xyzt = calloc(4, sizeof(double));

	t_zephyros_config *cfg = NULL;

	//load configuration
	zephyros_config_read(cfg_filename, additional_output_filename, &cfg);

	for ( i = 0; i < n; i++ ) {
		printf("wrapwindfield, step %i/%i\n", i, n);
		xyzt[0] = x[i];
		xyzt[1] = y[i];
		xyzt[2] = z[i];
		xyzt[3] = t[i];
		
		util_windfield_fuvw(cfg->simulation->windfield, xyzt, 
			u + i, v + i, w + i,
			0,
			dummy, dummy, dummy,
			0,
			1.
			);	
	}

	if (cfg->general->additional_output->print_configuration) zephyros_config_print(cfg, cfg->fp_ao);

	zephyros_config_free(&cfg);
	free(xyzt);
}
