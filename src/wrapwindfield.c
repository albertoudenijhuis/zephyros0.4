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
	double *xyzt = malloc(4 * sizeof(double));

	t_zephyros_config *cfg = NULL;

	//load configuration
	zephyros_config_read(cfg_filename, additional_output_filename, &cfg);

	for ( i = 0; i < n; i++ ) {
		printf("wrapwindfield, step %i/%i\n", i, n);
		xyzt[0] = x[i];
		xyzt[1] = y[i];
		xyzt[2] = z[i];
		xyzt[3] = t[i];
		
		util_windfield_fuvw(cfg->simulation->windfield, xyzt, 0, u + i, 0, dummy);
		util_windfield_fuvw(cfg->simulation->windfield, xyzt, 1, v + i, 0, dummy);
		util_windfield_fuvw(cfg->simulation->windfield, xyzt, 2, w + i, 0, dummy);
	}

	zephyros_config_free(&cfg);
	free(xyzt);
}
