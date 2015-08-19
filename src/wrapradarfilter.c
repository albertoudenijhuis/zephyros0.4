#include <math.h>
#include <string.h> 
#include <stdlib.h>
#include <stdio.h>

#include "radarfilter.h"
#include "wrapradarfilter.h"
#include "zephyros_config.h"
#include "util.h"

void wrapradarfilter(char cfg_filename[8192], char additional_output_filename[8192], char measurements_file_name[8192])
{
	int i;
	int n;
	
	t_zephyros_config   		*cfg = NULL;
	t_radarmeasurement			**radarmeasurement;
	t_radarfilter_todolist		*todo;
	
	//load configuration
	zephyros_config_read(cfg_filename, additional_output_filename, &cfg);
	
	//initialize
	radarfilter_initialize_todolist(cfg, 0, &todo);
	radarfilter_read_measurements(&n, &radarmeasurement, measurements_file_name);
	radarfilter_exec(cfg, 0, todo, n, radarmeasurement);
	
	radarfilter_write_measurements(cfg, 0, n, radarmeasurement, cfg->fp_ao);
	if (cfg->general->additional_output->print_detailed_analysis) radarfilter_write_measurements_detailed_analysis(cfg, 0, n, radarmeasurement, cfg->fp_ao);
	
	radarfilter_free_todolist(&todo);
	radarfilter_free_radarmeasurement(n, &radarmeasurement);	
	zephyros_config_free(&cfg);	
}
