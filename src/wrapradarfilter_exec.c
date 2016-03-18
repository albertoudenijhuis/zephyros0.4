#include "wrapradarfilter.h"
#include <stdio.h>
#include <string.h> 

int main ( int argc, char *argv[] )
{
	char cfg_filename[8192];
	char additional_output_filename[8192];
	char measurements_file_name[8192];

    if ( argc != 4 ) {
		printf("Incorrect number of arguments. Exiting\n");
	} else {
		strcpy(cfg_filename, argv[1]);
		strcpy(additional_output_filename, argv[2]);
		strcpy(measurements_file_name, argv[3]);
		wrapradarfilter(cfg_filename, additional_output_filename, measurements_file_name);
	}

    return 0;
}
