#include "wrapretrieval.h"
#include <stdlib.h>
#include <stdio.h>
#include "func.h"

int main () {
	char line[8192];
	int i;
	
	char cfg_filename[8192];
	char additional_output_filename[8192];
	int n;
	double *o_enu_x0;
	double *o_enu_y0;
	double *o_enu_z0;
	double *o_azel_r1_m;
	double *o_azel_r2_m;
	double *o_azel_alpha_rad;
	double *o_azel_gamma_rad;
	double *o_beam_FWHM0_rad;
	double *o_beam_FWHM1_rad;
	double *o_t;
	double *o_dt;
	double *o_refl;
	double *o_srefl;
	double *o_vr;
	double *o_svr;
	double *o_spectralwidth;
	double *o_sspectralwidth;
	double *windvector_u;
	double *windvector_su;
	double *windvector_v;
	double *windvector_sv;
	double *windvector_w;
	double *windvector_sw;

	fpos_t pos_thisline;
	fpos_t pos_nextline;
	char identifier[8192];
	FILE *fp;

	fp = fopen("input.z","r"); // read mode
	while (1) {
		fgetpos(fp, &pos_thisline);		//store position of where this line starts
		//read full line
		if (fgets (line, sizeof(line), fp) == NULL) {
			break;
		}
				
		fgetpos(fp, &pos_nextline);		//store position of where next line starts
			
		//read identifier
		fsetpos(fp, &pos_thisline);
		fscanf(fp, "%s", identifier);			

		if ( strcmp(identifier,"cfg_filename") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %s", identifier, cfg_filename);			
		}
		if ( strcmp(identifier,"additional_output_filename") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %s", identifier, additional_output_filename);			
		}
		if ( strcmp(identifier,"o_enu_x0") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_enu_x0 = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_enu_x0 + i);
			}
		}
		if ( strcmp(identifier,"o_enu_y0") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_enu_y0 = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_enu_y0 + i);
			}
		}
		if ( strcmp(identifier,"o_enu_z0") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_enu_z0 = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_enu_z0 + i);
			}
		}
		if ( strcmp(identifier,"o_azel_r1_m") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_azel_r1_m = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_azel_r1_m + i);
			}
		}
		if ( strcmp(identifier,"o_azel_r2_m") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_azel_r2_m = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_azel_r2_m + i);
			}
		}
		if ( strcmp(identifier,"o_azel_alpha_rad") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_azel_alpha_rad = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_azel_alpha_rad + i);
			}
		}
		if ( strcmp(identifier,"o_azel_gamma_rad") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_azel_gamma_rad = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_azel_gamma_rad + i);
			}
		}
		if ( strcmp(identifier,"o_beam_FWHM0_rad") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_beam_FWHM0_rad = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_beam_FWHM0_rad + i);
			}
		}
		if ( strcmp(identifier,"o_beam_FWHM1_rad") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_beam_FWHM1_rad = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_beam_FWHM1_rad + i);
			}
		}
		if ( strcmp(identifier,"o_t") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_t = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_t + i);
			}
		}
		if ( strcmp(identifier,"o_dt") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_dt = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_dt + i);
			}
		}
		if ( strcmp(identifier,"o_refl") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_refl = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_refl + i);
			}
		}
		if ( strcmp(identifier,"o_srefl") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_srefl = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_srefl + i);
			}
		}
		if ( strcmp(identifier,"o_vr") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_vr = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_vr + i);
			}
		}
		if ( strcmp(identifier,"o_svr") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_svr = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_svr + i);
			}
		}
		if ( strcmp(identifier,"o_spectralwidth") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_spectralwidth = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_spectralwidth + i);
			}
		}
		if ( strcmp(identifier,"o_sspectralwidth") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			o_sspectralwidth = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", o_sspectralwidth + i);
			}
		}
		fsetpos(fp, &pos_nextline);		//set position to next line	
	}
	fclose(fp);

	windvector_u = malloc(n * sizeof(double));
	windvector_su = malloc(n * sizeof(double));
	windvector_v = malloc(n * sizeof(double));
	windvector_sv = malloc(n * sizeof(double));
	windvector_w = malloc(n * sizeof(double));
	windvector_sw = malloc(n * sizeof(double));

	//TBD
	/*
	wrapwindvectors(cfg_filename,
					additional_output_filename,
					n,
					o_enu_x0,
					o_enu_y0,
					o_enu_z0,
					o_azel_r1_m,
					o_azel_r2_m,
					o_azel_alpha_rad,
					o_azel_gamma_rad,
					o_beam_FWHM0_rad,
					o_beam_FWHM1_rad,
					o_t,
					o_dt,
					o_refl,
					o_srefl,
					o_vr,
					o_svr,
					o_spectralwidth,
					o_sspectralwidth,
					windvector_u,
					windvector_su,
					windvector_v,
					windvector_sv,
					windvector_w,
					windvector_sw	
					);
	*/
	
    fp = fopen("output.z", "w");
	
	fprintf(fp, "%-30s %-15i", 
				"windvector_u",
				n);
	for (i=0; i<n; i++) {
		fprintf(fp, " %-15.3e", pnan(windvector_u + i));
	}
	fprintf(fp, "\n");	
	
	fprintf(fp, "%-30s %-15i", 
				"windvector_su",
				n);
	for (i=0; i<n; i++) {
		fprintf(fp, " %-15.3e", pnan(windvector_su + i));
	}
	fprintf(fp, "\n");

	
	fprintf(fp, "%-30s %-15i", 
				"windvector_v",
				n);
	for (i=0; i<n; i++) {
		fprintf(fp, " %-15.3e", pnan(windvector_v + i));
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "%-30s %-15i", 
				"windvector_sv",
				n);
	for (i=0; i<n; i++) {
		fprintf(fp, " %-15.3e", pnan(windvector_sv + i));
	}
	fprintf(fp, "\n");

	
	fprintf(fp, "%-30s %-15i", 
				"windvector_w",
				n);
	for (i=0; i<n; i++) {
		fprintf(fp, " %-15.3e", pnan(windvector_w + i));
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "%-30s %-15i", 
				"windvector_sw",
				n);
	for (i=0; i<n; i++) {
		fprintf(fp, " %-15.3e", pnan(windvector_sw + i));
	}
	fprintf(fp, "\n");
	
	fclose(fp);
	
	
    return 0;
}
