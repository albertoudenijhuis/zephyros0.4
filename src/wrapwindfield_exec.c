#include "wrapwindfield.h"
#include <stdlib.h>
#include <stdio.h>
#include "func.h"

int main () {
	char line[8192];
	char cfg_filename[8192];
	char additional_output_filename[8192];
	int n;
	int i;
	double *x;
	double *y;
	double *z;
	double *t;
	double *u;
	double *v;
	double *w;
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
		if ( strcmp(identifier,"x") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			x = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", x + i);
			}
		}
		if ( strcmp(identifier,"y") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			y = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", y + i);
			}
		}
		if ( strcmp(identifier,"z") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			z = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", z + i);
			}
		}
		if ( strcmp(identifier,"t") == 0 ) {
			fsetpos(fp, &pos_thisline);
			fscanf(fp, "%s %i", identifier, &n );
			//allocate
			t = malloc(n * sizeof(double));
			for ( i = 0; i < n; i++ ) {
				fscanf(fp, "%lf", t + i);
			}
		}
		fsetpos(fp, &pos_nextline);		//set position to next line	
	}
	fclose(fp);

	u = malloc(n * sizeof(double));
	v = malloc(n * sizeof(double));
	w = malloc(n * sizeof(double));

	wrapwindfield(cfg_filename, additional_output_filename, n, x, y, z, t, u, v, w);

    fp = fopen("output.z", "w");

	fprintf(fp, "%-30s %-15i", 
				"u",
				n);
	for (i=0; i<n; i++) {
		fprintf(fp, " %-15.3e", pnan(u + i));
	}
	fprintf(fp, "\n");

	fprintf(fp, "%-30s %-15i", 
				"v",
				n);
	for (i=0; i<n; i++) {
		fprintf(fp, " %-15.3e", pnan(v + i));
	}
	fprintf(fp, "\n");

	fprintf(fp, "%-30s %-15i", 
				"w",
				n);
	for (i=0; i<n; i++) {
		fprintf(fp, " %-15.3e", pnan(w + i));
	}
	fprintf(fp, "\n");
	
	fclose(fp);
	
	
    return 0;
}