#include "particles_mishchenko2000.h"

#include <stdio.h>



void mishchenko2000_calc_tmatrix(t_mishchenko2000_widget *wid)
{
	int i,j,k;

	//update fortran io
	mishchenko_io_.axi 	= wid->axi;
	mishchenko_io_.rat 	= wid->rat;
	mishchenko_io_.lam 	= wid->lam;
	mishchenko_io_.mrr 	= wid->mrr;
	mishchenko_io_.mri 	= wid->mri;
	mishchenko_io_.eps 	= wid->eps;
	mishchenko_io_.np 	= wid->np;
	mishchenko_io_.ndgs 	= wid->ndgs;
	mishchenko_io_.ddelt = wid->ddelt;

	//execute fortran code to calculate t-matrix
	calc_tmatrix_();
	
	//save results to widget
	wid->nmax = mishchenko_io_.nmax;		
	
	for (i=0; i < wid->nmax; i++) {
		for (j=0; j < wid->nmax; j++) {
			for (k=0; k < (wid->nmax + 1); k++) {			
				wid->rt11[i][j][k] = tmat_.rt11[i][j][k];
				wid->rt12[i][j][k] = tmat_.rt12[i][j][k];
				wid->rt21[i][j][k] = tmat_.rt21[i][j][k];
				wid->rt22[i][j][k] = tmat_.rt22[i][j][k];
				wid->it11[i][j][k] = tmat_.it11[i][j][k];
				wid->it12[i][j][k] = tmat_.it12[i][j][k];
				wid->it21[i][j][k] = tmat_.it21[i][j][k];
				wid->it22[i][j][k] = tmat_.it22[i][j][k];
			}
		}
	}
	mishchenko_io_.tmatrixnr = rand();
	wid->tmatrixnr = mishchenko_io_.tmatrixnr;
}

void mishchenko2000_amp_phase_matrices(t_mishchenko2000_widget *wid)
{
	int i,j,k;

	//update fortran io
	mishchenko_io_.axi 	= wid->axi;
	mishchenko_io_.rat 	= wid->rat;
	mishchenko_io_.lam 	= wid->lam;
	mishchenko_io_.mrr 	= wid->mrr;
	mishchenko_io_.mri 	= wid->mri;
	mishchenko_io_.eps 	= wid->eps;
	mishchenko_io_.np 	= wid->np;
	mishchenko_io_.ndgs 	= wid->ndgs;
	mishchenko_io_.ddelt = wid->ddelt;
	
	//update fortran io
	mishchenko_io_.alpha = wid->alpha;
	mishchenko_io_.beta 	= wid->beta;
	mishchenko_io_.thet0 = wid->thet0;
	mishchenko_io_.thet 	= wid->thet;
	mishchenko_io_.phi0 	= wid->phi0;
	mishchenko_io_.phi 	= wid->phi;

	//update fortran tmatrix
	if (mishchenko_io_.tmatrixnr != wid->tmatrixnr) {
		mishchenko_io_.nmax 	= wid->nmax;
		for (i=0; i < wid->nmax; i++) {
			for (j=0; j < wid->nmax; j++) {
				for (k=0; k < (wid->nmax + 1); k++) {				
					tmat_.rt11[i][j][k] = wid->rt11[i][j][k];
					tmat_.rt12[i][j][k] = wid->rt12[i][j][k];
					tmat_.rt21[i][j][k] = wid->rt21[i][j][k];
					tmat_.rt22[i][j][k] = wid->rt22[i][j][k];
					tmat_.it11[i][j][k] = wid->it11[i][j][k];
					tmat_.it12[i][j][k] = wid->it12[i][j][k];
					tmat_.it21[i][j][k] = wid->it21[i][j][k];
					tmat_.it22[i][j][k] = wid->it22[i][j][k];
				}
			}
		}
		mishchenko_io_.tmatrixnr = wid->tmatrixnr;
	}
      
	//execute fortran code to calculate matrices
	calc_amp_phase_matrices_();
    	
	//save results to widget
	wid->s11 = mishchenko_io_.s11r + 1.j * mishchenko_io_.s11i;
	wid->s12 = mishchenko_io_.s12r + 1.j * mishchenko_io_.s12i;
	wid->s21 = mishchenko_io_.s21r + 1.j * mishchenko_io_.s21i;
	wid->s22 = mishchenko_io_.s22r + 1.j * mishchenko_io_.s22i;
}
