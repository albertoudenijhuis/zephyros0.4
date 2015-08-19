#include "particles_mischenko2000.h"

#include <stdio.h>



void mischenko2000_calc_tmatrix(t_mischenko2000_widget *wid)
{
	int i,j,k;

	//update fortran io
	mischenko_io_.axi 	= wid->axi;
	mischenko_io_.rat 	= wid->rat;
	mischenko_io_.lam 	= wid->lam;
	mischenko_io_.mrr 	= wid->mrr;
	mischenko_io_.mri 	= wid->mri;
	mischenko_io_.eps 	= wid->eps;
	mischenko_io_.np 	= wid->np;
	mischenko_io_.ndgs 	= wid->ndgs;
	mischenko_io_.ddelt = wid->ddelt;

	//execute fortran code to calculate t-matrix
	calc_tmatrix_();
	
	//save results to widget
	wid->nmax = mischenko_io_.nmax;		
	
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
	mischenko_io_.tmatrixnr = rand();
	wid->tmatrixnr = mischenko_io_.tmatrixnr;
}

void mischenko2000_amp_phase_matrices(t_mischenko2000_widget *wid)
{
	int i,j,k;

	//update fortran io
	mischenko_io_.axi 	= wid->axi;
	mischenko_io_.rat 	= wid->rat;
	mischenko_io_.lam 	= wid->lam;
	mischenko_io_.mrr 	= wid->mrr;
	mischenko_io_.mri 	= wid->mri;
	mischenko_io_.eps 	= wid->eps;
	mischenko_io_.np 	= wid->np;
	mischenko_io_.ndgs 	= wid->ndgs;
	mischenko_io_.ddelt = wid->ddelt;
	
	//update fortran io
	mischenko_io_.alpha = wid->alpha;
	mischenko_io_.beta 	= wid->beta;
	mischenko_io_.thet0 = wid->thet0;
	mischenko_io_.thet 	= wid->thet;
	mischenko_io_.phi0 	= wid->phi0;
	mischenko_io_.phi 	= wid->phi;

	//update fortran tmatrix
	if (mischenko_io_.tmatrixnr != wid->tmatrixnr) {
		mischenko_io_.nmax 	= wid->nmax;
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
		mischenko_io_.tmatrixnr = wid->tmatrixnr;
	}
      
	//execute fortran code to calculate matrices
	calc_amp_phase_matrices_();
    	
	//save results to widget
	wid->s11 = mischenko_io_.s11r + 1.j * mischenko_io_.s11i;
	wid->s12 = mischenko_io_.s12r + 1.j * mischenko_io_.s12i;
	wid->s21 = mischenko_io_.s21r + 1.j * mischenko_io_.s21i;
	wid->s22 = mischenko_io_.s22r + 1.j * mischenko_io_.s22i;
}
