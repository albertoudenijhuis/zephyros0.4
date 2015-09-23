/*
Description: 
	bilinear interpolation of look-up tables

Revision History:
	2014

Functions:
	bilinear interpolation of look-up tables

Author:
	Albert Oude Nijhuis <albertoudenijhuis@gmail.com>

Institute:
	Delft University of Technology
	
Zephyros version:
	0.4

Project:
	EU FP7 program, the UFO project

Dissemination:
	Confidential, only for members of the UFO project. Potentially public in the future.

Acknowledgement and citation:
	Whenever this code is used for publication of scientific results,
	the code writer should be informed, acknowledged and referenced.

Note:
	If you have any suggestions for improvements or amendments, please inform the author of this code.

*/

/* System include files */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "specialfunctions.h"

#include "interpolation.h"

//uncomment next statement for debug mode
//#define _ZEPHYROS_INTERPOLATION_DEBUG

void interpolation_bilint3D2lut(int n1, double *vec1, int n2, double *vec2, int n3, double *vec3, double *variable, int special, t_zephyros_interpolation_bilint_lut **plut)
{
	t_zephyros_interpolation_bilint_lut *lut;
	int i;

	interpolation_initialize_lut(&lut);

	//LUT interpolation stuff
	lut->n_dim 	= 3;
	lut->n		= n1 * n2 * n3;
	lut->shape 	= malloc(lut->n_dim * sizeof(int));	
	lut->shape[0] = n1;
	lut->shape[1] = n2;
	lut->shape[2] = n3;
	
	lut->ax_values = malloc(lut->n_dim * sizeof(double*));
	lut->ax_values[0] = malloc(n1 * sizeof(double));
	lut->ax_values[1] = malloc(n2 * sizeof(double));
	lut->ax_values[2] = malloc(n3 * sizeof(double));
	
	for ( i = 0; i < n1; i++ ) {
		lut->ax_values[0][i] = vec1[i];
	}
	for ( i = 0; i < n2; i++ ) {
		lut->ax_values[1][i] = vec2[i];
	}
	for ( i = 0; i < n3; i++ ) {
		lut->ax_values[2][i] = vec3[i];
	}
	
	lut->special = special;
	
	if (lut->special == 2) {
		lut->mesh_y = malloc(lut->n * sizeof(double));
		lut->mesy_y_allocated = 1;
		for ( i = 0; i < lut->n ; i++ ) {
			lut->mesh_y[i] = fabs(variable[i]);
		}		
		lut->special = 0;
	}
	if (lut->special == 3) {
		lut->mesh_y = malloc(lut->n * sizeof(double));
		lut->mesy_y_allocated = 1;
		for ( i = 0; i < lut->n ; i++ ) {
			lut->mesh_y[i] = log(variable[i]);
		}		
		lut->special = 0;		
	} else {
		lut->mesh_y = variable;
	}

	*plut = lut;
}



void interpolation_linearinterpolation2lut(int n, double *x, double *variable, int special, t_zephyros_interpolation_bilint_lut **plut)
{
	t_zephyros_interpolation_bilint_lut *lut;
	int i;
	
	interpolation_initialize_lut(&lut);
	
	//LUT interpolation stuff
	lut->n_dim 	= 1;
	lut->n		= n;
	lut->shape 	= malloc(1 * sizeof(int));	
	
	lut->shape[0] = n;
	
	lut->ax_values = malloc(lut->n_dim * sizeof(double*));
	lut->ax_values[0] = malloc(lut->shape[0] * sizeof(double));
	
	for ( i = 0; i < n; i++ ) {
		lut->ax_values[0][i] = x[i];
	}
	
	lut->special = special;
	
	if (lut->special == 2) {
		lut->mesh_y = malloc(lut->n * sizeof(double));
		lut->mesy_y_allocated = 1;
		for ( i = 0; i < lut->n ; i++ ) {
			lut->mesh_y[i] = fabs(variable[i]);
		}		
		lut->special = 0;
	}
	if (lut->special == 3) {
		lut->mesh_y = malloc(lut->n * sizeof(double));
		lut->mesy_y_allocated = 1;
		for ( i = 0; i < lut->n ; i++ ) {
			lut->mesh_y[i] = log(variable[i]);
		}		
		lut->special = 0;		
	} else {
		lut->mesh_y = variable;
	}

	*plut = lut;
}



/*
special function for bilinear interpolation.
It calculates the grid dependence for yout.
calculates a term of the cost function for minimization.
calculates:
sum i (dPrK_i / dK1, dPrK_i / dK2, ..., dPrK_i / dKn) * alpha_i

can e.g. be used to calculate:
sum i (dPrK_i / dK1, dPrK_i / dK2, ..., dPrK_i / dKn) * 2 (v_r_i - PrK_i) / si_v_r_i**2
*/

void interpolation_initialize_lut(t_zephyros_interpolation_bilint_lut **plut)
{
	t_zephyros_interpolation_bilint_lut *lut = malloc(sizeof(t_zephyros_interpolation_bilint_lut));
	lut->shape = NULL;
	lut->ax_values = NULL;
	lut->mesh_y = NULL;
	lut->periodic_L = NULL;
	lut->n_dim = 0;
	lut->n = 0;
	lut->special = 0;
	lut->mesy_y_allocated = 0;
	lut->periodic = 0;

	*plut = lut;
}


void interpolation_free_lut(t_zephyros_interpolation_bilint_lut **plut)
{
	int i_axis;
	t_zephyros_interpolation_bilint_lut *lut = *plut;
	
	if (lut != NULL) {
		if (lut->ax_values != NULL) {
			for ( i_axis = 0; i_axis < lut->n_dim; i_axis++ ) {			
				free(lut->ax_values[i_axis]);
			}
			free(lut->ax_values);
		}
		if (lut->shape != NULL) free(lut->shape);
		if ((lut->mesy_y_allocated == 1) & (lut->mesh_y != NULL)) {free(lut->mesh_y);}		
		free(lut);
		lut = NULL;
	}
}

void interpolation_bilint(
	t_zephyros_interpolation_bilint_lut	*lut,
	double	*outx,				//output x vector		length: n_dim
	double	*outy,				//output y				length: 1
	int		calc_derivatives,	//0 = no, 1 = yes.
	double	*derivatives		//derivatives			length: n_dim
	)
{
	int 	i, j, k;
	int 	n_edges;
	int 	nr_edge;

	int		*myslice_i0;
	int		*myslice_i1;
	double	*myslice_di;
	double	*myslice_dval;
	
	double 	w_edge;
	double	*w_edge_deriv;
	
	int		*tmp_tup 		= malloc(lut->n_dim * sizeof(int));	
	int		*tmp_tup_edge	= malloc(lut->n_dim * sizeof(int));	
	int		*tmp_tup_edge2	= malloc(lut->n_dim * sizeof(int));	
	
	double 	outysin;
	double 	outycos;

	double 	*derivativespluscos;
	double 	*derivativesplussin;
	double 	*derivativesmincos;
	double 	*derivativesminsin;
	
	double 	tmpplus, tmpmin;

	if (lut == NULL) {
		printf("Error with interpolation"); fflush(stdout);
		exit(0);
	}

	#ifdef _ZEPHYROS_INTERPOLATION_DEBUG
		printf("interpolation_bilint\n"); fflush(stdout);
	#endif 
	
	
	//memory allocation
	myslice_i0 										= malloc(lut->n_dim * sizeof(int));
	myslice_i1 										= malloc(lut->n_dim * sizeof(int));
	myslice_di 										= malloc(lut->n_dim * sizeof(double));
	if (calc_derivatives == 1) {myslice_dval 		= malloc(lut->n_dim * sizeof(double));}
	if (calc_derivatives == 1) {w_edge_deriv 		= malloc(lut->n_dim * sizeof(double));}
	
	if (lut->special == 1) {
		if (calc_derivatives == 1) {
			derivativespluscos	= malloc(lut->n_dim * sizeof(double));	
			derivativesplussin	= malloc(lut->n_dim * sizeof(double));	
			derivativesmincos	= malloc(lut->n_dim * sizeof(double));	
			derivativesminsin	= malloc(lut->n_dim * sizeof(double));	
		}
	}
	
	// clean up
	*outy = 0.;
	if (lut->special == 1) {
		outysin = 0.;
		outycos = 0.;
		if (calc_derivatives == 1) {
			for ( i = 0; i < lut->n_dim; i++ ) {
				derivativespluscos[i] 	= 0.;
				derivativesplussin[i] 	= 0.;
				derivativesmincos[i] 	= 0.;
				derivativesminsin[i] 	= 0.;
			}
		}
	}
	if (calc_derivatives == 1) {
		for ( i = 0; i < lut->n_dim; i++ ) {
			derivatives[i] = 0.;
		}
	}
	
	//loop over parameters and calculate indices
	j = 0;
	for ( i = 0; i < lut->n_dim; i++ ) {
		interpolation_calcindices(
						lut,
						i,
						outx[i],
						myslice_i0 + i,
						myslice_i1 + i,
						myslice_di + i
						);

		//debugging
		//~ if ((myslice_di[i] > 1.) | (myslice_di[i] < 0.)) {
			//~ printf("dim %i, i0 = %i, i1 = %i, di = %.2e\n", i, myslice_i0[i], myslice_i1[i], myslice_di[i]);
		//~ }

		if (calc_derivatives == 1) {
			//if derivatives are requested, also calculate this one.
			myslice_dval[i] = lut->ax_values[i][myslice_i1[i]] - lut->ax_values[i][myslice_i0[i]];				
		}
	}

	//make a template tuple, represting the edges of the interpolation lut
	n_edges = pow(2,lut->n_dim);
	for ( i = 0; i < lut->n_dim; i++ ) {
		tmp_tup[i] = 2;
	}
	
	//loop over edges.
	for ( i = 0; i < n_edges; i++ ) {
		interpolation_nr2tup(&lut->n_dim, tmp_tup, &i, tmp_tup_edge);
		
		//for *yout, calculate weight for this edge
		w_edge = 1.;
		for ( k = 0; k < lut->n_dim; k++ ) {
			if (tmp_tup_edge[k] == 0) {
				w_edge = w_edge * (1. - myslice_di[k]);
			} else {
				w_edge = w_edge * myslice_di[k];
			}
		}

		if (calc_derivatives == 1) {
			//for derivatives, calculate weight for this edge
			for ( j = 0; j < lut->n_dim; j++ ) {		
				w_edge_deriv[j] = 1.;
				for ( k = 0; k < lut->n_dim; k++ ) {
					if (j==k) {
						if (myslice_dval[k] == 0.) {
							w_edge_deriv[j] = 0.;
						} else {
							if (tmp_tup_edge[k] == 0) {
								w_edge_deriv[j] =  -1. * w_edge_deriv[j] / myslice_dval[k];					
							} else {
								w_edge_deriv[j] = 1. * w_edge_deriv[j] / myslice_dval[k];					
							}
						}
					} else {
						if (tmp_tup_edge[k] == 0) {
							w_edge_deriv[j] = w_edge_deriv[j] * (1. - myslice_di[k]);
						} else {
							w_edge_deriv[j] = w_edge_deriv[j] * myslice_di[k];
						}
					}
				}				
			}
		}
		
		// Calculate position in mesh
		for ( k = 0; k < lut->n_dim; k++ ) {
			tmp_tup_edge2[k] = myslice_i0[k] + tmp_tup_edge[k];
		}
		interpolation_tup2nr(&lut->n_dim, lut->shape, tmp_tup_edge2, &nr_edge);

		if (nr_edge < lut->n) {// solution to one-valued ax bug
			// update *outy	
			if (lut->special == 2) {
				//abs interpolation
				*outy 	+= w_edge * fabs(lut->mesh_y[nr_edge]);
			} else if (lut->special == 3) {
				//log interpolation
				*outy 	+= w_edge * log(lut->mesh_y[nr_edge]);
			} else if (lut->special == 1) {
				//angular interpolation
				outysin += w_edge * sin(lut->mesh_y[nr_edge]);
				outycos += w_edge * cos(lut->mesh_y[nr_edge]);
			} else {
				//normal interpolation
				*outy 	+=  w_edge * lut->mesh_y[nr_edge];			
			}
			if (calc_derivatives == 1) {
				// update derivatives, 
				for ( j = 0; j < lut->n_dim; j++ ) {
					if (lut->special == 2) {
						//abs interpolation
						derivatives[j] 		+= w_edge_deriv[j] * fabs(lut->mesh_y[nr_edge]) ;
					} else if (lut->special == 3) {
						//log interpolation
						derivatives[j] 		+= w_edge_deriv[j] * log(lut->mesh_y[nr_edge]) ;
					} else if (lut->special == 1) {
						//angular interpolation
						if (w_edge_deriv[j] > 0.) {
							derivativesplussin[j] 	+= w_edge_deriv[j] * sin(lut->mesh_y[nr_edge]);
							derivativespluscos[j] 	+= w_edge_deriv[j] * cos(lut->mesh_y[nr_edge]);
						} else {
							derivativesminsin[j] 	+= -w_edge_deriv[j] * sin(lut->mesh_y[nr_edge]);
							derivativesmincos[j] 	+= -w_edge_deriv[j] * cos(lut->mesh_y[nr_edge]);
						}
					} else {
						//normal
						derivatives[j] 		+= w_edge_deriv[j] * lut->mesh_y[nr_edge] ;
					}
				}
			}
		}
	}
	
	if (lut->special == 1) {
		*outy = atan2(outysin, outycos);
		for ( j = 0; j < lut->n_dim; j++ ) {
			tmpplus = atan2(derivativesplussin[j], derivativespluscos[j]);
			tmpmin	= atan2(derivativesminsin[j], derivativesmincos[j]);
			derivatives[j] = (fmod(M_PI + tmpplus - tmpmin, 2. * M_PI) - M_PI) ;
		}
	}
	
	free(myslice_i0);
	free(myslice_i1);
	free(myslice_di);
	if (calc_derivatives == 1) {free(myslice_dval);}
	if (calc_derivatives == 1) {free(w_edge_deriv);}
		
	free(tmp_tup);
	free(tmp_tup_edge);
	free(tmp_tup_edge2);

	if (lut->special == 1) {
		if (calc_derivatives == 1) {
			free(derivativespluscos);
			free(derivativesplussin);
			free(derivativesmincos);
			free(derivativesminsin);
		}
	}

	#ifdef _ZEPHYROS_INTERPOLATION_DEBUG
		printf("interpolation_bilint: end\n"); fflush(stdout);
	#endif 
}

void interpolation_bilint_griddep(
	t_zephyros_interpolation_bilint_lut	*lut,
	double		*outx,					//output x vector		length: n_dim
	double		*griddep				//griddep				length: n_val
	)
{
	int i, j, k;
	int n_edges;
	int nr_edge;
	
	int 	*myslice_i0 = malloc(lut->n_dim * sizeof(int)); 		
	int 	*myslice_i1 = malloc(lut->n_dim * sizeof(int)); 		
	double 	*myslice_di = malloc(lut->n_dim * sizeof(double)); 		
	
	double 	w_edge;
	
	int 	*tmp_tup 		= malloc(lut->n_dim * sizeof(int)); 		
	int 	*tmp_tup_edge 	= malloc(lut->n_dim * sizeof(int));	
	int 	*tmp_tup_edge2 	= malloc(lut->n_dim * sizeof(int));	
			
	//check memory allocation
	if (((myslice_i0 == NULL) | (myslice_i1 == NULL) | (myslice_di == NULL) | (tmp_tup == NULL)
		| (tmp_tup_edge == NULL) | (tmp_tup_edge2 == NULL))) {
		printf("Memory allocation failed...\n");
		exit(0);
	}
		
	/* clean up */
	for ( i = 0; i < lut->n; i++ ) {
		griddep[i] = 0.;
	}
	
	//loop over parameters and calculate indices
	j = 0;
	for ( i = 0; i < lut->n_dim; i++ ) {
		interpolation_calcindices(
						lut,
						i,
						outx[i],
						myslice_i0 + i,
						myslice_i1 + i,
						myslice_di + i
					);
	}
	
	//make a template tuple, represting the edges of the interpolation lut
	n_edges = pow(2,lut->n_dim);
	for ( i = 0; i < lut->n_dim; i++ ) {
		tmp_tup[i] = 2;
	}
	
	//loop over edges.
	for ( i = 0; i < n_edges; i++ ) {
		interpolation_nr2tup(&lut->n_dim, tmp_tup, &i, tmp_tup_edge);
		
		//for *yout, calculate weight for this edge
		w_edge = 1.;
		for ( k = 0; k < lut->n_dim; k++ ) {
			if (tmp_tup_edge[k] == 0) {
				w_edge = w_edge * (1. - myslice_di[k]);
			} else {
				w_edge = w_edge * myslice_di[k];
			}
		}

		/* Calculate position in mesh */
		for ( k = 0; k < lut->n_dim; k++ ) {
			tmp_tup_edge2[k] = myslice_i0[k] + tmp_tup_edge[k];
		}
		interpolation_tup2nr(&lut->n_dim, lut->shape, tmp_tup_edge2, &nr_edge);

		if (nr_edge < lut->n) {// solution to one-valued ax bug
			//update *outy	
			griddep[nr_edge] = griddep[nr_edge] + w_edge;
		}
	}

	free(myslice_i0);
	free(myslice_i1);
	free(myslice_di);
		
	free(tmp_tup);
	free(tmp_tup_edge);
	free(tmp_tup_edge2);
}

void interpolation_calcindices(
	t_zephyros_interpolation_bilint_lut	*lut,
	int i_axis,	
	double		val,
	int			*i0,
	int			*i1,
	double		*di
	)
{
	int i, foundvalue;
	double val2;
	
	#ifdef _ZEPHYROS_INTERPOLATION_DEBUG
		printf("interpolation_calcindices\n"); fflush(stdout);
	#endif 
	
	if (lut->shape[i_axis] == 1) {
		/* only one value on the axis */
		*i0 = 0;
		*i1 = 0;
		*di = 0.;
	} else {
		/* multiple values on the axis */
		if (lut->periodic) {
			//assumption, ax_values are modulo periodic_L

			//scale of the domain
			val2 = interpolation_modulo(val, lut->periodic_L[i_axis]);
			
			if (val2 < 0.) {
				printf("val = %.2e\n; val2 = %.2e; L = %.2e\n", val, val2, lut->periodic_L[i_axis]);
				printf("i0 = %i; i1 = %i\n", *i0, *i1);
				printf("di = %.2e\n", *di);
			}
			
			foundvalue = 0;
			/* walk through axis */
			for ( i = 0; i < (lut->shape[i_axis] - 1); i++ ) {
				if ((lut->ax_values[i_axis][i] <= val2) & (val2 <= lut->ax_values[i_axis][i+1])) {
					*i0 = i;
					*i1 = i + 1;
					foundvalue = 1;
					break;
				}
			}

			if (foundvalue) {
				*di = (val2 - lut->ax_values[i_axis][*i0]) / (lut->ax_values[i_axis][*i1] - lut->ax_values[i_axis][*i0]);
			} else {
				/* value in the last interval */
				*i0 = lut->shape[i_axis] - 1;
				*i1 = 0;
				*di = (val2 - lut->ax_values[i_axis][*i0]) / (lut->periodic_L[i_axis] - lut->ax_values[i_axis][*i0]);
			}

		} else {
			foundvalue = 0;
			/* walk through axis */
			for ( i = 0; i < (lut->shape[i_axis] - 1); i++ ) {
				if ((lut->ax_values[i_axis][i] <= val) & (val <= lut->ax_values[i_axis][i+1])) {
					*i0 = i;
					*i1 = i + 1;
					foundvalue = 1;
					break;
				}
			}
			
			if (foundvalue == 0) {
				/* out of range, extrapolate */
				if (fabs(lut->ax_values[i_axis][0] - val) < fabs(lut->ax_values[i_axis][lut->shape[i_axis]-1] - val)) {
					*i0 = 0;
					*i1 = 1;
				} else {
					*i0 = lut->shape[i_axis]-2;
					*i1 = lut->shape[i_axis]-1;
				}
			}
			
			*di = (val - lut->ax_values[i_axis][*i0]) / (lut->ax_values[i_axis][*i1] - lut->ax_values[i_axis][*i0]);
		}
	}

	#ifdef _ZEPHYROS_CONFIG_DEBUG
		printf("interpolation_calcindices:end\n"); fflush(stdout);
	#endif 
}

// Calculate nr in 1d-array, corresponding to index tuple tup for an n-dimensional array with the shape arr_shape
// in, 	arr_shape: 		shape of the n-dimensional array
// in, 	tup:			indices of the point under consideration.
// out,	nr:				number in 1d-array.
void interpolation_tup2nr(
	int	*n,
	int	*arr_shape,
	int	*tup,
	int	*nr)
{
	int i, multiplier;
	
	//initialize
	*nr = 0;
	multiplier = 1;
	for ( i = 0; i < *n; i++ ) {
		*nr = *nr + tup[i] * multiplier;
		multiplier = multiplier * arr_shape[i];
	}
}

/*
Calculate index tuple tup for an n-dimensional array with the shape arr_shape, corresponding to nr in 1d-array.
in, 	arr_shape: 		shape of the n-dimensional array
in, 	nr:				number in 1d-array.
out, 	tup:			indices of the tuple
*/
void interpolation_nr2tup(
	int	*n,
	int	*arr_shape,
	int	*nr,	
	int	*tup)
{
	int i, j, tmp, tmp2;
	
	tmp = *nr;
	interpolation_intprod(n, arr_shape, &j);

	for ( i = *n - 1; i >= 0; i-- ) {
		j = j / arr_shape[i];
		tmp2 = tmp % j;
		tup[i] = (tmp - tmp2) / j;
		tmp = tmp2;
	}
}

/*
void interpolation_piecewise_int(
	int		*n,
	double	*xin,
	double	*yin,
	double	*xout,
	double	*yout)
{
	int	i0, i1;
	double di;

	interpolation_calcindices(n, xin, xout, &i0, &i1, &di, 0);
	*yout = (1. - di) * yin[i0] + di * yin[i1];
}
*/

/* product of an integer array */
void interpolation_intprod(
	int		*n,
	int		*arr,
	int		*result)
{
	int i;
	
	*result = 1.;
	for ( i = 0; i < *n; i++ ) {
		*result *= arr[i];
	}
}

//calculate nearest integer
int interpolation_nint(
	double *xin)
	{
	int ans = (int) floor(*xin);

	if (fabs(ans - *xin) < 0.5) {
		return ans;
	} else {
		return ans + 1;
	}
}

void interpolation_nq_sh_interp(
		int		*N,
		double 	*u,
		double 	*xin,
		double 	*q,
		double 	*ans)
{
		int		i;
		double 	x = fmod(*xin, (double) *N);
					
		*ans = 0.;
		for ( i = 0; i < *N; i++ ) {
			*ans += 
				u[i] * interpolation_nq_sh_f_xi(N, &x, &i, q);
		}
}

void interpolation_nq_sh_interp3D(
		int		*Nx,
		int		*Ny,
		int		*Nz,
		double 	*u,
		double 	*xin,
		double 	*yin,
		double 	*zin,
		double 	*q,
		double 	*ans)
{
	int i,j,k;
	double 	x = fmod(*xin, (double) *Nx);
	double 	y = fmod(*yin, (double) *Ny);
	double 	z = fmod(*zin, (double) *Nz);
	int n_dim = 3;
	int shape[3] = {*Nx, *Ny, *Nz};
	int tmp_tup[3];
	int nr;
		
	*ans = 0.;
	for ( i = 0; i < *Nx; i++ ) {
		for ( j = 0; j < *Ny; j++ ) {
			for ( k = 0; k < *Nz; k++ ) {
				tmp_tup[0] = i;
				tmp_tup[1] = j;
				tmp_tup[2] = k;
				interpolation_tup2nr(&n_dim, shape, tmp_tup, &nr);
				*ans += 
					u[nr] 
					* interpolation_nq_sh_f_xi(Nx, &x, &i, q)
					* interpolation_nq_sh_f_xi(Ny, &y, &j, q)
					* interpolation_nq_sh_f_xi(Nz, &z, &k, q);
			}
		}
	}
}

double interpolation_nq_sh_f_xi(
			int		*N,
			double 	*xin,
			int		*i,
			double 	*q)
{
	double x  = fmod(*xin, (double) *N);
	double a, b, ka, kb;
	
	//check input
	/*
	if (*N != *N) {
		printf("N = %.2e", *N);
		exit(0);
	}
	if (*xin != *xin) {
		printf("xin = %.2e", *xin);
		exit(0);
	}
	if (*i != *i) {
		printf("i = %.2e", *i);
		exit(0);
	}
	*/
	
	if (fabs(interpolation_nint(&x) - x) < 1.e-50) {
		//special limits of this function.
		if (interpolation_nint(&x) == *i) {
			return 1.;
		} else {
			return 0.;
		}
	} else {
		a 	= (x - *i) / (2. * *N);
		ka 	= -1. * (floor(a));
		
		b 	= ((x - *i) / (2. * *N)) - 0.5;
		kb 	= -1. * (floor(b));

		//~ printf("\n\n");
		//~ printf("*q: %.2e\n", *q);
		//~ printf("\n");
		//~ printf("(-1.0 * a) + *q + 1.0:    %.2e\n", (-1.0 * a) + *q + 1.0);
		//~ printf("(-1.0 * a) + 1.0 - ka:    %.2e\n", (-1.0 * a) + 1.0 - ka);
		//~ printf("a + *q + 1.0:           %.2e\n", a + *q + 1.0);
		//~ printf("a + ka:                 %.2e\n", a + ka);
		//~ printf("\n");
		//~ printf("(-1.0 * b) + *q + 1.0:    %.2e\n", (-1.0 * b) + *q + 1.0);
		//~ printf("(-1.0 * b) + 1.0 - kb:    %.2e\n", (-1.0 * b) + 1.0 - kb);
		//~ printf("b + *q + 1.0:           %.2e\n", b + *q + 1.0);
		//~ printf("b + kb:                 %.2e\n", b + kb);
		
		return (1.0 / (2.0 * M_PI * *N )) *
				( sin(M_PI * (x - *i )) *
					( 		-1.0 * specialfunctions_digamma((-1.0 * a) + *q + 1.0)
								+ specialfunctions_digamma((-1.0 * a) + 1.0 - ka)
								+ specialfunctions_digamma(a + *q + 1.0)
								- specialfunctions_digamma(a + ka))
				 + sin(M_PI * (x - *i - *N)) *
						( 	-1.0 * specialfunctions_digamma((-1.0 * b) + *q + 1.0)
								+ specialfunctions_digamma((-1.0 * b) + 1.0 - kb)
								+ specialfunctions_digamma(b + *q + 1.0)
								- specialfunctions_digamma(b + kb) )
				);
	}
}

double interpolation_modulo(double x, double y)
{
	int i;
	
	if (x < 0) {
		i = abs(x / y);
		return fmod(x + ((1 + i) * y), y);
	} else {
		return fmod(x,y);
	}
}
