#ifndef _ZEPHYROS_INTERPOLATION
#define _ZEPHYROS_INTERPOLATION

typedef struct st_zephyros_interpolation_bilint_lut
{
	int 		n_dim;
	int 		n;
	int 		*shape;
	double 		*ax_values;
	int 		special;			//special: 0 = nothing, 1 = apply angular interpolation [rad], 2 = interpolate abs value, 3 = interpolate log value
	double		*mesh_y;			//actual lut; y-values on the mesh; length: n.
	int			mesy_y_allocated;	//0 = not allocated, 1 = allocated, will be freed.
} t_zephyros_interpolation_bilint_lut;

void interpolation_bilint3D2lut(int n1, double *vec1, int n2, double *vec2, int n3, double *vec3, double *variable, int special, t_zephyros_interpolation_bilint_lut **plut);

void interpolation_linearinterpolation2lut(int n, double *x, double *variable, int special, t_zephyros_interpolation_bilint_lut **plut);

void interpolation_free_lut(t_zephyros_interpolation_bilint_lut **plut);

void interpolation_bilint(
	t_zephyros_interpolation_bilint_lut	*lut,
	double			*outx,				//output x vector		length: n_dim
	double			*outy,				//output y				length: 1
	int				calc_derivatives,	//0 = no, 1 = yes.
	double			*derivatives		//derivatives			length: n_dim
	);
	
void interpolation_bilint_griddep(
	t_zephyros_interpolation_bilint_lut	*lut,
	double		*outx,					//output x vector		length: n_dim
	double		*griddep				//griddep				length: n_val
	);

void interpolation_calcindices(
	int			*n,
	double		*ax_values,
	double		*val,
	int			*i0,
	int			*i1,
	double		*di
	);

void interpolation_tup2nr(
	int	*n,
	int	*arr_shape,
	int	*tup,
	int	*nr);

void interpolation_nr2tup(
	int	*n,
	int	*arr_shape,
	int	*nr,
	int	*tup);

void interpolation_piecewise_int(int *n, double	*xin, double *yin, double *xout, double	*yout);

void interpolation_intprod(int *n, int *arr, int *result);
	
int interpolation_nint(double *xin);

void interpolation_nq_sh_interp(
		int		*N,
		double 	*u,
		double 	*xin,
		double 	*q,
		double 	*ans);
		
void interpolation_nq_sh_interp3D(
		int		*Nx,
		int		*Ny,
		int		*Nz,
		double 	*u,
		double 	*xin,
		double 	*yin,
		double 	*zin,
		double 	*q,
		double 	*ans);
		
double interpolation_nq_sh_f_xi(
			int		*N,
			double 	*xin,
			int		*i,
			double 	*q);
		
#endif
