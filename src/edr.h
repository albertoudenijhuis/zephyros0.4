void retr_edr_variance(
	int 	*n,				//length of variables v,t,l,v2
	double 	*v,
	double 	*t,
	double 	*l,
	double 	*v2,
	int 	*retr_domain,	//0 = space, 1 = time and U0 from v2
	double 	*retr_C,		//Kolmogorov constant that is used.
	int		*retr_n,		//Number of points used for average and variance calculation
	double	*retr_av,		//output: retrieved average velocity
	double 	*retr_var,		//output: retrieved variance of velocity
	double	*retr_edr,		//output: retrieved edr
	double	*retr_edr13_err	//output: error in retrieved edr^(1/3) (std)
	);
	
void retr_edr_powerspectrum(
	int 	*n,				//length of variables v,t,l,v2
	double 	*v,
	double 	*t,
	double 	*l,
	double 	*v2,
	int 	*retr_domain,		//0 = space, 1 = time and U0 from v2
	double 	*retr_C,			//Kolmogorov constant that is used.
	int		*retr_n,			//Number of points used for average and variance calculation
	int		*periodic,			//1 = periodic, 0 = non-periodic
	int		*retr_nintervals,	//number of intervals used
	double	*retr_edr,			//output: retrieved edr
	double	*retr_edr13_err		//error in retrieved edr^(1/3) (std)
	);
	
void retr_edr_2nd_order_structure_function(
	int 	*n,				//length of variables v,t,l,v2
	double 	*v,
	double 	*t,
	double 	*l,
	double 	*v2,
	int 	*retr_domain,		//0 = space, 1 = time and U0 from v2
	double 	*retr_C,			//Kolmogorov constant that is used.
	int		*retr_n,			//Number of points used for average and variance calculation
	int		*periodic,			//1 = periodic, 0 = non-periodic
	double	*retr_edr,			//output: retrieved edr
	double	*retr_edr13_err		//error in retrieved edr^(1/3) (std)
	);
	
void retr_edr_3rd_order_structure_function(
	int 	*n,				//length of variables v,t,l,v2
	double 	*v,
	double 	*t,
	double 	*l,
	double 	*v2,
	int 	*retr_domain,		//0 = space, 1 = time and U0 from v2
	double 	*retr_C,			//Kolmogorov constant that is used.
	int		*retr_n,			//Number of points used for average and variance calculation
	int		*periodic,			//1 = periodic, 0 = non-periodic
	double	*retr_edr,			//output: retrieved edr
	double	*retr_edr13_err		//error in retrieved edr^(1/3) (std)
	);
	
void kolmogorov_skewness(
	int 	*n,
	double 	*x,
	double 	*y,
	int		*periodic
	);
	
void retr_skewness_structure_functions(
	int		*n,
	double 	*v,
	int		*retr_n,
	int		*periodic,
	double	*retr_skewness,
	double	*retr_constant
	);
		
void retr_factor_wind_components(
	int		*n,
	double	*u,
	double	*v,
	double 	*w,
	int		*retr_n,
	int		*periodic,
	double	*retr_frachor,	//retrieved factor between components
	double 	*retr_fracver	//retrieved factor between components
	);
	
void kolmogorov_constants(
	char   *choice,
	double *power_c,
	double *struc2_c,
	double *struc3_c
	);
	
void radial_kolmogorov_constants(
	double *azimuth_rad,
	double *azimuth0_rad,
	double *elevation_rad,
	double *power_c,
	double *struc2_c,
	double *struc3_c
	);
	
