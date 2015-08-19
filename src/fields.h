#ifndef _ZEPHYROS_FIELDS
#define _ZEPHYROS_FIELDS

/* field_x: array, field model, position x [m] */
/* field_y: array, field model, position y [m] */
/* field_z: array, field model, position z [m] */
/* field_t: array, field model, time     t [s] */
/* field_u: array: field model, speed in x direction [m/s] */
/* field_v: array: field model, speed in y direction [m/s] */
/* field_w: array: field model, speed in z direction [m/s] */

#include <complex.h>

typedef struct st_zephyros_field
{
    int 		n;					/* n = n_x * n_y * n_z * n_t */
	int			n_x;
	int			n_y;
	int			n_z;
	int			n_t;
    double		*vec_x;
    double		*vec_y;
    double		*vec_z;
    double		*vec_t;
	char		name[8192];	//for error reporting
	
	double(*x)(void *vdfield, int i);
	double(*y)(void *vdfield, int i);
	double(*z)(void *vdfield, int i);
	double(*t)(void *vdfield, int i);
} t_zephyros_field;

void fields_initialize(t_zephyros_field **fieldpt);

void fields_prepare(t_zephyros_field *field);

void fields_prepare2(
	t_zephyros_field **fieldpt,
	double xmin, double xmax, double xdif,
	double ymin, double ymax, double ydif,
	double zmin, double zmax, double zdif,
	double tmin, double tmax, double tdif
	);

void fields_copy(t_zephyros_field **pdst, t_zephyros_field *src);

void fields_free(t_zephyros_field **pfield);

double fields_give_xyzt(void *vdfield, int i, int i_xyzt);
double fields_give_x(void *vdfield, int i);
double fields_give_y(void *vdfield, int i);
double fields_give_z(void *vdfield, int i);
double fields_give_t(void *vdfield, int i);

#endif
