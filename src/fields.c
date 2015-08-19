/*
Description: 
	Field functions

Revision History:
	2014

Functions:
	Field functions

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
	Whenever this code used for publication of scientific results,
	the code writer should be informed, acknowledged and referenced.

Note:
	If you have any suggestions for improvements or amendments, please inform the author of this code.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "fields.h"
#include "interpolation.h"

void fields_initialize(t_zephyros_field **fieldpt)
{	
	t_zephyros_field *field = malloc(sizeof(t_zephyros_field));

	//set functions
	field->x = fields_give_x;
	field->y = fields_give_y;
	field->z = fields_give_z;
	field->t = fields_give_t;
	
	field->vec_x = NULL;
	field->vec_y = NULL;
	field->vec_z = NULL;
	field->vec_t = NULL;
	strcpy(field->name, "Untitled");
	
	*fieldpt = field;
}

void fields_prepare(t_zephyros_field *field)
{
	//it is expected that the following variables are set:
	//field->n_x, field->vec_x and for other y, z and t.
	if (field->vec_x == NULL) 	{printf("vec_x variable was not set for field with name ``%s''\n", field->name); exit(0);}
	if (field->vec_y == NULL) 	{printf("vec_y variable was not set for field with name ``%s''\n", field->name); exit(0);}
	if (field->vec_z == NULL) 	{printf("vec_z variable was not set for field with name ``%s''\n", field->name); exit(0);}
	if (field->vec_t == NULL) 	{printf("vec_t variable was not set for field with name ``%s''\n", field->name); exit(0);}
	if (field->n_x == 0) 		{printf("Size for vec_x was 0 for field with name ``%s''\n", field->name); exit(0);}
	if (field->n_y == 0) 		{printf("Size for vec_y was 0 for field with name ``%s''\n", field->name); exit(0);}
	if (field->n_z == 0) 		{printf("Size for vec_z was 0 for field with name ``%s''\n", field->name); exit(0);}
	if (field->n_t == 0) 		{printf("Size for vec_t was 0 for field with name ``%s''\n", field->name); exit(0);}
	
	field->n = field->n_x * field->n_y * field->n_z * field->n_t;
}

void fields_prepare2(
	t_zephyros_field **fieldpt,
	double xmin, double xmax, double xdif,
	double ymin, double ymax, double ydif,
	double zmin, double zmax, double zdif,
	double tmin, double tmax, double tdif
	)
{
	t_zephyros_field *field = malloc(sizeof(t_zephyros_field));
	double dummy_delta = 1.e-100;

	int i, j;
	int in;
	int ix, iy, iz, it;
	int tmpsum;
	
	*fieldpt = field;
		
	field->n_x = (int) (1. + (dummy_delta + xmax - xmin) / xdif);
	func_dbl_arr_calloc(field->n_x, &field->vec_x);
	for ( i = 0; i < field->n_x; i++ ) {
		field->vec_x[i] = xmin + (i * xdif);
	}

	field->n_y = (int) (1. + (dummy_delta + ymax - ymin) / ydif);
	func_dbl_arr_calloc(field->n_y, &field->vec_y);
	
	for ( i = 0; i < field->n_y; i++ ) {
		field->vec_y[i] = ymin + (i * ydif);
	}

	field->n_z = (int) (1. + (dummy_delta + zmax - zmin) / zdif);
	func_dbl_arr_calloc(field->n_z, &field->vec_z);

	for ( i = 0; i < field->n_z; i++ ) {
		field->vec_z[i] = zmin + (i * zdif);
	}

	field->n_t = (int) (1. + (dummy_delta + tmax - tmin) / tdif);
	func_dbl_arr_calloc(field->n_t, &field->vec_t);
	
	for ( i = 0; i < field->n_t; i++ ) {
		field->vec_t[i] = tmin + (i * tdif);
	}
	
	field->n = field->n_x * field->n_y * field->n_z * field->n_t;
}

void fields_copy(t_zephyros_field **pdst, t_zephyros_field *src)
{
	int i;
	t_zephyros_field *dst;
	
	if (src == NULL) {
		dst = NULL;
	} else {
		dst = malloc(sizeof(t_zephyros_field));

		memcpy(dst, src, sizeof(t_zephyros_field));
		
		dst->vec_x = malloc(src->n_x * sizeof(double));
		memcpy(dst->vec_x, src->vec_x, src->n_x * sizeof(double));
		dst->vec_y = malloc(src->n_y * sizeof(double));
		memcpy(dst->vec_y, src->vec_y, src->n_y * sizeof(double));
		dst->vec_z = malloc(src->n_z * sizeof(double));
		memcpy(dst->vec_z, src->vec_z, src->n_z * sizeof(double));
		dst->vec_t = malloc(src->n_t * sizeof(double));
		memcpy(dst->vec_t, src->vec_t, src->n_t * sizeof(double));	
	}
	*pdst = dst;
}


void fields_free(t_zephyros_field **pfield)
{
	t_zephyros_field *field = *pfield;
	
	if (field != NULL) {
		if (field->vec_x != NULL) {free(field->vec_x); field->vec_x = NULL;}
		if (field->vec_y != NULL) {free(field->vec_y); field->vec_y = NULL;}
		if (field->vec_z != NULL) {free(field->vec_z); field->vec_z = NULL;}
		if (field->vec_t != NULL) {free(field->vec_t); field->vec_t = NULL;}
			
		free(field);
		field = NULL;
	}
}


double fields_give_xyzt(void *vdfield, int i, int i_xyzt)
{
	t_zephyros_field *field =  (t_zephyros_field*) vdfield ;
	int ndim = 4;
	int *arr_shape = malloc(ndim * sizeof(double));
	int *tup = malloc(ndim * sizeof(double));
	
	arr_shape[0] = field->n_x;
	arr_shape[1] = field->n_y;
	arr_shape[2] = field->n_z;
	arr_shape[3] = field->n_t;
	
	interpolation_nr2tup(&ndim, arr_shape, &i, tup);

	if (i_xyzt == 0) return field->vec_x[tup[0]];
	if (i_xyzt == 1) return field->vec_y[tup[1]];
	if (i_xyzt == 2) return field->vec_z[tup[2]];
	if (i_xyzt == 3) return field->vec_t[tup[3]];
}


double fields_give_x(void *vdfield, int i)
{
	fields_give_xyzt(vdfield, i, 0);
}
double fields_give_y(void *vdfield, int i)
{
	fields_give_xyzt(vdfield, i, 1);
}
double fields_give_z(void *vdfield, int i)
{
	fields_give_xyzt(vdfield, i, 2);
}
double fields_give_t(void *vdfield, int i)
{
	fields_give_xyzt(vdfield, i, 3);
}

