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

//uncomment next statement for debug mode
#define _ZEPHYROS_FIELDS_DEBUG

void fields_initialize(t_zephyros_field **pfield)
{	
	t_zephyros_field *field = calloc(1, sizeof(t_zephyros_field));

	//set functions
	field->x = fields_give_x;
	field->y = fields_give_y;
	field->z = fields_give_z;
	field->t = fields_give_t;
	
	//default settings
	field->n_x 	= 1;
	field->vec_x = calloc(1, sizeof(double));
	field->n_y 	= 1;
	field->vec_y = calloc(1, sizeof(double));
	field->n_z 	= 1;
	field->vec_z = calloc(1, sizeof(double));
	field->n_t 	= 1;
	field->vec_t = calloc(1, sizeof(double));
	
	strcpy(field->name, "Untitled");
	field->prepared 	= 0;
	field->initialized 	= 1;
	
	*pfield = field;
}

void fields_assert_initialized(t_zephyros_field *field)
{
	if (field == NULL) {
		printf("Field was NULL pointer. Exiting.\n");
		fflush(stdout); exit(0);
	}
	if (field->initialized != 1) {
		printf("Field was not initialized. Exiting.\n");
		fflush(stdout); exit(0);
	}
}

void fields_assert_prepared(t_zephyros_field *field)
{
	fields_assert_initialized(field);
	if (field->prepared != 1) {
		printf("Field `%s' was not prepared. Exiting.\n", field->name);
		fflush(stdout); exit(0);
	}
}

void fields_free(t_zephyros_field **pfield)
{
	t_zephyros_field *field = *pfield;
	
	if (field != NULL) {
		fields_assert_initialized(field);
		if (field->vec_x != NULL) {free(field->vec_x);}
		if (field->vec_y != NULL) {free(field->vec_y);}
		if (field->vec_z != NULL) {free(field->vec_z);}
		if (field->vec_t != NULL) {free(field->vec_t);}
			
		free(field);
		*pfield = NULL;
	}
}


void fields_prepare(t_zephyros_field *field)
{
	fields_assert_initialized(field);
	
	//it is expected that the following variables are set:
	//field->n_x, field->vec_x and for other y, z and t.
	if (field->vec_x == NULL) 	{printf("vec_x variable was not set for field with name ``%s''\n", field->name); fflush(stdout); exit(0);}
	if (field->vec_y == NULL) 	{printf("vec_y variable was not set for field with name ``%s''\n", field->name); fflush(stdout); exit(0);}
	if (field->vec_z == NULL) 	{printf("vec_z variable was not set for field with name ``%s''\n", field->name); fflush(stdout); exit(0);}
	if (field->vec_t == NULL) 	{printf("vec_t variable was not set for field with name ``%s''\n", field->name); fflush(stdout); exit(0);}
	if (field->n_x == 0) 		{printf("Size for vec_x was 0 for field with name ``%s''\n", field->name); fflush(stdout); exit(0);}
	if (field->n_y == 0) 		{printf("Size for vec_y was 0 for field with name ``%s''\n", field->name); fflush(stdout); exit(0);}
	if (field->n_z == 0) 		{printf("Size for vec_z was 0 for field with name ``%s''\n", field->name); fflush(stdout); exit(0);}
	if (field->n_t == 0) 		{printf("Size for vec_t was 0 for field with name ``%s''\n", field->name); fflush(stdout); exit(0);}
	
	field->n = field->n_x * field->n_y * field->n_z * field->n_t;
	
	field->prepared = 1;	
}

void fields_prepare2(
	t_zephyros_field **fieldpt,
	double xmin, double xmax, double xdif,
	double ymin, double ymax, double ydif,
	double zmin, double zmax, double zdif,
	double tmin, double tmax, double tdif
	)
{
	t_zephyros_field *field;	
	double dummy_delta = 1.e-100;

	int i, j;
	int in;
	int ix, iy, iz, it;
	int tmpsum;
	
	fields_initialize(&field);
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
	
	field->prepared = 1;
}


void fields_copy(t_zephyros_field **pdst, t_zephyros_field *src)
{
	int i;
	t_zephyros_field *dst;
		
	if (src == NULL) {
		dst = NULL;
	} else {
		fields_assert_initialized(src);

		dst = calloc(1, sizeof(t_zephyros_field));

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




double fields_give_xyzt(t_zephyros_field *field, int i, int i_xyzt)
{
	int ndim = 4;
	int *arr_shape = calloc(ndim, sizeof(double));
	int *tup = calloc(ndim, sizeof(double));
	double res = -999.9;
	
	fields_assert_prepared(field);

	arr_shape[0] = field->n_x;
	arr_shape[1] = field->n_y;
	arr_shape[2] = field->n_z;
	arr_shape[3] = field->n_t;
	
	interpolation_nr2tup(&ndim, arr_shape, &i, tup);
	
	if (i_xyzt == 0) res = field->vec_x[tup[0]];
	if (i_xyzt == 1) res = field->vec_y[tup[1]];
	if (i_xyzt == 2) res = field->vec_z[tup[2]];
	if (i_xyzt == 3) res = field->vec_t[tup[3]];
	
	free(arr_shape);
	free(tup);
	
	return res;
}


double fields_give_x(t_zephyros_field *field, int i)
{
	return fields_give_xyzt(field, i, 0);
}
double fields_give_y(t_zephyros_field *field, int i)
{
	return fields_give_xyzt(field, i, 1);
}
double fields_give_z(t_zephyros_field *field, int i)
{
	return fields_give_xyzt(field, i, 2);
}
double fields_give_t(t_zephyros_field *field, int i)
{
	return fields_give_xyzt(field, i, 3);
}

