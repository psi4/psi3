#include <ruby.h>
#include "psirb.h"

VALUE create_array(unsigned int count, int *array);
VALUE create_array(unsigned int count, double *array);
int create_array(VALUE arr, int **array);
int create_array(VALUE arr, double **array);

VALUE create_array(unsigned int count, int *array)
{
	VALUE arr;
	int idx;
	
	// Create a new Ruby array
	arr = rb_ary_new();
	
	// Push the values of the values of the C-array to Ruby
	for (idx = 0; idx < count; ++idx)
		rb_ary_push(arr, INT2FIX(array[idx]));
		
	return arr;
}

VALUE create_array(unsigned int count, double *array)
{
	VALUE arr;
	int idx;
	
	// Create a new Ruby array
	arr = rb_ary_new();
	
	// Push the values
	for (idx = 0; idx < count; ++idx)
		rb_ary_push(arr, rb_float_new(array[idx]));
		
	return arr;
}

int create_array(VALUE arr, int **array)
{
	unsigned int count = 0;
	int idx;
	
	// Make sure we are sent an array
	if (TYPE(arr) == T_ARRAY) {
		count = RARRAY(arr)->len;
		*array = (int*)malloc(sizeof(int) * count);

		for (idx = 0; idx < count; ++idx)
			(*array)[idx] = NUM2INT(RARRAY(arr)->ptr[idx]);
	}
	else {
		rb_raise(rb_eTypeError, "expected an array");
	}
	return count;
}

int create_array(VALUE arr, double **array)
{
	unsigned int count = 0;
	int idx;
	
	// Make sure we are sent an array
	if (TYPE(arr) == T_ARRAY) {
		count = RARRAY(arr)->len;
		*array = (double*)malloc(sizeof(double) * count);

		for (idx = 0; idx < count; ++idx)
			(*array)[idx] = NUM2DBL(RARRAY(arr)->ptr[idx]);
	}
	else {
		rb_raise(rb_eTypeError, "expected an array");
	}
	return count;
}

extern "C" {
	int * init_int_array(int size)
	{
	   int *array;

	   if ((array = (int *) malloc(sizeof(int)*size))==NULL) {
	      fprintf(stderr,"init_array:  trouble allocating memory \n");
	      fprintf(stderr,"size = %d\n",size);
	      exit(0);
	      }
	   bzero(array,sizeof(int)*size);
	   return(array);
	}
}