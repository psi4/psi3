/*
 \file init_array.c
*/

/* $Log$
 * Revision 1.2  2002/03/25 02:43:45  sherrill
 * Update documentation
 *
/* Revision 1.1.1.1  2000/02/04 22:53:19  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.2  1999/11/01 20:10:56  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.1  1991/06/15 18:29:19  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

extern void resource_command(void);

/*
** init_array(): This function initializes an array of doubles of 
** length 'size' and returns a pointer to the first element
*/
double * init_array(int size)
   {
    double *array;

    if ((array = (double *) malloc(sizeof(double)*size))==NULL) {
       fprintf(stderr,"init_array:  trouble allocating memory \n");
       fprintf(stderr,"size = %d\n",size);
       resource_command();
       exit(2);
       }
    bzero(array,sizeof(double)*size);
    return(array);
   }

