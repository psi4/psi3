/*!
  \file init_matrix.c
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
/* Revision 2.1  1991/06/15 18:29:20  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

extern void resource_command(void);

/* allocates memory for an n x m matrix */
/* returns pointer to pointer to 1st element */

/*!
** init_matrix(): Initialize a matrix of doubles and return a pointer to
** the first row.  Note that this does not form a matrix which is 
** necessarily contiguous in memory.  Use block_matrix() for that.
**
** \param n = number of rows
** \param m = number of columns
*/

double ** init_matrix(int n, int m)
   {
    double **array=NULL;
    int i;

    if ((array = (double **) malloc(sizeof(double *)*n))==NULL) {
         fprintf(stderr,"init_matrix: trouble allocating memory \n");
         fprintf(stderr,"n = %d\n",n);
         resource_command();
         exit(2);
         }

    for (i = 0; i < n; i++) {
        if ((array[i] = (double *) malloc(sizeof(double)*m))==NULL) {
           fprintf(stderr,"init_matrix: trouble allocating memory \n");
           fprintf(stderr,"i = %d m = %d\n",i,m);
           resource_command();
           exit(3);
           }
        bzero(array[i],sizeof(double)*m);
        }
    return(array);
    }

/*!
** free_matrix(): Free a 2D matrix allocated with init_matrix().
**
** \param array = matrix to free
** \param size = number of rows
*/
void free_matrix(double **array, int size)
   {
      int i;

      for (i=0; i < size ; i++) {
         free(array[i]);
         }

      free(array);
      }
