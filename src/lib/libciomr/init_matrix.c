
/* $Log$
 * Revision 1.1  2000/02/04 22:53:19  evaleev
 * Initial revision
 *
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

double ** init_matrix(n,m)
   int n,m;

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

void free_matrix(array,size)
   double **array;
   int size;

   {
      int i;

      for (i=0; i < size ; i++) {
         free(array[i]);
         }

      free(array);
      }
