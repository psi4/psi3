
/* $Log$
 * Revision 1.2  2002/03/25 02:43:45  sherrill
 * Update documentation
 *
/* Revision 1.1.1.1  2000/02/04 22:53:18  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.2  1994/06/02 02:28:51  seidl
/* define ALLOC_GLOBALS
/*
 * Revision 2.1  1991/06/15  18:28:45  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#define ALLOC_GLOBALS
#include "includes.h"
#undef ALLOC_GLOBALS
#include "common.h"

/*!
** dot_mat():
** Takes the dot product between two 2D matrices a and b with dimensions
** n x n and returns the value
*/
void dot_mat(double **a, double **b, int n, double *value)
   {
      register int i,j;
      double *ta, *tb, tval;

      tval = 0.0;
      for (i=0; i < n; i++) {
         ta = a[i];
         tb = b[i];
         for (j=0; j < n; j++,ta++,tb++) {
            tval += (*ta) * (*tb);
            }
         }
      *value = tval;
      }

/*!
** dot_arr():
** Take the dot product of the first n elements of two arrays a and b
** and put the result in variable value.
*/
void dot_arr(double *a, double *b, int n, double *value)
   {
      register int i;
      double tval;

      tval = 0.0;
      for (i=0; i < n; i++) {
         tval += a[i]*b[i];
         }
      *value = tval;
      }
