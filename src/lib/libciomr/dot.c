
/* $Log$
 * Revision 1.1  2000/02/04 22:53:18  evaleev
 * Initial revision
 *
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

/* does dot product of 2d matrix */

void dot_mat(a,b,n,value)
   double **a, **b, *value;
   int n;

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

/* does dot product of matrix stored as linear array */

void dot_arr(a,b,n,value)
   double *a, *b, *value;
   int n;

   {
      register int i;
      double tval;

      tval = 0.0;
      for (i=0; i < n; i++) {
         tval += a[i]*b[i];
         }
      *value = tval;
      }
