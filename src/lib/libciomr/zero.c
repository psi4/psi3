
/* $Log$
 * Revision 1.1  2000/02/04 22:53:24  evaleev
 * Initial revision
 *
/* Revision 2.1  1991/06/15 18:30:19  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

void zero_arr(a,size)
   double *a;
   int size;

   {
      bzero(a,sizeof(double)*size);
      }

void zero_mat(a,n,m)
   double **a;
   int n,m;

   {
      register int i;

      for (i=0; i < n ; i++) {
         bzero(a[i],sizeof(double)*m);
         }
      }
