
/* $Log$
 * Revision 1.1  2000/02/04 22:53:24  evaleev
 * Initial revision
 *
/* Revision 2.2  1991/08/21 05:41:05  psi
/* pass in ioff
/*
 * Revision 2.1  1991/06/15  18:30:14  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"

/* converts matrix in lower triangle form to array stored in separate blocks */

void tri_to_block(a,b,num_ir,num_so,ioff)
   double *a,**b;
   int num_ir,num_so[],ioff[];

   {
      int i,j,k;
      int joff=0;
      int blk=0;

      for (k=0; k < num_ir ; k++) {
         if (num_so[k] > 0) {
            for (i=0; i < num_so[k] ; i++) {
               for (j=0; j <= i ; j++) {
                  b[blk][ioff[i]+j] = a[ioff[i+joff]+j+joff];
                  }
               }
            joff += num_so[k];
            blk++;
            }
         }
      }

void block_to_tri(a,b,num_ir,num_so,ioff)
   double *a,**b;
   int num_ir,num_so[],ioff[];

   {
      int i,j,k;
      int joff=0;
      int blk=0;

      for (k=0; k < num_ir ; k++) {
         if (num_so[k] > 0) {
            for (i=0; i < num_so[k] ; i++) {
               for (j=0; j <= i ; j++) {
                  a[ioff[i+joff]+j+joff] = b[blk][ioff[i]+j];
                  }
               }
            joff += num_so[k];
            blk++;
            }
         }
      }
