
/* $Log$
 * Revision 1.1  2000/02/04 22:53:24  evaleev
 * Initial revision
 *
/* Revision 2.1  1991/06/15 18:30:15  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

/* converts lower triangle to square matrix */

void tri_to_sq(amat,bmat,size)
     double *amat, **bmat;
     int size;

     {
        int i, j, ij;

        ij=0;
        for(i = 0 ; i < size ; i++) {
            for(j = 0 ; j <= i ; j++) {
               bmat[i][j] = amat[ij];
               bmat[j][i] = amat[ij++];
               }
            }
     }
