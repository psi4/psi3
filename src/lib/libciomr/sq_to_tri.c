
/* $Log$
 * Revision 1.1  2000/02/04 22:53:23  evaleev
 * Initial revision
 *
/* Revision 2.1  1991/06/15 18:30:08  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

/* converts square matrix to lower triangle */

void sq_to_tri(bmat,amat,size)
     double *amat, **bmat;
     int size;

     {
        int i, j, ij;

        ij=0;
        for(i = 0 ; i < size ; i++) {
            for(j = 0 ; j <= i ; j++) {
               amat[ij++] = bmat[i][j];
               }
            }
     }

