
/* $Log$
 * Revision 1.1  2000/02/04 22:53:21  evaleev
 * Initial revision
 *
/* Revision 2.2  1999/11/01 20:10:57  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.1  1991/06/15 18:29:31  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

extern void mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
	   int nr, int nl, int nc, int add);
extern void mxmbol();

/* multiplies two square matrices together */
/* if in (n=a,b,c) = 1, then normal multiply */
/* if jn (n=a,b,c) = 1, then multiply by transpose of n */

void mxmb(a,ia,ja, b,ib,jb, c,ic,jc, nrow, nlnk, ncol)
   double **a, **b, **c;
   int ia, ja, ib, jb, ic, jc, nrow, nlnk, ncol;

   {
      if (ic == 1) {
         if (ia == 1) {
            if (ib == 1) {
               mmult(a,0,b,0,c,0,nrow,nlnk,ncol,0);
               }
            else {
               if (jb == 1) {
                  mmult(a,0,b,1,c,0,nrow,nlnk,ncol,0);
                  }
               else {
                  mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                  }
               }
            }
         else {
            if (ja == 1) {
               if (ib == 1) {
                  mmult(a,1,b,0,c,0,nrow,nlnk,ncol,0);
                  }
               else {
                  if (jb == 1) {
                     mmult(a,1,b,1,c,0,nrow,nlnk,ncol,0);
                     }
                  else {
                     mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                     }
                  }
               }
            else {
               mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
               }
            }
         }
      else {
         if (jc == 1) {
            if (ia == 1) {
               if (ib == 1) {
                  mmult(a,0,b,0,c,1,nrow,nlnk,ncol,0);
                  }
               else {
                  if (jb == 1) {
                     mmult(a,0,b,1,c,1,nrow,nlnk,ncol,0);
                     }
                  else {
                     mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                     }
                  }
               }
            else {
               if (ja == 1) {
                  if (ib == 1) {
                     mmult(a,1,b,0,c,1,nrow,nlnk,ncol,0);
                     }
                  else {
                     if (jb == 1) {
                        mmult(a,1,b,1,c,1,nrow,nlnk,ncol,0);
                        }
                     else {
                        mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                        }
                     }
                  }
               else {
                  mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
                  }
               }
            }
         else {
            mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol);
            }
         }
      }
