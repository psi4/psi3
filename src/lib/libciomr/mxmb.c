/*!
** \file mxmb.c
** \ingroup (CIOMR)
*/
 
/* $Log$
 * Revision 1.3  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.2  2002/04/19 21:48:06  sherrill
/* Remove some unused functions and do doxygen markup of libciomr.
/*
/* Revision 1.1.1.1  2000/02/04 22:53:21  evaleev
/* Started PSI 3 repository
/*
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

/*!
** mxmb: multiplies two rectangular matrices together.  If in=1 (n=a,b,c),
** then normal multiply.  If jn=1 (n=a,b,c) then multiply the transpose
** of matrix n.
** \ingroup (CIOMR)
*/
void mxmb(double **a, int ia, int ja, double **b, int ib, int jb, 
          double **c, int ic, int jc, int nrow, int nlnk, int ncol)
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
