/*!
** \file sq_rsp.c
** \ingroup (CIOMR)
*/

/* $Log$
 * Revision 1.4  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.3  2002/04/19 21:48:06  sherrill
/* Remove some unused functions and do doxygen markup of libciomr.
/*
/* Revision 1.2  2002/04/04 22:24:38  evaleev
/* Converted allocation functions (init_array, etc.) to take unsigned long ints
/* to be able to allocate properly 2GB+ chunks). Some declarations cleaned up.
/*
/* Revision 1.1.1.1  2000/02/04 22:53:23  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.3  1999/11/01 20:11:00  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.2  1998/02/03 19:34:11  evaleev
/* Modified eigsort(), rsp(), and sq_rsp() to sort eigenvalues and
/* eigenvectors in either ascending OR descending order.
/*
 * Revision 2.1  1991/06/15  18:30:07  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

#include "includes.h"

extern void tred2(int, double **, double *, double *, int);
extern int tqli(int, double *, double **, double *, int, double);
extern void eigsort(double *, double **, int);

/* translation into c of a translation into FORTRAN77 of the EISPACK */
/* matrix diagonalization routines */

/*!
** sq_rsp: diagomalize a square matrix ('array').
**
**   \param nm     = rows of matrix
**   \param n      = columns of matrix
**   \param nv     = number of elements in lower triangle (n*(n+1)/2)
**   \param array  = matrix to diagonalize
**   \param e_vals = array to hold eigenvalues
**   \param matz   = 0 (no eigenvectors, eigenvals in ascending order)
**                 = 1 (eigenvectors and eigenvalues in ascending order)
**                 = 2 (no eigenvectors, eigenvalues in descending order)
**                 = 3 (eigenvectors and eigenvalues in descending order)
**   \param e_vecs = matrix of eigenvectors (one column for each eigvector)
**   \param toler  = tolerance for eigenvalues?  Often 1.0E-14.
**
** \ingroup (CIOMR)
*/
void sq_rsp(int nm, int n, double **array, double *e_vals, int matz,
            double **e_vecs, double toler)
   {
      int i, j, ii, ij, ierr;
      int ascend_order;
      double *fv1, **temp;
      double zero = 0.0;
      double one = 1.0;

/* Modified by Ed - matz can have the values 0 through 3 */
      
      if ((matz > 3) || (matz < 0)) {
        matz = 0;
        ascend_order = 1;
        }
      else
        if (matz < 2)
          ascend_order = 1;	/* Eigenvalues in ascending order */
        else {
          matz -= 2;
          ascend_order = 0;	/* Eigenvalues in descending order */
          }

      fv1 = (double *) init_array(n);
      temp = (double **) init_matrix(n,n);

      if (n > nm) {
         ierr = 10*n;
         fprintf(stderr,"n = %d is greater than nm = %d in rsp\n",n,nm);
         exit(ierr);
         }

      for (i=0; i < n; i++) {
         for (j=0; j < n; j++) {
            e_vecs[i][j] = array[i][j];
            }
          }

      tred2(n,e_vecs,e_vals,fv1,matz);

      for (i=0; i < n; i++)
         for (j=0; j < n; j++)
            temp[i][j]=e_vecs[j][i];

      tqli(n,e_vals,temp,fv1,matz,toler);

      for (i=0; i < n; i++)
         for (j=0; j < n; j++)
            e_vecs[i][j]=temp[j][i];

      if (ascend_order)
        eigsort(e_vals,e_vecs,n);
      else
        eigsort(e_vals,e_vecs,(-1)*n);

      free(fv1);
      free_matrix(temp,n);

      }

