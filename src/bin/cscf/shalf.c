/* $Log$
 * Revision 1.1  2000/02/04 22:52:33  evaleev
 * Initial revision
 *
/* Revision 1.3  1999/08/17 19:04:18  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.2  1999/08/11 18:39:03  evaleev
/* Added some checks on the lowest eigenvalue of the overlap matrix.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:28  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id$";

/* construct S-1/2 matrix 'sahalf' using CANONICAL orthogonalization  */

#define EXTERN
#include "includes.h"
#include "common.h"

void shalf()

  {
     int i,nn,num_mo;
     int ii,jj,kk;
     double *eig_vals, **eig_vecs;
     double tol = 1.0e-20;
     double min_eval = 100000.0;
     struct symm *s;
     
     eig_vals = (double *) init_array(nsfmax);
     eig_vecs = (double **) init_matrix(nsfmax,nsfmax);

     mxcoef = 0;
     nmo = 0;
     
/*  diagonalize smat to get eigenvalues and eigenvectors  */

     for (i=0; i < num_ir ; i++) {
        s = &scf_info[i];
        if (nn=s->num_so) {

           rsp(nn,nn,ioff[nn],s->smat,eig_vals,1,eig_vecs,tol);

	   /*--- Find the lowest eigenvalue ---*/
	   for(ii=0; ii < nn ; ii++)
	     if (eig_vals[ii] < min_eval)
	       min_eval = eig_vals[ii];

           if(print & 64) {
             fprintf(outfile,"\noverlap eigenstuff\n");
             eivout(eig_vecs,eig_vals,nn,nn,outfile);
             }

	   /*--- Go through the eigenvalues and "throw away"
	     the dangerously small ones ---*/
	   num_mo = 0;
	   for(ii=0; ii < nn ; ii++)
	     if (eig_vals[ii] >= SEVAL_CUTOFF) {
	       eig_vals[ii] = 1.0/sqrt(eig_vals[ii]);
	       num_mo++;
	     }
	   s->num_mo = num_mo;
	   mxcoef += num_mo * nn;
	   nmo += num_mo;
	   if (num_mo < nn)
	     fprintf(outfile,"\n  In symblk %d %d eigenvectors of S with eigenvalues < %lf are thrown away",
		     i,nn-num_mo,SEVAL_CUTOFF);

/* form 'sahalf' matrix sahalf = U*(s-1/2)  */

           for (ii=0; ii < nn ; ii++)
              for (jj=0; jj < num_mo ; jj++)
		s->sahalf[ii][jj] += eig_vecs[ii][jj+nn-num_mo]*eig_vals[jj+nn-num_mo];

	}
     }

     free(eig_vals);
     free(eig_vecs);

     fprintf(outfile,"\n  The lowest eigenvalue of the overlap matrix was %e\n\n",
	     min_eval);
   
   }
