/* $Log$
 * Revision 1.1  2000/02/04 22:52:30  evaleev
 * Initial revision
 *
/* Revision 1.2  1999/08/17 19:04:15  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:26  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"
#include "libciomr.h"

void form_vec()
{
   int i,nn,num_mo;
   int k,l;
   double **ctrans;
   double **temp;
   double **sqhmat;
   double **sqhmat2;
   double tol=1.0e-20;
   struct symm *s;

   ctrans = (double **) init_matrix(nsfmax,nsfmax); 
   temp = (double **) init_matrix(nsfmax,nsfmax);
   sqhmat = (double **) init_matrix(nsfmax,nsfmax);
   sqhmat2 = (double **) init_matrix(nsfmax,nsfmax);
   
   for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
	 num_mo = s->num_mo;
         tri_to_sq(s->hmat,sqhmat,nn);
         mmult(s->sahalf,1,sqhmat,0,temp,0,num_mo,nn,nn,0);
         mmult(temp,0,s->sahalf,0,sqhmat2,0,num_mo,nn,num_mo,0);
         sq_rsp(num_mo,num_mo,sqhmat2,s->hevals,1,ctrans,tol);
         mxmb(s->sahalf,1,0,ctrans,1,0,s->cmat,1,0,nn,num_mo,num_mo);
         for(k=0;k < nn; k++)
            for(l=0;l < num_mo; l++) s->sahalf[k][l]=s->cmat[k][l];

         if(print & 2) {
            fprintf(outfile,"\nguess vector for irrep %s\n",s->irrep_label);
            print_mat(s->cmat,nn,num_mo,outfile);
            }
         }
      }
   inflg = 0;
   free_matrix(ctrans,nsfmax);
   free_matrix(temp,nsfmax);
   free_matrix(sqhmat,nsfmax);
   free_matrix(sqhmat2,nsfmax);
   }
