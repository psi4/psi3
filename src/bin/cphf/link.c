
/* $Log$
 * Revision 1.1  2000/02/04 22:50:49  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:45:28  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

link()
{
   int i,j,ij,ii,it,ix;
   int m,n,mm,nn;
   int p1,p2;
   double f;
   double *temp,*hf,*bf,*t2;
   double **dipda,**dipdb,**epf;
   double ***hm;

   rfile(itap43);

   hf = (double *) init_array(nbatri);
   dipda = (double **) init_matrix(3,natom3);
   dipdb = (double **) init_matrix(3,natom3);

   hm = (double ***) malloc(sizeof(double **)*9);
   for(i=0; i < 9 ; i++) hm[i] = (double **) init_matrix(natom,nbatri);

   dipole_derivs(dipda,dipdb,hf,hm);
   }
