
/* $Log$
 * Revision 1.1  2000/02/04 22:50:48  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:46  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

/* function to calculate ha matrices (eqns 33-35 in ref 2) */

hamat()

{
   int i,j,ij,ii,mm,nn;
   int kt,ntril;
   int p1,p2,iabc;
   double elec1,elec2,valu;
   double *hm,**tm,*sm,*eps;

   ntril=(2*nbstri-1)/1024+1;
   hm = (double *) init_array(nbstri);
   sm = (double *) init_array(nbstri);
   eps = (double *) init_array(nbstri);
   tm = (double **) init_matrix(6,nbstri);

   mread(eps,24);

   srew(work);

   mm=ioff[nc]+nc;
   nn=ioff[nc+1]+nc+1;

   for(iabc=0; iabc < natom3 ; iabc++) {
      rread(itap44,(char *) sm,sizeof(double)*nbstri,sa_loc[iabc]);
      for(i=ij=0,valu=0.0; i < nocc ; i++) {
         for(j=0; j < i ; j++,ij++) {
            valu -= 2.0*sm[ij]*eps[ij];
            }
         valu -= sm[ij]*eps[ij];
         ij++;
         }
      e1t[iabc]=2.0*valu;
         
      rread(itap44,(char *) hm,sizeof(double)*nbstri,ha_loc[iabc]);
      for(ii=0; ii < 6 ; ii++) {
         kt=ii*natom3+iabc;
         kt=kt*ntril+1;
         rread(work,(char *) tm[ii],sizeof(double)*nbstri,kt);
         }
      elec1=elec2=0.0;
      h11=2.0*hm[mm]+tm[2][mm];
      h22=2.0*hm[nn]+tm[4][nn];
      for(i=0; i < nc ; i++) {
         ii=ioff[i]+i;
         elec1 += 2.0*hm[ii];
         elec2 += 2.0*tm[0][ii]-tm[1][ii];
         h11 += 4.0*tm[2][ii]-2.0*tm[3][ii];
         h22 += 4.0*tm[4][ii]-2.0*tm[5][ii];
         }
      ha11[iabc]=elec1+elec2+h11;
      ha22[iabc]=elec1+elec2+h22;
      ha12[iabc]=tm[5][mm];
      }

   c1sq=0.5*focc[1];
   c2sq=0.5*focc[2];
   c12p=2.0*beta[1][2];
   c1=sqrt(c1sq);
   c2=sqrt(c2sq);
   c2 = (c12p < 0) ? -c2 : c2;
   for(iabc=0; iabc < natom3 ; iabc++)
      e1t[iabc] += ha11[iabc]*c1sq+ha22[iabc]*c2sq+ha12[iabc]*c12p;

   if(print & 512) {
      fprintf(outfile,"\n    ha11          ha22          ha12      e1t\n");
      for(iabc=0; iabc < natom3 ; iabc++)
         fprintf(outfile,"%20.10f %20.10f %20.10f %20.10f\n",
                          ha11[iabc],ha22[iabc],ha12[iabc],e1t[iabc]);
      }

   srew(work);
   free(hm);
   free(sm);
   free(eps);
   free_matrix(tm,6);
   }
