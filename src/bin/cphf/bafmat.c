
/* $Log$
 * Revision 1.1  2000/02/04 22:50:46  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:32  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

/* function to calculate h matrix (eqns 5-7 in ref 2) */

bafmat()

{
   int i,j,ij;
   int iabc,ixyz;
   double fac1,fac2,fac4;
   double b11,b12,b21,b22,ba1,ba2;
   double *sa;
   double **e11,**e22,**e12;
   double **scr1,**scr2,*temp;

   sa = (double *) init_array(nbstri);
   temp = (double *) init_array(nbfso*nbfso);
   scr1 = (double **) init_matrix(nbfso,nbfso);
   scr2 = (double **) init_matrix(nbfso,nbfso);
   e11 = (double **) init_matrix(nbfso,nbfso);
   e22 = (double **) init_matrix(nbfso,nbfso);
   e12 = (double **) init_matrix(nbfso,nbfso);

   srew(work2);
   sread(work2,(char *) temp,sizeof(double)*nbfso*nbfso);
   for(i=ij=0; i < nbfso ; i++)
      for(j=0; j < nbfso ; j++,ij++)
         e11[j][i]=temp[ij];
   sread(work2,(char *) temp,sizeof(double)*nbfso*nbfso);
   for(i=ij=0; i < nbfso ; i++)
      for(j=0; j < nbfso ; j++,ij++)
         e22[j][i]=temp[ij];
   sread(work2,(char *) temp,sizeof(double)*nbfso*nbfso);
   for(i=ij=0; i < nbfso ; i++)
      for(j=0; j < nbfso ; j++,ij++)
         e12[j][i]=temp[ij];


   fac1= -c12p*c2;
   fac2=  c12p*c1;
   for(iabc=0; iabc < natom3x ; iabc++) {
      if(iabc < natom3) {
         rread(itap44,(char *) sa,sizeof(double)*nbstri,sa_loc[iabc]);
         b11 = -ha11[iabc]+e1t[iabc];
         b22 = -ha22[iabc]+e1t[iabc];
         b12 = -ha12[iabc];
         for(i=ij=0; i < nocc ; i++) {
            for(j=0; j <= i ; j++,ij++) {
               fac4 = (i==j) ? 2.0 : 4.0;
               b11 += sa[ij]*e11[j][i]*fac4;
               b22 += sa[ij]*e22[j][i]*fac4;
               b12 += sa[ij]*e12[j][i]*fac4;
               }
            }
         ba1=0.5*(b11*c1+b12*c2);
         ba2=0.5*(b12*c1+b22*c2);
         baf[iabc].one = ba1;
         baf[iabc].two = ba2;
         }
      else {
         ixyz=iabc-natom3;
         baf[iabc].one = (dpm[ixyz]-dpn[ixyz])*fac1;
         baf[iabc].two = (dpm[ixyz]-dpn[ixyz])*fac2;
         }
      }
   
   a22[0][0]=0.5*(h11-elec+c1*c1);
   a22[1][1]=0.5*(h22-elec+c2*c2);
   a22[1][0]=a22[0][1]=0.5*(h12+c1*c2);

   srew(work2);

   free(sa);
   free(temp);
   free_matrix(scr1,nbfso);
   free_matrix(scr2,nbfso);
   free_matrix(e11,nbfso);
   free_matrix(e12,nbfso);
   free_matrix(e22,nbfso);
   }
