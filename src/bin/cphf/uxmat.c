
/* $Log$
 * Revision 1.1  2000/02/04 22:50:50  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:57  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

uxmat()

{
   int i,j,ij,iabc;
   double *sm,*b,**u,*temp;

   sm = (double *) init_array(nbstri);
   b = (double *) init_array(nbstri);
   u = (double **) init_matrix(nbfso,nbfso);
   temp = (double *) init_array(nbfso*nbfso);

   for(iabc=0; iabc < natom3x ; iabc++) {
      rread(itap44,(char *) b,sizeof(double)*nbstri,ba_loc[iabc]);
      if(iabc < natom3) {
         rread(itap44,(char *) sm,sizeof(double)*nbstri,sa_loc[iabc]);
         for(i=ij=0; i < nbfso ; i++) {
            for(j=0; j < i ; j++,ij++) {
               u[i][j]=b[ij];
               u[j][i]= -b[ij]-sm[ij];
               }
            u[i][i] = -0.5*sm[ij];
            ij++;
            }
         }
      else {
         for(i=ij=0; i < nbfso ; i++,ij++) {
            u[i][i]=0.0;
            for(j=0; j < i ; j++,ij++) {
               u[i][j] =  b[ij];
               u[j][i] = -b[ij];
               }
            }
         }

      for(i=ij=0; i < nbfso ; i++)
         for(j=0; j < nbfso ; j++,ij++)
            temp[ij]=u[j][i];
      rwrit(itap44,(char *) temp,sizeof(double)*nbfso*nbfso,ua_loc[iabc]);

      if(print & 128) {
         fprintf(outfile,"\nu matrix in uxmat iabc = %5d\n",iabc);
         print_mat(u,nbfso,nbfso,outfile);
         }
      }
   free_matrix(u,nbfso);
   free(temp);
   free(sm);
   free(b);
   }
