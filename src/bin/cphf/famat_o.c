
/* $Log$
 * Revision 1.1  2000/02/04 22:50:48  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:45  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

/* function to calculate lagrangian derivative matrix (eqn 7) */
/* (eqn 32 in ref 2 for tcscf */

famat_o()

{
   int i,j,k,ij,it,kt;
   int p1,p2,iabc;
   int how_many,ntril;
   double fac,valu;
   double *hm,**tm,**epa;
   double *temp;

   how_many = (twocon) ? 6 : ntypes;
   ntril = (2*nbstri-1)/1024+1;

   hm = (double *) init_array(nbstri);
   epa = (double **) init_matrix(nbfso,nbfso);
   temp = (double *) init_array(nbfso*nbfso);
   tm = (double **) init_matrix(how_many,nbstri);

   srew(work);

   /*       fi*(dh/da)       */

   for(iabc=0; iabc < natom3 ; iabc++) {
      rread(itap44,(char *) hm,sizeof(double)*nbstri,ha_loc[iabc]);
      for(it=0; it < how_many ; it++) {
         kt=it*natom3+iabc;
         kt=kt*ntril+1;
         rread(work,(char *) tm[it],sizeof(double)*nbstri,kt);
         }
      zero_mat(epa,nbfso,nbfso);
      if(!twocon) {
         for(it=0; it < how_many ; it++) {
            fac = focc[it]*0.5;
            for(i=nstart[it]; i < nend[it] ; i++) {
               for(j=0; j < nbfso ; j++) {
                  p1=MAX0(i,j);
                  p2=MIN0(i,j);
                  ij=ioff[p1]+p2;
                  epa[i][j] = hm[ij]*fac+tm[it][ij];
                  }
               }
            }
         }
      else {
         for(i=0; i < nbfso ; i++) {
            it=motyp[i];
            fac=focc[it]*0.5;
            for(j=0; j < nbfso ; j++) {
               p1=MAX0(i,j);
               p2=MIN0(i,j);
               ij=ioff[p1]+p2;
               valu = hm[ij]*fac;
               for(k=kt=0; kt < ntypes ; kt++,k+=2) {
                  valu += alpa[it][kt]*tm[k][ij] +
                          beta[it][kt]*tm[k+1][ij];
                  }
               epa[i][j]=valu;
               }
            }
         }

      for(i=ij=0; i < nbfso ; i++)
         for(j=0; j < nbfso ; j++,ij++)
            temp[ij] = epa[j][i];

      rwrit(itap44,(char *) temp,sizeof(double)*nbfso*nbfso,ea_loc[iabc]);

      if(print & 8) {
         fprintf(outfile,"\nepa matrix iabc = %5d\n",iabc);
         print_mat(epa,nbfso,nbfso,outfile);
         }
      }

   srew(work);
   free(temp);
   free(hm);
   free_matrix(epa,nbfso);
   free_matrix(tm,how_many);
   }
