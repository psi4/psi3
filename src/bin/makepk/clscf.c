
/* $Log$
 * Revision 1.1  2000/02/04 22:51:32  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:06:29  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

clscf()

   {
      int i,j,k,ij;
      double *lag_ao,*lag_mo;

      lag_ao = (double *) init_array(nbatri);
      lag_mo = (double *) init_array(nbstri);

      for(i=0; i < nbfso ; i++) 
         if(occ_num[i]) lag_mo[ioff[i]+i] = 2.0*e_vals[i];

      mwrit(lag_mo,24);

      for (i=0,ij=0; i < nbfao ; i++)
         for (j=0; j <= i ; j++) {
            double val = 0.0;
            for (k=0; k < nbfso ; k++)
               val += lag_mo[ioff[k]+k]*e_vecs_ao[i][k]*e_vecs_ao[j][k];
            lag_ao[ij] = val;
            ij++;
            }

      mwrit(lag_ao,23);

      if(print) {
         fprintf(outfile,"\n lagrangian in mo basis\n");
         print_array(lag_mo,nbfso,outfile);
         fprintf(outfile,"\n lagrangian in ao basis\n");
         print_array(lag_ao,nbfao,outfile);
         }
      }
