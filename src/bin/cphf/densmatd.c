
/* $Log$
 * Revision 1.1  2000/02/04 22:50:47  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:45:28  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

densmatd(dp,bb,nvec)
   double **dp,**bb;
   int nvec;

{
   int i,j,ij,ii,iabc;
   double f,fb;
   double *esi1,*esi2,*esj1,*esj2;
   double *dp_ij,*bb_ii,**dq_ii;

   for(i=ij=0; i < nbfso ; i++) {
      for(j=0; j <= i ; j++,ij++) {
         esi1=so_vecs_k[i];
         esi2=so_vecs_l[i];
         esj1=so_vecs_l[j];
         esj2=so_vecs_k[j];
         dq_ii=bb;
         for(ii=nind; ii ; ii--,esi1++,esi2++,esj1++,esj2++,dq_ii++) {
            f = *esi1 * *esj1;
            fb = *esi2 * *esj2;
            f += fb;
            if(f) {
               f *= 2.0;
               dp_ij=dp[ij];
               bb_ii= *dq_ii;
               for(iabc=nvec; iabc ; iabc--,dp_ij++,bb_ii++)
                  *dp_ij += *bb_ii*f;
               }
            }
         }
      }

   for(i=ij=0; i < nbfso ; i++)
      for(j=0; j <= i ; j++,ij++)
         for(iabc=0; iabc < nvec ; iabc++)
            dp[ij][iabc] = (i == j) ? dp[ij][iabc] : 2.0*dp[ij][iabc];

   }
