
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

densmatd_o(dp,dq,bb,nvec)
   double ***dp,***dq,**bb;
   int nvec;

{
   int i,j,ij,it,iabc,kt,lt;
   int ii;
   double vala,valb,fb,f;
   double *esi1,*esi2,*esj1,*esj2;
   double *dp_ij_it,*dq_ij_it;
   double *alpa_kt,*alpa_lt;
   double *beta_kt,*beta_lt;
   double *bb_ii;

      for(i=ij=0; i < nbfso ; i++) {
         esi1=so_vecs_k[i];
         esi2=so_vecs_l[i];
         for(j=0; j <= i ; j++,ij++) {
            esj1=so_vecs_l[j];
            esj2=so_vecs_k[j];
            for(ii=0; ii < nind ; ii++) {
               f = esi1[ii]*esj1[ii];
               fb = esi2[ii]*esj2[ii];
               f += fb;
               if(f) {
                  kt = indep[ii].it;
                  lt = indep[ii].jt;
                  alpa_kt=alpa[kt];
                  alpa_lt=alpa[lt];
                  beta_kt=beta[kt];
                  beta_lt=beta[lt];
                  for(it=0; it < ntypes ; it++) {
                     vala = f*(alpa_kt[it]-alpa_lt[it]);
                     valb = f*(beta_kt[it]-beta_lt[it]);
                     dp_ij_it=dp[ij][it];
                     dq_ij_it=dq[ij][it];
                     bb_ii=bb[ii];
                     for(iabc=nvec; iabc ; 
                                    iabc--,dp_ij_it++,dq_ij_it++,bb_ii++) {
                        *dp_ij_it += *bb_ii*vala; 
                        *dq_ij_it += *bb_ii*valb; 
                        }
                     }
                  }
               }
            }
         }
         

     for(i=ij=0; i < nbfso ; i++)
       for(j=0; j <= i ; j++,ij++) 
        for(it=0; it < ntypes ; it++) 
          for(iabc=0; iabc < nvec ; iabc++) {
            dp[ij][it][iabc] = (i == j) ? dp[ij][it][iabc] : 2.0*dp[ij][it][iabc];
            dq[ij][it][iabc] = (i == j) ? dq[ij][it][iabc] : 2.0*dq[ij][it][iabc];
            }
   }
