#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <file30.h>
#include <qt.h>
#include "globaldefs.h"
#include "globals.h"

/*
** form_diag_mo_hess
**
** Calculates the approximate diagonal MO Hessian
**
** I am assuming that the pairs (p,q) are always given such that p>=q
**
** C. David Sherrill
** April 1998
*/
void form_diag_mo_hess(int npairs, int *ppair, int *qpair, double *F_core, 
                       double *tei, double **opdm, double *tpdm, double *F_act,
                       int firstact, int lastact, double *hess)
{

  int pair, p, q, pq, pp, qq;
  int i,ii,a,aa,t,tt,u,tu,v,w,vw,tuvw;
  double value;

  /* loop over the independent pairs */
  for (pair=0; pair<npairs; pair++) {
    p = ppair[pair];
    q = qpair[pair];
    pq = ioff[p] + q;
    pp = ioff[p] + p;
    qq = ioff[q] + q;
  
    /* H_{ai,ai}, i.e., inactive virt/inactive occ */
    if (p >= lastact && q < firstact) {
      hess[pair] = 4.0 * (F_core[pp] + F_act[pp] - F_core[qq] - F_act[qq]);
    }

    /* H_{at,at}, i.e., inactive virt with active orb */
    else if (p >= lastact && q >= firstact) {
      a = p;  t = q;
      aa = ioff[a] + a;

      hess[pair] = 2.0 * opdm[t][t] * (F_core[aa] + F_act[aa]);

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        value += opdm[t][u] * F_core[tu];
      }
      hess[pair] -= 2.0 * value;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
 
        /* loop over active v,w ... later restrict to symmetry allowed */
        for (v=firstact; v<lastact; v++) {
          for (w=firstact; w<lastact; w++) {
            vw = INDEX(v,w);
            tuvw = INDEX(tu,vw); 
            value += tpdm[tuvw] * tei[tuvw];
          }
        }
         
      }
      hess[pair] -= 2.0 * value; 
       
    } 

    /* H_{ti,ti}, i.e., active orb with inactive occ */
    else if (p >= firstact && q < firstact) {
      t = p;  i = q;
      tt = ioff[t] + t;
      ii = ioff[i] + i; 
      
      hess[pair] = 2.0 * opdm[t][t] * (F_core[ii] + F_act[ii]);
      hess[pair] += 4.0 * (F_core[tt] + F_act[tt] - F_core[ii] - F_act[ii]);

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
        value += opdm[t][u] * F_core[tu];
      }
      hess[pair] -= 2.0 * value;

      value = 0.0;
      for (u=firstact; u<lastact; u++) {
        tu = INDEX(t,u);
 
        /* loop over active v,w ... later restrict to symmetry allowed */
        for (v=firstact; v<lastact; v++) {
          for (w=firstact; w<lastact; w++) {
            vw = INDEX(v,w);
            tuvw = INDEX(tu,vw); 
            value += tpdm[tuvw] * tei[tuvw];
          }
        }
         
      }
      hess[pair] -= 2.0 * value; 
       
    }

    else {
      fprintf(outfile, 
             "(form_diag_mo_hess): Error, unrecognized class of indep pair\n");
    }

  } /* end loop over pairs */

}


