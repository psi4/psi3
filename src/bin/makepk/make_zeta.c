
/* $Log$
 * Revision 1.1  2000/02/04 22:51:33  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:06:29  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

make_zeta(ii,jj,kk,ll,value)
   unsigned int ii,jj,kk,ll;
   double value;
{

   int ityp,iiii,iijj,iikk,iill,kkkk,jjjj,kkll,jjkk,jjll,kkjj;

   if(ii==jj) {
     if(ii==kk) {
       if(ii==ll) {           /* type (ii/ii) */
         iiii=ioff[ii]+ii;
         for(ityp=0; ityp < ntypes ; ityp++) {
           zeta_so[ityp][iiii] += (densa[ityp][iiii]+densb[ityp][iiii])*value;
           }
         }
       else {                /* type (ii/il) */
         iiii=ioff[ii]+ii;
         iill=ioff[ii]+ll;
         for(ityp=0; ityp < ntypes ; ityp++) {
           zeta_so[ityp][iiii] += (densa[ityp][iill]+densb[ityp][iill])*value*2.0;
           zeta_so[ityp][iill] += (densa[ityp][iiii]+densb[ityp][iiii])*value;
           }
         }
       }
     else if (kk==ll) {       /* type (ii/kk) */
       iiii=ioff[ii]+ii;
       kkkk=ioff[kk]+kk;
       iikk=ioff[ii]+kk;
       for(ityp=0; ityp < ntypes ; ityp++) {
         zeta_so[ityp][iiii] += densa[ityp][kkkk]*value;
         zeta_so[ityp][kkkk] += densa[ityp][iiii]*value;
         zeta_so[ityp][iikk] += densb[ityp][iikk]*value;
         }
       }
     else {               /* type (ii/kl)  */
       iiii=ioff[ii]+ii;
       kkll=ioff[kk]+ll;
       iikk=ioff[ii]+kk;
       iill=ioff[ii]+ll;
       for(ityp=0; ityp < ntypes ; ityp++) {
         zeta_so[ityp][iiii] += densa[ityp][kkll]*value*2.0;
         zeta_so[ityp][kkll] += densa[ityp][iiii]*value;
         zeta_so[ityp][iikk] += densb[ityp][iill]*value;
         zeta_so[ityp][iill] += densb[ityp][iikk]*value;
         }
       }
     }
   else if(ii==kk) {
     if(jj==ll) {       /* type (ij/ij) */
       iijj=ioff[ii]+jj;
       iiii=ioff[ii]+ii;
       jjjj=ioff[jj]+jj;
       for(ityp=0; ityp < ntypes ; ityp++) {
         zeta_so[ityp][iijj] += (2.0*densa[ityp][iijj]+densb[ityp][iijj])*value;
         zeta_so[ityp][iiii] += densb[ityp][jjjj]*value;
         zeta_so[ityp][jjjj] += densb[ityp][iiii]*value;
         }
       }
     else {          /* type (ij/il) */
       iijj=ioff[ii]+jj;
       iill=ioff[ii]+ll;
       iiii=ioff[ii]+ii;
       jjll=ioff[jj]+ll;
       for(ityp=0; ityp < ntypes ; ityp++) {
         zeta_so[ityp][iijj] += (2.0*densa[ityp][iill]+densb[ityp][iill])*value;
         zeta_so[ityp][iill] += (2.0*densa[ityp][iijj]+densb[ityp][iijj])*value;
         zeta_so[ityp][iiii] += densb[ityp][jjll]*value*2.0;
         zeta_so[ityp][jjll] += densb[ityp][iiii]*value;
         }
       }
     }
   else if(jj==kk) {
     if(jj==ll) {     /* type (ij/jj) */
       iijj=ioff[ii]+jj;
       jjjj=ioff[jj]+jj;
       for(ityp=0; ityp < ntypes ; ityp++) {
         zeta_so[ityp][iijj] += (densa[ityp][jjjj]+densb[ityp][jjjj])*value;
         zeta_so[ityp][jjjj] += (densa[ityp][iijj]+densb[ityp][iijj])*value*2.0;
         }
       }
     else {         /* type (ij/jl) */
       iijj=ioff[ii]+jj;
       jjll=ioff[jj]+ll;
       iill=ioff[ii]+ll;
       jjjj=ioff[jj]+jj;
       for(ityp=0; ityp < ntypes ; ityp++) {
         zeta_so[ityp][iijj] += (2.0*densa[ityp][jjll]+densb[ityp][jjll])*value;
         zeta_so[ityp][jjll] += (2.0*densa[ityp][iijj]+densb[ityp][iijj])*value;
         zeta_so[ityp][jjjj] += densb[ityp][iill]*value*2.0;
         zeta_so[ityp][iill] += densb[ityp][jjjj]*value;
         }
       }
     }
   else if(jj==ll) {        /* type (ij/kj) */
     iijj=ioff[ii]+jj;
     kkjj=ioff[kk]+jj;
     iikk=ioff[ii]+kk;
     jjjj=ioff[jj]+jj;
     for(ityp=0; ityp < ntypes ; ityp++) {
       zeta_so[ityp][iijj] += (2.0*densa[ityp][kkjj]+densb[ityp][kkjj])*value;
       zeta_so[ityp][kkjj] += (2.0*densa[ityp][iijj]+densb[ityp][iijj])*value;
       zeta_so[ityp][jjjj] += densb[ityp][iikk]*value*2.0;
       zeta_so[ityp][iikk] += densb[ityp][jjjj]*value;
       }
     }
   else if(kk==ll) {        /* type (ij/kk) */
     iijj=ioff[ii]+jj;
     iikk=ioff[ii]+kk;
     kkkk=ioff[kk]+kk;
     jjkk=(jj > kk) ? ioff[jj]+kk : ioff[kk]+jj;
     for(ityp=0; ityp < ntypes ; ityp++) {
       zeta_so[ityp][iijj] += densa[ityp][kkkk]*value;
       zeta_so[ityp][kkkk] += densa[ityp][iijj]*value*2.0;
       zeta_so[ityp][iikk] += densb[ityp][jjkk]*value;
       zeta_so[ityp][jjkk] += densb[ityp][iikk]*value;
       }
     }
   else {               /* type (ij/kl) */
     iijj=ioff[ii]+jj;
     iikk=ioff[ii]+kk;
     iill=ioff[ii]+ll;
     kkll=ioff[kk]+ll;
     jjkk=(jj > kk) ? ioff[jj]+kk : ioff[kk]+jj;
     jjll=(jj > ll) ? ioff[jj]+ll : ioff[ll]+jj;
     for(ityp=0; ityp < ntypes ; ityp++) {
       zeta_so[ityp][iijj] += densa[ityp][kkll]*value*2.0;
       zeta_so[ityp][kkll] += densa[ityp][iijj]*value*2.0;
       zeta_so[ityp][iikk] += densb[ityp][jjll]*value;
       zeta_so[ityp][jjll] += densb[ityp][iikk]*value;
       zeta_so[ityp][iill] += densb[ityp][jjkk]*value;
       zeta_so[ityp][jjkk] += densb[ityp][iill]*value;
       }
     }
   }
     
