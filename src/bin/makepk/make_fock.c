
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

make_fock(ii,jj,kk,ll,value,fock_av,den_av)
   unsigned int ii,jj,kk,ll;
   double value,*fock_av,*den_av;
{

   int iiii,iijj,iikk,iill,kkkk,jjjj,kkll,jjkk,jjll,kkjj;

   if(ii==jj) {
     if(ii==kk) {
       if(ii==ll) {           /* type (ii/ii) */
         iiii=ioff[ii]+ii;
         fock_av[iiii] += den_av[iiii]*value;
         }
       else {                /* type (ii/il) */
         iiii=ioff[ii]+ii;
         iill=ioff[ii]+ll;
         fock_av[iiii] += den_av[iill]*value*2.0;
         fock_av[iill] += den_av[iiii]*value;
         }
       }
     else if (kk==ll) {       /* type (ii/kk) */
       iiii=ioff[ii]+ii;
       kkkk=ioff[kk]+kk;
       iikk=ioff[ii]+kk;
       fock_av[iiii] += den_av[kkkk]*value*2.0;
       fock_av[kkkk] += den_av[iiii]*value*2.0;
       fock_av[iikk] -= den_av[iikk]*value;
       }
     else {               /* type (ii/kl)  */
       iiii=ioff[ii]+ii;
       kkll=ioff[kk]+ll;
       iikk=ioff[ii]+kk;
       iill=ioff[ii]+ll;
       fock_av[iiii] += den_av[kkll]*value*4.0;
       fock_av[kkll] += den_av[iiii]*value*2.0;
       fock_av[iikk] -= den_av[iill]*value;
       fock_av[iill] -= den_av[iikk]*value;
       }
     }
   else if(ii==kk) {
     if(jj==ll) {       /* type (ij/ij) */
       iijj=ioff[ii]+jj;
       iiii=ioff[ii]+ii;
       jjjj=ioff[jj]+jj;
       fock_av[iijj] += 3.0*den_av[iijj]*value;
       fock_av[iiii] -= den_av[jjjj]*value;
       fock_av[jjjj] -= den_av[iiii]*value;
       }
     else {          /* type (ij/il) */
       iijj=ioff[ii]+jj;
       iill=ioff[ii]+ll;
       iiii=ioff[ii]+ii;
       jjll=ioff[jj]+ll;
       fock_av[iijj] += den_av[iill]*value*3.0;
       fock_av[iill] += den_av[iijj]*value*3.0;
       fock_av[iiii] -= den_av[jjll]*value*2.0;
       fock_av[jjll] -= den_av[iiii]*value;
       }
     }
   else if(jj==kk) {
     if(jj==ll) {     /* type (ij/jj) */
       iijj=ioff[ii]+jj;
       jjjj=ioff[jj]+jj;
       fock_av[iijj] += den_av[jjjj]*value;
       fock_av[jjjj] += den_av[iijj]*value*2.0;
       }
     else {         /* type (ij/jl) */
       iijj=ioff[ii]+jj;
       jjll=ioff[jj]+ll;
       iill=ioff[ii]+ll;
       jjjj=ioff[jj]+jj;
       fock_av[iijj] += den_av[jjll]*value*3.0;
       fock_av[jjll] += den_av[iijj]*value*3.0;
       fock_av[jjjj] -= den_av[iill]*value*2.0;
       fock_av[iill] -= den_av[jjjj]*value;
       }
     }
   else if(jj==ll) {        /* type (ij/kj) */
     iijj=ioff[ii]+jj;
     kkjj=ioff[kk]+jj;
     iikk=ioff[ii]+kk;
     jjjj=ioff[jj]+jj;
     fock_av[iijj] += den_av[kkjj]*value*3.0;
     fock_av[kkjj] += den_av[iijj]*value*3.0;
     fock_av[jjjj] -= den_av[iikk]*value*2.0;
     fock_av[iikk] -= den_av[jjjj]*value;
     }
   else if(kk==ll) {        /* type (ij/kk) */
     iijj=ioff[ii]+jj;
     iikk=ioff[ii]+kk;
     kkkk=ioff[kk]+kk;
     jjkk=(jj > kk) ? ioff[jj]+kk : ioff[kk]+jj;
     fock_av[iijj] += den_av[kkkk]*value*2.0;
     fock_av[kkkk] += den_av[iijj]*value*4.0;
     fock_av[iikk] -= den_av[jjkk]*value;
     fock_av[jjkk] -= den_av[iikk]*value;
     }
   else {               /* type (ij/kl) */
     iijj=ioff[ii]+jj;
     iikk=ioff[ii]+kk;
     iill=ioff[ii]+ll;
     kkll=ioff[kk]+ll;
     jjkk=(jj > kk) ? ioff[jj]+kk : ioff[kk]+jj;
     jjll=(jj > ll) ? ioff[jj]+ll : ioff[ll]+jj;
     fock_av[iijj] += den_av[kkll]*value*4.0;
     fock_av[kkll] += den_av[iijj]*value*4.0;
     fock_av[iikk] -= den_av[jjll]*value;
     fock_av[jjll] -= den_av[iikk]*value;
     fock_av[iill] -= den_av[jjkk]*value;
     fock_av[jjkk] -= den_av[iill]*value;
     }
   }
     
