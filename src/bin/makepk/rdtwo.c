
/* $Log$
 * Revision 1.1  2000/02/04 22:51:34  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:52:32  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:06:29  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

double *pa, *pb;
unsigned int *lbij,*lbkl;
int intmx=25920;

void rdtwo()
{
   int ilsti, nbuf;
   int ibufsz = 8942;
   int ibufqt = 2236;
   int i;
   int iii;
   int ior, ism, jor, jsm;
   int kor, ksm, lor, lsm;
   unsigned int ii,jj,kk,ll;
   int p1,p2;
   double pki_int;
   int *tmp;
   int num_read=0;

   union psi_buffer inbuf;

   if(nbfso > 150) intmx *= 4;
   nint=0;

   inbuf.pki = (double *) init_array(ibufsz/2);
   pa = (double *) init_array(intmx);
   pb = (double *) init_array(intmx);
   lbij = (unsigned int *) init_array(intmx/2);
   lbkl = (unsigned int *) init_array(intmx/2);

   do {
      sread(itap34,(char *) inbuf.lbli,sizeof(int)*ibufsz);
      ilsti=inbuf.lbli[0];
      nbuf=inbuf.lbli[1];
      num_read += nbuf;

      if (print > 6) fprintf(outfile,"%5d\n",nbuf);

      tmp = &inbuf.lbli[2];

      for (i=0 ; i < nbuf ; i++,tmp += 2) {
         lor = *(tmp+1) >> 2;
         lsm = lor >> 8;
         kor = lsm >> 3;
         ksm = kor >> 8;
         kor = (kor & 255) - 1;
         lsm = lsm & 7;
         lor = (lor & 255) - 1;
         iii = *(tmp+1) & 3;
         jsm = *tmp >> 8;
         ior = jsm >> 3;
         ism = ior >> 8;
         ior = (ior & 255) - 1;
         jsm = jsm & 7;
         jor = (*tmp & 255) - 1;
         ii = ior+ideg[ism];
         jj = jor+ideg[jsm];
         kk = kor+ideg[ksm];
         ll = lor+ideg[lsm];
         pki_int = inbuf.pki[i+ibufqt];

         if(!ci_calc) {
            if (ii == jj && ii == kk || jj == kk && jj == ll) {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,5);
               }
            else if (ii == kk || jj == ll) {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,3);

               p1 = MAX0(jj,ll);
               p2 = MIN0(jj,ll);

               findit(ii,kk,p1,p2,ism,ksm,pki_int,4);
               }
            else if (jj == kk) {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,3);

               p1 = MAX0(jj,kk);
               p2 = MIN0(jj,kk);

               findit(ii,ll,jj,kk,ism,ksm,pki_int,4);
               }
            else if (ii == jj || kk == ll) {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,1);

               p1 = MAX0(jj,ll);
               p2 = MIN0(jj,ll);

               findit(ii,kk,p1,p2,ism,ksm,pki_int,2);
               }
            else {
               findit(ii,jj,kk,ll,ism,ksm,pki_int,1);

               p1 = MAX0(jj,ll);
               p2 = MIN0(jj,ll);

               findit(ii,kk,p1,p2,ism,ksm,pki_int,2);

               p1 = MAX0(jj,kk);
               p2 = MIN0(jj,kk);

               findit(ii,ll,p1,p2,ism,ksm,pki_int,2);
               }
            }
         if(iopen || ci_calc) make_zeta(ii,jj,kk,ll,pki_int);

         if(iii && nint) packit(lbij,lbkl,0);
         }
      } while(!ilsti);

   fprintf(outfile,"\n%8cread  %10d integrals from file34\n",' ',num_read);
   packit(lbij,lbkl,1);

   free(inbuf.pki);
   free(pa);
   free(pb);
   free(lbij);
   free(lbkl);
   fflush(outfile);
   }
