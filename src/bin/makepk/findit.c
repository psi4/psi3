
/* $Log$
 * Revision 1.2  2002/05/10 05:44:06  crawdad
 * Changed to variable name "nint" to "nnint" to avoid conflict with typedef
 * in Tru64's math.h.
 * -TDC
 *
/* Revision 1.1.1.1  2000/02/04 22:51:33  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.1  1991/06/15 22:06:29  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

extern double *pa, *pb;
extern unsigned int *lbij,*lbkl;
static int *inext;

static int old_nint=25920;
extern int intmx;

void findit(ii,jj,kk,ll,ism,ksm,value,iab)
   double value;
   int ii,jj,kk,ll;
   int ism, ksm, iab;
{
   register int i,j;
   unsigned int *ijtmp, *kltmp;
   int *nxtmp;
   int next,start;
   int lij,lkl,ityp;
   int keep=128;
   int keep2=127;
   int d2i = sizeof(double)/sizeof(int);
   double value2;
   double *patmp, *pbtmp;

   if(nbfso > 150) {
     keep = 1024;
     keep2 = 1023;
     }

   if(!inext) inext = (int *) init_array((keep+intmx)/2);

   lij = ii*(ii+1)/2 + jj;
   lkl = kk*(kk+1)/2 + ll;

   if(!nnint) {
      bzero(inext,sizeof(int)*old_nint);
      bzero(&inext[intmx],sizeof(int)*keep);
      }

   start = 2*lij + lkl;
   start = (start & keep2) + intmx;

L1:
   next=inext[start];
   if(next) {
      if (lbij[next-1] == lij && lbkl[next-1] == lkl) i=next-1;
      else {
         start = next;
         goto L1;
         }
      }
   else {
      i=nnint;
      inext[start] = ++nnint;
      if(nnint >= intmx) {
        fprintf(outfile,"\n  increasing size of buffers in findit\n");
        fprintf(outfile,"  intmx was %d, is %d\n",intmx,intmx*2);
        fflush(outfile);
        intmx*=2;

 /* i don't use realloc because strange things were happening */

        nxtmp = (int *) init_array((int) (keep+intmx)/d2i);
        bcopy(inext,nxtmp,(int)sizeof(int)*(intmx/2));
        for(j=0; j < keep ; j++) nxtmp[j+intmx]=inext[j+intmx/2];
        free(inext);
        inext=nxtmp;

        ijtmp = (unsigned int *) init_array(intmx/d2i);
        bcopy(lbij,ijtmp,sizeof(int)*(intmx/2));
        free(lbij);
        lbij = ijtmp;

        kltmp = (unsigned int *) init_array(intmx/d2i);
        bcopy(lbkl,kltmp,sizeof(int)*(intmx/2));
        free(lbkl);
        lbkl = kltmp;

        patmp = (double *) init_array(intmx);
        bcopy(pa,patmp,sizeof(double)*(intmx/2));
        free(pa);
        pa = patmp;

        pbtmp = (double *) init_array(intmx);
        bcopy(pb,pbtmp,sizeof(double)*(intmx/2));
        free(pb);
        pb = pbtmp;

        if(inext==NULL || lbij==NULL || lbkl==NULL) {
          fprintf(outfile,"\n pathological problems with realloc in findit\n");
          fprintf(outfile," try upping intmx to %d\n",intmx);
          fflush(outfile);
          exit(1);
          }
        }
      lbij[i] = lij;
      lbkl[i] = lkl;
      pa[i] = pb[i] = 0.0;
      }

   value2 = (lij == lkl) ? value*0.5 : value;

   switch(iab) {
      case 1:
         pa[i] += value2;
         break;
      case 2:
         pa[i] -= 0.25*value2;
         pb[i] += 0.25*value2;
         break;
      case 3:
         pa[i] += 0.75*value2;
         pb[i] += 0.25*value2;
         break;
      case 4:
         pa[i] -= 0.5*value2;
         pb[i] = 0.5*value2;
         break;
      case 5:
         pa[i] = 0.5*value2;
         pb[i] = 0.5*value2;
      }
   old_nint=nnint;
   }
