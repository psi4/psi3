
/* $Log$
 * Revision 1.1  2000/02/04 22:51:34  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:52:30  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:06:29  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

extern double *pa, *pb;
static int ibl=0,num_bufs=0,num_ints=0;
static struct o_pkints {
         unsigned int ij;
         unsigned int kl;
         double pval;
         double kval;
         } *o_outbuf;

static struct c_pkints {
         unsigned int ij;
         unsigned int kl;
         double pval;
         } *c_outbuf;


void packit(lbij,lbkl,endflg)
   unsigned int *lbij, *lbkl;
   int endflg;

   {
      register int i;
      unsigned int ij,kl;
      double pval, qval;
      double tol = pow(10.0,(double) -toler);

      if(!c_outbuf && !o_outbuf) {
         if(iopen) {
            o_outbuf = 
               (struct o_pkints *) malloc(sizeof(struct o_pkints)*maxbuf);
            if(!o_outbuf) {
               fprintf(stderr,"cannot allocate memory for outbuf in packit\n");
               exit(4);
               }
            }
         else {
            c_outbuf = 
               (struct c_pkints *) malloc(sizeof(struct c_pkints)*maxbuf);
            if(!c_outbuf) {
               fprintf(stderr,"cannot allocate memory for outbuf in packit\n");
               exit(4);
               }
            }
         }

      if(!endflg) {
         for(i=0; i < nint ; i++) {
            pval=pa[i];
            qval=pb[i];
            ij = lbij[i];
            kl = lbkl[i];
            if(print > 5) fprintf(outfile,"%5d%5d%9.5f%9.5f\n",ij,kl,pval,qval);
            if (fabs(pval) >= tol || fabs(qval) >= tol) {
               if(iopen) {
                  o_outbuf[ibl].ij = ij;
                  o_outbuf[ibl].kl = kl;
                  o_outbuf[ibl].pval = pval+qval;
                  o_outbuf[ibl].kval = qval+qval;
                  }
               else {
                  c_outbuf[ibl].ij = ij;
                  c_outbuf[ibl].kl = kl;
                  c_outbuf[ibl].pval = pval;
                  }

               ibl++;

               if (ibl >= maxbuf) {
                  if(iopen) {
                     o_outbuf[0].ij = o_outbuf[0].ij*(maxbuf+1)+maxbuf;
                     o_outbuf[0].kl *= (maxbuf+1); 
                     swrit(itap37,(char *) o_outbuf,sizeof(struct o_pkints)*maxbuf);
                     }
                  else {
                     c_outbuf[0].ij = c_outbuf[0].ij*(maxbuf+1)+maxbuf;
                     c_outbuf[0].kl *= (maxbuf+1); 
                     swrit(itap37,(char *) c_outbuf,sizeof(struct c_pkints)*maxbuf);
                     }
                  num_ints += ibl;
                  num_bufs++;
                  ibl=0;
                  }
               }
            }
         nint=0;
         }
      else {
         num_ints += ibl;
         num_bufs++;
         if(iopen) {
            o_outbuf[0].ij = o_outbuf[0].ij*(maxbuf+1)+ibl;
            o_outbuf[0].kl = o_outbuf[0].kl*(maxbuf+1)+1;
            swrit(itap37,(char *) o_outbuf,sizeof(struct o_pkints)*maxbuf);
            }
         else {
            c_outbuf[0].ij = c_outbuf[0].ij*(maxbuf+1)+ibl;
            c_outbuf[0].kl = c_outbuf[0].kl*(maxbuf+1)+1;
            swrit(itap37,(char *) c_outbuf,sizeof(struct c_pkints)*maxbuf);
            }

         fprintf(outfile,"%8cwrote %10d integrals to file37\n",' ',num_ints);
         if(iopen) free(o_outbuf);
         else free(c_outbuf);
         }
      }
