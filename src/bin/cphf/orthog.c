
/* $Log$
 * Revision 1.1  2000/02/04 22:50:50  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:45:28  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#include "includes.h"

orthog(b,nind,nabc,nvec)
   double **b;
   int nind,nabc,*nvec;

   {
      int i,k,ii,ix;
      int iflag;
      int ntemp = nabc;
      int n = nind;
      int ndum[360];
      double *w,eu;

      w = (double *) init_array(nind);

      bzero(ndum,sizeof(int)*ntemp);

      for(i=0; i < n ; i++) w[i] = b[i][0];

      normalize(w,nind,n,&iflag);

      for(i=0; i < n ; i++) b[i][0] = w[i];

      *nvec = 1;
      if(ntemp == 1) return;

      for(ii=1; ii < ntemp ; ii++) {
         for(i=0; i < n ; i++) w[i] = b[i][ii];
         for(k=0; k < ii ; k++) {
            if(!ndum[k]) {
               eu=0.0;
               for(i=0; i < n ; i++) eu += b[i][k]*b[i][ii];
               for(i=0; i < n ; i++) w[i] -= eu*b[i][k];
               }
            }
         normalize(w,nind,n,&iflag);
         if(iflag) ndum[ii]=1;
         if(!iflag) {
            (*nvec)++;
            for(i=0; i < n ; i++) b[i][ii]=w[i];
            }
         }

      ix=0;
      for(ii=0; ii < ntemp ; ii++) {
         if(!ndum[ii]) {
            for(i=0; i < n ; i++) w[i]=b[i][ii];
            for(i=0; i < n ; i++) b[i][ix]=w[i];
            ix++;
            }
         }
      free(w);
      }

normalize(w,nind,n,flag)
   double *w;
   int nind,n,*flag;

   {
      int i;
      double sql,r,rlim=1.0e-38,tol=1.0e-9;

      *flag=0;
      sql = 0.0;

      for(i=0; i < n ; i++) sql += w[i]*w[i];
      r = sqrt(sql);

      if(r < rlim) {
         fprintf(stderr,"r too small in normalize %11.10e\n",r);
         exit(1);
         }

      if(r <= tol) {
         *flag=1;
         return;
         }

      for(i=0; i < n ; i++) w[i] /= r;
      }
