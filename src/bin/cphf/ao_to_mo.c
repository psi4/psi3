
/* $Log$
 * Revision 1.1  2000/02/04 22:50:46  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:45:28  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#include "includes.h"

ao_to_mo(a,b,c,d,nao,nso)
   double *a,**b,**c,**d;
   int nao,nso;

{
   int i,j,ij,ii;
   double t,*dt,*st;

   for(i=ij=0; i < nao ; i++)
      for(j=0; j <= i ; j++,ij++)
         b[i][j]=b[j][i]=a[ij];

   mmult(c,1,b,0,d,0,nso,nao,nao,0);

   for(i=0; i < nso ; i++)
      for(j=0; j < nao ; j++)
         b[i][j]=c[j][i];

   for(i=ij=0; i < nso ; i++) {
      for(j=0; j <= i ; j++,ij++) {
         dt=d[i];
         st=b[j];
         for(ii=nao,t=0.0; ii ; ii--,dt++,st++)
            t += *dt * *st;
         a[ij] = t;
         }
      }
   }
