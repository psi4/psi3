
/* $Log$
 * Revision 1.1  2000/02/04 22:51:34  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:06:29  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

zetmat()
   {
   int i,j,ij,ityp,nn,nstr,nend;
   int na,ns;
   double *lag_ao,*lag_mo;
   double **zeta_mo;
   double **scr1,**scr2;

   zeta_mo = (double **) init_matrix(ntypes,nbstri);
   scr1 = (double **) init_matrix(nbfao,nbfao);
   scr2 = (double **) init_matrix(nbfao,nbfao);
   lag_ao = (double *) init_array(nbatri);
   lag_mo = (double *) init_array(nbstri);

/* transform zeta matrices from so to mo basis and write to master file */

   if(print) {
      for (ityp=0; ityp < ntypes ; ityp++) {
         fprintf(outfile,"\nzeta_so ityp = %d\n",ityp);
         print_array(zeta_so[ityp],nbfso,outfile);
         }
      }

   nn=nbfso;
   for(i=0; i < ntypes ; i++) {
      tri_to_sq(zeta_so[i],scr1,nbfso);

      mmult(e_vecs_so,1,scr1,0,scr2,0,nn,nn,nn,0);
      mmult(scr2,0,e_vecs_so,0,scr1,0,nn,nn,nn,0);

      sq_to_tri(scr1,zeta_mo[i],nbfso);
      mwrit(zeta_mo[i],26+i);

      if(print) {
         fprintf(outfile,"\nzeta_mo ityp = %d\n",i);
         print_array(zeta_mo[i],nbfso,outfile);
         }
      }

/* then form lagrangian in mo and ao bases */

   nstr=nend=0;
   for(ityp=0; ityp < ntypes ; ityp++) {
      if(nn=nsorb[ityp]) {
         nstr=nend;
         nend += nn;
         for(i=nstr; i < nend ; i++) {
            for(j=0; j <= i ; j++) {
               ij=ioff[i]+j;
               lag_mo[ij]=zeta_mo[ityp][ij];
               }
            }
         }
      }
   mwrit(lag_mo,24);

   na=nbfao;
   ns = nbfso;
   tri_to_sq(lag_mo,scr1,nbfso);
   mmult(e_vecs_ao,0,scr1,0,scr2,0,na,ns,ns,0);
   mmult(scr2,0,e_vecs_ao,1,scr1,0,na,ns,na,0);
   sq_to_tri(scr1,lag_ao,nbfao);
   for(i=0; i < nbatri ; i++) lag_ao[i] *= 2.0;

   mwrit(lag_ao,23);

   if(print > 1) {
      fprintf(outfile,"\n lagrangian matrix mo basis\n");
      print_array(lag_mo,nbfso,outfile);
      fprintf(outfile,"\n lagrangian matrix ao basis\n");
      print_array(lag_ao,nbfao,outfile);
      }
   }
