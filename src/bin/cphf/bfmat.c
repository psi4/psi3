
/* $Log$
 * Revision 1.1  2000/02/04 22:50:46  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:35  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

bfmat()

   {
      int i,j,ij,ii,it,ix;
      int m,n,mm,nn;
      int p1,p2;
      double f;
      double *temp,*hf,*bf,*t2;
      double **dipda,**dipdb,**epf;
      double ***hm;

      rfile(itap43);

      hf = (double *) init_array(nbatri);
      dipda = (double **) init_matrix(3,natom3);
      dipdb = (double **) init_matrix(3,natom3);

      hm = (double ***) malloc(sizeof(double **)*9);
      for(i=0; i < 9 ; i++) hm[i] = (double **) init_matrix(natom,nbatri);

      dipole_derivs(dipda,dipdb,hf,hm);

      for(i=0; i < 9 ; i++) free_matrix(hm[i],natom);

      temp = (double *) init_array(natom3*3);

      if(iopen) {
         for(i=0; i < 3 ; i++) hm[i] = (double **) init_matrix(ntypes,nbstri);
         t2 = (double *) init_array(nbfso*nbfso);
         bf = (double *) init_array(nbstri);
         epf = (double **) init_matrix(nbfso,nbfso);
         }

      srew(itap43);
      sread(itap43,(char *) temp,sizeof(double)*natom3*3);
      sread(itap43,(char *) temp,sizeof(double)*natom3*3);

      if(print & 1) {
         fprintf(outfile,"\ndipda matrix\n");
         print_mat(dipda,3,natom3,outfile);
         fprintf(outfile,"\ndipdb matrix\n");
         print_mat(dipdb,3,natom3,outfile);
         }

      m=nc;n=nc+1;
      mm=ioff[m]+m;
      nn=ioff[n]+n;
      for(ij=0; ij < 3 ; ij++) {
         sread(itap43,(char *) hf,sizeof(double)*nbstri);
         rwrit(itap44,(char *) hf,sizeof(double)*nbstri,ha_loc[ij+natom3]);
         if(!iopen)
            rwrit(itap44,(char *) hf,sizeof(double)*nbstri,ba_loc[ij+natom3]);
         else {
            for(j=0; j < ntypes ; j++) {
               f = 0.5*focc[j];
               for(i=0; i < nbstri ; i++) hm[ij][j][i] = hf[i]*f;
               }
            }
         if(twocon) {
            dpm[ij]=hf[mm];
            dpn[ij]=hf[nn];
            for(i=0,f=0.0; i < nc ; i++) {
               ii=ioff[i]+i;
               f += hf[ii];
               }
            h11f[ij]=2.0*(f+hf[mm]);
            h22f[ij]=2.0*(f+hf[nn]);
            }

         if(print & 32) {
            fprintf(outfile,"\nhf matrix ij = %d\n",ij);
            print_array(hf,nbfso,outfile);
            }
         }

/* store derivative lagrangian matrices */

      if(iopen) {
         for(ix=0; ix < 3 ; ix++) {
            rread(itap44,(char *) hf,sizeof(double)*nbstri,ha_loc[ix+natom3]);
            bzero(bf,sizeof(double)*nbstri);
            zero_mat(epf,nbfso,nbfso);

            for(it=0; it < ntypes ; it++)
               for(i=nstart[it]; i < nend[it] ; i++)
                  for(j=0; j < nbfso ; j++) {
                     p1=MAX0(i,j);
                     p2=MIN0(i,j);
                     ij = ioff[p1]+p2;
                     epf[i][j] = hm[ix][it][ij];
                     }

            for(it=0; it < nind ; it++) {
               i=indep[it].ii;
               j=indep[it].jj;
               ij = indep[it].ij;
               bf[ij] = epf[i][j]-epf[j][i];
               }

            rwrit(itap44,(char *) bf,sizeof(double)*nbstri,ba_loc[ix+natom3]);

            for(i=ij=0; i < nbfso ; i++)
               for(j=0; j < nbfso ; j++,ij++)
                  t2[ij] = epf[j][i];
   
            rwrit(itap44,(char *) t2,sizeof(double)*nbstri,ea_loc[ix+natom3]);
            if(print & 32) {
               fprintf(outfile,"\nepf matrix ix = %d\n",ix);
               print_mat(epf,nbfso,nbfso,outfile);
               }
            }
         }

      srew(itap43);

      free(hf);
      free(temp);
      free_matrix(dipda,3);
      free_matrix(dipdb,3);
      if(iopen) {
         free_matrix(epf,nbfso);
         free(t2);
         free(bf);
         for(i=0; i < 3 ; i++) free_matrix(hm[i],ntypes);
         }
      free(hm);
      }
