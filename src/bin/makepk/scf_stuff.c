
/* $Log$
 * Revision 1.1  2000/02/04 22:51:34  evaleev
 * Initial revision
 *
/* Revision 1.4  1997/09/12 13:54:51  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 1.3  1997/08/25  21:52:33  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.2  1993/05/31  12:50:46  seidl
 * small bug fix for same symmetry singlets
 *
 * Revision 1.1  1991/06/15  22:06:29  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

scf_stuff()

   {
      int i,j,ij;
      PSI_FPTR junk,loccal,locvec;
      int junkpoint;
      int nn,to;
      int mconst,mpoint,ncalcs,mcalcs;
      double *coeffs;
      double *temp;
      int i30[400];
      double r30[200];

      e_vals = (double *) init_array(nbfso);
      occ_num = (double *) init_array(nbfso);
      coeffs = (double *) init_array(mxcoef);
      charges = (double *) init_array(natom);
      e_vecs_so = (double **) init_matrix(nbfso,nbfso);
      e_vecs_ao = (double **) init_matrix(nbfao,nbfao);
      ao_to_so = (double **) init_matrix(nbfao,nbfso);
      temp = (double *) init_array(nbfao*nbfao);
      coords = (double **) init_matrix(natom,3);
      if(iopen) {
         alpha = (double *) init_array(iopen);
         beta = (double *) init_array(iopen);
         }

      junk = (PSI_FPTR) sizeof(int)*100;
      wreadw(itap30,(char *) i30,sizeof(int)*200,junk,&junk);

      mpoint = i30[1];
      mconst = i30[2];
      mcalcs = i30[3];
      ncalcs = i30[4];

      junk = (PSI_FPTR) sizeof(int)*(100+mconst);
      wreadw(itap30,(char *) i30,sizeof(int)*mpoint,junk,&junk);
      
/* get nuclear charges */

      junk = (PSI_FPTR) sizeof(int)*(i30[0]-1);
      wreadw(itap30,(char *) charges,sizeof(double)*natom,junk,&junk);

      mwrit(charges,7);

/* get ao to so transformation matrix */

      junk = (PSI_FPTR) sizeof(int)*(i30[28]-1);
      wreadw(itap30,(char *) temp,sizeof(double)*nbfao*nbfso,junk,&junk);

      for (i=0,ij=0; i < nbfso ; i++)
         for (j=0; j < nbfao ; j++,ij++)
            ao_to_so[j][i] = temp[ij];

      mwrit(temp,8);

/* read in geometry */

      junk = (PSI_FPTR) sizeof(int)*(100+mconst+mpoint+ncalcs-1);
      wreadw(itap30,(char *) &junkpoint,sizeof(int),junk,&junk);
      loccal = (PSI_FPTR) junkpoint;
      wreadw(itap30,(char *) i30,sizeof(int)*60,
             (PSI_FPTR) sizeof(int)*(loccal-1),&loccal);
      wreadw(itap30,(char *) i30,sizeof(int)*20,loccal,&loccal);

      locvec = (PSI_FPTR) sizeof(int)*(i30[0]-1);

      for (i=0; i < natom ; i++)
         wreadw(itap30,(char *) coords[i],sizeof(double)*3,loccal,&loccal);

/* read in scf energy */

      wreadw(itap30,(char *) r30,sizeof(double)*10,loccal,&loccal);
      enuc = r30[0];
      escf = r30[1];

/* get mo coefficients in the so basis   */

      wreadw(itap30,(char *) coeffs,sizeof(double)*mxcoef,locvec,&locvec);

      wreadw(itap30,(char *) e_vals,sizeof(double)*nbfso,locvec,&locvec);

      mwrit(e_vals,9);

/* get occupation numbers and alpha and beta arrays */

      locvec += (PSI_FPTR) sizeof(int)*n_so_typs;
      wreadw(itap30,(char *) nlamda,sizeof(int)*n_so_typs,locvec,&locvec);
      wreadw(itap30,(char *) nclosd,sizeof(int)*n_so_typs,locvec,&locvec);

      if(iopen) {
         wreadw(itap30,(char *) nopen,sizeof(int)*n_so_typs,locvec,&locvec);
         wreadw(itap30,(char *) alpha,sizeof(double)*iopen,locvec,&locvec);
         wreadw(itap30,(char *) beta,sizeof(double)*iopen,locvec,&locvec);
         if(twocon) {
            two_occ[0] = 2.0/(1.0-alpha[0]);
            two_occ[1] = 2.0/(1.0-alpha[2]);
            }
         }

/* set up occupations */
         
      for (i=j=nn=to=0; i < num_ir ; i++) {
         if(num_so[i]) {
            ij = block_num[i];
            for ( ; j < nn+nclosd[ij] ; j++) occ_num[j] = 2.0;
            if(nopen[ij]) {
              if(twocon) {
                int k;
                for(k=0; k < nopen[ij]; k++,j++,to++)
                  occ_num[j] = two_occ[to];
                }
              else for ( ; j < nn+nclosd[ij]+nopen[ij] ; j++) occ_num[j] = 1.0;
              }
            for ( ; j < nn+num_so[i] ; j++) occ_num[j] = 0.0;
            nn += num_so[i];
            }
         }

      mwrit(occ_num,10);

/* form so eigenvector matrix */

      for (i=j=ij=0; ij < n_so_typs ; ij++) {
         int k,l,n=nlamda[ij];
         for (k=j; k < j+n ; k++)
            for (l=j; l < j+n ; l++,i++)
               e_vecs_so[l][k] = coeffs[i];
         j += n;
         }

      for (i=0,ij=0; i < nbfso ; i++)
         for (j=0; j < nbfso ; j++,ij++)
            temp[ij] = e_vecs_so[j][i];

      mwrit(temp,11);

/* transform eigenvectors from so to ao basis */

      mmult(ao_to_so,0,e_vecs_so,0,e_vecs_ao,0,nbfao,nbfso,nbfso,0);

      for (i=0,ij=0; i < nbfso ; i++)
         for (j=0; j < nbfao ; j++,ij++)
            temp[ij] = e_vecs_ao[j][i];

      mwrit(temp,12);

      if(print > 2) {
         fprintf(outfile,"\n e_vecs_ao matrix\n");
         eigout(e_vecs_ao,e_vals,occ_num,nbfao,nbfao,outfile);
         }

      free(temp);
      free(coeffs);
      }
