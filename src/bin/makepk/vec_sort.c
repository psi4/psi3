
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

vec_sort()

   {
      int i,j,ij,col;
      int nc,no,nv;
      int colc,colo,colv;
      int otest,where;
      int params[1024];
      int ia47[192];
      int itap47=47;
      int *order;
      double last;
      double **sort_so,*sort_evals,*sort_occ;
      double *temp;
      double real_par[1024];

      hsos=1;

      bzero(real_par,sizeof(double)*1024);
      bzero(params,sizeof(int)*1024);
      bzero(nsorb,sizeof(int)*10);

      sort_so = (double **) init_matrix(nbfso,nbfso);
      sort_evals = (double *) init_array(nbfso);
      sort_occ = (double *) init_array(nbfso);
      temp = (double *) init_array(nbfao*nbfao);
      order = (int *) init_array(nbfso);

      mread(params,1);

      for (i=nc=0; i < 10 ; i++) nc += nclosd[i];
      for (i=no=0; i < 10 ; i++) no += nopen[i];
      nv = nbfso-(nc+no);

      colc=colo=colv=0;
      for (i=0; i < nbfso ; i++) {
         otest = (int) occ_num[i];
         if (twocon && occ_num[i] != 2.0 && occ_num[i] != 0.0) otest = 1;
         switch(otest) {
            case 0:
               where = colv+nc+no;
               sort_occ[where]=occ_num[i];
               sort_evals[where] = e_vals[i];
               order[where] = i+1;
               for (j=0; j < nbfso ; j++) sort_so[j][where]=e_vecs_so[j][i];
               colv++;
               break;
            case 1:
               where = colo+nc;
               sort_occ[where]=occ_num[i];
               sort_evals[where] = e_vals[i];
               order[where] = i+1;
               for (j=0; j < nbfso ; j++) sort_so[j][where]=e_vecs_so[j][i];
               colo++;
               break;
            case 2:
               sort_occ[colc]=occ_num[i];
               sort_evals[colc] = e_vals[i];
               order[colc] = i+1;
               for (j=0; j < nbfso ; j++) sort_so[j][colc]=e_vecs_so[j][i];
               colc++;
               break;
            default:
               fprintf(stderr,"something's amiss in vec_sort\n");
               exit(otest);
            }
         }

/* sort closed, open, and virtual parts of eigenvector matrix */
/* don't sort open part of tcscf */

      for (i=0; i < nc-1 ; i++) {
         colc=i;
         last=sort_evals[i];
         for (j=i+1; j < nc ; j++) {
            if (sort_evals[j] < last) {
               colc=j;
               last = sort_evals[j];
               }
            }
         if (colc != i) {
            sort_evals[colc]=sort_evals[i];
            sort_evals[i]=last;
            colo = order[i];
            order[i]=order[colc];
            order[colc] = colo;
            for (j=0; j < nbfso ; j++) {
               last = sort_so[j][i];
               sort_so[j][i] = sort_so[j][colc];
               sort_so[j][colc] = last;
               }
            }
         }

      if(!twocon) {
         for (i=nc; i < nc+no-1 ; i++) {
            colc=i;
            last=sort_evals[i];
            for (j=i+1; j < nc+no ; j++) {
               if (sort_evals[j] < last) {
                  colc=j;
                  last = sort_evals[j];
                  }
               }
            if (colc != i) {
               sort_evals[colc]=sort_evals[i];
               sort_evals[i]=last;
               colo = order[i];
               order[i]=order[colc];
               order[colc] = colo;
               for (j=0; j < nbfso ; j++) {
                  last = sort_so[j][i];
                  sort_so[j][i] = sort_so[j][colc];
                  sort_so[j][colc] = last;
                  }
               }
            }
         }

      for (i=nc+no; i < nbfso-1 ; i++) {
         colv=i;
         last=sort_evals[i];
         for (j=i+1; j < nbfso ; j++) {
            if (sort_evals[j] < last) {
               colv=j;
               last = sort_evals[j];
               }
            }
         if (colv != i) {
            sort_evals[colv]=sort_evals[i];
            sort_evals[i]=last;
            colo = order[i];
            order[i]=order[colv];
            order[colv] = colo;
            for (j=0; j < nbfso ; j++) {
               last = sort_so[j][i];
               sort_so[j][i] = sort_so[j][colv];
               sort_so[j][colv] = last;
               }
            }
         }

      for (i=0,ij=0; i < nbfso ; i++)
         for (j=0; j < nbfso ; j++,ij++)
            temp[ij] = sort_so[j][i];

      mwrit(temp,18);

/* transform eigenvectors from so to ao basis */

      mmult(ao_to_so,0,sort_so,0,e_vecs_ao,0,nbfao,nbfso,nbfso,0);

      for (i=0,ij=0; i < nbfso ; i++)
         for (j=0; j < nbfao ; j++,ij++)
            temp[ij] = e_vecs_ao[j][i];

      mwrit(temp,19);
      mwrit(sort_occ,17);
      mwrit(sort_evals,16);
      mwrit(order,3);

/* create a few more parameters */

      if(ci_calc && iopen) {
         for(i=0; i < iopen ; i++)
            if(alpha[i] != 0.0 || beta[i] != 1.0) hsos=0;
         }
      else if(iopen) {
         for(i=0; i < iopen ; i++)
            if(alpha[i] != 0.0 || beta[i] != -1.0) hsos=0;
         }

      if(iopen)
         if(hsos) ntypes=2;
         else ntypes=3;
      else ntypes=1;

      if(!iopen) {
         focc[0] = 2.0;
         focc[1] = 0.0;
         alpc[0]=alpc[1]=alpc[2]=0.0;
         betc[0]=betc[1]=betc[2]=0.0;
         nsorb[0]=nc;
         nsorb[1]=nv;
         }
      else if (!twocon) {
         int nnp=ioff[ntypes+1];
         for (i=0; i < nnp ; i++) alpc[i]=betc[i]=0.0;
         focc[0]=2.0;
         nsorb[0]=nc;
         if(hsos) {
            focc[1]=1.0;
            nsorb[1]=no;
            nsorb[2]=nv;
            for(i=0;i < 6; i++) alpc[i]=betc[i]=0.0;
            betc[2] = (ci_calc) ? 1.0 : -1.0;
            }
         else {  /* for the moment hard-wired for open-shell singlet */
            focc[1]=focc[2]=1.0;
            nsorb[1]=1;
            nsorb[2]=1;
            nsorb[3]=nv;
            for(i=0;i < 10; i++) alpc[i]=betc[i]=0.0;
            if(!ci_calc) {
               betc[2] = -1.0;
               betc[4] = 3.0;
               betc[5] = -1.0;
               }
            else {
               betc[2] = 1.0;
               betc[4] = -3.0;
               betc[5] = 1.0;
               }
            }
         fprintf(outfile,"\n   alpc matrix\n");
         print_array(alpc,ntypes+1,outfile);
         fprintf(outfile,"\n   betc matrix\n");
         print_array(betc,ntypes+1,outfile);
         }
      else {
         for (i=0; i < 15 ; i++) alpc[i]=betc[i]=0.0;
         focc[0]=2.0;
         focc[1]=sort_occ[nc];
         focc[2]=sort_occ[nc+1];
         focc[3]=0.0;
         nsorb[0]=nc;
         nsorb[1]=1;
         nsorb[2]=1;
         nsorb[3]=nv;
         alpc[2]=1.0 - 1.0/focc[1];
         alpc[4]=1.0;
         alpc[5]=1.0 - 1.0/focc[2];
         betc[2]=1.0;
         betc[4]=beta[1];
         betc[5]=1.0;
         print_array(alpc,3,outfile);
         print_array(betc,3,outfile);
         }
         
      real_par[0] = enuc;
      real_par[1] = escf;
      for (i=0; i < 15; i++) {
         real_par[i+10]=alpc[i];
         real_par[i+25]=betc[i];
         }
      for (i=0; i < 10 ; i++) real_par[i+40]=focc[i];

      params[5] = (twocon) ? -iopen : iopen;
      params[6] = nc;
      params[7] = no;
      params[8] = no+nc;
      params[17] = ntypes;
      if(iopen) 
         if (twocon) params[18] = 4;
         else params[18] = 3;
      else params[18] = 1;
   
/* for now hardwired for scf 2nd  or ci 1st */
      params[19]=(ci_calc) ? 2 : 1;
      params[20]=(ci_calc) ? 1 : 2;

      for(i=0; i < 10; i++) {
         params[i+30]=num_so[i];
         params[i+40]=nsorb[i];
         }

      mwrit(params,1);
      mwrit(real_par,2);


      if(print > 2) {
         fprintf(outfile,"\n sorted e_vecs_ao matrix\n");
         eigout(e_vecs_ao,sort_evals,sort_occ,nbfao,nbfao,outfile);
         }

      for(i=0; i < nbfso ; i++) {
         occ_num[i]=sort_occ[i];
         e_vals[i]=sort_evals[i];
         for (j=0; j < nbfso ; j++) e_vecs_so[i][j]=sort_so[i][j];
         }

      free(order);
      free(sort_occ);
      free(sort_evals);
      free_matrix(sort_so,nbfso);
      free(temp);
      }

