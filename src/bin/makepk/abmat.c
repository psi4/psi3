
/* $Log$
 * Revision 1.1  2000/02/04 22:51:32  evaleev
 * Initial revision
 *
/* Revision 1.3  1997/09/12 13:54:47  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 1.2  1997/08/25  21:52:23  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.1  1991/06/15  22:06:29  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

abmat()

{
   int i,j,ij,k;
   int mm,nn,blk;
   PSI_FPTR junk;
   int iijj,ii,jj;
   int nstri,nendi,nstrj,nendj;
   double fi,fj;

   alpb = (double *) init_array(nbatri);
   betb = (double *) init_array(nbatri);

   nstri=nendi=0;
   for(i=0; i < ntypes ; i++) {
      if(mm=nsorb[i]) {
         fi=focc[i];
         nstri = nendi;
         nendi = nstri+mm;
         nstrj=nendj=0;
         for (j=0; j <= i ; j++) {
            if(nn=nsorb[j]) {
               fj=focc[j];
               nstrj = nendj;
               nendj = nstrj+nn;
               for (ii=nstri; ii < nendi ; ii++) {
                  for (jj=nstrj; jj < nendj ; jj++) {
                     ij = MAX0(ii,jj);
                     iijj = MIN0(ii,jj);
                     iijj = ioff[ij]+iijj;
                     ij = ioff[i]+j;
                     alpb[iijj] = (1.0-alpc[ij])*fi*fj*0.5;
                     betb[iijj] = -(1.0-betc[ij])*fi*fj*0.25;
                     }
                  }
               }
            }
         }
      }

   if(print > 2) {
      fprintf(outfile,"\n alpb matrix\n");
      print_array(alpb,nbfso,outfile);
      fprintf(outfile,"\n betb matrix\n");
      print_array(betb,nbfso,outfile);
      }

/* transform alpd and betd to scf ordering */

   if(ci_calc) {
      int ncl,nop;
      int itap47=47,ia47[192];
      int *order,*lorder;
      double valc,valo,valo2;
      double *alp,*bet,*alpd,*betd;
      double **scr1,**scr2,**scr3,*temp;
      double *denc,*deno,*deno2;

      denc = (double *) init_array(nbatri);
      deno = (double *) init_array(nbatri);
      deno2 = (double *) init_array(nbatri);
      temp = (double *) init_array(nbfao*nbfao);
      alp = (double *) init_array(nbatri);
      bet = (double *) init_array(nbatri);
      alpd = (double *) init_array(nbatri);
      betd = (double *) init_array(nbatri);
      order = (int *) init_array(nbfso);
      lorder = (int *) init_array(nbfso);
      scr2 = (double **) init_matrix(nbfao,nbfao);

      rfile(itap47);
      rfile(49);
      wreadw(itap47,(char *) ia47,(PSI_FPTR) sizeof(int)*192,0,&junk);

      mread(order,3);
      wreadw(itap47,(char *) alpd,sizeof(double)*nbstri,
             (PSI_FPTR) sizeof(int)*(ia47[105]-1),&junk);
      wreadw(itap47,(char *) betd,sizeof(double)*nbstri,
             (PSI_FPTR) sizeof(int)*(ia47[106]-1),&junk);

      for(i=0; i < nbfso ; i++) {
         mm=order[i]-1;
         for (j=0; j <= i ; j++) {
            nn=order[j]-1;
            ij = ioff[i]+j;
            ii=MAX0(mm,nn);
            jj=MIN0(mm,nn);
            iijj = ioff[ii]+jj;
            alp[iijj] = alpd[ij];
            bet[iijj] = betd[ij];
            }
         }

      mwrit(alp,39);
      mwrit(bet,40);
      wwritw(itap47,(char *) alp,sizeof(double)*nbstri,
             (PSI_FPTR) sizeof(int)*(ia47[105]-1),&junk);
      wwritw(itap47,(char *) bet,sizeof(double)*nbstri,
             (PSI_FPTR) sizeof(int)*(ia47[106]-1),&junk);
            
      if(print > 2) {
         fprintf(outfile,"\n alpd matrix\n");
         print_array(alpd,nbfso,outfile);
         fprintf(outfile,"\n betd matrix\n");
         print_array(betd,nbfso,outfile);
   
         fprintf(outfile,"\n alp matrix\n");
         print_array(alp,nbfso,outfile);
         fprintf(outfile,"\n bet matrix\n");
         print_array(bet,nbfso,outfile);
         }

 /* write the rest of file47 and file49 */

      mread(temp,8);
      wwritw(itap47,(char *) temp,sizeof(double)*nbfso*nbfso,
                     (PSI_FPTR)  sizeof(int)*(ia47[107]-1),&junk);
      mread(temp,11);
      wwritw(itap47,(char *) temp,sizeof(double)*nbfso*nbfso,
                     (PSI_FPTR)  sizeof(int)*(ia47[108]-1),&junk);
      mread(temp,12);
      wwritw(itap47,(char *) temp,sizeof(double)*nbfao*nbfso,
                     (PSI_FPTR)  sizeof(int)*(ia47[109]-1),&junk);

      wwritw(itap47,(char *) order,sizeof(int)*nbfso,
             (PSI_FPTR) sizeof(int)*(ia47[116]-1),&junk);

      mread(temp,18);
      wwritw(itap47,(char *) temp,sizeof(double)*nbfso*nbfso,
                     (PSI_FPTR)  sizeof(int)*(ia47[117]-1),&junk);
      mread(temp,19);
      wwritw(itap47,(char *) temp,sizeof(double)*nbfao*nbfso,
                     (PSI_FPTR)  sizeof(int)*(ia47[118]-1),&junk);

      for(i=0; i < nbfso ; i++) {
         j=order[i];
         lorder[j-1]=i+1;
         }
      wwritw(itap47,(char *) lorder,sizeof(int)*nbfso,
             (PSI_FPTR) sizeof(int)*(ia47[101]-1),&junk);
                                       

 /* make density matrices for file49 */

      for(i=ncl=nop=0; i < n_so_typs ; i++) {
         ncl+=nclosd[i];
         nop+=nopen[i];
         }

      for(i=ij=0; i < nbfao ; i++) {
         for(j=0; j <= i ; j++,ij++) {
            for(k=0,valc=0.0; k < ncl ; k++)
               valc += e_vecs_ao[i][k]*e_vecs_ao[j][k];
            denc[ij]=valc;
            if(hsos) {
               for(k=ncl,valo=0.0; k < ncl+nop ; k++)
                  valo += e_vecs_ao[i][k]*e_vecs_ao[j][k];
               deno[ij]=valo;
               }
            else if(iopen) {
               deno[ij]=e_vecs_ao[i][ncl]*e_vecs_ao[j][ncl];
               deno2[ij]=e_vecs_ao[i][ncl+1]*e_vecs_ao[j][ncl+1];
               }
            }
         }

      sread(49,(char *) temp,sizeof(int)*4);
      swrit(49,(char *) denc,sizeof(double)*nbatri);
      if(iopen) swrit(49,(char *) deno,sizeof(double)*nbatri);
      if(ntypes > 2) swrit(49,(char *) deno2,sizeof(double)*nbatri);

      wwritw(itap47,(char *) denc,sizeof(double)*nbatri,
             (PSI_FPTR) sizeof(int)*(ia47[113]-1),&junk);
      wwritw(itap47,(char *) deno,sizeof(double)*nbatri,
             (PSI_FPTR) sizeof(int)*(ia47[114]-1),&junk);
      wwritw(itap47,(char *) deno2,sizeof(double)*nbatri,
             (PSI_FPTR) sizeof(int)*(ia47[115]-1),&junk);

      if(print > 2) {
         fprintf(outfile,"\n closed density matrix\n");
         print_array(denc,nbfao,outfile);
         fprintf(outfile,"\n open density matrix\n");
         print_array(deno,nbfao,outfile);
         fprintf(outfile,"\n open density matrix 2\n");
         print_array(deno2,nbfao,outfile);
         }

/* transform lagrangian from mo to ao basis */

      mread(temp,12);
      for(i=ij=0; i < nbfso ; i++)
         for(j=0; j < nbfao ; j++,ij++)
            scr2[j][i]=temp[ij];
      wreadw(itap47,(char *) temp,sizeof(double)*nbstri,
             (PSI_FPTR) sizeof(int)*(ia47[110]-1),&junk);
      for(ii=iijj=0; ii < nbfao ; ii++) {
         for(jj=0; jj <= ii ; jj++,iijj++) {
            for(i=0,valc=0.0; i < nbfso ; i++) {
               for(j=0; j < nbfso ; j++) {
                  ij=(i > j) ? ioff[i]+j : ioff[j]+i;
                  valc += temp[ij]*scr2[ii][i]*scr2[jj][j];
                  }
               }
            denc[iijj]=valc;
            }
         }

      swrit(49,(char *) denc,sizeof(double)*nbatri);
      wwritw(itap47,(char *) denc,sizeof(double)*nbatri,
             (PSI_FPTR) sizeof(int)*(ia47[111]-1),&junk);
   
      rclose(itap47,3);
      rclose(49,3);

      if(print > 2) {
         fprintf(outfile,"\n lagrangian in ao basis\n");
         print_array(denc,nbfao,outfile);
         }

      free(temp);
      free(lorder);
      free(order);
      free(alp);
      free(bet);
      free(alpd);
      free(betd);
      free(denc);
      free(deno);
      free(deno2);
      free_matrix(scr2,nbfso);
      }
   }
