
/* $Log$
 * Revision 1.1  2000/02/04 22:50:48  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:47  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

/* function to calculate h matrix (eqns 5-7 in ref 2) */

hemat()

{
   int i,j,ij,k,kl,ik,ii,m,n,im,in,mm,nn,mn,it;
   int p1,p2,iabc,num,ilast;
   double sumc,summ,sumn,esim,esin;
   double elec1c,elec1m,elec1n;
   double elec2c,elec2m,elec2n;
   double elec3c,elec3m,elec3n;
   double pval,qval;
   double a1,a2;
   double *hmo;
   double **denp,**denq,**denm;
   double **e11,**e22,**e12;
   double **scr1,**scr2,*temp;
   struct o_pkints *o_pkpt;

   hmo = (double *) init_array(nbstri);
   temp = (double *) init_array(nbfso*nbfso);
   scr1 = (double **) init_matrix(nbfso,nbfso);
   scr2 = (double **) init_matrix(nbfso,nbfso);
   e11 = (double **) init_matrix(nbfso,nbfso);
   e22 = (double **) init_matrix(nbfso,nbfso);
   e12 = (double **) init_matrix(nbfso,nbfso);
   denp = (double **) init_matrix(5,nbstri);
   denq = (double **) init_matrix(5,nbstri);
   denm = (double **) init_matrix(5,nbstri);

   a12 = (struct off_diags *) malloc(sizeof(struct off_diags)*nind);

   work2=92;
   rfile(work2);

   mread(hmo,14);
   ao_to_mo(hmo,scr1,e_vecs_so,scr2,nbfso,nbfso);

   for(i=0,elec1c=0.0; i < nc ; i++) {
      ii=ioff[i]+i;
      elec1c += 2.0*hmo[ii];
      }

   m=nc;n=nc+1;
   mm=ioff[nc]+nc;
   nn=ioff[nc+1]+nc+1;
   elec1m = 2.0*hmo[mm];
   elec1n = 2.0*hmo[nn];

/* form density matrices */

   for(i=ij=0; i < nbfso ; i++) {
      esim = e_vecs_so[i][nc];
      esin = e_vecs_so[i][nc+1];
      for(j=0; j <= i ; j++,ij++) {
         for(k=0,sumc=0.0; k < nc ; k++)
            sumc += e_vecs_so[i][k]*e_vecs_so[j][k];
         denp[0][ij] = 2.0*sumc;
         denq[0][ij] = -sumc;
         summ = esim*e_vecs_so[j][nc];
         sumn = esin*e_vecs_so[j][nc+1];
         denp[1][ij] = 2.0*summ;
         denq[1][ij] = -summ;
         denp[2][ij] = 2.0*sumn;
         denq[2][ij] = -sumn;
         denp[3][ij] = 0.0;
         denq[3][ij] = summ;
         denp[4][ij] = 0.0;
         denq[4][ij] = sumn;
         }
      }

   for(it=0; it < 5 ; it++) {
      for(i=ij=0; i < nbfso ; i++) {
         for(j=0; j <= i ; j++,ij++) {
            denp[it][ij] = (i==j) ? denp[it][ij] : 2.0*denp[it][ij];
            denq[it][ij] = (i==j) ? denq[it][ij] : 2.0*denq[it][ij];
            }
         }
      }

/* for half-transformed density matrix */

   srew(itap37);

   do {
      sread(itap37,(char *) o_pkbuf,sizeof(struct o_pkints)*maxbuf);
      num = (o_pkbuf[0].ij % (maxbuf+1));
      ilast = o_pkbuf[0].kl % (maxbuf+1);
      o_pkbuf[0].ij /= (maxbuf+1);
      o_pkbuf[0].kl /= (maxbuf+1);

      o_pkpt=o_pkbuf;
      for (i=num; i ; i--,o_pkpt++) {
         ij = (*o_pkpt).ij;
         kl = (*o_pkpt).kl;
         pval = (*o_pkpt).pval;
         qval = (*o_pkpt).kval;

         for(it=0; it < 5 ; it++) {
            denm[it][ij] += denp[it][kl]*pval + denq[it][kl]*qval;
            denm[it][kl] += denp[it][ij]*pval + denq[it][ij]*qval;
            }
         }
      } while(!ilast);

   for(it=0; it < 5 ; it++) ao_to_mo(denm[it],scr1,e_vecs_so,scr2,nbfso,nbfso);

/* calculate bare lagrangian matrices (eqns 18-22 ref 2) */

   elec2c=elec2m=elec2n=0.0;
   for(i=0; i < nc ; i++) {
      ii=ioff[i]+i;
      elec2c += denm[0][ii];
      elec2m += 2.0*denm[1][ii];
      elec2n += 2.0*denm[2][ii];
      }
   elec2m += denm[3][mm];
   elec2n += denm[4][nn];

   h11=elec1c+elec1m+elec2c+elec2m;
   h22=elec1c+elec1n+elec2c+elec2n;
   h12=denm[3][nn];
   h21=denm[4][mm];

   elec=h11*c1sq+h22*c2sq+h12*c12p;
   fprintf(outfile,"\n\telec = %20.10f\n",elec);

   mn=ioff[n]+m;
   for(i=0; i < nbfso ; i++) {
      im = (i > m) ? ioff[i]+m : ioff[m]+i;
      in = (i > n) ? ioff[i]+n : ioff[n]+i;
      e12[i][m]=0.5*denm[4][im];
      e12[i][n]=0.5*denm[3][in];
      for(k=0; k < nocc ; k++) {
         ik = (i > k) ? ioff[i]+k : ioff[k]+i;
         if(k!=n) e11[i][k]=hmo[ik]+denm[0][ik]+denm[1][ik];
         if(k!=m) e22[i][k]=hmo[ik]+denm[0][ik]+denm[2][ik];
         }
      }

   for(i=0; i < nbfso ; i++)
      for(j=0; j < nbfso ; j++)
         scr1[i][j]=e11[i][j]*c1sq+e22[i][j]*c2sq+e12[i][j]*c12p;

   if(print & 2048) {
      fprintf(outfile,"\n e11\n");
      print_mat(e11,nbfso,nbfso,outfile);
      fprintf(outfile,"\n e22\n");
      print_mat(e22,nbfso,nbfso,outfile);
      fprintf(outfile,"\n e12\n");
      print_mat(e12,nbfso,nbfso,outfile);
      fprintf(outfile,"\n ee\n");
      print_mat(scr1,nbfso,nbfso,outfile);
      }

   for(ii=0; ii < nind ; ii++) {
      i=indep[ii].ii;
      j=indep[ii].jj;
      a1=(e11[i][j]-e11[j][i])*c1 + (e12[i][j]-e12[j][i])*c2;
      a2=(e12[i][j]-e12[j][i])*c1 + (e22[i][j]-e22[j][i])*c2;
      a12[ii].one = 2.0*a1;
      a12[ii].two = 2.0*a2;
      }

   for(i=ij=0; i < nbfso ; i++)
      for(j=0; j < nbfso ; j++,ij++)
         temp[ij]=e11[j][i];
   swrit(work2,(char *) temp,sizeof(double)*nbfso*nbfso);

   for(i=ij=0; i < nbfso ; i++)
      for(j=0; j < nbfso ; j++,ij++)
         temp[ij]=e22[j][i];
   swrit(work2,(char *) temp,sizeof(double)*nbfso*nbfso);

   for(i=ij=0; i < nbfso ; i++)
      for(j=0; j < nbfso ; j++,ij++)
         temp[ij]=e12[j][i];
   swrit(work2,(char *) temp,sizeof(double)*nbfso*nbfso);

   srew(work2);

   free(hmo);
   free(temp);
   free_matrix(scr1,nbfso);
   free_matrix(scr2,nbfso);
   free_matrix(e11,nbfso);
   free_matrix(e12,nbfso);
   free_matrix(e22,nbfso);
   free_matrix(denp,5);
   free_matrix(denq,5);
   free_matrix(denm,5);
   }
