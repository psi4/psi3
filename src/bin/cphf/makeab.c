
/* $Log$
 * Revision 1.1  2000/02/04 22:50:49  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:51  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

makeab(bb,bb1,nvec)
   double **bb,**bb1;
   int nvec;

{
   register unsigned i;
   int j,ii,ij,kl;
   int iabc,num,ilast;
   double pval,eji;
   double **scr1,**scr2;
   double **denp,**denmt,*denm;
   double *denp_ij,*denp_kl;
   double *denm_ij,*denm_kl;
   struct c_pkints *c_pkpt;

   scr1 = (double **) init_matrix(nbfso,nbfso);
   scr2 = (double **) init_matrix(nbfso,nbfso);
   denm = (double *) init_array(nbstri);

   denp = (double **) init_matrix(nbstri,nvec);
   denmt = (double **) init_matrix(nbstri,nvec);

   densmatd(denp,bb,nvec);

/* form half-transformed density matrix */

   srew(itap37);

   do {
      sread(itap37,(char *) c_pkbuf,sizeof(struct c_pkints)*maxbuf);
      num = (c_pkbuf[0].ij % (maxbuf+1));
      ilast = c_pkbuf[0].kl % (maxbuf+1);
      c_pkbuf[0].ij /= (maxbuf+1);
      c_pkbuf[0].kl /= (maxbuf+1);

      c_pkpt=c_pkbuf;
      for (i=num; i ; i--,c_pkpt++) {
         ij = (*c_pkpt).ij;
         kl = (*c_pkpt).kl;
         pval = (*c_pkpt).pval;
         denm_ij=denmt[ij];
         denm_kl=denmt[kl];
         denp_ij=denp[ij];
         denp_kl=denp[kl];
         for(iabc=nvec; iabc ; iabc--,denm_ij++,denm_kl++,denp_ij++,denp_kl++) {
            *denm_ij += *denp_kl*pval;
            *denm_kl += *denp_ij*pval;
            }
         }
      } while(!ilast);

   for(iabc=0; iabc < nvec ; iabc++) {
      for(ii=0; ii < nbstri ; ii++) denm[ii]=denmt[ii][iabc];
      ao_to_mo(denm,scr1,e_vecs_so,scr2,nbfso,nbfso);

      for(ii=0; ii < nind ; ii++) {
         i=indep[ii].ii;
         j=indep[ii].jj;
         ij=indep[ii].ij;
         eji = e_vals[j]-e_vals[i];
         bb1[ii][iabc] = eji*bb[ii][iabc]-denm[ij];
         }
      }

   if(print & 64) {
      fprintf(outfile,"\nb1 matrix\n");
      print_mat(bb1,nind,nvec,outfile);
      }


   free(denm);
   free_matrix(scr1,nbfso);
   free_matrix(scr2,nbfso);
   free_matrix(denp,nbstri);
   free_matrix(denmt,nbstri);
   }
