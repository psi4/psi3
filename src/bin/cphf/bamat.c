
/* $Log$
 * Revision 1.1  2000/02/04 22:50:46  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:33  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

bamat(max_words)
   int max_words;
{
   int i,j,k,l,ij,kl,ii;
   int iabc,p1,p2,num,ilast;
   int first,last,num_vecs,vec;
   int core_left,core_needed,num_passes;
   double pval,esi_k,f;
   double *b0,*esi,*esj;
   double **denp,*denm,**denmt,**sm;
   double *dp_ij,*dp_kl,*dm_ij,*dm_kl;
   double **scr1,**scr2;
   struct c_pkints *c_pkpt;

   b0 = (double *) init_array(nbstri);
   scr1 = (double **) init_matrix(nbfso,nbfso);
   scr2 = (double **) init_matrix(nbfso,nbfso);
   denm = (double *) init_array(nbstri);

   core_left = max_words-2*nbfso*nind-2*nbstri-4*nbfso*nbfso;
   core_needed = 3*nbstri*natom3;
   if(core_needed < core_left) num_vecs=natom3;
   else {
      num_passes = core_needed/core_left+1;
      num_vecs = natom3/num_passes+1;
      }

   if(num_vecs <= 0) num_vecs=1;

   num_vecs = MIN0(num_vecs,natom3);
   fprintf(outfile,"\n\tallocating %d vectors with %d words in bamat\n",
          num_vecs,3*nbstri*num_vecs);
   fflush(outfile);

   first= -num_vecs;
   last=0;
   do {
      first += num_vecs;
      last += num_vecs;
      last = MIN0(last,natom3);

      sm = (double **) init_matrix(num_vecs,nbstri);

      for(vec=first,iabc=0; vec < last ; vec++,iabc++) {
         rread(itap44,(char *) sm[iabc],sizeof(double)*nbstri,sa_loc[vec]);
         rread(itap44,(char *) b0,sizeof(double)*nbstri,fa_loc[vec]);

         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j <= i ; j++,ij++)
               b0[ij] -= e_vals[j]*sm[iabc][ij];

         if(print & 16) {
            fprintf(outfile,"\nb0 in bamat iabc = %5d\n",iabc);
            print_array(b0,nbfso,outfile);
            }
         rwrit(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[vec]);
         }

/* form density-like matrix in so basis */

      denp = (double **) init_matrix(nbstri,num_vecs);

      for(i=ij=0; i < nbfso ; i++) {
         esi = e_vecs_so[i];
         for(j=0; j <= i ; j++,ij++) {
            esj = e_vecs_so[j];
            for(k=0; k < nc ; k++) {
               esi_k = 2.0*esi[k];
               for(l=0; l < nc ; l++) {
                  p1=MAX0(k,l);
                  p2=MIN0(k,l);
                  kl=ioff[p1]+p2;
                  f = esi_k*esj[l];
                  if(f) 
                     for(vec=first,iabc=0; vec < last ; vec++,iabc++)
                        denp[ij][iabc] -= sm[iabc][kl]*f;
                  }
               }
            }
         }

      free_matrix(sm,num_vecs);

      for(i=ij=0; i < nbfso ; i++)
         for(j=0; j <= i ; j++,ij++)
            for(iabc=0; iabc < num_vecs ; iabc++)
               denp[ij][iabc] = (i == j) ? denp[ij][iabc] : 2.0*denp[ij][iabc];


      if(print & 16) {
         fprintf(outfile,"denp in bamat\n");
         print_mat(denp,nbstri,num_vecs,outfile);
         }

/* for half-transformed density matrix */

      denmt = (double **) init_matrix(nbstri,num_vecs);

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
            dp_ij=denp[ij];
            dp_kl=denp[kl];
            dm_ij=denmt[ij];
            dm_kl=denmt[kl];

            for(vec=first;vec < last;vec++,dm_ij++,dm_kl++,dp_ij++,dp_kl++) {
               *dm_ij += *dp_kl*pval;
               *dm_kl += *dp_ij*pval;
               }
            }
         } while(!ilast);

      free_matrix(denp,nbstri);

      if(print & 16) {
         fprintf(outfile,"denmt in bamat\n");
         print_mat(denmt,nbstri,num_vecs,outfile);
         }

      for(vec=first,iabc=0; vec < last ; vec++,iabc++) {
         rread(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[vec]);
         for(ii=0; ii < nbstri ; ii++) denm[ii]=denmt[ii][iabc];
         ao_to_mo(denm,scr1,e_vecs_so,scr2,nbfso,nbfso);
         for(ii=0; ii < nbstri ; ii++)
            b0[ii] += denm[ii];
         rwrit(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[vec]);

         if(print & 16) {
            fprintf(outfile,"\nb0 in bamat iabc = %5d\n",iabc);
            print_array(b0,nbfso,outfile);
            }
         }

      free_matrix(denmt,nbstri);

      } while(last != natom3);

   free_matrix(scr1,nbfso);
   free_matrix(scr2,nbfso);
   free(b0);
   free(denm);
   }
