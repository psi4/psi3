
/* $Log$
 * Revision 1.1  2000/02/04 22:50:50  evaleev
 * Initial revision
 *
/* Revision 1.3  1997/08/25 21:53:58  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.2  1995/01/19  19:48:41  seidl
 * replace some nulls with spaces
 *
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

wamat_c(max_words)
   int max_words;
{
   int i,j,ij,k,l,kl;
   int iabc,p1,p2,ii;
   int num,ilast;
   int first,last,vec,num_vecs;
   int core_left,core_needed,num_passes;
   double pval,f;
   double **denp,**denmt,**sa,***ua;
   double *dp_ij,*dp_kl;
   double *dm_ij,*dm_kl;
   double *esi,*esj,esik,esjk;
   double *temp,**scr1,**scr2;
   struct c_pkints *c_pkpt;

   temp = (double *) init_array(nbfso*nbfso);
   scr1 = (double **) init_matrix(nbfso,nbfso);
   scr2 = (double **) init_matrix(nbfso,nbfso);

   core_left = max_words-2*nbfso*nind-3*nbfso*nbfso;
   core_needed = 3*nbstri*natom3+nbfso*nbfso*natom3;
   if(core_needed < core_left) num_vecs=natom3;
   else {
      num_passes = core_needed/core_left+1;
      num_vecs = natom3/num_passes+1;
      }
   if(num_vecs <= 0) num_vecs=1;
   num_vecs = MIN0(num_vecs,natom3);

   fprintf(outfile,"\n%8callocating %d vectors with %d words in wamat\n",
          ' ',num_vecs,(nbfso*nbfso+nbstri*3)*num_vecs);
   fflush(outfile);

   ua = (double ***) malloc(sizeof(double **)*num_vecs);
   for(i=0; i < num_vecs ; i++) ua[i]=(double **) init_matrix(nbfso,nbfso);

   srew(work);

   first = -num_vecs;
   last=0;
   do {
      first += num_vecs;
      last += num_vecs;
      last = MIN0(last,natom3);

      sa = (double **) init_matrix(num_vecs,nbstri);

      for(vec=first,iabc=0; vec < last ; vec++,iabc++) {
         rread(itap44,(char *) sa[iabc],sizeof(double)*nbstri,sa_loc[vec]);
         rread(itap44,(char *) temp,sizeof(double)*nbfso*nbfso,ua_loc[vec]);
         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++) 
               ua[iabc][j][i]=temp[ij];
         }

/* form density-like matrix in so basis */

      denp = (double **) init_matrix(nbstri,num_vecs);

      for(i=ij=0; i < nbfso ; i++) {
         esi=e_vecs_so[i];
         for(j=0; j <= i ; j++,ij++) {
            esj=e_vecs_so[j];
            for(k=0; k < nocc ; k++) {
               esik=esi[k];
               for(l=0; l < nocc ; l++) {
                  f = esik*esj[l];
                  if(f) {
                     p1=MAX0(k,l);
                     p2=MIN0(k,l);
                     kl = ioff[p1]+p2;
                     for(vec=first,iabc=0; vec < last ; vec++,iabc++)
                        denp[ij][iabc] += sa[iabc][kl]*f;
                     }
                  }
               }
            for(k=nocc; k < nbfso ; k++) {
               esik=esi[k];
               esjk=esj[k];
               for(l=0; l < nocc ; l++) {
                  f=esik*esj[l]+esi[l]*esjk;
                  if(f) {
                     for(vec=first,iabc=0; vec < last ; vec++,iabc++)
                        denp[ij][iabc] -= ua[iabc][k][l]*f;
                     }
                  }
               }
            }
         }

      free_matrix(sa,num_vecs);

      for(i=ij=0; i < nbfso ; i++)
         for(j=0; j <= i ; j++,ij++)
            for(vec=first,iabc=0; vec < last ; vec++,iabc++)
               denp[ij][iabc] = (i==j) ? 2.0*denp[ij][iabc] : 4.0*denp[ij][iabc];
                                             
      if(print & 256) {
         fprintf(outfile,"\ndenp iabc = %5d\n",vec);
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

            dm_ij = denmt[ij];
            dm_kl = denmt[kl];
            dp_ij = denp[ij];
            dp_kl = denp[kl];
            for(vec=first;vec < last;vec++,dm_ij++,dm_kl++,
                                              dp_ij++,dp_kl++) {
               *dm_ij += *dp_kl*pval;
               *dm_kl += *dp_ij*pval;
               }
            }
         } while(!ilast);

      free_matrix(denp,nbstri);

      for(vec=first,iabc=0; vec < last ; vec++,iabc++) {
         for(ii=0; ii < nbstri ; ii++)
            temp[ii]=denmt[ii][iabc];
         ao_to_mo(temp,scr1,e_vecs_so,scr2,nbfso,nbfso);

         swrit(work,(char *) temp,sizeof(double)*nbstri);

         if(print & 256) {
            fprintf(outfile,"\nwa iabc = %5d\n",vec);
            print_array(temp,nbfso,outfile);
            }
         }

      free_matrix(denmt,nbstri);

      } while(last != natom3);

   free_matrix(scr1,nbfso);
   free_matrix(scr2,nbfso);
   free(temp);
   for(i=0; i < num_vecs ; i++) {
      free_matrix(ua[i],nbfso);
      }
   free(ua);
   }
