
/* $Log$
 * Revision 1.1  2000/02/04 22:50:50  evaleev
 * Initial revision
 *
/* Revision 1.3  1997/08/25 21:53:59  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.2  1995/01/19  19:48:42  seidl
 * replace some nulls with spaces
 *
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

wamat_o(max_words)
   int max_words;
{
   int i,j,ij,k,l,kl;
   int ii,it,jt,ik,jk,lt;
   int iabc,p1,p2;
   int num,ilast;
   int first,last,vec,num_vecs;
   int core_left,core_needed,num_passes;
   double f1,pval,qval;
   double val1,val2,val3,f,f2,vala,valb;
   double **zeta,**epa,**wa,*b0;
   double ***denp,***denq,**denm,***denmt,**sm;
   double *dp_ij_it,*dp_kl_it;
   double *dq_ij_it,*dq_kl_it;
   double *dm_ij_it,*dm_kl_it;
   double *esi,*esj,esik,esjk;
   double *temp,**scr1,**scr2;
   struct o_pkints *o_pkpt;

   b0 = (double *) init_array(nbstri);
   zeta = (double **) init_matrix(ntype1,nbstri);
   denm = (double **) init_matrix(ntype1,nbstri);
   epa = (double **) init_matrix(nbfso,nbfso);
   wa = (double **) init_matrix(nbfso,nbfso);
   temp = (double *) init_array(nbfso*nbfso);
   scr1 = (double **) init_matrix(nbfso,nbfso);
   scr2 = (double **) init_matrix(nbfso,nbfso);

   core_left = max_words-2*nbfso*nind-5*nbfso*nbfso-nbstri-2*ntype1*nbstri;
   core_needed = nbstri*(3*ntypes+1)*natom3;
   if(core_needed < core_left) num_vecs=natom3;
   else {
      num_passes = core_needed/core_left+1;
      num_vecs = natom3/num_passes+1;
      }
   if(num_vecs <= 0) num_vecs=1;
   num_vecs = MIN0(num_vecs,natom3);

   fprintf(outfile,"\n%8callocating %d vectors with %d words in wamat_o\n",
          ' ',num_vecs,nbstri*(3*ntypes+1)*num_vecs);
   fflush(outfile);

   denp = (double ***) malloc(sizeof(double **)*nbstri);
   denq = (double ***) malloc(sizeof(double **)*nbstri);
   denmt = (double ***) malloc(sizeof(double **)*nbstri);
   sm = (double **) init_matrix(num_vecs,nbstri);

   for(i=0; i < nbstri ; i++) {
      denp[i] = (double **) init_matrix(ntypes,num_vecs);
      denq[i] = (double **) init_matrix(ntypes,num_vecs);
      denmt[i] = (double **) init_matrix(ntypes,num_vecs);
      }

   for(i=0; i < ntypes ; i++)
      mread(zeta[i],26+i);

   srew(work);

   first = -num_vecs;
   last=0;
   do {
      first += num_vecs;
      last += num_vecs;
      last = MIN0(last,natom3);

      for(i=0; i < nbstri ; i++) {
         zero_mat(denp[i],ntypes,num_vecs);
         zero_mat(denq[i],ntypes,num_vecs);
         zero_mat(denmt[i],ntypes,num_vecs);
         }

      for(vec=first,iabc=0; vec < last ; vec++,iabc++)
         rread(itap44,(char *) sm[iabc],sizeof(double)*nbstri,sa_loc[vec]);

/* form density-like matrix in so basis as in appendix */

      for(i=ij=0; i < nbfso ; i++) {
         esi=e_vecs_so[i];
         for(j=0; j <= i ; j++,ij++) {
            esj = e_vecs_so[j];
            for(k=0; k < nocc ; k++) {
               esik = esi[k];
               esjk = esj[k];
               for(l=0; l < nocc ; l++) {
                  f = esik*esj[l]+esi[l]*esjk;
                  if(f) {
                     p1=MAX0(k,l);
                     p2=MIN0(k,l);
                     kl = ioff[p1]+p2;
                     lt = motyp[l];
                     for(it=0; it < ntypes ; it++) {
                        vala=f*alpa[it][lt];
                        valb=f*beta[it][lt];
                        dp_ij_it=denp[ij][it];
                        dq_ij_it=denq[ij][it];
                        for(vec=first,iabc=0; vec < last ;
                                       iabc++,vec++,dp_ij_it++,dq_ij_it++) {
                           *dp_ij_it += sm[iabc][kl]*vala;
                           *dq_ij_it += sm[iabc][kl]*valb;
                           }
                        }
                     }
                  }
               }
            }
         }

      for(i=ij=0; i < nbfso ; i++)
         for(j=0; j <= i ; j++,ij++)
            for(it=0; it < ntypes ; it++)
               for(vec=first,iabc=0; vec < last ; vec++,iabc++) {
                  denp[ij][it][iabc] = (i==j) ? denp[ij][it][iabc] :
                                                2.0*denp[ij][it][iabc];
                  denq[ij][it][iabc] = (i==j) ? denq[ij][it][iabc] :
                                                2.0*denq[ij][it][iabc];
                  }

/* for half-transformed density matrix as in appendix */

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

            for(it=0; it < ntypes ; it++) {
               dm_ij_it=denmt[ij][it];
               dm_kl_it=denmt[kl][it];
               dp_ij_it=denp[ij][it];
               dp_kl_it=denp[kl][it];
               dq_ij_it=denq[ij][it];
               dq_kl_it=denq[kl][it];
               for(vec=first;vec < last;vec++,dm_ij_it++,dm_kl_it++,
                                              dp_ij_it++,dp_kl_it++,
                                              dq_ij_it++,dq_kl_it++) {
                  *dm_ij_it += *dp_kl_it*pval + *dq_kl_it*qval;
                  *dm_kl_it += *dp_ij_it*pval + *dq_ij_it*qval;
                  }
               }
            }
         } while(!ilast);

  /* transform to mo basis */

      for(vec=first,iabc=0; vec < last ; vec++,iabc++) {
         for(it=0; it < ntypes ; it++) {
            for(ii=0; ii < nbstri ; ii++)
               denm[it][ii]=denmt[ii][it][iabc];
            ao_to_mo(denm[it],scr1,e_vecs_so,scr2,nbfso,nbfso);
            }

         rread(itap44,(char *) sm[iabc],sizeof(double)*nbstri,sa_loc[vec]);
         rread(itap44,(char *) temp,sizeof(double)*nbfso*nbfso,ea_loc[vec]);
         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++)
               epa[j][i]=temp[ij];

       /* now construct wa as in eqn 6 */

         for(i=ij=0; i < nbfso ; i++) {
            it=motyp[i];
            for(j=0; j < nbfso ; j++) {
               jt=motyp[j];
               p1=MAX0(i,j);
               p2=MIN0(i,j);
               ij = ioff[p1]+p2;
               val1=2*epa[j][i];
               val2=denm[jt][ij];
               val3=0.0;
               for(k=0; k < nocc ; k++) {
                  p1=MAX0(i,k);
                  p2=MIN0(i,k);
                  ik=ioff[p1]+p2;
                  p1=MAX0(j,k);
                  p2=MIN0(j,k);
                  jk=ioff[p1]+p2;
                  f1=zeta[it][ik]+zeta[jt][ik];
                  f2=2.0*zeta[jt][jk];
                  val3 += sm[iabc][jk]*f1+sm[iabc][ik]*f2;
                  }
               wa[i][j]=val1-val2-val3;
               }
            }
         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++)
               temp[ij]=wa[j][i];
         swrit(work,(char *) temp,sizeof(double)*nbfso*nbfso);

         if(print & 256) {
            fprintf(outfile,"\nwa iabc = %5d\n",vec);
            print_mat(wa,nbfso,nbfso,outfile);
            }
         }
      } while(last != natom3);

      free_matrix(zeta,ntype1);
      free_matrix(denm,ntype1);
      free_matrix(epa,nbfso);
      free_matrix(wa,nbfso);
      free_matrix(scr1,nbfso);
      free_matrix(scr2,nbfso);
      free_matrix(sm,num_vecs);
      free(b0);
      free(temp);
      for(i=0; i < nbstri ; i++) {
         free_matrix(denp[i],ntypes);
         free_matrix(denq[i],ntypes);
         free_matrix(denmt[i],ntypes);
         }
      free(denp);
      free(denq);
      free(denmt);
      }
