
/* $Log$
 * Revision 1.1  2000/02/04 22:50:46  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:34  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

/* calculate ba0 for independent pairs (eqn 13) */

bamat_o(max_words)
   int max_words;
{
   int i,j,k,l,ij;
   int ii,it,jt,kt,kl,il,jl;
   int iabc,lim,p1,p2,num,ilast;
   int first,last,num_vecs,vec;
   int core_left,core_needed,num_passes;
   double pval,qval;
   double val1,val2,valt,fac2,vala,valb,f;
   double esik,esjk;
   double *esi,*esj;
   double **zeta,**epa,*b0;
   double *zetait,*zetajt;
   double ***denp,***denq,**denm,***denmt,**sm;
   double *temp,**scr1,**scr2;
   double *dp_ij_it,*dp_kl_it;
   double *dq_ij_it,*dq_kl_it;
   double *dm_ij_it,*dm_kl_it;
   double *sm_iabc;
   struct o_pkints *o_pkpt;

   zeta = (double **) init_matrix(ntype1,nbstri);
   epa = (double **) init_matrix(nbfso,nbfso);
   temp = (double *) init_array(nbfso*nbfso);
   b0 = (double *) init_array(nbstri);
   scr1 = (double **) init_matrix(nbfso,nbfso);
   scr2 = (double **) init_matrix(nbfso,nbfso);
   denm = (double **) init_matrix(ntype1,nbstri);

   core_left = max_words-2*nbfso*nind-nbstri*(2*ntype1+1)-6*nbfso*nbfso;
   core_needed = nbstri*(3*ntypes+1)*natom3;
   if(core_needed < core_left) num_vecs=natom3;
   else {
      num_passes = core_needed/core_left+1;
      num_vecs = natom3/num_passes+1;
      }
   if(num_vecs <= 0) num_vecs=1;

   num_vecs = MIN0(num_vecs,natom3);
   fprintf(outfile,"\n\tallocating %d vectors with %d words in bamat_o\n",
          num_vecs,(3*ntypes+1)*nbstri*num_vecs);
   fflush(outfile);

   sm = (double **) init_matrix(num_vecs,nbstri);
   denp = (double ***) malloc(sizeof(double **)*nbstri);
   denq = (double ***) malloc(sizeof(double **)*nbstri);
   denmt = (double ***) malloc(sizeof(double **)*nbstri);
   for(i=0; i < nbstri ; i++) {
      denp[i] = (double **) init_matrix(ntypes,num_vecs);
      denq[i] = (double **) init_matrix(ntypes,num_vecs);
      denmt[i] = (double **) init_matrix(ntypes,num_vecs);
      }

   for(i=0; i < ntypes ; i++) {
      mread(zeta[i],26+i);
      if(print & 16) {
         fprintf(outfile,"\nzeta i = %5d\n",i);
         print_array(zeta[i],nbfso,outfile);
         }
      }

   first = -num_vecs;
   last=0;
   do {
      first += num_vecs;
      last += num_vecs;
      last = MIN0(last,natom3);

   /* eps_ij-eps_ji + sum sa_kl*(del_ki(e_jl-zet_jl)-del_kj(e_il-zet_il)) */

      for(vec=first,iabc=0; vec < last ; vec++,iabc++) {
         rread(itap44,(char *) sm[iabc],sizeof(double)*nbstri,sa_loc[vec]);
         rread(itap44,(char *) temp,sizeof(double)*nbfso*nbfso,ea_loc[vec]);
         sm_iabc=sm[iabc];

         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++)
               epa[j][i] = temp[ij];

         bzero(b0,sizeof(double)*nbstri);
         
         for(ii=0; ii < nind ; ii++) {
            i=indep[ii].ii;
            j=indep[ii].jj;
            ij=indep[ii].ij;
            it=indep[ii].it;
            jt=indep[ii].jt;
            zetait=zeta[it];
            zetajt=zeta[jt];

            if(i!=j) {
               val1 = epa[i][j]-epa[j][i];
               val2 = 0.0;
   
               lim = MIN0(nocc-1,j-1);
               kl = ioff[j];
               for (l=0; l <= lim ; l++,kl++) {
                  fac2 = 0.0;
                  p1 = MAX0(i,l);
                  p2 = MIN0(i,l);
                  il = ioff[p1]+p2;
                  fac2 += zetait[il]-zetajt[il];
                  if(fac2) val2 -= sm_iabc[kl]*fac2;
                  }
               lim = MIN0(nocc-1,i-1);
               kl = ioff[i];
               for (l=0; l <= lim ; l++,kl++) {
                  fac2 = 0.0;
                  p1 = MAX0(j,l);
                  p2 = MIN0(j,l);
                  jl = ioff[p1]+p2;
                  fac2 += zetait[jl]-zetajt[jl];
                  if(fac2) val2 -= sm_iabc[kl]*fac2;
                  }
               valt = val1+val2;
               }
            else {
               val2 = 0.0;

               lim = MIN0(nocc-1,i-1);
               kl = ioff[i];
               for (l=0; l <= lim ; l++,kl++) {
                  fac2 = 0.0;
                  p1 = MAX0(j,l);
                  p2 = MIN0(j,l);
                  jl = ioff[p1]+p2;
                  fac2 += zetait[jl]-zetajt[jl];
                  if(fac2) val2 -= 2.0*sm_iabc[kl]*fac2;
                  }
               valt = val2;
               }
            b0[ij] = valt;
            }

         if(print & 16) {
            fprintf(outfile,"\nb0 in bamat iabc = %5d\n",iabc);
            print_array(b0,nbfso,outfile);
            }
         rwrit(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[vec]);
         }

      for(i=0; i < nbstri ; i++) {
         zero_mat(denp[i],ntypes,num_vecs);
         zero_mat(denq[i],ntypes,num_vecs);
         zero_mat(denmt[i],ntypes,num_vecs);
         }

/* form density-like matrix in so basis  (eqn A-2) */

      for(i=ij=0; i < nbfso ; i++) {
         esi=e_vecs_so[i];
         for(j=0; j <= i ; j++,ij++) {
            esj = e_vecs_so[j];
            for(k=0; k < nocc ; k++) {
               esik = esi[k];
               esjk = esj[k];
               kt = motyp[k];
               for(l=0; l <= k ; l++) {
                  f = esik*esj[l]+esi[l]*esjk;
                  if(f) {
                     kl = ioff[k]+l;
                     if(k==l) f/=2.0;
                     for(it=0; it < ntypes ; it++) {
                        vala=f*alpa[it][kt];
                        valb=f*beta[it][kt];
                        dp_ij_it=denp[ij][it];
                        dq_ij_it=denq[ij][it];
                        for(vec=first,iabc=0; vec < last ;
                                       iabc++,vec++,dp_ij_it++,dq_ij_it++) {
                           *dp_ij_it -= sm[iabc][kl]*vala;
                           *dq_ij_it -= sm[iabc][kl]*valb;
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


/* for half-transformed density matrix  (eqn A-3) */

      srew(itap37);

      do {
         sread(itap37,(char *) o_pkbuf,sizeof(struct o_pkints)*maxbuf);
         num = (o_pkbuf[0].ij % (maxbuf+1));
         ilast = o_pkbuf[0].kl % (maxbuf+1);
         o_pkbuf[0].ij /= (maxbuf+1);
         o_pkbuf[0].kl /= (maxbuf+1);
         o_pkpt = o_pkbuf;

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

   /* transform to mo basis (eqn A-4) */

      for(vec=first,iabc=0; vec < last ; vec++,iabc++) {
         rread(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[vec]);
         for(it=0; it < ntypes ; it++) {
            for(ii=0; ii < nbstri ; ii++)
               denm[it][ii]=denmt[ii][it][iabc];
               ao_to_mo(denm[it],scr1,e_vecs_so,scr2,nbfso,nbfso);
            }
 
   /* and add in last term of eqn 13 */

         for(ii=0; ii < nind ; ii++) {
            it=indep[ii].it;
            jt=indep[ii].jt;
            ij = indep[ii].ij;
            b0[ij] += denm[it][ij]-denm[jt][ij];
            }
         rwrit(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[vec]);

         if(print & 16) {
            fprintf(outfile,"\nb0 in bamat iabc = %5d\n",iabc);
            print_array(b0,nbfso,outfile);
            }
         }
      } while (last != natom3);

   free_matrix(zeta,ntype1);
   free_matrix(epa,nbfso);
   free(temp);
   free(b0);
   free_matrix(scr1,nbfso);
   free_matrix(scr2,nbfso);
   free_matrix(denm,ntype1);
   free_matrix(sm,num_vecs);

   for(i=0; i < nbstri ; i++) {
      free_matrix(denp[i],ntypes);
      free_matrix(denq[i],ntypes);
      free_matrix(denmt[i],ntypes);
      }
   free(denp);
   free(denq);
   free(denmt);
   }
