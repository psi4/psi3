
/* $Log$
 * Revision 1.1  2000/02/04 22:50:49  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:52  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

/* form A*B' in eqn 20 */

makeab_o(bb,bb1,nvec,zeta)
   double **bb,**bb1,**zeta;
   int nvec;

   {
      register unsigned i;
      int j,k,l,ij;
      int iabc;
      int ii,jj,it,jt,kl,il,jl;
      int ik,jk;
      int p1,p2;
      int num,ilast;
      double pval,qval;
      double z1,z2,z3,z4,zt;
      double **scr1,**scr2;
      double ***denp,***denq,***denmt,**denm;
      double *dm_ij_it,*dm_kl_it;
      double *dp_ij_it,*dp_kl_it;
      double *dq_ij_it,*dq_kl_it;
      double *bt,*btt;
      double *zeta_it,*zeta_jt;
      struct o_pkints *o_pkpt;

      scr1 = (double **) init_matrix(nbfso,nbfso);
      scr2 = (double **) init_matrix(nbfso,nbfso);
      denm = (double **) init_matrix(ntype1,nbstri);

      denp = (double ***) malloc(sizeof(double **)*nbstri);
      denq = (double ***) malloc(sizeof(double **)*nbstri);
      denmt = (double ***) malloc(sizeof(double **)*nbstri);

      for(i=0; i < nbstri ; i++) {
          denp[i] = (double **) init_matrix(ntypes,nvec);
          denq[i] = (double **) init_matrix(ntypes,nvec);
          denmt[i] = (double **) init_matrix(ntypes,nvec);
          }

 /* eqns 22 & 23 */

#if TIME
      fprintf(outfile,"\nbefore denmat\n");
      resource_command();
#endif
      densmatd_o(denp,denq,bb,nvec);
#if TIME
      fprintf(outfile,"\nafter denmat\n");
      resource_command();
#endif

 /* form half-transformed density matrix as in appendix */

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
               for(iabc=nvec; iabc ; iabc--,dm_ij_it++,dm_kl_it++,
                                            dp_ij_it++,dp_kl_it++,
                                            dq_ij_it++,dq_kl_it++) {
                  *dm_ij_it += *dp_kl_it*pval + *dq_kl_it*qval;
                  *dm_kl_it += *dp_ij_it*pval + *dq_ij_it*qval;
                  }
               }
            }
         } while(!ilast);

#if TIME
      fprintf(outfile,"\nafter supermatrix\n");
      resource_command();
#endif

  /* transform to mo basis */

      if(twocon) zero_mat(bb1,nind+2,nvec);

      for(iabc=0; iabc < nvec ; iabc++) {
         for(it=0; it < ntypes ; it++) {
            for(ii=0; ii < nbstri ; ii++)
               denm[it][ii]=denmt[ii][it][iabc];
            ao_to_mo(denm[it],scr1,e_vecs_so,scr2,nbfso,nbfso);
            }

    /*  eqn 21 */

         for(ii=0; ii < nind ; ii++) {
            it=indep[ii].it;
            jt=indep[ii].jt;
            ij=indep[ii].ij;
            bb1[ii][iabc] = denm[it][ij]-denm[jt][ij];
            }
         }

#if TIME
      fprintf(outfile,"\nafter moconv\n");
      resource_command();
#endif

      if(print & 64) {
         int p2 = (twocon) ? nind+2 : nind;
         fprintf(outfile,"\nb1 matrix\n");
         print_mat(bb1,p2,nvec,outfile);
         }

 /* final terms in eqn 12 */

   for(ii=0; ii < nind ; ii++) {
      i=indep[ii].ii;
      j=indep[ii].jj;
      it=indep[ii].it;
      jt=indep[ii].jt;
      zeta_it=zeta[it];
      zeta_jt=zeta[jt];
      for(jj=0; jj < nind ; jj++) {
         k=indep[jj].ii;
         l=indep[jj].jj;
         z1=z2=z3=z4=0.0;

         if(j == k) {
            p1 = MAX0(i,l);
            p2 = MIN0(i,l);
            il = ioff[p1]+p2;
            z1 = zeta_it[il]-zeta_jt[il];
            }
         if(i == k) {
            p1 = MAX0(j,l);
            p2 = MIN0(j,l);
            jl = ioff[p1]+p2;
            z2 = zeta_it[jl]-zeta_jt[jl];
            }
         if(j == l) {
            p1 = MAX0(i,k);
            p2 = MIN0(i,k);
            ik = ioff[p1]+p2;
            z3 = zeta_jt[ik]-zeta_it[ik];
            }
         if(i == l) {
            p1 = MAX0(j,k);
            p2 = MIN0(j,k);
            jk = ioff[p1]+p2;
            z4 = zeta_jt[jk]-zeta_it[jk];
            }
         zt = z1+z2+z3+z4;

         bt = bb1[ii];
         btt = bb[jj];
         for(iabc=nvec; iabc ; iabc--,bt++,btt++)
            *bt += *btt*zt;
         }
      if(twocon) {
         for(iabc=0; iabc < nvec ; iabc++) {
            bb1[ii][iabc] += bb[nind][iabc] * a12[ii].one;
            bb1[ii][iabc] += bb[nind+1][iabc]*a12[ii].two;
            }
         }
      }

   if(twocon) {
      for(jj=0; jj < nind ; jj++) {
         for(iabc=0; iabc < nvec ; iabc++) {
            bb1[nind][iabc]   += bb[jj][iabc]*a12[jj].one;
            bb1[nind+1][iabc] += bb[jj][iabc]*a12[jj].two;
            }
         }
      for(ii=0; ii < 2 ; ii++) {
         for(jj=0; jj < 2 ; jj++) {
            for(iabc=0; iabc < nvec ; iabc++) {
               bb1[nind+ii][iabc] += bb[nind+jj][iabc]*a22[ii][jj];
               }
            }
         }
      }

   if(print & 64) {
      int nindp2 = (twocon) ? nind+2 : nind;
      fprintf(outfile,"\nfinal b1 matrix\n");
      print_mat(bb1,nindp2,nvec,outfile);
      fprintf(outfile,"\n%f %f %f %f\n",a22[0][0],a22[0][1],a22[1][0],a22[1][1]);
      }

#if TIME
      fprintf(outfile,"\nend of makeab\n");
      resource_command();
#endif

      free_matrix(scr1,nbfso);
      free_matrix(scr2,nbfso);
      free_matrix(denm,ntype1);
      for(i=0; i < nbstri ; i++) {
         free_matrix(denp[i],ntypes);
         free_matrix(denq[i],ntypes);
         free_matrix(denmt[i],ntypes);
         }
      free(denp);
      free(denq);
      free(denmt);
      }
