
/* $Log$
 * Revision 1.1  2000/02/04 22:51:33  evaleev
 * Initial revision
 *
/* Revision 1.1  1991/06/15 22:06:29  seidl
/* Initial revision
/* */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

/* calculate density matrices and one-electron part of zeta */

density()

   {
      int i,j,ij,k,kl;
      int ityp,it,itk;
      int num,ilast,nn,nocc;
      double f,vala,valb;
      double *hmat;

      hmat = (double *) init_array(nbstri);
      densa = (double **) init_matrix(ntypes,nbstri);
      densb = (double **) init_matrix(ntypes,nbstri);
      zeta_so = (double **) init_matrix(ntypes,nbstri);

      mread(hmat,14);

      for(i=nocc=0; i < 10 ; i++) nocc += nclosd[i]+nopen[i];

      for (ityp=it=0; ityp < ntypes ; ityp++) {
         if(ityp) it += nsorb[ityp-1];
         f = 0.5*occ_num[it];
         for(i=ij=0; i < nbfso ; i++) {
            for(j=0; j <= i ; j++,ij++) {
               vala=valb=0.0;
               for(k=0; k < nocc ; k++) {
                  itk = (it < k) ? ioff[k]+it : ioff[it]+k;
                  vala += alpb[itk]*e_vecs_so[i][k]*e_vecs_so[j][k];
                  valb += betb[itk]*e_vecs_so[i][k]*e_vecs_so[j][k];
                  }
               densa[ityp][ij]=vala;
               densb[ityp][ij]=valb;
               }
            }
         for(i=0; i < nbstri ; i++)
            zeta_so[ityp][i]=hmat[i]*f;

         if(print > 2) {
            fprintf(outfile,"\ndensa ityp = %d\n",ityp);
            print_array(densa[ityp],nbfso,outfile);

            fprintf(outfile,"\ndensb ityp = %d\n",ityp);
            print_array(densb[ityp],nbfso,outfile);
            }
         if(print > 3) {
            fprintf(outfile,"\nzeta_so ityp = %d\n",ityp);
            print_array(zeta_so[ityp],nbfso,outfile);
            }
         }
      }
