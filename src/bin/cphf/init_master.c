
/* $Log$
 * Revision 1.1  2000/02/04 22:50:48  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:50  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

init_master()
   {
      int i,j,ij;
      int ii,jj;
      int ntri2,ntril,nbasql,n3trl,n3qrl;
      int ipara[1024];
      double rpara[1024];

      nsect=1024;
      maxbuf=8192;

      itap40=40;
      rfile(itap40);

      block_locs = (int *) init_array(1024/2);
      rread(itap40,(char *) block_locs,sizeof(int)*nsect,1);

      mread(ipara,1);
      mread(rpara,2);

      nbstri = ipara[3];
      nbfao = ipara[10];
      nbfso = ipara[11];
      nc = ipara[6];
      no = ipara[7];
      nocc = ipara[8];
      iopen = ipara[5];
      natom = ipara[9];
      natom3 = ipara[16];
      natom3x = (dipole) ? natom3+3 : natom3;
      n_so_typs = ipara[12];
      ntypes = ipara[17];
      ntype1 = ntypes+1;

      nbasis = nbfso;
      ntri = nbstri;
      nbatri = nbfao*(nbfao+1)/2;
      ntril = (2*ntri-1)/nsect + 1;
      nbasql = (2*nbfso*nbfso-1)/nsect + 1;

      n3trl = natom3x*ntril;
      n3qrl = natom3x*nbasql;

      enuc = rpara[0];
      escf = rpara[1];

      if(iopen) 
         o_pkbuf = (struct o_pkints *) malloc(sizeof(struct o_pkints)*maxbuf);
      else
         c_pkbuf = (struct c_pkints *) malloc(sizeof(struct c_pkints)*maxbuf);

      if(twocon) {
         ha11 = (double *) init_array(natom3);
         ha12 = (double *) init_array(natom3);
         ha22 = (double *) init_array(natom3);
         e1t = (double *) init_array(natom3);
         baf = (struct off_diags *) malloc(sizeof(struct off_diags)*natom3x);
         }

/* calculate some pointers */

      for (i=0; i < natom3x ; i++) {
         sa_loc[i] = i*ntril+1;
         ha_loc[i] = i*ntril+n3trl+1;
         fa_loc[i] = i*ntril+n3trl+n3trl+1;
         ea_loc[i] = i*nbasql+n3trl+n3trl+n3trl+1;
         ba_loc[i] = i*ntril+n3trl+n3trl+n3trl+n3qrl+1;
         ua_loc[i] = i*nbasql+n3trl+n3trl+n3trl+n3qrl+n3trl+1;
         }

      for(i=0; i < 10 ; i++) {
         focc[i] = rpara[i+40];
         nsorb[i] = ipara[40+i];
         }

      nstart[0] = 0;
      nend[0] = nsorb[0];
      for(i=1; i <= ntypes ; i++) {
         nstart[i] = nend[i-1];
         nend[i] = nstart[i]+nsorb[i];
         }

      for(i=0; i <= ntypes ; i++)
         for(j=nstart[i]; j < nend[i] ; j++)
            motyp[j] = i;

      for(i=1,ij=0; i <= ntypes ; i++)
         for(j=0; j < i ; j++)
            for(ii=nstart[i]; ii < nend[i] ; ii++)
               for(jj=nstart[j]; jj < nend[j] ; jj++,ij++) 
      nind = ij+1;

      indep = (struct ind_pairs *) malloc(sizeof(struct ind_pairs)*nind);
      for(i=1,ij=0; i <= ntypes ; i++)
         for(j=0; j < i ; j++)
            for(ii=nstart[i]; ii < nend[i] ; ii++)
               for(jj=nstart[j]; jj < nend[j] ; jj++,ij++) {
                  indep[ij].ii = ii;
                  indep[ij].jj = jj;
                  indep[ij].ij = ioff[ii]+jj;
                  indep[ij].it = motyp[ii];
                  indep[ij].jt = motyp[jj];
                  }

      for(ii=ij=0; ii <= ntypes ; ii++)
         for(i=nstart[ii]; i < nend[ii] ; i++)
            for(j=nstart[ii]; j <= i ; j++,ij++)
      ndep = ij+1;

      if(print & 2) {
         fprintf(outfile,"\nnsorb nstart nend \n");
         for(i=0; i < 10 ; i++)
            fprintf(outfile,"%5d %5d %5d\n",nsorb[i],nstart[i],nend[i]);

         fprintf(outfile,"\nmotyp \n");
         for(i=0; i < nbasis ; i++)
            fprintf(outfile,"%5d \n",motyp[i]);
         
         fprintf(outfile,"\nnij1 nij2 kij1 kij2\n");
         for(i=0; i < nind ; i++)
            fprintf(outfile,"%5d %5d %5d %5d %5d\n",
               indep[i].ii,indep[i].jj,indep[i].ij,indep[i].it,indep[i].jt);
         }

/* calculate alpa and beta matrices */

      if(iopen) {
         alpa = (double **) init_matrix(ntype1,ntype1);
         beta = (double **) init_matrix(ntype1,ntype1);

         for(i=ij=0; i < ntypes ; i++)
            for(j=0; j <= i ; j++) {
               alpa[j][i]=alpa[i][j] = (1.0-rpara[10+ij])*focc[i]*focc[j]*0.5;
               beta[j][i]=beta[i][j] = -(1.0-rpara[25+ij])*focc[i]*focc[j]*0.25;
               ij++;
               }
         }
      }
