
/* $Log$
 * Revision 1.1  2000/02/04 22:50:50  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:53  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

polar()

   {
      int ix,jx,i,j,ij,p1,p2;
      double bohr = 0.52917706;
      double a33 = bohr*bohr*bohr;
      double vp,f,c1g,c2g;
      double **hf,**uf,*eig,**eiv,*a;
      double **pol,**polc,**poc;
      double *temp;

      hf = (double **) init_matrix(3,nbstri);
      uf = (double **) init_matrix(nbfso,nbfso);
      eiv = (double **) init_matrix(3,3);
      eig = (double *) init_array(3);
      a = (double *) init_array(6);
      pol = (double **) init_matrix(3,3);
      poc = (double **) init_matrix(3,3);
      polc = (double **) init_matrix(3,3);
      temp = (double *) init_array(nbfso*nbfso);

      for(ix=0; ix < 3 ; ix++)
         rread(itap44,(char *) hf[ix],sizeof(double)*nbstri,ha_loc[ix+natom3]);

      for(ix=0; ix < 3 ; ix++) {
         if(twocon) {
            c1g=c1a[natom3+ix];
            c2g=c2a[natom3+ix];
            }
         rread(itap44,(char *) temp,sizeof(double)*nbfso*nbfso,ua_loc[ix+natom3]);
         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++)
               uf[j][i] = temp[ij];

         for(jx=0; jx < 3 ; jx++) {
            vp=0.0;
            for(i=0; i < nocc ; i++) {
               f=occ_num[i]*0.5;
               for(j=0; j < nbfso ; j++) {
                  p1=MAX0(i,j);
                  p2=MIN0(i,j);
                  ij=ioff[p1]+p2;
                  vp -= uf[j][i]*hf[jx][ij]*f;
                  }
               }
            pol[ix][jx]=vp*4.0;
            if(twocon) {
               poc[ix][jx] = -2.0*(c1g*c1*h11f[jx]+c2g*c2*h22f[jx]);
               pol[ix][jx] += poc[ix][jx];
               }
            }
         }

      if(twocon) {
         fprintf(outfile,"\n\n      poc matrix (a.u.)\n");
         print_mat(pol,3,3,outfile);
         }
      fprintf(outfile,"\n\n      polarizability tensor matrix (a.u.)\n");
      print_mat(pol,3,3,outfile);

      for(i=0; i < 3 ; i++)
         for(j=0; j < 3 ; j++)
            polc[i][j] = pol[i][j]*a33;

      fprintf(outfile,"\n      polarizability tensor matrix (a**3)\n");
      print_mat(polc,3,3,outfile);
               
      for(i=ij=0; i < 3 ; i++)
         for(j=0; j <= i ; j++,ij++)
            a[ij] = pol[i][j];

      rsp(3,3,6,a,eig,1,eiv,10.0e-15);

      zero_mat(pol,3,3);
      for(i=0; i < 3 ; i++) pol[i][i]=eig[i];

      fprintf(outfile,"\n principal polarizability tensor matrix (a.u.)\n");
      print_mat(pol,3,3,outfile);

      for(i=0; i < 3 ; i++)
         for(j=0; j < 3 ; j++)
            polc[i][j] = pol[i][j]*a33;

      fprintf(outfile,"\n principal polarizability tensor matrix (a**3)\n");
      print_mat(polc,3,3,outfile);
         
      if(print & 1) {
         fprintf(outfile,"\n eigenvector matrix\n");
         print_mat(eiv,3,3,outfile);
         }
      }
