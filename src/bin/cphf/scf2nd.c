
/* $Log$
 * Revision 1.1  2000/02/04 22:50:50  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:55  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

scf2nd()
   {
      int i,iabc,jabc,ij,p1,p2,j;
      double vp,vm,vq,vr;
      double v11,v22,v12;
      double f;
      double **sa,**fa,*temp,**e2p,**e2m,**e2a,**e2c;
      double **e11,**e22,**e12;
      double *ba,**ua;
      FILE *itap15;

      ffile(&itap15,"file15.dat",0);

      sa = (double **) init_matrix(natom3,nbstri);
      fa = (double **) init_matrix(natom3,nbstri);
      ba = (double *) init_array(nbstri);
      e2p = (double **) init_matrix(natom3,natom3);
      e2m = (double **) init_matrix(natom3,natom3);
      e2a = (double **) init_matrix(natom3,natom3);
      e2c = (double **) init_matrix(natom3,natom3);
      if(twocon) {
         e11 = (double **) init_matrix(nbfso,nbfso);
         e22 = (double **) init_matrix(nbfso,nbfso);
         e12 = (double **) init_matrix(nbfso,nbfso);
         }
      ua = (double **) init_matrix(nbfso,nbfso);
      temp = (double *) init_array(nbfso*nbfso*natom);

      srew(itap42);
      srew(work);

      sread(itap42,(char *) temp,sizeof(double)*natom3);
      sread(itap42,(char *) temp,sizeof(double)*natom3*natom3);
      for(i=ij=0; i < natom3 ; i++)
         for(j=0; j < natom3 ; j++,ij++)
            e2a[j][i] = temp[ij];

      for(iabc=0; iabc < natom3 ; iabc++) {
         rread(itap44,(char *) sa[iabc],sizeof(double)*nbstri,sa_loc[iabc]);
         rread(itap44,(char *) fa[iabc],sizeof(double)*nbstri,fa_loc[iabc]);
         }

      for(iabc=0; iabc < natom3 ; iabc++) {
         rread(itap44,(char *) temp,sizeof(double)*nbfso*nbfso,ua_loc[iabc]);
         sread(work,(char *) ba,sizeof(double)*nbstri);
         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++)
               ua[j][i]=temp[ij];

         for(jabc=0; jabc < natom3 ; jabc++) {
            vp=vm=vr=vq=0.0;
            for(i=0; i < nocc ; i++) {
               for(j=0; j < nocc ; j++) {
                  p1=MAX0(i,j);
                  p2=MIN0(i,j);
                  ij=ioff[p1]+p2;
                  vm += ba[ij]*sa[jabc][ij];
                  vp -= (sa[iabc][ij]*fa[jabc][ij]+sa[jabc][ij]*fa[iabc][ij]);
                  vq += (sa[iabc][ij]*sa[jabc][ij])*(e_vals[i]*3.0+e_vals[j]);
                  }
               for(j=nocc; j < nbfso ; j++) {
                  p1=MAX0(i,j);
                  p2=MIN0(i,j);
                  ij=ioff[p1]+p2;
                  f = fa[jabc][ij]-e_vals[i]*sa[jabc][ij];
                  vr += f*ua[j][i];
                  }
               }
            e2m[iabc][jabc] = vm*2.0+vq;
            e2p[iabc][jabc] = vp*2.0+4.0*vr;
            }
         }

      if(print & 2) {
         fprintf(outfile,"\n e2p matrix\n");
         print_mat(e2p,natom3,natom3,outfile);
         fprintf(outfile,"\n e2m matrix\n");
         print_mat(e2m,natom3,natom3,outfile);
         }

      for(i=0; i < natom3 ; i++)
         for(j=0; j < natom3 ; j++)
            e2m[i][j] += e2p[i][j];

 /* contribution from ci coefficient change */

      if(twocon) {
         srew(work2);
         sread(work2,(char *) temp,sizeof(double)*nbfso*nbfso);
         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++)
               e11[j][i]=temp[ij];
         sread(work2,(char *) temp,sizeof(double)*nbfso*nbfso);
         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++)
               e22[j][i]=temp[ij];
         sread(work2,(char *) temp,sizeof(double)*nbfso*nbfso);
         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++)
               e12[j][i]=temp[ij];

         for(iabc=0; iabc < natom3 ; iabc++) {
            v11=ha11[iabc];
            v22=ha22[iabc];
            v12=ha12[iabc];
            for(i=0; i < nocc ; i++) {
               for(j=0; j < nocc ; j++) {
                  ij=(i > j) ? ioff[i]+j : ioff[j]+i;
                  v11 -= sa[iabc][ij]*e11[i][j]*2.0;
                  v22 -= sa[iabc][ij]*e22[i][j]*2.0;
                  v12 -= sa[iabc][ij]*e12[i][j]*2.0;
                  }
               }
            for(jabc=0; jabc < natom3 ; jabc++) 
               e2c[jabc][iabc]= 2.0*(c1a[jabc]*c1*v11 + c1a[jabc]*c2*v12
                                   + c2a[jabc]*c1*v12 + c2a[jabc]*c2*v22);
            }
         if(print & 2) {
            fprintf(outfile,"\n e2c matrix\n");
            print_mat(e2c,natom3,natom3,outfile);
            }
         }

      if(print & 2) {
         fprintf(outfile,"\n\tsummary of first order cphf calculation\n");

         fprintf(outfile,"\n e2a matrix\n");
         print_mat(e2a,natom3,natom3,outfile);
         fprintf(outfile,"\n e2m matrix\n");
         print_mat(e2m,natom3,natom3,outfile);
         }

      for(i=0; i < natom3 ; i++)
         for(j=0; j < natom3 ; j++)
            e2m[i][j] += e2a[i][j]+e2c[i][j];

      if(print & 1) {
         fprintf(outfile,"\n scf second derivative matrix\n");
         print_mat(e2m,natom3,natom3,outfile);
         }

      fprintf(itap15,"%5d%5d\n",natom,natom3*2);
      for(i=ij=0; i < natom3 ; i++)
         for(j=0; j < natom3 ; j++,ij++)
            temp[ij] = e2m[j][i];
      for(ij=0; ij < natom3*natom3 ; ij += 3)
         fprintf(itap15,"%20.10f%20.10f%20.10f\n",temp[ij],temp[ij+1],temp[ij+2]);

      rewind(itap15);
      fclose(itap15);

      free(temp);
      free(ba);
      free_matrix(sa,natom3);
      free_matrix(fa,natom3);
      free_matrix(e2p,natom3);
      free_matrix(e2m,natom3);
      free_matrix(e2a,natom3);
      free_matrix(e2c,natom3);
      free_matrix(ua,nbfso);
      if(twocon) {
         free_matrix(e11,nbfso);
         free_matrix(e12,nbfso);
         free_matrix(e22,nbfso);
         }
      }
