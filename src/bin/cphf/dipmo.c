
/* $Log$
 * Revision 1.1  2000/02/04 22:50:47  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:41  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

dipmo()

{
   int i,j,ij,ix;
   int i00,ixx,iyy,izz,kabc,iabc;
   int p1,p2;
   double debye=2.541765480;
   double bohr=0.52917706;
   double deb4 = 4.0*debye;
   double deb2 = 2.0*debye;
   double rbohr = 1.0/bohr;
   double f;
   double **dip,***ua,**sa;
   double **dipda,**dipdm,**dipdt,**dipde,**dipdn,**dipdc;
   double *temp,*tmp2;
   FILE *itap17;

   ffile(&itap17,"file17.dat",0);

   dip = (double **) init_matrix(3,nbstri);
   sa = (double **) init_matrix(3,nbstri);
   dipda = (double **) init_matrix(3,natom3);
   dipdm = (double **) init_matrix(3,natom3);
   dipdn = (double **) init_matrix(3,natom3);
   dipde = (double **) init_matrix(3,natom3);
   dipdt = (double **) init_matrix(3,natom3);
   dipdc = (double **) init_matrix(3,natom3);
   ua = (double ***) malloc(sizeof(double **)*3);
   for(i=0; i < 3 ; i++) ua[i] = (double **) init_matrix(nbfso,nbfso);

   temp = (double *) init_array(natom3*3);
   tmp2 = (double *) init_array(nbfso*nbfso*natom);

   srew(itap43);

   sread(itap43,(char *) temp,sizeof(double)*natom3*3);
   for(i=ij=0; i < natom3 ; i++)
      for(j=0; j < 3 ; j++,ij++)
         dipda[j][i] = temp[ij];

   sread(itap43,(char *) temp,sizeof(double)*natom3*3);
   for(i=ij=0; i < natom3 ; i++)
      for(j=0; j < 3 ; j++,ij++)
         dipdn[j][i] = temp[ij];

   for(ix=0; ix < 3 ; ix++)
      rread(itap44,(char *) dip[ix],sizeof(double)*nbstri,ha_loc[natom3+ix]);

   for(iabc=kabc=0; iabc < natom ; iabc++) {
      i00=iabc*3;
      ixx=i00;
      iyy=i00+1;
      izz=i00+2;     

      for(ix=0; ix < 3 ; ix++,kabc++) {
         rread(itap44,(char *) sa[ix],sizeof(double)*nbstri,sa_loc[kabc]);
         rread(itap44,(char *) tmp2,sizeof(double)*nbfso*nbfso,ua_loc[kabc]);
         for(i=ij=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++,ij++)
               ua[ix][j][i] = tmp2[ij];
         }

      for(i=0; i < nocc ; i++) {
         f = (iopen) ? 0.5*occ_num[i] : 0.5;
         if(!iopen) {
            for(j=0; j < nc ; j++) {
               p1=MAX0(i,j);
               p2=MIN0(i,j);
               ij=ioff[p1]+p2;
               dipdm[0][ixx] -= sa[0][ij]*dip[0][ij]*f;
               dipdm[0][iyy] -= sa[1][ij]*dip[0][ij]*f;
               dipdm[0][izz] -= sa[2][ij]*dip[0][ij]*f;
               dipdm[1][ixx] -= sa[0][ij]*dip[1][ij]*f;
               dipdm[1][iyy] -= sa[1][ij]*dip[1][ij]*f;
               dipdm[1][izz] -= sa[2][ij]*dip[1][ij]*f;
               dipdm[2][ixx] -= sa[0][ij]*dip[2][ij]*f;
               dipdm[2][iyy] -= sa[1][ij]*dip[2][ij]*f;
               dipdm[2][izz] -= sa[2][ij]*dip[2][ij]*f;
               }
            for(j=nc; j < nbfso ; j++) {
               p1=MAX0(i,j);
               p2=MIN0(i,j);
               ij=ioff[p1]+p2;
               dipdm[0][ixx] += ua[0][j][i]*dip[0][ij];
               dipdm[0][iyy] += ua[1][j][i]*dip[0][ij];
               dipdm[0][izz] += ua[2][j][i]*dip[0][ij];
               dipdm[1][ixx] += ua[0][j][i]*dip[1][ij];
               dipdm[1][iyy] += ua[1][j][i]*dip[1][ij];
               dipdm[1][izz] += ua[2][j][i]*dip[1][ij];
               dipdm[2][ixx] += ua[0][j][i]*dip[2][ij];
               dipdm[2][iyy] += ua[1][j][i]*dip[2][ij];
               dipdm[2][izz] += ua[2][j][i]*dip[2][ij];
               }
            }
         else {
            for(j=0; j < nbfso ; j++) {
               p1=MAX0(i,j);
               p2=MIN0(i,j);
               ij=ioff[p1]+p2;
               dipdm[0][ixx] += ua[0][j][i]*dip[0][ij]*f;
               dipdm[0][iyy] += ua[1][j][i]*dip[0][ij]*f;
               dipdm[0][izz] += ua[2][j][i]*dip[0][ij]*f;
               dipdm[1][ixx] += ua[0][j][i]*dip[1][ij]*f;
               dipdm[1][iyy] += ua[1][j][i]*dip[1][ij]*f;
               dipdm[1][izz] += ua[2][j][i]*dip[1][ij]*f;
               dipdm[2][ixx] += ua[0][j][i]*dip[2][ij]*f;
               dipdm[2][iyy] += ua[1][j][i]*dip[2][ij]*f;
               dipdm[2][izz] += ua[2][j][i]*dip[2][ij]*f;
               }
            }
         }
      if(twocon) {
         dipdc[0][ixx] = c1*c1a[ixx]*h11f[0] + c2*c2a[ixx]*h22f[0];
         dipdc[0][iyy] = c1*c1a[iyy]*h11f[0] + c2*c2a[iyy]*h22f[0];
         dipdc[0][izz] = c1*c1a[izz]*h11f[0] + c2*c2a[izz]*h22f[0];
         dipdc[1][ixx] = c1*c1a[ixx]*h11f[1] + c2*c2a[ixx]*h22f[1];
         dipdc[1][iyy] = c1*c1a[iyy]*h11f[1] + c2*c2a[iyy]*h22f[1];
         dipdc[1][izz] = c1*c1a[izz]*h11f[1] + c2*c2a[izz]*h22f[1];
         dipdc[2][ixx] = c1*c1a[ixx]*h11f[2] + c2*c2a[ixx]*h22f[2];
         dipdc[2][iyy] = c1*c1a[iyy]*h11f[2] + c2*c2a[iyy]*h22f[2];
         dipdc[2][izz] = c1*c1a[izz]*h11f[2] + c2*c2a[izz]*h22f[2];
         }
      }

   for(i=0; i < 3 ; i++)
      for(j=0; j < natom3 ; j++) {
         dipdm[i][j] *= -deb4;
         if(twocon) dipdc[i][j] *= -deb2;
         }

   if(print & 128) {
      fprintf(outfile,"\n dipdm matrix\n");
      print_mat(dipdm,3,natom3,outfile);
      if(twocon) {
         fprintf(outfile,"\n dipdc matrix\n");
         print_mat(dipdc,3,natom3,outfile);
         }
      }

   for(i=0; i < 3 ; i++)
      for(j=0; j < natom3 ; j++)
         dipde[i][j] = dipda[i][j]+dipdm[i][j]+dipdc[i][j];

   if(print & 128) {
      fprintf(outfile,"\n dipde matrix\n");
      print_mat(dipde,3,natom3,outfile);
      }

   for(i=0; i < 3 ; i++)
      for(j=0; j < natom3 ; j++)
         dipdt[i][j] = dipde[i][j]+dipdn[i][j];
     
   if(print & 1) {
      fprintf(outfile,
                  "\n first derivatives of dipole moments (debye/bohr)\n");
      print_mat(dipdt,3,natom3,outfile);
      }
     
   for(i=0; i < 3 ; i++)
      for(j=0; j < natom3 ; j++)
         dipdt[i][j] *= rbohr;
     
   if(print & 1) {
      fprintf(outfile,"\n first derivatives of dipole moments (debye/a)\n");
      print_mat(dipdt,3,natom3,outfile);
      }

   fprintf(itap17,"%5d%5d\n",natom,natom3);
   for(i=ij=0; i < 3 ; i++)
      for(j=0; j < natom3 ; j++,ij++)
         tmp2[ij] = dipdt[i][j];
   for(ij=0; ij < natom3*3 ; ij += 3)
      fprintf(itap17,"%20.10f%20.10f%20.10f\n",tmp2[ij],tmp2[ij+1],tmp2[ij+2]);

   rewind(itap17);
   fclose(itap17);
   rclose(itap43,3);
   }
