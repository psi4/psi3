
/* $Log$
 * Revision 1.1  2000/02/04 22:50:47  evaleev
 * Initial revision
 *
/* Revision 1.5  1997/09/12 13:54:45  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 1.4  1997/08/25  21:53:42  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.3  1991/07/30  04:51:08  seidl
 * fix for f and g functions
 *
 * Revision 1.2  1991/07/13  08:51:31  seidl
 * change dimensions of contr and ex, ey
 * and ez to accomodate f and g functions
 *
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"
#include "dipder.h"

/* calculates dipole derivatives                                    */
/* just a translation of the HFS group code DIPDER, which is mostly */
/* ripped off from HONDO (I believe)                                */

dipole_derivs(dipda,dipdn,p,dder)
   double **dipda,**dipdn,*p,***dder;
{
   int i,j,ij,mm,nn;
   int mconst,mcalcs,mpoint,ncalcs;
   PSI_FPTR loccal,junk,locvec;
   int junkpoint;
   int nprim,nshell;
   int p1,p2;
   int iatom,jatom,ish,jsh,igmin,igmax,jgmin,jgmax;
   int lit,ljt,mini,minj,maxi,maxj,jmax,loci,locj;
   int ig,jg,ipqr,iabc;
   int nx,ny,nz,ii,jj,idxi,idxj,idxs;
   int in,ni,nj,jn;
   int ixx,iyy,izz,i00;
   int iandj,maxij,equal,mij;
   int itap30=30;
   int itap31=31;
   int itol=10,nqlim=7;
   int i30[200];
   double value,tol=23.0258;
   double x,y,z,xfac,yfac,zfac,yzfac,xzfac,xyfac;
   double ax,ay,az,bx,by,bz,xab,yab,zab,rab;
   double axi,ayi,azi,ai,csi;
   double bxj,byj,bzj,bj,csj;
   double pp,tpp,tp,tmp,dij,pfac,ptwo,px,py,pz,pax,pay,paz,pbx,pby,pbz;
   double cx,cy,cz,pcx,pcy,pcz;
   double sx11,sy11,sz11;
   double dx1a,dy1a,dz1a,dx1b,dy1b,dz1b;
   double dxy1,dxz1,dyx1,dyz1,dzx1,dzy1,t;
   double dex,dey,dez,fac,dnx,dny,dnz,dtx,dty,dtz,dipol,znuc;
   double a12,b12,a11,b11;
   double pi212=1.1283791670955;
   double pi32=5.56832799683170;
   double debye=2.54176548;
   double *dipx,*dipy,*dipz;
   double **uu,**tt,**so_to_ao,*temp;

   double charges[MAX_BASIS],coords[MAX_BASIS][3];
   double zetas[MAX_BASIS],contr_coeffs[MAX_BASIS],contr[MAX_BASIS*5];
   double ex[7][7][20],ey[7][7][20],ez[7][7][20];
   double sx0[36],sy0[36],sz0[36],dx0[36],dy0[36],dz0[36];
   double sx1[36],sy1[36],sz1[36],dx1[36][2],dy1[36][2],dz1[36][2];
   
   int ijx[MAX_BASIS],ijy[MAX_BASIS],ijz[MAX_BASIS],ijpos[MAX_BASIS];
   int knuc[MAX_BASIS],ktype[MAX_BASIS],kstart[MAX_BASIS],kloc[MAX_BASIS];
   int kmin[MAX_BASIS],kmax[MAX_BASIS],kng[MAX_BASIS];

   dipx = (double *) init_array(nbatri);
   dipy = (double *) init_array(nbatri);
   dipz = (double *) init_array(nbatri);
   temp = (double *) init_array(nbfao*nbfao);
   uu = (double **) init_matrix(nbfao,nbfao);
   tt = (double **) init_matrix(nbfao,nbfao);
   so_to_ao = (double **) init_matrix(nbfao,nbfso);

   rfile(itap30);
   rfile(itap31);

   mread(temp,8);
   for(i=ij=0; i < nbfso ; i++)
      for(j=0; j < nbfao ; j++,ij++)
         so_to_ao[j][i] = temp[ij];

 /* get density matrix from file31 */

   sread(itap31,(char *) temp,sizeof(double)*nbstri);

   for(i=ij=0; i < nbfso ; i++)
      for(j=0; j <= i ; j++,ij++)
         tt[i][j]=tt[j][i] = (i==j) ? temp[ij] : 0.5*temp[ij];

   mmult(so_to_ao,0,tt,0,uu,0,nbfao,nbfso,nbfso,0);
   mmult(uu,0,so_to_ao,1,tt,0,nbfao,nbfso,nbfao,0);

   for(i=ij=0; i < nbfao ; i++)
      for(j=0; j <= i ; j++,ij++)
         p[ij]=tt[i][j];

   if(print & 32) {
      fprintf(outfile,"\ndensity matrix in dipder\n");
      print_array(p,nbfao,outfile);
      }

   rclose(itap31,3);

 /* now get stuff from file30 */

   junk = sizeof(int)*100;
   wreadw(itap30,(char *) i30,sizeof(int)*200,junk,&junk);

   mpoint = i30[1];
   mconst = i30[2];
   mcalcs = i30[3];
   ncalcs = i30[4];
   nshell = i30[26];
   nprim = i30[31];

/* read in geometry */

   junk = sizeof(int)*(100+mconst+mpoint+ncalcs-1);
   wreadw(itap30,(char *) &junkpoint,sizeof(int),junk,&junk);
   loccal = (PSI_FPTR) junkpoint;
   wreadw(itap30,(char *) i30,sizeof(int)*60,(PSI_FPTR) sizeof(int)*(loccal-1),&loccal);
   wreadw(itap30,(char *) i30,sizeof(int)*20,loccal,&loccal);

   locvec = sizeof(int)*(i30[0]-1);

   for (i=0; i < natom ; i++)
      wreadw(itap30,(char *) coords[i],sizeof(double)*3,loccal,&loccal);

   junk = sizeof(int)*(100+mconst);
   wreadw(itap30,(char *) i30,sizeof(int)*mpoint,junk,&junk);

/* get nuclear charges, etc. */

   junk = sizeof(int)*(i30[0]-1);
   wreadw(itap30,(char *) charges,sizeof(double)*natom,junk,&junk);
   junk = sizeof(int)*(i30[4]-1);
   wreadw(itap30,(char *) zetas,sizeof(double)*nprim,junk,&junk);
   junk = sizeof(int)*(i30[5]-1);
   wreadw(itap30,(char *) contr,sizeof(double)*nprim*5,junk,&junk);
   junk = sizeof(int)*(i30[6]-1);
   wreadw(itap30,(char *) kstart,sizeof(double)*nshell,junk,&junk);
   junk = sizeof(int)*(i30[7]-1);
   wreadw(itap30,(char *) knuc,sizeof(double)*nshell,junk,&junk);
   junk = sizeof(int)*(i30[8]-1);
   wreadw(itap30,(char *) ktype,sizeof(double)*nshell,junk,&junk);
   junk = sizeof(int)*(i30[9]-1);
   wreadw(itap30,(char *) kng,sizeof(double)*nshell,junk,&junk);
   junk = sizeof(int)*(i30[10]-1);
   wreadw(itap30,(char *) kloc,sizeof(double)*nshell,junk,&junk);
   junk = sizeof(int)*(i30[11]-1);
   wreadw(itap30,(char *) kmin,sizeof(double)*nshell,junk,&junk);
   junk = sizeof(int)*(i30[12]-1);
   wreadw(itap30,(char *) kmax,sizeof(double)*nshell,junk,&junk);

   for(i=ij=0; i < 5 ; i++)
      for(j=0; j < nprim ; j++,ij++)
         if(contr[ij]) contr_coeffs[j]=contr[ij];

   if(print & 32) {
      fprintf(outfile,"\n i   zeta coeff\n");
      for(i=0; i < nprim ; i++)
         fprintf(outfile,"%5d %20.10f %20.10f\n",i,zetas[i],contr_coeffs[i]);
      fprintf(outfile,"\n i    knuc ktype kstart kloc kmin kmax kng \n");
      for(i=0; i < nshell ; i++)
         fprintf(outfile,"%5d %5d %5d %5d %5d %5d %5d %5d\n",i,knuc[i],ktype[i],
                kstart[i],kloc[i],kmin[i],kmax[i],kng[i]);
      fprintf(outfile,"\n i   x   y   z   charge\n");
      for(i=0; i < natom ; i++)
         fprintf(outfile,"%5d %20.10f %20.10f %20.10f %20.10f\n",i,
           coords[i][0],coords[i][1],coords[i][2],charges[i]);
      }

   for(i=0; i < nshell ; i++) {
      knuc[i]--;
      ktype[i]--;
      kstart[i]--;
      kloc[i]--;
      kmin[i]--;
      kmax[i]--;
      kng[i]--;
      }
      
 /* well, let's get to it and calculate the dipole derivatives */

 /* i shell */

   for(ish=0; ish < nshell ; ish++) {
      iatom=knuc[ish];
      ax=coords[iatom][0];
      ay=coords[iatom][1];
      az=coords[iatom][2];
      igmin=kstart[ish];
      igmax=igmin+kng[ish];
      lit=ktype[ish];
      mini=kmin[ish];
      maxi=kmax[ish];
      loci=kloc[ish];

     /*  j shell  */

      for(jsh=0; jsh <= ish ; jsh++) {
         jatom=knuc[jsh];
         bx=coords[jatom][0];
         by=coords[jatom][1];
         bz=coords[jatom][2];
         jgmin=kstart[jsh];
         jgmax=jgmin+kng[jsh];
         ljt=ktype[jsh];
         minj=kmin[jsh];
         maxj=kmax[jsh];
         locj=kloc[jsh];
         iandj=(ish==jsh);

     /* work out indices for combining 2d integrals */

         mij=0; jmax=maxj;
         for(i=mini,ii=0; i <= maxi ; i++,ii++) {
            nx=ix[i];
            ny=iy[i];
            nz=iz[i];
            if(iandj) jmax=i;
            for(j=minj,jj=0; j <= jmax ; j++,jj++,mij++) {
               ijx[mij]=nx+jx[j]-1;
               ijy[mij]=ny+jy[j]-1;
               ijz[mij]=nz+jz[j]-1;
               idxi=loci+ii;
               idxj=locj+jj;
               p1=MAX0(idxi,idxj);
               p2=MIN0(idxi,idxj);
               ijpos[mij]=ioff[p1]+p2;
               }
            }

         xab=ax-bx;
         yab=ay-by;
         zab=az-bz;
         rab=xab*xab+yab*yab+zab*zab;
         maxij=MAX0(lit+2,ljt+2);

         for(ig=igmin; ig <= igmax ; ig++) {
            csi=contr_coeffs[ig]*pi32;
            ai=zetas[ig];
            axi=ai*ax;
            ayi=ai*ay;
            azi=ai*az;

            if(iandj) jgmax=ig;
            for(jg=jgmin; jg <= jgmax ; jg++) {
               csj=contr_coeffs[jg];
               bj=zetas[jg];
               bxj=bj*bx;
               byj=bj*by;
               bzj=bj*bz;
            
               pp=ai+bj;
               tpp=1.0/pp;
               tmp=(ai*bj*rab)*tpp;

               if(tmp > tol) continue;

               tp=sqrt(tpp);
               dij=csi*csj*exp(-tmp)*tpp;
               if(iandj && (ig != jg)) dij += dij;
               pfac=dij*tp;
               ptwo=tpp*0.5;

               px=(axi+bxj)*tpp;
               py=(ayi+byj)*tpp;
               pz=(azi+bzj)*tpp;
               pax=px-ax;
               pay=py-ay;
               paz=pz-az;
               pbx=px-bx;
               pby=py-by;
               pbz=pz-bz;

               ecal(pax,pay,paz,pbx,pby,pbz,ptwo,maxij,ex,ey,ez);

               cx=cy=cz=0.0;
               pcx=px-cx;
               pcy=py-cy;
               pcz=pz-cz;

               a12=ai+ai;
               b12=bj+bj;

               in= -nqlim;
               for(ni=0; ni <= lit ; ni++) {
                  in += nqlim;
                  a11 = (double) ni;
                  for(nj=0; nj <= ljt ; nj++) {
                     jn=in+nj;
                     b11 = (double) nj;

                /* normal overlap integrals */
                     sx0[jn]=ex[ni][nj][0];
                     sy0[jn]=ey[ni][nj][0];
                     sz0[jn]=ez[ni][nj][0];

                /* first derivatives of overlap integrals */
                     sx11 = ex[ni+1][nj][0]*a12;
                     sy11 = ey[ni+1][nj][0]*a12;
                     sz11 = ez[ni+1][nj][0]*a12;
                     if(ni) {
                        sx11 -= ex[ni-1][nj][0]*a11;
                        sy11 -= ey[ni-1][nj][0]*a11;
                        sz11 -= ez[ni-1][nj][0]*a11;
                        }
                     sx1[jn]=sx11;
                     sy1[jn]=sy11;
                     sz1[jn]=sz11;

                 /* normal dipole integrals */
                     dx0[jn]=ex[ni][nj][1]+ex[ni][nj][0]*pcx;
                     dy0[jn]=ey[ni][nj][1]+ey[ni][nj][0]*pcy;
                     dz0[jn]=ez[ni][nj][1]+ez[ni][nj][0]*pcz;

                 /* first derivatives of dipole integrals */
                     dx1a=(ex[ni+1][nj][1]+ex[ni+1][nj][0]*pcx)*a12;
                     dy1a=(ey[ni+1][nj][1]+ey[ni+1][nj][0]*pcy)*a12;
                     dz1a=(ez[ni+1][nj][1]+ez[ni+1][nj][0]*pcz)*a12;
                     if(ni) {
                        dx1a -= (ex[ni-1][nj][1]+ex[ni-1][nj][0]*pcx)*a11;
                        dy1a -= (ey[ni-1][nj][1]+ey[ni-1][nj][0]*pcy)*a11;
                        dz1a -= (ez[ni-1][nj][1]+ez[ni-1][nj][0]*pcz)*a11;
                        }
                     dx1b=(ex[ni][nj+1][1]+ex[ni][nj+1][0]*pcx)*b12;
                     dy1b=(ey[ni][nj+1][1]+ey[ni][nj+1][0]*pcy)*b12;
                     dz1b=(ez[ni][nj+1][1]+ez[ni][nj+1][0]*pcz)*b12;
                     if(nj) {
                        dx1b -= (ex[ni][nj-1][1]+ex[ni][nj-1][0]*pcx)*b11;
                        dy1b -= (ey[ni][nj-1][1]+ey[ni][nj-1][0]*pcy)*b11;
                        dz1b -= (ez[ni][nj-1][1]+ez[ni][nj-1][0]*pcz)*b11;
                        }
                     dx1[jn][0]=dx1a;
                     dy1[jn][0]=dy1a;
                     dz1[jn][0]=dz1a;
                     dx1[jn][1]=dx1b;
                     dy1[jn][1]=dy1b;
                     dz1[jn][1]=dz1b;
                     }
                  }

               for(i=0; i < mij ; i++) {
                  nx=ijx[i];
                  ny=ijy[i];
                  nz=ijz[i];
                  x=sx0[nx];
                  y=sy0[ny];
                  z=sz0[nz];
                  ij=ijpos[i];
                  xfac=x*pfac;
                  yfac=y*pfac;
                  zfac=z*pfac;
                  yzfac=y*zfac;
                  xzfac=x*zfac;
                  xyfac=x*yfac;
                  dipx[ij] += dx0[nx]*yzfac;
                  dipy[ij] += dy0[ny]*xzfac;
                  dipz[ij] += dz0[nz]*xyfac;
            
                  dder[0][iatom][ij] += dx1[nx][0]*yzfac;
                  dder[0][jatom][ij] += dx1[nx][1]*yzfac;
                  dxy1 = dx0[nx]*sy1[ny]*zfac;
                  dder[1][iatom][ij] += dxy1;
                  dder[1][jatom][ij] -= dxy1;
                  dxz1 = dx0[nx]*sz1[nz]*yfac;
                  dder[2][iatom][ij] += dxz1;
                  dder[2][jatom][ij] -= dxz1;
                  dyx1 = sx1[nx]*dy0[ny]*zfac;
                  dder[3][iatom][ij] += dyx1;
                  dder[3][jatom][ij] -= dyx1;
                  dder[4][iatom][ij] += dy1[ny][0]*xzfac;
                  dder[4][jatom][ij] += dy1[ny][1]*xzfac;
                  dyz1 = dy0[ny]*sz1[nz]*xfac;
                  dder[5][iatom][ij] += dyz1;
                  dder[5][jatom][ij] -= dyz1;
                  dzx1 = sx1[nx]*dz0[nz]*yfac;
                  dder[6][iatom][ij] += dzx1;
                  dder[6][jatom][ij] -= dzx1;
                  dzy1 = sy1[ny]*dz0[nz]*xfac;
                  dder[7][iatom][ij] += dzy1;
                  dder[7][jatom][ij] -= dzy1;
                  dder[8][iatom][ij] += dz1[nz][0]*xyfac;
                  dder[8][jatom][ij] += dz1[nz][1]*xyfac;
                  }
               }
            }
         }
      }
                  
   if(print & 32) {
      fprintf(outfile,"\n dipx matrix\n");
      print_array(dipx,nbfao,outfile);
      fprintf(outfile,"\n dipy matrix\n");
      print_array(dipy,nbfao,outfile);
      fprintf(outfile,"\n dipz matrix\n");
      print_array(dipz,nbfao,outfile);
      }

 /* now calculate dipole moments */

   dex=dey=dez=0.0;
   for(i=ij=0; i < nbfao ; i++) {
      for(j=0; j <= i ; j++,ij++) {
         fac = (i==j) ? p[ij] : 2.0*p[ij];
         dex += dipx[ij]*fac;
         dey += dipy[ij]*fac;
         dez += dipz[ij]*fac;
         }
      }
   dex *= -debye;
   dey *= -debye;
   dez *= -debye;

   dnx=dny=dnz=0.0;
   for(i=0; i < natom ; i++) {
      dnx += charges[i]*coords[i][0];
      dny += charges[i]*coords[i][1];
      dnz += charges[i]*coords[i][2];
      }
   dnx *= debye;
   dny *= debye;
   dnz *= debye;

   dtx = dnx+dex;
   dty = dny+dey;
   dtz = dnz+dez;

   dipol = sqrt(dtx*dtx+dty*dty+dtz*dtz);
   fprintf(outfile,"\n\ttotal dipole moment = %15.10f debye\n",dipol);
   fflush(outfile);

   if(print & 1) {
      fprintf(outfile,"\n\tcomponents of dipole moment\n");
      fprintf(outfile,"\n\t      de         dn         dt\n");
      fprintf(outfile,"\t x %10.7f %10.7f %10.7f\n",dex,dey,dez);
      fprintf(outfile,"\t y %10.7f %10.7f %10.7f\n",dnx,dny,dnz);
      fprintf(outfile,"\t z %10.7f %10.7f %10.7f\n",dtx,dty,dtz);
      fflush(outfile);
      }

 /* transform dipole moment integrals from ao to mo basis */
   ao_to_mo(dipx,tt,e_vecs_ao,uu,nbfao,nbfso);
   ao_to_mo(dipy,tt,e_vecs_ao,uu,nbfao,nbfso);
   ao_to_mo(dipz,tt,e_vecs_ao,uu,nbfao,nbfso);

   if(print & 32) {
      fprintf(outfile,"\n dipx matrix (mo basis)\n");
      print_array(dipx,nbfso,outfile);
      fprintf(outfile,"\n dipy matrix (mo basis)\n");
      print_array(dipy,nbfso,outfile);
      fprintf(outfile,"\n dipz matrix (mo basis)\n");
      print_array(dipz,nbfso,outfile);
      }

   for(iabc=0; iabc < natom ; iabc++) {
      znuc = charges[iabc];
      i00=iabc*3;
      ixx=i00;
      iyy=i00+1;
      izz=i00+2;

      for(i=ij=0; i < nbfao ; i++) {
         for(j=0; j <= i ; j++,ij++) {
            fac = (i==j) ? p[ij] : 2.0*p[ij];
            dipda[0][ixx] += dder[0][iabc][ij]*fac;
            dipda[0][iyy] += dder[1][iabc][ij]*fac;
            dipda[0][izz] += dder[2][iabc][ij]*fac;
            dipda[1][ixx] += dder[3][iabc][ij]*fac;
            dipda[1][iyy] += dder[4][iabc][ij]*fac;
            dipda[1][izz] += dder[5][iabc][ij]*fac;
            dipda[2][ixx] += dder[6][iabc][ij]*fac;
            dipda[2][iyy] += dder[7][iabc][ij]*fac;
            dipda[2][izz] += dder[8][iabc][ij]*fac;
            }
         }
      dipdn[0][ixx] = znuc;
      dipdn[1][iyy] = znuc;
      dipdn[2][izz] = znuc;
      }
   
   for(i=0 ; i < 3 ; i++)
      for(j=0; j < natom3 ; j++) {
         dipda[i][j] *= -debye;
         dipdn[i][j] *= debye;
         }

   for(i=ij=0; i < natom3 ; i++)
      for(j=0; j < 3 ; j++,ij++)
         temp[ij]=dipda[j][i];

   srew(itap43);
   swrit(itap43,(char *) temp,sizeof(double)*natom3*3);

   for(i=ij=0; i < natom3 ; i++)
      for(j=0; j < 3 ; j++,ij++)
         temp[ij]=dipdn[j][i];

   swrit(itap43,(char *) temp,sizeof(double)*natom3*3);
   swrit(itap43,(char *) dipx,sizeof(double)*nbstri);
   swrit(itap43,(char *) dipy,sizeof(double)*nbstri);
   swrit(itap43,(char *) dipz,sizeof(double)*nbstri);

   for(iabc=0; iabc < natom ; iabc++) {
      for(i=0; i < 3 ; i++) {
         for(j=0; j < 3 ; j++) {
            ij=i+j*3;
            ao_to_mo(dder[ij][iabc],tt,e_vecs_ao,uu,nbfao,nbfso);
            swrit(itap43,(char *) dder[ij][iabc],sizeof(double)*nbstri);

            if(print & 32) {
               fprintf(outfile,"\n dder matrix %d %d %d %d\n",iabc,i,j,ij);
               print_array(dder[ij][iabc],nbfso,outfile);
               }
            }
         }
      }

   rclose(itap30,3);
   free_matrix(tt,nbfao);
   free_matrix(uu,nbfao);
   free_matrix(so_to_ao,nbfao);
   free(dipx);
   free(dipy);
   free(dipz);
   free(temp);
   }


ecal(pax,pay,paz,pbx,pby,pbz,ptwo,maxij,ex,ey,ez)
   double pax,pay,paz,pbx,pby,pbz,ptwo;
   int maxij;
   double ex[7][7][20],ey[7][7][20],ez[7][7][20];

{
   int ii,jj,ijtot,equal,kmax,km1,kk,kp1;
   double pxx,pyy,pzz,qxx,qyy,qzz,t;

   ex[0][0][0]=ey[0][0][0]=ez[0][0][0]=1.0;

   for(ii=1; ii <= maxij ; ii++) {
      for(jj=0; jj <= ii ; jj++) {
         equal = (ii==jj);
         ijtot = ii-1+jj;
         kmax = ii+jj;
         for(kk=0; kk <= kmax ; kk++) {
            pxx=pyy=pzz=qxx=qyy=qzz=0.0;

            km1=kk-1;
            if(km1 > -1) {
               pxx += ex[ii-1][jj][km1]*ptwo;
               pyy += ey[ii-1][jj][km1]*ptwo;
               pzz += ez[ii-1][jj][km1]*ptwo;
               if(!equal) {
                  qxx += ex[jj][ii-1][km1]*ptwo;
                  qyy += ey[jj][ii-1][km1]*ptwo;
                  qzz += ez[jj][ii-1][km1]*ptwo;
                  }
               }
            if(km1 < ijtot) {
               pxx += ex[ii-1][jj][kk]*pax;
               pyy += ey[ii-1][jj][kk]*pay;
               pzz += ez[ii-1][jj][kk]*paz;
               if(!equal) {
                  qxx += ex[jj][ii-1][kk]*pbx;
                  qyy += ey[jj][ii-1][kk]*pby;
                  qzz += ez[jj][ii-1][kk]*pbz;
                  }
               }
            if(kk < ijtot) {
               t = (double) kk+1;
               pxx += ex[ii-1][jj][kk+1]*t;
               pyy += ey[ii-1][jj][kk+1]*t;
               pzz += ez[ii-1][jj][kk+1]*t;
               if(!equal) {
                  qxx += ex[jj][ii-1][kk+1]*t;
                  qyy += ey[jj][ii-1][kk+1]*t;
                  qzz += ez[jj][ii-1][kk+1]*t;
                  }
               }
            ex[ii][jj][kk]=pxx;
            ey[ii][jj][kk]=pyy;
            ez[ii][jj][kk]=pzz;
            if(!equal) {
               ex[jj][ii][kk]=qxx;
               ey[jj][ii][kk]=qyy;
               ez[jj][ii][kk]=qzz;
               }
            }
         }
      }
   }
