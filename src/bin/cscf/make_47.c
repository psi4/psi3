/* $Log$
 * Revision 1.1  2000/02/04 22:52:31  evaleev
 * Initial revision
 *
/* Revision 1.1.1.1  1999/04/12 16:59:27  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
/*
 * Revision 1.3  1997/09/12  13:54:29  crawdad
 * Changing marco name from ULL to PSI_FPTR.
 *
 * Revision 1.2  1997/08/25  21:51:30  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.1  1991/06/15  20:22:33  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

make_47(lag_mo,temp,ci_type,natoms,ntri,nbfao)
   double *lag_mo,*temp;
   char *ci_type;
   int natoms,ntri,nbfao;

{
   int i,j,k,ij,nn;
   int os1,ob1,os2,ob2,n1,n2;
   int ncs,nos,nvs;
   int nbsq,nind,nbatri,nbasq;
   PSI_FPTR junk;
   int data[4],nsorb[10],ia47[192];
   int size_of_nada;
   int first=1;
   int itap49=49;
   int itap47=47;
   char *nada;
   double focc[10],*occs;
   double *alp,*bet;
   struct symm *s;
   union some_people_are_idiots {
      int dumb[2];
      double dumber;
      } stupid;

   size_of_nada = (ntri > 4096) ? ntri : 4096;
   nada = (char *) init_array(size_of_nada);
   bzero(ia47,sizeof(int)*192);
   occs = (double *) init_array(nbasis*2);

   ntri=ioff[nbasis];
   nbsq = nbasis*nbasis;
   nbasq = nbfao*nbfao;
   nbatri = ioff[nbfao];

   alp = (double *) init_array(ntri);
   bet = (double *) init_array(ntri);

   rfile(itap47);
   srew(itap47);
   srew(itap49);

   bzero(nsorb,sizeof(int)*10);
   bzero(focc,sizeof(double)*10);

/* first let's get file47 out of the way */

   ncs=nos=0;
   for (i=ij=k=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
         for(j=0; j < nn ; j++) temp[k++] = s->fock_evals[j];
         for(j=0; j < s->nclosed ; j++,ij+=2) occs[ij]=2.0;
         for(j=s->nclosed;j < s->nclosed+s->nopen ;j++,ij+=2) 
            occs[ij+1]=s->occ_num[j];
         for(j=s->nclosed+s->nopen; j < nn ; j++,ij+=2) occs[ij]=0.0;

         ncs += s->nclosed;
         nos += s->nopen;
         }
      }
   nvs=nbasis-ncs-nos;

   if(!iopen) nind=ncs*nvs;
   else if(hsos) nind=ncs*nvs+ncs*nos+nos*nvs;
   else nind=2*ncs+ncs*nvs+2*nvs+1;
   if(twocon) nind += 2;

/* data[0] = ntypes, data[1] = 3 for grscf, 4 for tcscf */
/* data[2] = dertyp, data[3] = 2 for ci, 4 for gvbci */

   data[0]=data[2]=1;
   data[1]=3;
   data[3]=2;
   if(hsos) {
      data[0]=2;
      }
   else if(singlet) {
      data[0]=3;
      }
   else if(twocon) {
      data[0]=3;
      data[1]=4;
      if(!strcmp(ci_type,"GVBCI")) data[3]=4;
      }
   swrit(itap49,(char *) data,sizeof(int)*4);

   if(!twocon) save_ci1=save_ci2=0.0;

   ia47[0]=nbasis;
   ia47[1]=ncs;
   ia47[2]=nos;
   ia47[3]=ncs+nos;
   ia47[4]=ntri;
   ia47[5]=natoms;
   ia47[6]=data[0];
   ia47[7]=data[1];
   ia47[8]=data[2];
   ia47[9]=data[3];

   stupid.dumber=repnuc;
   ia47[11]=stupid.dumb[0];
   ia47[12]=stupid.dumb[1];
   stupid.dumber=save_ci1;
   ia47[13]=stupid.dumb[0];
   ia47[14]=stupid.dumb[1];
   stupid.dumber=save_ci2;
   ia47[15]=stupid.dumb[0];
   ia47[16]=stupid.dumb[1];

   ia47[100]=1;
   ia47[101]=ia47[100]+192;
   ia47[102]=ia47[101]+nbasis;
   ia47[103]=ia47[102]+2*nbasis;
   ia47[104]=ia47[103]+10;
   ia47[105]=ia47[104]+5;
   ia47[106]=ia47[105]+2*ntri;
   ia47[107]=ia47[106]+2*ntri;
   ia47[108]=ia47[107]+2*nbsq;
   ia47[109]=ia47[108]+2*nbsq;
   ia47[110]=ia47[109]+2*nbasis*nbfao;
   ia47[111]=ia47[110]+2*ntri;
   ia47[112]=ia47[111]+2*nbatri;
   ia47[113]=ia47[112]+4*nbasis;
   ia47[114]=ia47[113]+2*nbatri;
   ia47[115]=ia47[114]+2*nbatri;
   ia47[116]=ia47[115]+2*nbatri;
   ia47[117]=ia47[116]+nbasis;
   ia47[118]=ia47[117]+2*nbsq;
   ia47[119]=ia47[118]+2*nbasis*nbfao;
   ia47[120]=ia47[119]+2*nbsq;
   ia47[121]=ia47[120]+2*nbatri;
   ia47[122]=ia47[121]+2*nbatri;
   ia47[123]=ia47[122]+2*ntri;
   ia47[124]=ia47[123]+2*nind;
   ia47[125]=ia47[124]+2*nbsq*3*natoms;
   ia47[126]=ia47[125]+2*natoms*natoms*9;
   ia47[127]=ia47[126]+nbatri;
   ia47[128]=ia47[127]+2*nbsq;
   ia47[129]=ia47[128]+2*nbsq;
   ia47[130]=ia47[129]+2*3*natoms;
   ia47[131]=ia47[130]+2*nbatri;
   ia47[150]=192;
   ia47[151]=nbasis;
   ia47[152]=2*nbasis;
   ia47[153]=10;
   ia47[154]=5;
   ia47[155]=2*ntri;
   ia47[156]=2*ntri;
   ia47[157]=2*nbsq;
   ia47[158]=2*nbsq;
   ia47[159]=2*nbasis*nbfao;
   ia47[160]=2*ntri;
   ia47[161]=2*nbatri;
   ia47[162]=4*nbasis;
   ia47[163]=2*nbatri;
   ia47[164]=2*nbatri;
   ia47[165]=2*nbatri;
   ia47[166]=nbasis;
   ia47[167]=2*nbsq;
   ia47[168]=2*nbasis*nbfao;
   ia47[169]=2*nbsq;
   ia47[170]=2*nbatri;
   ia47[171]=2*nbatri;
   ia47[172]=2*ntri;
   ia47[173]=2*nind;
   ia47[174]=2*nbsq*3*natoms;
   ia47[175]=2*natoms*natoms*9;
   ia47[176]=nbatri;
   ia47[177]=2*nbsq;
   ia47[178]=2*nbsq;
   ia47[179]=2*3*natoms;
   ia47[180]=2*nbatri;
   ia47[181]=2*nbatri*nbasis;

/* zero out file47      */
   junk = sizeof(int)*(ia47[131]-1);
   n1 = i2sec(junk);
   bzero(nada,4096);
   for(i=0,junk=1; i < n1 ; i++) 
      wwritw(itap47,(char *) nada,4096,junk,&junk);
   wwritw(itap47,(char *) ia47,sizeof(int)*192,
          (PSI_FPTR) sizeof(int)*(ia47[100]-1),&junk);

 /* write out mo lagrangian */

   wwritw(itap47,(char *) lag_mo,sizeof(double)*ntri,
          (PSI_FPTR) sizeof(int)*(ia47[110]-1),&junk);
   
 /* eigenvalues */

   wwritw(itap47,(char *) temp,sizeof(double)*nbasis,
          (PSI_FPTR) sizeof(int)*(ia47[102]-1),&junk);

/* make space in file49 for some density matrices to be added later */

   swrit(itap49,(char *) nada,sizeof(double)*ntri);
   if(iopen) swrit(itap49,(char *) nada,sizeof(double)*ntri);
   if(twocon || singlet) swrit(itap49,(char *) nada,sizeof(double)*ntri);
   swrit(itap49,(char *) lag_mo,sizeof(double)*ntri);
   
/* focc and nsorb */

   if(!iopen) {
      focc[0]=2.0;
      focc[1]=0.0;
      nsorb[0]=ncs;
      nsorb[1]=nvs;
      swrit(itap49,(char *) focc,sizeof(double));
      swrit(itap49,(char *) &focc[1],sizeof(double));
      swrit(itap49,(char *) &focc[1],sizeof(double));
      }
   else if(hsos) {
      focc[0]=2.0;
      focc[1]=1.0;
      focc[2]=0.0;
      nsorb[0]=ncs;
      nsorb[1]=nos;
      nsorb[2]=nvs;
      swrit(itap49,(char *) focc,sizeof(double)*2);
      swrit(itap49,(char *) &focc[2],sizeof(double)*3);
      temp[0]=temp[1]=0.0;
      temp[2]= -1.0;
      swrit(itap49,(char *) temp,sizeof(double)*3);
      }
   else if(singlet) {
      focc[0]=2.0;
      focc[1]=1.0;
      focc[2]=1.0;
      focc[3]=0.0;
      nsorb[0]=ncs;
      nsorb[1]=1;
      nsorb[2]=1;
      nsorb[3]=nvs;
      swrit(itap49,(char *) focc,sizeof(double)*3);
      swrit(itap49,(char *) &focc[3],sizeof(double)*6);
      temp[0]=temp[1]=temp[3]=0.0;
      temp[2]=temp[5]= -1.0;
      temp[4]=3.0;
      swrit(itap49,(char *) temp,sizeof(double)*6);
      }
   else if(twocon) {
      focc[0]=2.0;
      focc[1]=scf_info[opblk1].occ_num[opshl1];
      focc[2]=scf_info[opblk2].occ_num[opshl2];
      focc[3]=0.0;
      nsorb[0]=ncs;
      nsorb[1]=1;
      nsorb[2]=1;
      nsorb[3]=nvs;
      swrit(itap49,(char *) focc,sizeof(double)*3);
      temp[0]=temp[1]=temp[3]=0.0;
      temp[4]=1.0;
      temp[2]=1.0-(1.0/focc[1]);
      temp[5]=1.0-(1.0/focc[2]);
      swrit(itap49,(char *) temp,sizeof(double)*6);
      temp[0]=temp[1]=temp[3]=0.0;
      temp[2]=temp[5]= 1.0;
      temp[4]=beta[1];
      swrit(itap49,(char *) temp,sizeof(double)*6);
      }

   wwritw(itap47,(char *) focc,sizeof(double)*5,
          (PSI_FPTR) sizeof(int)*(ia47[103]-1),&junk);
   wwritw(itap47,(char *) nsorb,sizeof(int)*5,
          (PSI_FPTR) sizeof(int)*(ia47[104]-1),&junk);
   wwritw(itap47,(char *) occs,sizeof(double)*2*nbasis,
          (PSI_FPTR) sizeof(int)*(ia47[112]-1),&junk);

/* now make alpha and beta matrices for ngrcphf */

   bzero(alp,sizeof(double)*ntri);
   bzero(bet,sizeof(double)*ntri);
   for(i=ij=0; i < ncs ; i++) {
      for(j=0; j <= i ; j++,ij++) {
         alp[ij]=2.0;
         bet[ij]= -1.0;
         }
      }
   if(hsos) {
      for(i=ncs; i < ncs+nos ; i++) {
         for(j=0; j < ncs ; j++,ij++) {
            alp[ij]=1.0;
            bet[ij]= -0.5;
            }
         for(j=ncs; j <= i ; j++,ij++) {
            alp[ij]=0.5;
            bet[ij]= -0.5;
            }
         }
      }
   else if(singlet) {
      for(j=0; j < ncs ; j++,ij++) {
         alp[ij]=1.0;
         bet[ij] = -0.5;
         }
      ij++;
      for(j=0; j < ncs ; j++,ij++) {
         alp[ij]=1.0;
         bet[ij] = -0.5;
         }
      ij=ioff[ncs+1]+ncs;
      alp[ij]=0.5;
      bet[ij]=0.5;
      }
   else if(twocon) {
      for(j=0; j < ncs ; j++,ij++) {
         alp[ij]=focc[1];
         bet[ij] = -focc[1]*0.5;
         }
      ij++;
      for(j=0; j < ncs ; j++,ij++) {
         alp[ij]=focc[2];
         bet[ij] = -focc[2]*0.5;
         }
      ij=ioff[ncs]+ncs;
      alp[ij]=0.5*focc[1];
      bet[ij]=0.0;
      ij=ioff[ncs+1]+ncs;
      alp[ij]=0.0;
      bet[ij] = -0.5*(sqrt(focc[1]*focc[2]));
      ij++;
      alp[ij]=0.5*focc[2];
      bet[ij]=0.0;
      }

   wwritw(itap47,(char *) alp,sizeof(double)*ntri,
          (PSI_FPTR) sizeof(int)*(ia47[105]-1),&junk);
   wwritw(itap47,(char *) bet,sizeof(double)*ntri,
          (PSI_FPTR) sizeof(int)*(ia47[106]-1),&junk);
   }
