
/* $Log$
 * Revision 1.1  2000/02/04 22:50:47  evaleev
 * Initial revision
 *
/* Revision 1.3  1997/08/25 21:53:38  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.2  1995/01/19  19:48:39  seidl
 * replace some nulls with spaces
 *
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

cphf_iter(max_words)
   int max_words;

   {
      int i,j,ij,iii;
      int k,l;
      int iabc,it,jt,ii,jj;
      int nabc,nvec,nvtot;
      int npvec=0;
      int core_left,core_needed,num_passes;
      int ilast,ifirst,igroup,niter;
      int nadd=natom3x;
      int naddv=0;
      int maxvc,nindp2;
      double bnorm[1024],f,xbb,xtt,sum;
      double aval,bval,determ;
      double **zeta,**a,**bb,**w;
      double **bt,**b1t;
      double *aa,*tt,*b0;
      double *temp;

      if(twocon) {
         c1a = (double *) init_array(natom3x);
         c2a = (double *) init_array(natom3x);
         }
      nindp2 = (twocon) ? nind+2 : nind;
      maxvc = MIN0(600,nadd*20);
      core_left = max_words-2*nbfso*nindp2-4*nbfso*nbfso-5*nbstri;
      if(iopen) core_left -= (nbfso*nbfso+2*ntype1*nbstri);
      core_needed = maxvc*(maxvc+nadd)+2*maxvc*nbstri;
      if(iopen) 
         core_needed += 3*ntypes*nbstri*nadd; /* denp,denq,denmt */
      else
         core_needed += 2*nbstri*nadd; /* denp,denm */

      if(core_needed < core_left) nadd=natom3x;
      else {
         num_passes = core_needed/core_left+1;
         nadd = natom3x/num_passes+2;
         }
      maxvc = MIN0(600,nadd*20);
      core_needed = maxvc*(maxvc+nadd)+2*maxvc*nbstri;
      if(iopen) 
         core_needed += 3*ntypes*nbstri*nadd; /* denp,denq,denmt */
      else
         core_needed += 2*nbstri*nadd; /* denp,denm */

      fprintf(outfile,"\n%8callocating approximately %d real words of memory\n",
                       ' ',core_needed);
      fflush(outfile);

      if(iopen) {
         zeta = (double **) init_matrix(ntype1,nbstri);
         a = (double **) init_matrix(nbfso,nbfso);
         }
      bb = (double **) init_matrix(nindp2,2*maxvc);
      bt = (double **) malloc(sizeof(double *)*nindp2);
      b1t = (double **) malloc(sizeof(double *)*nindp2);
      w = (double **) init_matrix(maxvc,maxvc+nadd);
      aa = (double *) init_array(nbstri);
      tt = (double *) init_array(nbstri);
      b0 = (double *) init_array(nbstri);
      temp = (double *) init_array(nbstri);

      if(iopen) {
         for(i=0; i < ntypes ; i++) 
            mread(zeta[i],26+i);

   /* form a matrix (eqn A-1) */

         for(i=0; i < nbfso ; i++)
            for(j=0; j < nbfso ; j++) {
               it = motyp[i];
               jt = motyp[j];
               ii = ioff[i]+i;
               jj = ioff[j]+j;
               a[i][j] += zeta[jt][ii]-zeta[it][ii]+zeta[it][jj]-zeta[jt][jj];
               }

         if(print & 64) {
            fprintf(outfile,"\na matrix\n");
            print_mat(a,nbfso,nbfso,outfile);
            }
         }

/* scaling matrix (T matrix in eqn 24) */

      for(ii=0; ii < nind ; ii++) {
         i=indep[ii].ii;
         j=indep[ii].jj;
         if (iopen) 
            aa[ii] = 1.0/sqrt(fabs(a[i][j]));
         else
            aa[ii] = 1.0/sqrt(fabs(e_vals[i]-e_vals[j]));
         }
      if(twocon) {
         aa[nind]=1.0/sqrt(fabs(a22[0][0]));
         aa[nind+1]=1.0/sqrt(fabs(a22[1][1]));
         }

/* start grouping the b vectors */

      igroup=ilast=0;
      do {
         igroup++;
         ifirst = ilast+1;
         ilast = ifirst+nadd-1;
         if(ilast > natom3x) ilast=natom3x;
         nabc = ilast-ifirst+1;

         fprintf(outfile,"\n\tigroup = %5d   ifirst = %5d   nabc = %5d\n\n",
                              igroup,ifirst,nabc);
         fflush(outfile);
         srew(work);

/* read in b0 matrices */

         for(iabc=0; iabc < nabc ; iabc++) {
            iii=iabc+ifirst-1;
            rread(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[iii]);
            for(ii=0; ii < nind ; ii++) {
               ij=indep[ii].ij;
               bb[ii][iabc]=b0[ij];
               }
            if(twocon) {
               bb[nind][iabc]=baf[iii].one;
               bb[nind+1][iabc]=baf[iii].two;
               }
            }

         if(print & 64) {
            fprintf(outfile,"\nbb matrix\n");
            print_mat(bb,nindp2,nabc,outfile);
            }

/* scale b0 matrices */

         for(ii=0; ii < nindp2 ; ii++)
            for(iabc=0; iabc < nabc ; iabc++)
               bb[ii][iabc] *= aa[ii];

         if(print & 64) {
            fprintf(outfile,"\nscaled bb matrix\n");
            print_mat(bb,nindp2,nabc,outfile);
            }

/* orthonormalize b0 matrix */

         orthog(bb,nindp2,nabc,&nvec);
         for(iabc=0; iabc < nvec ; iabc++) bnorm[iabc]=1.0;
         if(print & 64) {
            fprintf(outfile,"\northonormalized bb matrix\n");
            print_mat(bb,nindp2,nvec,outfile);
            }

         if(print & 64) {
            double val;
            for(i=0; i < nvec ; i++)
               for(j=0; j < nvec ; j++) {
                  for(k=0,val=0.0; k < nindp2 ; k++) val += bb[k][i]*bb[k][j];
                  w[i][j]=val;
                  }
            fprintf(outfile,"\northonormality check matrix\n");
            print_mat(w,nvec,nvec,outfile);
            }

/* iteration starts here */

         npvec=0;
         niter=0;
         do {
            niter++;
            for(ii=0; ii < nindp2 ; ii++) {
               bt[ii] = bb[ii]+npvec;
               b1t[ii] = bt[ii]+nvec;
               for(iabc=0; iabc < nvec ; iabc++)
                  bt[ii][iabc] *= aa[ii];
               }

      /* form a*b matrix (eqns 12,21-23) */

            if(iopen) makeab_o(bt,b1t,nvec,zeta);
            else makeab(bt,b1t,nvec);

            for(ii=0; ii < nindp2 ; ii++)
               for(iabc=0; iabc < nvec ; iabc++) {
                  bt[ii][iabc] /= aa[ii];
                  b1t[ii][iabc] *= aa[ii];
                  }

      /* store a*b matrix */

            for(iabc=0; iabc < nvec ; iabc++) {
               for(ii=0; ii < nindp2 ; ii++) temp[ii]=b1t[ii][iabc];
               swrit(work,(char *) temp,sizeof(double)*nindp2);
               }

      /* form (1-a)*b matrix (eqn 20) */

            if(iopen) {
               for(ii=0; ii < nindp2 ; ii++)
                  for(iabc=0; iabc < nvec ; iabc++)
                     b1t[ii][iabc] = bt[ii][iabc]-b1t[ii][iabc];
               }
            else {
               for(ii=0; ii < nind ; ii++)
                  for(iabc=0; iabc < nvec ; iabc++)
                     b1t[ii][iabc] = -bt[ii][iabc]-b1t[ii][iabc];
               }

            if(print & 64) {
               fprintf(outfile,"\n(1-a)*b matrix\n");
               print_mat(b1t,nindp2,nvec,outfile);
               }

       /* form additional B' matrices (eqn 17) */

            naddv=0;
            k=npvec+nvec;
            for(iabc=0; iabc < nvec ; iabc++) {
               for(ii=0; ii < nindp2 ; ii++) tt[ii]=b1t[ii][iabc];
               for(l=0; l < k ; l++) {
                  for(ii=0,xbb=0.0; ii < nindp2 ; ii++)
                      xbb += bb[ii][l]*b1t[ii][iabc];
                  f = xbb/bnorm[l];
                  for(ii=0; ii < nindp2 ; ii++) tt[ii] -= bb[ii][l]*f;
                  }
               for(ii=0,xtt=0.0; ii < nindp2 ; ii++) xtt += tt[ii]*tt[ii];
               if(fabs(xtt) > pow(10.0,(double) -(conv))) {
                  naddv++;
                  bnorm[k] = xtt;
                  if(iopen)
                     for(ii=0; ii < nindp2 ; ii++) bb[ii][k] = tt[ii];
                  else
                     for(ii=0; ii < nindp2 ; ii++) bb[ii][k] = -tt[ii];
                  k++;
                  }
               }

            fprintf(outfile,"\titeration = %3d   nvec = %3d   naddv = %3d\n",
                                                               niter,nvec,naddv);
            fflush(outfile);
            npvec += nvec;
            nvec = naddv;

            if((print & 64) && naddv) {
               fprintf(outfile,"\nnorms\n");
               for(ii=npvec; ii < npvec+nvec ;ii++)
                          fprintf(outfile,"%5d %20.10f\n",ii,bnorm[ii]);

               fprintf(outfile,"\nbb matrix\n");
               print_mat(b1t,nindp2,nvec,outfile);
               }
            } while(naddv);

         nvtot = npvec;
         fprintf(outfile,"\n\tend of iteration:  nvtot = %5d\n",nvtot);
         if(nvtot > maxvc) {
            fprintf(outfile,"\ntoo many vectors in cphf\n");
            fprintf(outfile,"\nmaxvc = %5d   nvtot = %5d\n",maxvc,nvtot);
            }
         
         srew(work);
 
      /* read back a*b matrix */

         for(iabc=nvtot; iabc < 2*nvtot ; iabc++) {
            sread(work,(char *) temp,sizeof(double)*nindp2);
            for(ii=0; ii < nindp2 ; ii++) bb[ii][iabc]=temp[ii];
            }

         if(print & 64) {
            fprintf(outfile,"\nexpanded b*a matrix\n");
            print_mat(bb,nindp2,2*nvtot,outfile);
            }

         for(ii=0; ii < nindp2 ; ii++) bt[ii]=bb[ii]+nvtot;

      /* lhs of eqn 19 */

         for(i=0; i < nvtot ; i++)
            for(j=0; j < nvtot ; j++) {
               for(ii=0,sum=0.0; ii < nindp2 ; ii++) sum+= bb[ii][i]*bt[ii][j];
               w[i][j]=sum;
               }

      /* read back b0 matrices */

         for(iabc=0; iabc < nabc ; iabc++) {
            iii=iabc+ifirst-1;
            rread(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[iii]);
            for(ii=0; ii < nind ; ii++) {
               ij=indep[ii].ij;
               bt[ii][iabc]=b0[ij]*aa[ii];
               }
            if(twocon) {
               bt[nind][iabc]=baf[iii].one*aa[nind];
               bt[nind+1][iabc]=baf[iii].two*aa[nind+1];
               }
            }

         if(print & 64) {
            fprintf(outfile,"\nscaled bb matrix\n");
            print_mat(bt,nindp2,nabc,outfile);
            }

      /* and put b0*b in w (rhs of eqn 19) */

         for(i=0; i < nvtot ; i++)
            for(j=0; j < nabc ; j++) {
               for(ii=0,sum=0.0; ii < nindp2 ; ii++) sum+= bb[ii][i]*bt[ii][j];
               w[i][nvtot+j]=sum;
               }

         if(print & 64) {
            fprintf(outfile,"\nw matrix\n");
            print_mat(w,nvtot,nvtot+nabc,outfile);
            }

         for(i=0; i < nvtot ; i++) {
            aval = 1.0/sqrt(bnorm[i]);
            for(j=0; j < nvtot ; j++) {
               bval = 1.0/sqrt(bnorm[j]);
               w[i][j] *= aval*bval;
               }
            }
         for(i=0; i < nvtot ; i++) {
            aval = 1.0/sqrt(bnorm[i]);
            for(j=0; j < nabc ; j++) w[i][nvtot+j] *= aval;
            }
               
         if(print & 64) {
            fprintf(outfile,"\nscaled w matrix\n");
            print_mat(w,nvtot,nvtot+nabc,outfile);
            }

         /* solve eqn 19 */

         linear_eq(w,nvtot,nabc,&determ);

         if(fabs(determ) < 10.0e-75) {
            fprintf(outfile,"\nlinear equation cannot be solved\n");
            fprintf(outfile,"determ = %10.7e\n",determ);
            exit(1);
            }
         fprintf(outfile,"\n\tlinear equation solved:  determ = %10.7e\n",determ);

         for(i=0; i < nvtot ; i++) {
            aval = 1.0/sqrt(bnorm[i]);
            for(j=0; j < nabc ; j++) w[i][nvtot+j] *= aval;
            }

/* form independent elements of ua matrix (eqn 18) */

         for(ii=0; ii < nindp2 ; ii++)
            for(iabc=0; iabc < nabc ; iabc++) {
               for(i=0,sum=0.0; i < nvtot ; i++) 
                  sum += bb[ii][i]*w[i][iabc+nvtot];
               bt[ii][iabc] = sum*aa[ii];
               }
               
         for(iabc=0; iabc < nabc ; iabc++) {
            iii=iabc+ifirst-1;
            rread(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[iii]);
            for(ii=0; ii < nind ; ii++) {
               ij=indep[ii].ij;
               b0[ij]=bt[ii][iabc];
               }
            rwrit(itap44,(char *) b0,sizeof(double)*nbstri,ba_loc[iii]);
            if(twocon) {
               c1a[iii]=bt[nind][iabc];
               c2a[iii]=bt[nind+1][iabc];
               }
            if(print & 64) {
               fprintf(outfile,"\nb0 iabc = %5d\n",iabc);
               print_array(b0,nbfso,outfile);
               if(twocon) 
                  fprintf(outfile,"\n %20.10f%20.10f\n",c1a[iii],c2a[iii]);
               }
            }
         fflush(outfile);
         
         } while(ilast < natom3x);
      srew(work);

      if(iopen) {
         free_matrix(zeta,ntype1);
         free_matrix(a,nbfso);
         }
      free_matrix(bb,nindp2);
      free_matrix(w,maxvc);
      free(bt);
      free(b1t);
      free(aa);
      free(tt);
      free(b0);
      free(temp);
      }
