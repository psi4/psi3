/* $Log$
 * Revision 1.1  2000/02/04 22:52:31  evaleev
 * Initial revision
 *
/* Revision 1.5  1999/11/04 19:24:30  localpsi
/* STB (11/4/99) - Added the orb_mix feature which is equivalent to guess = mix
/* in G94 and also fixed restarting so that if you have different wavefuntions,
/* everything works.  Also if you specify no DOCC and SOCC and restart, if the
/* wavefunctions are different, it will guess again.
/*
/* Revision 1.4  1999/11/02 23:55:58  localpsi
/* Shawn Brown - (11/2/99) Modified to the code in a few major ways.
/*
/* 1.  Added the capability to do UHF.  All of the features available with the
/* other refrences have been added for UHF.
/*
/* 2.  For UHF, I had to alter the structure of file30. (See cleanup.c for a
/* map)  This entailed adding a pointer array right after the header in the SCF
/* section of file30 that pointed to all of the data for the SCF caclulation.
/* Functions were added to libfile30 to account for this and they are
/* incorporated in this code.
/*
/* 3.  Updated and fixed all of the problems associated with my previous
/* guessing code.  The code no longer uses OPENTYPE to specify the type of
/* occupation.  The keword REFERENCE and MULTP can now be used to indicate any
/* type of calculation.  (e.g. ROHF with MULTP of 1 is an open shell singlet
/* ROHF calculation)  This code was moved to occ_fun.c.  The code can also
/* guess at any multplicity in a highspin case, provided enough electrons.
/*
/* Revision 1.3  1999/08/17 19:04:16  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.2  1999/07/24 18:13:52  crawdad
/* Renamed variable "nint" to "cscf_nint" to avoid DEC compiler type conflict.
/* -Daniel
/*
 * Revision 1.1.1.1  1999/04/12  16:59:27  evaleev
 * Added a version of CSCF that can work with CINTS.
 * -Ed
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

int last;

extern double *pa, *pb;
static int ibl=0;
static double *gtmp,*ptmp;
struct c_pkints {
         int ij;
         int kl;
         double pval;
         } *c_outbuf;

void packit_closed(lbij,lbkl,endflg)
   unsigned int *lbij, *lbkl;
   int endflg;

{
   int i,j,k,joff,ij,kl;
   int lmax,ijkl,l;
   int tmpsiz,nn;
   double pval;
   double tol = 10e-14;
   
   if(!c_outbuf) {
      if((c_outbuf=(struct c_pkints *) malloc(maxbuf*sizeof(struct c_pkints)))
                                                                   ==NULL) {
         fprintf(stderr,"cannot allocate memory for c_outbuf in packit\n");
         exit(4);
         }
      }

   tmpsiz = ioff[nbasis];

   if(gtmp==NULL) {
      gtmp = (double *) init_array(tmpsiz);
      ptmp = (double *) init_array(tmpsiz);
      /*testpk = (double *) init_array(ioff[tmpsiz]);*/
      
      for(k=joff=0; k < num_ir ; k++) {
         if(nn=scf_info[k].num_so) {
            for(i=0; i < nn ; i++)
               for(j=0; j <= i ; j++)
                  ptmp[ioff[i+joff]+j+joff] = scf_info[k].pmat[ioff[i]+j];
            joff += nn;
            }
         }
      }

   if(!endflg) {
      for(i=0; i < cscf_nint ; i++) {
         pval=pa[i];
         ij = lbij[i];
         kl = lbkl[i];
         if(print & 128) fprintf(outfile,"%5d%5d%9.5f\n",ij,kl,pval);
         if (fabs(pval) >= tol) {
            c_outbuf[ibl].ij = ij;
            c_outbuf[ibl].kl = kl;
            c_outbuf[ibl].pval = pval;

            gtmp[ij] += ptmp[kl]*pval;
            gtmp[kl] += ptmp[ij]*pval;

	    /*testpk[ioff[ij]+kl] = pval;*/

            ibl++;

            if (ibl >= maxbuf) {
               if(readflg)
                  swrit(itap92,(char *) c_outbuf,sizeof(struct c_pkints)*maxbuf);
               num_ints += ibl;
               if(print & 16)
                  fprintf(outfile,"buf %3d: ibl = %10d\n",num_bufs,ibl);
               fflush(outfile);
               num_bufs++;
               ibl=0;
               }
            }
         }
      cscf_nint=0;
      }
   /* testing stuff */
   else {
       /*for(i=0;i<nbasis;i++){
	   for(j=0;j<=i;j++){
	       for(k=0;k<=i;k++){
		   lmax = (k==i) ? j : k;
		   for(l=0;l<=lmax;l++){
		       ijkl = ioff[ioff[i]+j]+ioff[k]+l;
		       fprintf(JK,"\n%5d%5d%5d%5d%5d%5d%5d\t%10.10lf",
			       i,j,k,l,ioff[i]+j,ioff[k]+l
			       ,ijkl,testpk[ijkl]);
		   }
	       }
	   }
	   }*/
       
       num_ints += ibl;
       if(print & 16) fprintf(outfile,"buf %3d: ibl = %10d\n",num_bufs,ibl);
       fflush(outfile);
       num_bufs++;
       last = ibl;
       if(readflg) swrit(itap92,(char *) c_outbuf,sizeof(struct c_pkints)*ibl);
       
      for(k=joff=0; k < num_ir ; k++) {
         if(nn=scf_info[k].num_so) {
            for(i=0; i < nn ; i++)
               for(j=0; j <= i ; j++)
                  scf_info[k].gmat[ioff[i]+j] = gtmp[ioff[i+joff]+j+joff];
            joff += nn;
            }
         }
      free(gtmp);
      free(ptmp);
      }
   }
