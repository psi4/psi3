/* $Log$
 * Revision 1.5  2000/08/23 17:15:16  sbrown
 * Added portions to separate out the correlation and exchange energy at the
 * end the calculation as well as do the consistency check on the integrated
 * density.
 *
/* Revision 1.4  2000/06/26 19:04:09  sbrown
/* Added DFT capapbilities to interface with cints using direct scf
/*
/* Revision 1.3  2000/06/22 22:15:00  evaleev
/* Modifications for KS DFT. Reading in XC Fock matrices and XC energy in formg_direct need to be uncommented (at present those are not produced by CINTS yet).
/*
/* Revision 1.2  2000/06/02 13:32:16  kenny
/*
/*
/* Added dynamic integral accuracy cutoffs for direct scf.  Added a few global
/* variables.  Added keyword 'dyn_acc'; true--use dynamic cutoffs.  Use of
/* 'dconv' and 'delta' to keep track of density convergence somewhat awkward,
/* but avoids problems when accuracy is switched and we have to wipe out density
/* matrices.  Also added error message and exit if direct rohf singlet is
/* attempted since it doesn't work.
/* --Joe Kenny
/*
/* Revision 1.1.1.1  2000/02/04 22:52:30  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.3  1999/11/02 23:55:56  localpsi
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
/* Revision 1.2  1999/08/17 19:04:14  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:26  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

static double twocut=1.0;
static double eelec;       /*--- elec. energy from the previous iteration ---*/
double dconv;

int ecalc(incr)
   double incr;

{
   int i,j,k,ij,nn;
   double edif;
   double plimit = pow(10.0,(double) -iconv);
   double neelec = 0.0;
   double ir_energy, dtmp, dtmp1;
   double cinext;
   struct symm *s;

   delta=0.0;
   for (k=0; k < num_ir ; k++) {
      s = &scf_info[k];
      if (nn=s->num_so) {
         ir_energy = 0.0;
         for (i=ij=0; i < nn ; i++) {
            for (j = 0 ; j <= i ; j++,ij++) {
               if(uhf) {
		    ir_energy += 0.5*((s->pmat[ij]*s->hmat[ij])
				      +(spin_info[0].scf_spin[k].pmat[ij]
					*spin_info[0].scf_spin[k].fock_pac[ij])
				      +(spin_info[1].scf_spin[k].pmat[ij]
					*spin_info[1].scf_spin[k].fock_pac[ij]));
	       }
	       else if(!iopen) {
		   ir_energy += 0.5*s->pmat[ij]*(s->hmat[ij]+s->fock_pac[ij]);
	       }
	       else {
		   ir_energy += 0.5*s->pmat[ij]*(s->hmat[ij]+s->fock_pac[ij])
		       - 0.5*s->pmato[ij]*s->gmato[ij];
	       }
	    }
	 }
         neelec += ir_energy;
         if (iter) {
	     if(uhf){
		 for (i = 0; i < ioff[nn] ; i++) {
		     dtmp = spin_info[0].scf_spin[k].dpmat[i];
		     dtmp1 = spin_info[1].scf_spin[k].dpmat[i];
		     delta += dtmp*dtmp;
		     delta += dtmp1*dtmp1;
		 }
	     }
	     else {
		 for (i = 0; i < ioff[nn] ; i++) {
		     dtmp = s->dpmat[i];
		     delta += dtmp*dtmp;
		 }
	     }
	 }
      }
   }
 
   /*JPK(6/1/00) dynamic integral accuracy modifications*/
   dconv = sqrt(delta)/mxcoef2;
   delta = dconv;
   if(acc_switch==1 || iter==0) {
      delta=1.0;
      acc_switch=0;
   }                            

   if (ksdft){
       coulomb_energy = neelec;
       neelec += exc;
   /*printf("XC_energy = %10.10lf",exc);*/
   }
   
   etot = repnuc + neelec;
   edif =  eelec - neelec;
   ediff = edif;

   if (!iter) fprintf(outfile,"\n  iter       total energy        delta E         delta P          diiser\n");
   fprintf(outfile, "%5d %20.10f %15.6e %15.6e %15.6e\n", 
                                    ++iter, etot, edif, dconv, diiser);
   fflush(outfile);
   diiser=0.0;

   if ( delta < plimit && iter > 1) {
      converged=1;
      if(!iopen || iopen && fock_typ >= 2) cleanup();
      }

   eelec = neelec;

   cinext = pow(10.0,-twocut);
   if (delta < cinext && delta && !converged) {
      twocut += incr;
      return(1);
      }
   else return(0);
   }
