/* $Log$
 * Revision 1.1  2000/02/04 22:52:30  evaleev
 * Initial revision
 *
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
static double eelec;
double delta;

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
				      +(spin_info[0].scf_spin[k].pmato[ij]
					*spin_info[0].scf_spin[k].fock_pac[ij])
				      +(spin_info[1].scf_spin[k].pmato[ij]
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
   
   delta = sqrt(delta)/mxcoef2;
   etot = repnuc + neelec;
   edif =  eelec - neelec;
   ediff = edif;

   if (!iter) fprintf(outfile,"\n  iter       total energy        delta E         delta P          diiser\n");
   fprintf(outfile, "%5d %20.10f %15.6e %15.6e %15.6e\n", 
                                    ++iter, etot, edif, delta, diiser);
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
