/* $Log$
 * Revision 1.1  2000/02/04 22:52:34  evaleev
 * Initial revision
 *
/* Revision 1.3  1999/11/17 19:40:47  evaleev
/* Made all the adjustments necessary to have direct UHF working. Still doesn't work though..
/*
/* Revision 1.2  1999/11/04 19:24:31  localpsi
/* STB (11/4/99) - Added the orb_mix feature which is equivalent to guess = mix
/* in G94 and also fixed restarting so that if you have different wavefuntions,
/* everything works.  Also if you specify no DOCC and SOCC and restart, if the
/* wavefunctions are different, it will guess again.
/*
/* Revision 1.1  1999/11/02 23:56:00  localpsi
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
/* Revision 1.1.1.1  1999/04/12 16:59:28  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
/*
 * Revision 1.5  1998/06/30  14:11:12  sbrown
 * *************************************************************
 * *Program Modification                                       *
 * *By: Shawn Brown                                            *
 * *Date: June 30, 1998                                        *
 * *Altered program to make a guess at occupations from the    *
 * *diagonalized core hamiltonian matrix.  Program can now     *
 * *make a guess at the beginning of the calculation or at     *
 * *or at every iteration.  Use the latter at your own risk.   *
 * *See man pages for details on new keywords.                 *
 * *************************************************************
 *
 * Revision 1.4  1995/07/21  17:37:15  psi
 * Made Jan 1st 1995 cscf the current accepted version of cscf.  Some
 * unidentified changes made after that date were causing problems.
 *
 * Revision 1.1  1991/06/15  20:22:40  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

void uhf_iter()

{
   int i,j,l,m,t,ij;
   int nn,newci;
   double cimax;
   double **scr;
   double **fock_c;
   double **fock_ct;
   double **ctrans;
   double tol = 1.0e-14;
   struct symm *s;
   struct spin *sp;

   diiser = 0.0;
   scr = (double **) init_matrix(nsfmax,nsfmax);
   fock_c = (double **) init_matrix(nsfmax,nsfmax);
   fock_ct = (double **) init_matrix(nsfmax,nsfmax);
   ctrans = (double **) init_matrix(nsfmax,nsfmax);
   
/* and iterate */

   for (iter=0; iter < itmax ; ) {
       for(t = 0; t<2 ;t++){
	   sp = &spin_info[t];

	   if(print & 4) {
	     for(m=0; m < num_ir ; m++) {
	       if (nn=scf_info[m].num_so) {
		 fprintf(outfile,
			 "\n%s gmat for irrep %s",sp->spinlabel,scf_info[m].irrep_label);
		 print_array(sp->scf_spin[m].gmat,nn,outfile);
	       }
	     }
	   }
	   
	   for (m=0; m < num_ir ; m++) {
	       s = &scf_info[m];
	       
	       if (nn=s->num_so) {
		   
		   /*  form fock matrix = h+g */
		   add_arr(s->hmat,sp->scf_spin[m].gmat,
			   sp->scf_spin[m].fock_pac,ioff[nn]);
	       }
	   }
	   
	   if(t==1) 
	       ecalc(tol);
	   
	   
	   
	   /* create new fock matrix in fock_pac or fock_eff */
	   if(t==1)
	       if(!diisflg) diis_uhf(scr,fock_c,fock_ct);
       }
       
       for(t=0;t<2;t++){
	   sp = &spin_info[t];
	   for (m=0; m < num_ir ; m++) {
	       s = &scf_info[m];
	       if (nn=s->num_so) {
		   /* transform fock_pac to mo basis */
		   tri_to_sq(sp->scf_spin[m].fock_pac,fock_ct,nn);
		   mxmb(sp->scf_spin[m].cmat,nn,1
			,fock_ct,1,nn,scr,1,nn,nn,nn,nn);
		   mxmb(scr,1,nn,sp->scf_spin[m].cmat,1
			,nn,fock_c,1,nn,nn,nn,nn);
		   
		   /*  diagonalize fock_c to get ctrans */
		   sq_rsp(nn,nn,fock_c,sp->scf_spin[m].fock_evals
			  ,1,ctrans,tol);
		   
		   if(print & 4) {
		       fprintf(outfile,"\n %s eigenvector for irrep %s\n"
			       ,sp->spinlabel,s->irrep_label);
		       eivout(ctrans,sp->scf_spin[m].fock_evals,nn,nn,outfile);
		   }
		   
		   mxmb(sp->scf_spin[m].cmat,1,nn,
			ctrans,1,nn,scr,1,nn,nn,nn,nn);
				   
		   if(print & 4) {
		       fprintf(outfile,"\n %s eigenvector after irrep %s\n",
			       sp->spinlabel,s->irrep_label);
		       print_mat(scr,nn,nn,outfile);
		   }
		   
		   for (i=0; i < nn; i++)
		       for (j=0; j < nn; j++)
			   sp->scf_spin[m].cmat[i][j] = scr[i][j];
               }
	   }
	   
	   if(converged) {
	       free_matrix(scr,nsfmax);
	       free_matrix(fock_c,nsfmax);
	       free_matrix(fock_ct,nsfmax);
	       free_matrix(ctrans,nsfmax);
	       exit(1);
	       cleanup();
	   }
       }
       schmit_uhf(1);
       
       if(print & 4) {
	   for(j=0; j < 2; j++){
	       for(i=0; i < num_ir ; i++) {
		   s = &scf_info[i];
		   if (nn=s->num_so) {
		       fprintf(outfile,"\northogonalized mos irrep %s\n",
			       s->irrep_label);
		       print_mat(spin_info[j].scf_spin[i].cmat,nn,nn,outfile);
		   }
	       }
	   }
       }
       
       if(mixing && iter ==1)
	   orb_mix();

       /* form new density matrix */
       dmatuhf();
       
       /* and form new fock matrix */
       if(iter < itmax) {
	   if (!direct_scf)
	       formg_open();
	   else
	       formg_direct();
       }
   }
}





