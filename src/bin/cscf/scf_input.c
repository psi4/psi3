/* $Log$
 * Revision 1.10  2001/01/04 14:13:35  sbrown
 * Fixed the problem with iconv:  The new versions of linux had iconv already
 * assigned to something else so I changed all references of it to scf_conv.
 *
/* Revision 1.9  2000/10/13 19:51:22  evaleev
/* Cleaned up a lot of stuff in order to get CSCF working with the new "Mo-projection-capable" INPUT.
/*
/* Revision 1.8  2000/08/23 17:15:16  sbrown
/* Added portions to separate out the correlation and exchange energy at the
/* end the calculation as well as do the consistency check on the integrated
/* density.
/*
/* Revision 1.7  2000/08/21 00:28:58  sbrown
/* Included dft_inputs.c and now cscf has information about both the dft
/* functionals and grids.
/*
/* Revision 1.6  2000/07/10 18:03:33  sbrown
/* Enabling cscf to send over just the occupied SCF eigenvector for DFT
/* calculations.  Only done for the RHF case.
/*
/* Revision 1.5  2000/06/27 21:08:10  evaleev
/* Fixed a minor string manipulation problem in scf_input.c
/*
/* Revision 1.4  2000/06/26 19:04:11  sbrown
/* Added DFT capapbilities to interface with cints using direct scf
/*
/* Revision 1.3  2000/06/22 22:15:02  evaleev
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
/* Revision 1.1.1.1  2000/02/04 22:52:32  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.8  1999/11/17 19:40:47  evaleev
/* Made all the adjustments necessary to have direct UHF working. Still doesn't work though..
/*
/* Revision 1.7  1999/11/11 21:15:15  localpsi
/* Altered cscf to do some guess at the multiplicity from SOCC. -STB (11/11/99)
/*
/* OH and in case your wondering who localpsi is, it is the superuser on my pc
/* that contains my psi files.
/*
/* Revision 1.6  1999/11/04 19:24:30  localpsi
/* STB (11/4/99) - Added the orb_mix feature which is equivalent to guess = mix
/* in G94 and also fixed restarting so that if you have different wavefuntions,
/* everything works.  Also if you specify no DOCC and SOCC and restart, if the
/* wavefunctions are different, it will guess again.
/*
/* Revision 1.5  1999/11/02 23:55:59  localpsi
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
/* Revision 1.4  1999/11/02 18:10:14  evaleev
/* Direct SCF improved
/*
/* Revision 1.3  1999/10/22 19:47:19  evaleev
/* A direct SCF-enabled version (set DIRECT_SCF=TRUE in input.dat).
/*
/* Revision 1.2  1999/08/17 19:04:17  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:28  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <ip_libv1.h>
#include <file30.h>

void scf_input(ipvalue)
   ip_value_t *ipvalue;
{
   int i,j,k,ijk,m;
   double elast;     
   double **scr_mat;
   char *alabel,*bool="YES",*optyp,*wfn,*dertype,*guess;
   char *grid_str;
   char cjunk[80];
   int norder,*iorder,reordr;
   int nc,no,nh,nn,num_mo;
   int mpoint,mconst,mcalcs,ncalcs;
   PSI_FPTR junk,locvec,loccal;
   int optri,ierr,nat;
   int io_locate();
   int errcod;
   int size;
   int phase_chk;
   int mo_offset, so_offset;
   struct symm *s;
   reftype reftmp;
   int depth;

   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":SCF");
  
   if(ipvalue) ip_print_value(stdout,ipvalue);

   errcod = ip_string("LABEL",&alabel,0);
   if(errcod == IPE_OK) fprintf(outfile,"  label       = %s\n",alabel);

   direct_scf = 0;
   errcod = ip_boolean("DIRECT_SCF",&direct_scf,0);
   /* Can do KS DFT direct only */
   if (ksdft) direct_scf=1;

   mixing = 0;
   errcod = ip_boolean("ORB_MIX",&mixing,0);

   /*-----------------------------------------------------
     Which type of guess to use. Sets inflg to:
     0 (AUTO,default) - check if there's an old vector
                        in file30, set inflg to 1 on yes,
		        to 2 otherwise.
     1                - use old vector in file30
     2 (GUESS=CORE)   - use core guess
    -----------------------------------------------------*/
   guess = strdup("AUTO");
   errcod = ip_string("GUESS",&guess,0);
   if (!strcmp(guess,"AUTO"))
       inflg=0;
   else if (!strcmp(guess,"CORE"))
       inflg=2;

   reordr = 0;
   norder = 0;
   errcod = ip_boolean("REORDER",&reordr,0);
   if(reordr) {
     norder = 1;

     errcod = ip_count("MOORDER",&size,0);
     errchk(errcod,"MOORDER");
     if(errcod != IPE_OK) {
       fprintf(outfile,"\ncannot find MOORDER. calculation continuing\n");
       norder=0;
       reordr=0;
       }
     else {
       if(size != nbasis) {
         fprintf(outfile,"\n you have not given enough mos to MOORDER\n");
         exit(size);
         }
       iorder = (int *) malloc(sizeof(int)*size);
       for(i=0; i < size ; i++) {
         errcod = ip_data("MOORDER","%d",&iorder[i],1,i);
         errchk(errcod,"MOORDER");
         }
       }
     }
   
   /* Remove after debugging.  Stop cscf right before going to cints */
   exitflag = 0;
   errcod = ip_boolean("EXIT_CINTS",&exitflag,0);
   
   itmax = 40;
   errcod = ip_data("MAXITER","%d",&itmax,0);

   it_diis = 0;
   errcod = ip_data("DIISSTART","%d",&it_diis,0);
   
   print = 0;
   errcod = ip_data("IPRINT","%d",&print,0);

   fock_typ = 0;
   errcod = ip_data("FOCK_TYPE","%d",&fock_typ,0);

   second_root = 0;
   if (twocon) {
       errcod = ip_boolean("SECOND_ROOT",&second_root,0);
     }
   
   icheck_rot = 1;
   errcod = ip_boolean("CHECK_ROT",&icheck_rot,0);

   ndiis = (iopen) ? 4 : 6;
   if(twocon) ndiis = 3;
   errcod = ip_data("NDIIS","%d",&ndiis,0);

   if(ipvalue) ip_print_tree(stdout,NULL);

   scf_conv = 7;
   if(ipvalue) ip_print_value(stdout,ipvalue);
   errcod = ip_string("WFN",&wfn,0);
   if(ipvalue) ip_print_value(stdout,ipvalue);
   errcod = ip_string("DERTYPE",&dertype,0);
   if(errcod == IPE_KEY_NOT_FOUND) {
     dertype = (char *) malloc(sizeof(char)*5);
     strcpy(dertype,"NONE");
     }
   if(strcmp(wfn,"SCF")) scf_conv = 10;
   if(!strcmp(dertype,"SECOND")) scf_conv = 12;
   errcod = ip_data("CONVERGENCE","%d",&scf_conv,0);

   if (ksdft){
       functional = (char *)determine_functional();
       grid_str = (char *)determine_grid();
   }
   
   if(ipvalue) ip_print_value(stdout,ipvalue);
   fprintf(outfile,"  wfn          = %s\n",wfn);
   fprintf(outfile,"  reference    = %s\n",reference);
   if (ksdft) {
   fprintf(outfile,"  functional   = %s\n",functional);
   fprintf(outfile,"  DFT grid     = %s\n",grid_str);
   }
   fprintf(outfile,"  multiplicity = %d\n",multp);
   fprintf(outfile,"  charge       = %d\n",charge);
   fprintf(outfile,"  direct SCF   = %s\n",(direct_scf) ? "true" : "false");
   if(direct_scf)
   fprintf(outfile,"  dyn_acc      = %s\n",(dyn_acc) ? "true" : "false");
   fprintf(outfile,"  dertype      = %s\n",dertype);
   fprintf(outfile,"  convergence  = %d\n",scf_conv);
   fprintf(outfile,"  maxiter      = %d\n",itmax);
   fprintf(outfile,"  guess        = %s\n",guess);
   if(print) fprintf(outfile,"  iprint       = %d\n",print);
   if (second_root)
     fprintf(outfile,"  second_root = TRUE\n");


   diisflg = 0;
   errcod = ip_string("DIIS",&bool,0);
   if(errcod == IPE_OK) {
      if(!strcmp(bool,"NO") || !strcmp(bool,"FALSE") || !strcmp(bool,"0"))
        diisflg=1;
      else diisflg=0;
      }

   fprintf (outfile,"\n  nuclear repulsion energy %22.13f\n",repnuc);
   fflush(outfile);

   mpoint = MPOINT;
   mconst = MCONST;
   mcalcs = MCALCS;
   nat    = file30_rd_natom();
   ncalcs = file30_rd_ncalcs();

/* if inflg is 0 and this isn't the first calc, then get the old vector */
/* from file30.  if inflg is 2, just use core hamiltonian guess */
/* if inflg is 1, get old vector no matter what ncalcs is */

   if ((inflg==0 && ncalcs) || inflg == 1) {
       
       inflg = 1;
       mxcoef = file30_rd_mxcoef();
       optri = abs(file30_rd_iopen());
       reftmp = file30_rd_ref();
       
       fprintf(outfile,"\n  using old vector from file30 as initial guess\n");
       
/* get old energy from file30 */
       
       elast = file30_rd_escf();
       fprintf(outfile,"  energy from old vector: %14.8f\n",elast);
       
       so_offset = 0;
       mo_offset = 0;

/* ----------------------------------------------------
** This is the UHF part of the restarting algorithm
** STB (10/29/99)
**
**----------------------------------------------------*/
       
       if(uhf){
	   
	   /* if the reference is not UHF, then just read in the vector for the 
	      restricted calculation into both */
	   if(reftmp != ref_uhf && reftmp != ref_uks){
	       
	       for(k=0; k < num_ir ; k++) {
		   s = &scf_info[k];
		   if(nn=s->num_so) {
		       spin_info[0].scf_spin[k].cmat = file30_rd_blk_scf(k);
		       spin_info[1].scf_spin[k].cmat = file30_rd_blk_scf(k);
		   }
	       }
	   }
	   else{
	       for(k=0; k < num_ir ; k++) {
		   s = &scf_info[k];
		   if(nn=s->num_so) {
		       spin_info[0].scf_spin[k].cmat = file30_rd_alpha_blk_scf(k);
		       spin_info[1].scf_spin[k].cmat = file30_rd_beta_blk_scf(k);
		   }
	       }
	   }
	   
	   for(m=0;m<2;m++){
	       for(k=0; k < num_ir; k++) {
		   s = &scf_info[k];
		   if(nn=s->num_so) {
		       for(j=0; j < nn; j++)
			   for(i=0; i < nn; i++) 
			       spin_info[m].scf_spin[k].cmat_orig[i][j] 
				   = spin_info[m].scf_spin[k].cmat[i][j];
		   }
	       }
	   }
	   phase_check = 1;
	   
/* reorder vector if norder = 1 */
	   
	   if (norder) {
	       int loff = 0;
	       int jnew;
	       
	       /* TDC(6/19/96) - If the vector is re-ordered, don't allow
		  phase_checking */
	       phase_check = 0;
	       
	       fprintf(outfile,"\n  mo's will be reordered\n");
	       for (m=0;m<2;m++){
		   for (i=0; i < num_ir; i++) {
		       s = &scf_info[i];
		       if (nn=s->num_so) {
			   scr_mat = (double **) init_matrix(nn,nn);
			   for (j=0; j < nn; j++) {
			       jnew = iorder[j+loff]-1;
			       for (k=0; k < nn ; k++) {
				   scr_mat[k][j]
				       =spin_info[m].scf_spin[i].cmat[k][jnew];
			       }
			   }
			   for (j=0; j < nn ; j++)
			       for (k=0; k < nn ; k++) 
				   spin_info[m].scf_spin[i].cmat[j][k] = scr_mat[j][k];
			   
			   fprintf(outfile,"\n reordered %s mo's for irrep %s\n",
				   spin_info[m].spinlabel,s->irrep_label);
			   print_mat(spin_info[m].scf_spin[i].cmat,nn,nn,outfile);
			   loff += nn;
			   free_matrix(scr_mat,nn);
		       }
		   }
	       }
	       free(iorder);
	   }
       }
   
       else{
	   for(k=0; k < num_ir ; k++) {
	       s = &scf_info[k];
	       if(nn=s->num_so) {
		   s->cmat = file30_rd_blk_scf(k);
	       }
	   }
	   
/* TDC(6/19/96) - Make a copy of the vector for later MO phase
   checking and temporarily set the phase_check flag to true */
	   
	   for(k=0; k < num_ir; k++) {
	       s = &scf_info[k];
	       if(nn=s->num_so) {
		   for(j=0; j < nn; j++)
		       for(i=0; i < nn; i++) s->cmat_orig[i][j] = s->cmat[i][j];
	       }
	   }
	   
	   phase_check = 1;
	   
/* reorder vector if norder = 1 */
	   
	   if (norder) {
	       int loff = 0;
	       int jnew;
	       
	       /* TDC(6/19/96) - If the vector is re-ordered, don't allow
		  phase_checking */
	       phase_check = 0;
	       
	       fprintf(outfile,"\n  mo's will be reordered\n");
	       for (i=0; i < num_ir; i++) {
		   s = &scf_info[i];
		   if (nn=s->num_so) {
		       scr_mat = (double **) init_matrix(nn,nn);
		       for (j=0; j < nn; j++) {
			   jnew = iorder[j+loff]-1;
			   for (k=0; k < nn ; k++) {
			       scr_mat[k][j]=s->cmat[k][jnew];
			   }
		       }
		       for (j=0; j < nn ; j++)
			   for (k=0; k < nn ; k++) s->cmat[j][k] = scr_mat[j][k];
		       
		       fprintf(outfile,"\n reordered mo's for irrep %s\n",
			       s->irrep_label);
		       print_mat(s->cmat,nn,nn,outfile);
		       loff += nn;
		       free_matrix(scr_mat,nn);
		   }
	       }
	       free(iorder);
	   }
       }
   }
   else {
       inflg = 2;
       fprintf(outfile,"  first run, so defaulting to core-hamiltonian guess\n");
       /* TDC(6/19/96) - If not starting from old vector, don't allow
	  phase checking */
       phase_check = 0;
   }

/* TDC(6/20/96) - Check to see if the user will let us do phase
   correction.  The default has already been set above. */
   phase_chk = 1;
   errcod = ip_boolean("PHASE",&phase_chk,0);
   if(phase_check && phase_chk) phase_check = 1;

/* read in damping factor and level shift */

   dampsv= (iopen) ? 0.02 : 0.0;
   if(twocon) dampsv = 0.01;
   errcod = ip_data("DIISDAMP","%lf",&dampsv,0);

   lshift=1.0;
   errcod = ip_data("LEVELSHIFT","%lf",&lshift,0);
   if(!iopen && fabs(lshift) > 0.0) lshift = 0.1;

   dampd=1.0;
   errcod = ip_data("DAMPD","%lf",&dampd,0);

   dampo=1.0;
   errcod = ip_data("DAMPO","%lf",&dampo,0);

   fprintf(outfile,"\n  level shift                      = %f\n",lshift);
   if(!diisflg) {
      fprintf(outfile,"  diis scale factor                = %f\n",dampsv+1.0);
      fprintf(outfile,"  iterations before extrapolation  = %d\n",it_diis);
      fprintf(outfile,"  %d error matrices will be kept\n",ndiis);
      }
   else fprintf(outfile,"\n  diis turned off\n");

   switch (fock_typ) {
      case 0:
         break;
      case 1:
         fprintf(outfile,"\n  a fock matrix for high spin will be used\n");
         fprintf(outfile,"  this form may not work well with diis\n");
         break;
      default:
         fprintf(outfile,"\n  an experimental fock matrix will be used\n");
         fprintf(outfile,"  the management will not be held responsible for the results\n");
      }


   /* EFV 10/24/98 Check if delete integrals */
   delete_ints = 0;
   if(!strcmp(wfn,"SCF") && (!strcmp(dertype,"FIRST") || !strcmp(dertype,"NONE")))
     delete_ints = 1;
   errcod = ip_boolean("DELETE_INTS",&delete_ints,0);
     /* These keywords will work only with IWL format */
   if (use_iwl) {
     delete_1e = delete_ints;
     errcod = ip_boolean("DELETE_1E",&delete_1e,0);
     delete_2e = delete_ints;
     errcod = ip_boolean("DELETE_2E",&delete_2e,0);
   }
   
   
   fflush(outfile);
}
