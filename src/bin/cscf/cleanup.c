/* $Log$
 * Revision 1.1  2000/02/04 22:52:28  evaleev
 * Initial revision
 *
/* Revision 1.7  1999/11/11 21:04:36  evaleev
/* A very minor fix.
/*
 * Revision 1.6  1999/11/02  23:55:54  localpsi
 * Shawn Brown - (11/2/99) Modified to the code in a few major ways.
 *
 * 1.  Added the capability to do UHF.  All of the features available with the
 * other refrences have been added for UHF.
 *
 * 2.  For UHF, I had to alter the structure of file30. (See cleanup.c for a
 * map)  This entailed adding a pointer array right after the header in the SCF
 * section of file30 that pointed to all of the data for the SCF caclulation.
 * Functions were added to libfile30 to account for this and they are
 * incorporated in this code.
 *
 * 3.  Updated and fixed all of the problems associated with my previous
 * guessing code.  The code no longer uses OPENTYPE to specify the type of
 * occupation.  The keword REFERENCE and MULTP can now be used to indicate any
 * type of calculation.  (e.g. ROHF with MULTP of 1 is an open shell singlet
 * ROHF calculation)  This code was moved to occ_fun.c.  The code can also
 * guess at any multplicity in a highspin case, provided enough electrons.
 *
/* Revision 1.5  1999/10/22 19:47:17  evaleev
/* A direct SCF-enabled version (set DIRECT_SCF=TRUE in input.dat).
/*
/* Revision 1.4  1999/10/11 17:03:17  evaleev
/* Modified the location of nmo in mconst array in file 30.
/*
/* Revision 1.3  1999/08/19 14:44:16  evaleev
/* Tiny fix in cleanup().
/*
/* Revision 1.2  1999/08/17 19:04:13  evaleev
/* Changed the default symmetric orthogonalization to the canonical
/* orthogonalization. Now, if near-linear dependencies in the basis are found,
/* eigenvectors of the overlap matrix with eigenvalues less than 1E-6 will be
/* left out. This will lead to num_mo != num_so, i.e. SCF eigenvector is no
/* longer a square matrix. Had to rework some routines in libfile30, and add some.
/* The progrem prints out a warning if near-linear dependencies are found. TRANSQT
/* and a whole bunch of other codes has to be fixed to work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:25  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <ip_libv1.h>

/* TDC(6/20/96) - Prototype for phase() */
int phase(void);
double ssquare(void);

void cleanup()

{
   int i,j,k,ij,ijk,m,nn,num_mo;
   PSI_FPTR junk,junk1,junk2;
   int mpoint,mconst,mcalcs,loccal;
   PSI_FPTR locvec;
   int nx,ntri;
   int newvec;
   int nat,iend,ierr,ci_calc,irot,nbfao;
   int numso;
   int n_there[20],nc[10],no[10];
   int errcod;
   char *ci_type="SCF";
   char *der_type="FIRST";
   double occj,occk;
   double ekin,epot,enpot,ovlp,virial,num_elec,s2;
   double *scr_arr, *lagrangian;
   double **scr1, **scr2;
   double *scrtmp,**scrbtmp;
   double *temp;
   struct symm *s;
   int *pointers;  /* Array to hold pointers to the scf info */
   char **labs;
   
   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":SCF");

   ci_calc=irot=0;
   errcod = ip_string("WFN",&ci_type,0);
   errcod = ip_string("DERTYPE",&der_type,0);

   if(strcmp(ci_type,"SCF")) ci_calc=1;
   if(ci_calc && iopen) irot=1;
   if(strcmp(der_type,"FIRST") && strcmp(der_type,"NONE")) irot=0;

   errcod = ip_boolean("ROTATE",&irot,0);

 /* TDC(6/19/96) - If we're not rotating, check the phases on the MOs,
   and correct them, if possible. */

   if(!irot && !uhf) {
       if(phase_check) phase_check = phase();
       }
   
/* first print mo's, then rotate if this is a ci */

   if(uhf)
       print_mos_uhf();
   else
       print_mos(); 

   if(irot) rotate_vector();
   else 
     if(iopen)
       fprintf(outfile,"\n%8cWFN is %s so no rotation\n",' ',ci_type);

/* TDC(6/20/96) - If we've rotated the orbitals, check the phases on
   the MOs, and correct them, if possible. */

   if(irot) {
       if(phase_check) phase_check = phase();
     }

/* update ncalcs in file30 */
   wreadw(itap30,(char *) i10,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&junk);

   if(!i10[44]) i10[44]++;

   iend = i10[0];
   mpoint = i10[1];
   mconst = i10[2];
   mcalcs = i10[3];
   nat = i10[18];
   nx = nmo*(nmo+1)/2;
   nbfao = i10[21];
   i10[40] = n_so_typs;
   i10[41] = mxcoef;
   i10[42] = ioff[n_open];
   if(twocon) i10[42] = -i10[42];
   i10[45] = nmo;
   
   /* STB(10/28/99) - Flag to tell what reference is being used*/
   i10[46] = refnum;
   
   /* TDC(6/19/96) - Set the phase_check flag here */
   i10[50] = phase_check;
  
   wwritw(itap30,(char *) i10,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&junk);

   /* STB(10/28/99) - Changing the file30 structure to utilize the pointer
      space */
   
   /*  The structure of the SCF section is as follows:
   (listed first is the length in integer words, 
   variable names correspond to those given in CSCF3.0)


   | 20 (header) | 40 (don't know)  |  20 (pointers written    |
   |             |                  |      this function)      |           
   _____________________________________________________________


   |  6*natoms      | 10 ( 5 double      |  mxcoef*2           |
   |  (geometry)    |      contants [a]) | (Alpha SCF eigenvec)|  
   -------------------------------------------------------------

   | mxcoef*2           | nmo*2               | nmo*2              |
   | (Beta SCF eigenvec)| (Alpha SCF eigenval)| (Beta SCF eigenval)|
   -----------------------------------------------------------------

   | num_ir*4               | n_so_typs                  |
   | (Non-zero irrep labels)| (number of orbitals/irrep) |
   -------------------------------------------------------
   
   | n_so_typs                   | n_so_typs                 |
   | (number closed shells/irrep)| (number open shells/irrep)|
   -----------------------------------------------------------

   | ioff[n_open]*2 | ioff[n_open]*2 | nmo*(nmo+1)/2*2   |
   | (Alpha coefs)  | (Beta coefs)   | (Alpha Lagrangian)|
   -------------------------------------------------------

   | nmo*(nmo+1)/2*2  |
   | (Beta Lagrangian)|
   --------------------*/

   /*--------------------------------------------------------------*/
   /* Pointers will be an array that will contain the positions in */
   /*   file30 for the following informaiton                       */
   /* [0] = Alpha SCF eigenvector (the only one stored if not UHF) */
   /* [1] = Beta  SCF eigenvector                                  */
   /* [2] = Alpha SCF eigenvalues (the only one stored if not UHF) */
   /* [3] = Beta  SCF eigenvalues                                  */
   /* [4] = Irrep Labels for irreps with non-zero occupancies      */
   /* [5] = Number of orbitals per irrep for non-zero occupancies  */
   /* [6] = Number of closed shells per irrep with non-zero occ    */
   /* [7] = Number of open shells per irrep with non-zero occ      */
   /* [8] = Alpha Coupling Constants                               */
   /* [9] = Beta  Coupling Constants                               */
   /* [10] = Alpha Lagrangian (the only one stored if not UHF)     */
   /* [11] = Beta  Lagrangian                                      */
   /*--------------------------------------------------------------*/
  
   pointers = init_int_array(20);
   
   /* get pointers to calculations */

   junk = (PSI_FPTR) sizeof(int)*(100+mconst+mpoint);
   wreadw(itap30,(char *) i10,sizeof(int)*mcalcs,junk,&junk);
   
   /* loccal holds the value that points to the beginning of the 
      SCF section of file30 */
   loccal = (PSI_FPTR) sizeof(int)*(i10[0]-1);

   /* pointer to scf vector is first in calculation pointers */

   junk = loccal+sizeof(int)*60; /* skip the header */
   
   /* Let's skip the pointer section for now */
   
   wreadw(itap30,(char *) i10,sizeof(int)*20,junk,&junk1);
   
   newvec = 0;
   locvec = (PSI_FPTR) sizeof(int)*(i10[0]-1);
   if (i10[0] <= 0) {
      locvec = (PSI_FPTR) sizeof(int)*(iend-1);
      i10[0] = iend;
      wwritw(itap30,(char *) i10,sizeof(int)*20,junk,&junk1);
      newvec = 1;
      }
   
   junk = loccal + ((PSI_FPTR) sizeof(int)*(60+20+6*nat+2));
   
   /* Set up the pointer array now that we know where everything is */
   
   pointers[0] = locvec;
   
/* write new energy to file30 */

   wwritw(itap30,(char *) &etot,sizeof(double)*1,junk,&junk);

/* write new vector and eigenvalues to file30 and file49 */
   scr_arr = (double *) init_array(mxcoef);
   scrtmp = (double *) init_array(mxcoef);
   scrbtmp = (double **) init_matrix(scf_info[0].num_so,scf_info[0].num_so);

   if(uhf){
       for(m=0;m<2;m++){
	   for(k=ijk=0; k < num_ir ; k++) {
	       s= &scf_info[k];
	       if(nn=s->num_so) {
		   num_mo = s->num_mo;
		   for(j=0; j < num_mo ; j++)
		       for(i=0; i < nn ; i++,ijk++){
			   scr_arr[ijk]=spin_info[m].scf_spin[k].cmat[i][j]; 
		       }
	       }
	   }
	   wwritw(itap30,(char *) scr_arr,sizeof(double)*mxcoef,locvec,&locvec);
	   pointers[m+1] = locvec;
       }
   }
   else{
       for(k=ijk=0; k < num_ir ; k++) {
	   s= &scf_info[k];
	   if(nn=s->num_so) {
	       num_mo = s->num_mo;
	       for(j=0; j < num_mo ; j++)
		   for(i=0; i < nn ; i++,ijk++)
		       scr_arr[ijk]=s->cmat[i][j];
	   }
       }
       
       wwritw(itap30,(char *) scr_arr,sizeof(double)*mxcoef,locvec,&locvec);
       pointers[1] = 0;
       pointers[2] = locvec;
   }
   
   if(uhf){
       for(m=0;m<2;m++){
	   for (i=k=ijk=0; i < num_ir ; i++) {
	       s= &scf_info[i];
	       if (nn=s->num_so) {
		   num_mo = s->num_mo;
		   for (j=0; j < num_mo ; j++,k++){
		       scr_arr[k] = spin_info[m].scf_spin[i].fock_evals[j];
		   }
		   ijk++;
	       }
	   }
	   wwritw(itap30,(char *) scr_arr,sizeof(double)*nmo,locvec,&locvec);
	   pointers[m+3] = locvec;
       }  
   }
   else{
       for (i=k=ijk=0; i < num_ir ; i++) {
	   s= &scf_info[i];
	   if (nn=s->num_so) {
	       num_mo = s->num_mo;
	       for (j=0; j < num_mo ; j++,k++){
		   scr_arr[k] = s->fock_evals[j];
	       }
	       ijk++;
	   }
       }
       wwritw(itap30,(char *) scr_arr,sizeof(double)*nmo,locvec,&locvec);
       pointers[3] = 0;
       pointers[4] = locvec;
   }

   for (i=k=ijk=0; i < num_ir ; i++) {
       s= &scf_info[i];
       if (nn=s->num_so) {
	   num_mo = s->num_mo;
	   nc[ijk]=s->nclosed;
	   no[ijk]=s->nopen;
	   n_there[ijk]=num_mo;
	   ijk++;
       }
   }
   
   for (i=0 ; i < num_ir ; i++) 
      if (scf_info[i].num_so)
	  wwritw(itap30,(char *) scf_info[i].irrep_label,sizeof(char)*4,locvec,&locvec);
   
   pointers[5] = locvec;
   
   wwritw(itap30,(char *) n_there,sizeof(int)*n_so_typs,locvec,&locvec);
   
   pointers[6] = locvec;
   
   wwritw(itap30,(char *) nc,sizeof(int)*n_so_typs,locvec,&locvec);
   
   pointers[7] = locvec;
   
   if(iopen || uhf) {
      wwritw(itap30,(char *) no,sizeof(int)*n_so_typs,locvec,&locvec);
      pointers[8] = locvec;
   }
   if(iopen){
      if(twocon) {
         double c1i = scf_info[opblk1].occ_num[opshl1];
         double c1ii = scf_info[opblk2].occ_num[opshl2];
         alpha[0] = beta[0] = 1.0 - 2.0/c1i;
         alpha[1] = 1.0;
         beta[1] = 1.0 + 1.0/(save_ci1*save_ci2);
         alpha[2] = beta[2] = 1.0 - 2.0/c1ii;
         }
      if(!twocon) for (i=0; i < ioff[n_open] ; i++) beta[i] = -beta[i];
      if(ci_calc && !twocon)
         for (i=0; i < ioff[n_open] ; i++) beta[i] = -beta[i];
      wwritw(itap30,(char *) alpha,sizeof(double)*ioff[n_open],locvec,&locvec);
      pointers[9] = locvec;
      wwritw(itap30,(char *) beta,sizeof(double)*ioff[n_open],locvec,&locvec);
      pointers[10] = locvec;
      }

/* calculate mo lagrangian and write to file30 */
/* also write mo fock matrices to file49       */
      
   lagrangian = (double *) init_array(nx);
   if(uhf){
       for(m=0;m<2;m++){
	   for (i=0; i < num_ir ; i++) {
	       s = &scf_info[i];
	       if(nn=s->num_so) {
		   num_mo = s->num_mo;
		   bzero(spin_info[m].scf_spin[i].fock_pac,sizeof(double)*ioff[nn]);
		   for (j=0; j < num_mo ; j++)
		       if (spin_info[m].scf_spin[i].occ_num[j] == 1.0)
			   spin_info[m].scf_spin[i].fock_pac[ioff[j]+j] 
			       = 1.0*spin_info[m].scf_spin[i].fock_evals[j];
	       }
	   }
	   
	   for(k=ijk=0; k < num_ir ; k++) {
	       if(nn=scf_info[k].num_so) {
		   num_mo = scf_info[k].num_mo;
		   for(i=0; i < num_mo ; i++)
		       for(j=0; j <= i ; j++)
			   lagrangian[ioff[i+ijk]+j+ijk] 
			       = spin_info[m].scf_spin[k].fock_pac[ioff[i]+j];
		   ijk += num_mo;
	       }
	   }
	   wwritw(itap30,(char *) lagrangian,sizeof(double)*nx,locvec,&locvec); 
	   pointers[m+11] = locvec; 
       }
   }
   else{
       if (iopen) {
	   scr1 = (double **) init_matrix(nsfmax,nsfmax);
	   scr2 = (double **) init_matrix(nsfmax,nsfmax);
	   
	   for (i=0; i < num_ir ; i++) {
	       s = &scf_info[i];
	       if(nn=s->num_so) {
		   num_mo = s->num_mo;
		   
		   tri_to_sq(s->fock_pac,scr2,nn);
/*            mxmb(s->cmat,nn,1,scr2,1,nn,scr1,1,nn,nn,nn,nn);
	      mxmb(scr1,1,nn,s->cmat,1,nn,scr2,1,nn,nn,nn,nn);*/
		   mmult(s->cmat,1,scr2,0,scr1,0,num_mo,nn,nn,0);
		   mmult(scr1,0,s->cmat,0,scr2,0,num_mo,nn,num_mo,0);
		   sq_to_tri(scr2,s->gmat,num_mo);
		   
		   tri_to_sq(s->fock_open,scr2,nn);
/*            mxmb(s->cmat,nn,1,scr2,1,nn,scr1,1,nn,nn,nn,nn);
	      mxmb(scr1,1,nn,s->cmat,1,nn,scr2,1,nn,nn,nn,nn);*/
		   mmult(s->cmat,1,scr2,0,scr1,0,num_mo,nn,nn,0);
		   mmult(scr1,0,s->cmat,0,scr2,0,num_mo,nn,num_mo,0);
		   sq_to_tri(scr2,s->gmato,num_mo);
		   
		   bzero(s->fock_pac,sizeof(double)*ioff[nn]);
		   for (j=ij=0; j < num_mo ; j++) {
		       for (k=0; k <= j ; k++,ij++) {
			   occj = s->occ_num[j];
			   occk = s->occ_num[k];
			   if(!twocon) {
			       if (occj == 2.0 && occk == 2.0 || 
				   occj == 2.0 && occk == 1.0 ||
				   occj == 1.0 && occk == 2.0) 
				   s->fock_pac[ij] = 2.0*s->gmat[ij];
			       else if (occj == 1.0 && occk == 1.0)
				   s->fock_pac[ij] = s->gmato[ij];
			       else
				   s->fock_pac[ij] = 0.0;
			   }
			   else {
			       if (occj == 2.0 && occk || 
				   occk == 2.0 && occj)
				   s->fock_pac[ij] = 2.0*s->gmat[ij];
			       else if (occj && occk)
				   s->fock_pac[ij] = occj*s->gmato[ij];
			       else
				   s->fock_pac[ij] = 0.0;
			   }
		       }
		   }
	       }
	   }
       }
       
       else {
	   for (i=0; i < num_ir ; i++) {
	       s = &scf_info[i];
	       if(nn=s->num_so) {
		   num_mo = s->num_mo;
		   bzero(s->fock_pac,sizeof(double)*ioff[nn]);
		   for (j=0; j < num_mo ; j++)
		       if (s->occ_num[j] == 2.0)
			   s->fock_pac[ioff[j]+j] = 2.0*s->fock_evals[j];
	       }
	   }
       }
   
   for(k=ijk=0; k < num_ir ; k++) {
      if(nn=scf_info[k].num_so) {
	 num_mo = scf_info[k].num_mo;
         for(i=0; i < num_mo ; i++)
            for(j=0; j <= i ; j++)
              lagrangian[ioff[i+ijk]+j+ijk] = scf_info[k].fock_pac[ioff[i]+j];
         ijk += num_mo;
         }
      }
   wwritw(itap30,(char *) lagrangian,sizeof(double)*nx,locvec,&locvec);
   pointers[11]=0;
   pointers[12]=locvec;
   }
   
   if(newvec) {
      iend = (int) locvec/sizeof(int)+1;
      wwritw(itap30,(char *) &iend,sizeof(int)*1,(PSI_FPTR) sizeof(int)*100,&locvec);
      }
   
   /* STB(10/28/99) - Convert and write the pointers to file30*/
      for(i=0;i<20;i++){
	  pointers[i] = pointers[i]/sizeof(int) + 1;
      }
      
      junk = loccal + 60*sizeof(int);
      wwritw(itap30,(char *) pointers,sizeof(int)*20,junk,&junk);
      
      if(ci_calc && iopen && irot) {
	  fprintf(outfile,
		  "\n ci_typ is %s so mo vector will be rotated\n",ci_type);
	  print_mos();
      }
      num_elec = ekin = enpot = ovlp = 0.0;
      for (i=0; i < num_ir ; i++) {
	  s = &scf_info[i];
	  if (nn=s->num_so) {
	      num_mo = s->num_mo;
	      for (j=0; j < ioff[nn] ; j++) {
		  ekin += s->pmat[j]*s->tmat[j];
		  enpot += s->pmat[j]*s->hmat[j];
		  ovlp += s->pmat[j]*s->smat[j];
	      }
	      if(uhf)
		  num_elec += spin_info[1].scf_spin[i].noccup
		      +spin_info[0].scf_spin[i].noccup;
	      else
		  for (j=0; j < num_mo ; j++)
		      num_elec += s->occ_num[j];
	  }
      }
      
      ovlp /= num_elec;
      epot = etot-ekin;
      enpot -= ekin;
      virial = epot/etot;
      
      if(uhf)
	  s2 = ssquare();
/* Print MOs in a special format */
      
      /*if(print & 1){
	print_mos_new();
	}*/
      
      
      fprintf(outfile,"\n%8ctotal energy       = %20.12f\n",' ',etot);
      fprintf(outfile,"%8ckinetic energy     = %20.12f\n",' ',ekin);
      fprintf(outfile,"%8cnuc. attr. energy  = %20.12f\n",' ',enpot);
      fprintf(outfile,"%8celec. rep. energy  = %20.12f\n",' ',epot-enpot);
      fprintf(outfile,"%8cpotential energy   = %20.12f\n",' ',epot);
      fprintf(outfile,"%8cvirial theorem     = %20.12f\n",' ',virial);
      fprintf(outfile,"%8cwavefunction norm  = %20.12f\n",' ',ovlp);
      if(uhf)
	  fprintf(outfile,"%8c<S^2>              = %20.12f\n",' ',s2); 
/*--- Warn user if some basis functions were dropped ---*/
      if (nmo != nbasis) {
	  fprintf(outfile,"\n  Near-linear dependencies in the basis were eliminated.\n");
	  fprintf(outfile,  "  Proceed at your own risk!\n");
      }
      
      file30_close();
      if(!direct_scf){
	  rclose(itap92,4);
	  rclose(itap93,4);
      }
      
      if(!converged)
	  fprintf(outfile,"\n%8cCalculation has not converged!\n",' ');
      
      tstop(outfile);
      ip_done();
      
      if(!converged) exit(3);
      exit(0);

}

print_mos()
{
   int i,nn,num_mo;
   struct symm *s;
   

   for (i=0; i < num_ir ; i++) {
      s = &scf_info[i];
      if (nn=s->num_so) {
	 num_mo = s->num_mo;
         fprintf(outfile,"\n molecular orbitals for irrep %s\n",s->irrep_label);
         eigout(s->cmat,s->fock_evals,s->occ_num,nn,num_mo,outfile);
         }
      }
}

print_mos_uhf()
{
    int i,j,nn;
    struct symm *s;
    
    for(j=0; j<2;j++){
	for (i=0; i < num_ir ; i++) {
	    s = &scf_info[i];
	    if (nn=s->num_so) {
		fprintf(outfile,"\n %s molecular orbitals for irrep %s\n"
			,spin_info[j].spinlabel,s->irrep_label);
		eigout(spin_info[j].scf_spin[i].cmat
		       ,spin_info[j].scf_spin[i].fock_evals
		       ,spin_info[j].scf_spin[i].occ_num,nn,nn,outfile);
	    }
	}
    }
}

/* STB(11/2/99) - This function does not at the moment, hence why it is commented out above*/
print_mos_new() {
  int i;
  int ncl_tot = 0;
  int nop_tot = 0;
  int num_mo = 0;
  
  for(i=0;i < num_ir;i++){
      ncl_tot+= scf_info[i].nclosed;
      nop_tot+= scf_info[i].nopen;
      num_mo += scf_info[i].num_mo;
    }
  sortev();

  
  fprintf(outfile,"\n\n        Molecular Orbitals\n");

  fprintf(outfile,"  Symmetry   OCC         Energy\n");
  fprintf(outfile,"  --------   ---   ------------------\n");

  for(i = 0; i < ncl_tot; i++)
      fprintf(outfile,"     %4s     2      %14.7lf\n",scf_info[symm_tot[i]].irrep_label,ener_tot[i]);
  for(i = ncl_tot; i < ncl_tot + nop_tot; i++)
      fprintf(outfile,"     %4s     1      %14.7lf\n",scf_info[symm_tot[i]].irrep_label,ener_tot[i]);
  fprintf(outfile,"  -----------------------------------\n");
  for(i = ncl_tot + nop_tot; i < num_mo; i++)
      fprintf(outfile,"     %4s     0      %14.7lf\n",scf_info[symm_tot[i]].irrep_label,ener_tot[i]);
  fprintf(outfile,"  --------   ---   ------------------\n\n");

  fflush(outfile);
}

/* Function to calculate the S^2 value for UHF */
/* See Szabo and Ostlund pg. 107 */
/* for the formula               */

double ssquare(void){
    
    int i,j,k,nn,n;
    int na=0;
    int nb=0;
    int nh=0;
    int num_mo;
    double ss=0.0;
    int nm=0;
    double **scr1,**scr2,**S;
    struct symm *s;
    
    scr1 = (double **)init_matrix(nsfmax,nsfmax);
    scr2 = (double **)init_matrix(nsfmax,nsfmax);
    
    /* Calculate the overlap matrix elements */
    
    for(i = 0;i < num_ir;i++){
	
	na += spin_info[0].scf_spin[i].noccup;
	nb += spin_info[1].scf_spin[i].noccup;
	
	s = &scf_info[i];
	if(nn = s->num_so){
	    num_mo = s->num_mo;
	    tri_to_sq(s->smat,scr1,nn);
	    
	    /* Transform the Overlap matrix to the MO basis */
	    
	    mmult(spin_info[0].scf_spin[i].cmat,1
		  ,scr1,0,scr2,0,num_mo,nn,nn,0);
	    mmult(scr2,0,spin_info[1].scf_spin[i].cmat,0
		  ,scr1,0,num_mo,nn,num_mo,0);
	    
	    for(j = 0; j < spin_info[0].scf_spin[i].noccup; j++){
		for(k = 0;k < spin_info[1].scf_spin[i].noccup; k++){
		    ss -= scr1[j][k]*scr1[j][k];
		}
	    }
	}
    }
    
    
      /* Calculate the occupation part of the equation */
    
    nm = (na-nb)/2.0;
    nh = (nm*(nm+1))+nb;
        
    ss += (nm*(nm+1))+nb;
    
    free_matrix(scr1,nsfmax);
    free_matrix(scr2,nsfmax);
    
    return fabs(ss);
}














