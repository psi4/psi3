/* $Log$
 * Revision 1.25  2003/05/06 20:47:22  evaleev
 * CSCF can now find frzvpi from the eigenvalues.
 *
/* Revision 1.24  2003/05/02 15:39:23  evaleev
/* Added ability to figure out the number of frozen doubly occupied orbitals in each irrep.
/*
/* Revision 1.23  2003/04/14 17:25:47  sherrill
/* Change "total energy" to "SCF total energy" to make more explicit for
/* new users.  Yeah, this will probably break some test case perl scripts
/* temporarily :)
/*
/* Revision 1.22  2002/12/22 17:01:14  evaleev
/* Updated cints, cscf, psi3 (probably not complete) and transqt to use psi_start/psi_stop.
/*
/* Revision 1.21  2002/12/06 20:39:08  evaleev
/* Write total SCF energy as reference energy as well.
/*
/* Revision 1.20  2002/12/06 15:50:32  crawdad
/* Changed all exit values to PSI_RETURN_SUCCESS or PSI_RETURN_FAILURE as
/* necessary.  This is new for the PSI3 execution driver.
/* -TDC
/*
/* Revision 1.19  2002/11/24 22:52:17  crawdad
/* Merging the gbye-file30 branch into the main trunk.
/* -TDC
/*
/* Revision 1.18.2.4  2002/11/23 21:54:45  crawdad
/* Removal of mxcoef stuff for chkpt runs.
/* -TDC
/*
/* Revision 1.18.2.3  2002/11/23 19:35:14  sherrill
/* re-institute writing scf vector to file30 if !USE_LIBCHKPT
/*
/* Revision 1.18.2.2  2002/10/01 22:16:28  sherrill
/* Fix a few minor libchkpt/libfile30 things (one define condition was
/* backwards), turn off libfile30 unless we're really using it, and
/* clean up a stupid mistake for the orthog_only option ---CDS
/*
/* Revision 1.18.2.1  2002/07/29 23:08:30  evaleev
/* A major set of changes designed to convert all psi modules to use libchkpt.
/*
/* Revision 1.18  2002/05/30 20:16:49  crawdad
/* Accidentally left dmalloc calls in place.  Fixed.
/* -TDC
/*
/* Revision 1.17  2002/05/30 12:57:08  crawdad
/* Buf fix.  psio_done() was called before chkpt_close().
/* -TDC
/*
/* Revision 1.16  2002/05/07 22:40:56  sherrill
/* Fix missing s=&scf_info[k] in UHF case for evals, fix missing bracket
/* in RHF case.
/*
/* Revision 1.15  2002/04/28 04:34:10  crawdad
/* Finshed initial additions for mirroring old file30 with new PSIF_CHKPT.  I
/* believe that everything cscf wrote to file30 is also written to PSIF_CHKPT.
/* Now ready to start converting other codes to libchkpt.
/* -TDC
/*
/* Revision 1.14  2002/04/27 22:28:48  crawdad
/* More changes to cleanup in preparation for libchkpt conversion.
/* -TDC
/*
/* Revision 1.13  2002/04/27 18:33:20  crawdad
/* Working on changes for new libchkpt code.  Current version does no reading
/* from chkpt yet.
/* -TDC
/*
/* Revision 1.12  2002/03/25 03:16:51  sherrill
/* Changed name of mxcoef keyword to Mxcoef
/*
/* Revision 1.11  2002/03/25 03:05:45  crawdad
/* More additions for new chkpoint file.
/* -TDC
/*
/* Revision 1.10  2002/03/25 02:17:36  janssen
/* Get rid of tmpl.  Use new naming scheme for libipv1 includes.
/*
/* Revision 1.9  2002/03/25 01:07:59  crawdad
/* Some changes to cleanup et al. to write SCF-generated data to both old
/* file30 and new chkpt.
/* -TDC
/*
/* Revision 1.8  2002/03/25 00:02:00  sherrill
/* Add libpsio
/*
/* Revision 1.7  2002/01/04 18:03:24  crawdad
/* Minor change to set phase_check flag to true when starting from a core
/* guess.  This is to allow correlated calculations that might have stopped
/* due to slow convergence to restart.
/* -TDC
/*
/* Revision 1.6  2001/06/29 20:39:27  evaleev
/* Modified cscf to use libpsio to store supermatrix files.
/*
/* Revision 1.5  2001/05/31 01:12:25  sherrill
/* fix up printing orbital eigenvalues, now does TCSCF too!
/*
/* Revision 1.4  2001/04/11 19:31:46  sherrill
/* I removed printing all MO's by default, since this gets pretty big
/* and useless for many of the molecules of interest today.  I added
/* a nice subroutine to print out orbital eigenvalues.  It seems
/* to work for RHF/ROHF/UHF but is probably broken for TCSCF, which
/* I'll try to check later.
/*
/* Revision 1.3  2000/10/13 19:51:19  evaleev
/* Cleaned up a lot of stuff in order to get CSCF working with the new 
/* "Mo-projection-capable" INPUT.
/*
/* Revision 1.2  2000/08/23 17:15:15  sbrown
/* Added portions to separate out the correlation and exchange energy at the
/* end the calculation as well as do the consistency check on the integrated
/* density.
/*
/* Revision 1.1.1.1  2000/02/04 22:52:28  evaleev
/* Started PSI 3 repository
/*
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
/* longer a square matrix. Had to rework some routines in libfile30, and add 
/* some.  The progrem prints out a warning if near-linear dependencies are 
/* found. TRANSQT and a whole bunch of other codes has to be fixed to 
/* work with such basis sets.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:25  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"
#include <libipv1/ip_lib.h>
#if USE_LIBCHKPT
#  include <libchkpt/chkpt.h>
#else
#  include <libfile30/file30.h>
#endif


/* TDC(6/20/96) - Prototype for phase() */
int phase(void);
double ssquare(void);
static int* compute_frzcpi(int);
static int* compute_frzvpi(int);

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
  int mo_print;
  int errcod;
  char *ci_type="SCF";
  char *der_type="FIRST";
  double occj,occk;
  double ekin,epot,enpot,ovlp,virial,num_elec,s2;
  double *scr_arr, *lagrangian, **lagr, **ccvecs;
  double **scr1, **scr2;
  double *temp;
  struct symm *s;
  int *pointers;  /* Array to hold pointers to the scf info */
  char **labs;
  void print_mo_eigvals(void);
  psio_address chkptr;
  int tmp_iopen;
  int row, col;
  int nfzc, nfzv, *frzcpi, *frzvpi;

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
  mo_print = 0;
  errcod = ip_boolean("PRINT_MOS",&mo_print,0);

  /* TDC(6/19/96) - If we're not rotating, check the phases on the MOs,
     and correct them, if possible. */

  if(!irot && !uhf) {
    if(phase_check) phase_check = phase();
  }
   
  /* first print mo's, then rotate if this is a ci */

  if (mo_print) {
    if(uhf)
      print_mos_uhf();
    else
      print_mos(); 
  }

  if(irot) rotate_vector();
  else 
    if(iopen)
      fprintf(outfile,"\n%8cWFN is %s so no rotation\n",' ',ci_type);

  /* TDC(6/20/96) - If we've rotated the orbitals, check the phases on
     the MOs, and correct them, if possible. */

  if(irot) {
    if(phase_check) phase_check = phase();
  }

  /* TDC(1/4/02) If we've started from a core guess, allow a restart of 
     correlated calcs */
  if(inflg==2) phase_check = 1;

#if !USE_LIBCHKPT
  wreadw(itap30,(char *) i10,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&junk);

  if(!i10[44]) i10[44]++;

  iend = i10[0];
  mpoint = i10[1];
  mconst = i10[2];
  mcalcs = i10[3];
  nat = i10[18];
  i10[40] = n_so_typs;
  i10[41] = mxcoef;
  i10[42] = ioff[n_open];
  if(twocon) i10[42] = -i10[42];
  i10[45] = nmo;
#else
  /* psio_write calls for above */
  chkpt_wt_nsymhf(n_so_typs);
  chkpt_wt_nmo(nmo);

  tmp_iopen = ioff[n_open];
  if(twocon) tmp_iopen = -tmp_iopen;
  chkpt_wt_iopen(tmp_iopen);
#endif
  nx = nmo*(nmo+1)/2;
   
  /* STB(10/28/99) - Flag to tell what reference is being used*/
#if !USE_LIBCHKPT
  i10[46] = refnum;
#else
  chkpt_wt_ref(refnum);
#endif
   
  /* TDC(6/19/96) - Set the phase_check flag here */
#if !USE_LIBCHKPT
  i10[50] = phase_check;
#else
  chkpt_wt_phase_check(phase_check);
#endif

#if !USE_LIBCHKPT
  wwritw(itap30,(char *) i10,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&junk);
#endif

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

      | n_so_typs*4            | n_so_typs                  |
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

#if !USE_LIBCHKPT
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
#endif

#if USE_LIBCHKPT
  /* write to new checkpoint file */
  chkpt_wt_etot(etot);
  chkpt_wt_escf(etot);
  chkpt_wt_eref(etot);

  /* These new arrays for PSIF_CHKPT contain data for ALL irreps */
  for(i=0; i < num_ir; i++) {
    s=&scf_info[i];
    n_there[i]=s->num_mo;
    nc[i]=s->nclosed;
    no[i]=s->nopen;
  }
  chkpt_wt_orbspi(n_there);
  chkpt_wt_clsdpi(nc);
  chkpt_wt_openpi(no);
  for(i=0; i < num_ir; i++) n_there[i] = nc[i] = no[i] = 0;

  /* Figure out frozen core orbitals in each irrep and write them out*/
  nfzc = chkpt_rd_nfzc();
  nfzv = chkpt_rd_nfzv();
  frzcpi = compute_frzcpi(nfzc);
  frzvpi = compute_frzvpi(nfzv);
  chkpt_wt_frzcpi(frzcpi);
  chkpt_wt_frzvpi(frzvpi);

  /* Write eigenvectors and eigenvalues to new PSIF_CHKPT */
  scr_arr = init_array(nmo);

  if(uhf) { 
    for(m=0; m<2; m++) {
      for(k=0,i=0; k < num_ir; k++) {
        s=&scf_info[k];
	for(j=0; j < s->num_mo; j++,i++) 
	  scr_arr[i] = spin_info[m].scf_spin[k].fock_evals[j];
      }
      if(m==0) chkpt_wt_alpha_evals(scr_arr);
      else chkpt_wt_beta_evals(scr_arr);

    }
  }
  else {
    for(k=0,i=0; k < num_ir; k++) {
      s=&scf_info[k];
      for(j=0; j < s->num_mo; j++,i++) {
	scr_arr[i] = s->fock_evals[j];
      }
    }

    chkpt_wt_evals(scr_arr);

  }
  free(scr_arr);

  /* This will write the full SCF matrices (including zeroes) */
  scr1 = block_matrix(nbfso,nmo);
  if(uhf) {
    for(m=0; m < 2; m++) {
      zero_mat(scr1, nbfso, nmo);
      for(k=0,row=0,col=0; k < num_ir; k++) {
	s=&scf_info[k];
	for(i=0; i < s->num_so; i++) {
	  for(j=0; j < s->num_mo; j++) {
	    scr1[i+row][j+col] = spin_info[m].scf_spin[k].cmat[i][j];
	  }
	}
	row += s->num_so;
	col += s->num_mo;
      }
      if(m==0) chkpt_wt_alpha_scf(scr1);
      else chkpt_wt_beta_scf(scr1);
    }
  }
  else {
    for(k=0,row=0,col=0; k < num_ir; k++) {
      s=&scf_info[k];
      for(i=0; i < s->num_so; i++) {
	for(j=0; j < s->num_mo; j++) {
	  scr1[i+row][j+col] = s->cmat[i][j];
	}
      }
      row += s->num_so;
      col += s->num_mo;
    }
    chkpt_wt_scf(scr1);
  }
  free_block(scr1);
#else
  /* write new vector and eigenvalues to file30 and file49 */
  scr_arr = (double *) init_array(mxcoef);

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
    if (scf_info[i].num_so) {
      wwritw(itap30,(char *) scf_info[i].irrep_label,sizeof(char)*4,locvec,&locvec);
    }
   
  pointers[5] = locvec;
   
  wwritw(itap30,(char *) n_there,sizeof(int)*n_so_typs,locvec,&locvec);
  pointers[6] = locvec;
   
  wwritw(itap30,(char *) nc,sizeof(int)*n_so_typs,locvec,&locvec);
  pointers[7] = locvec;
   
  if(iopen || uhf) {
    wwritw(itap30,(char *) no,sizeof(int)*n_so_typs,locvec,&locvec);
    pointers[8] = locvec;
  }
#endif

  /* write open-shell coupling coefficients */
  if(iopen){
#if !USE_LIBCHKPT
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
#else
    ccvecs = block_matrix(2,ioff[n_open]);
    for(i=0; i < ioff[n_open]; i++) {
      ccvecs[0][i] = alpha[i];
      ccvecs[1][i] = beta[i];
    }
    chkpt_wt_ccvecs(ccvecs);
    free_block(ccvecs);
#endif
  }

  /* calculate mo lagrangian and write to file30 */
  /* also write mo fock matrices to file49       */
      
  lagrangian = (double *) init_array(nx);
  lagr = block_matrix(nmo,nmo);
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
	    for(j=0; j <= i ; j++) {
	      lagrangian[ioff[i+ijk]+j+ijk] 
		= spin_info[m].scf_spin[k].fock_pac[ioff[i]+j];
	      lagr[i][j] = lagr[j][i] = spin_info[m].scf_spin[k].fock_pac[ioff[i]+j];
	    }
	  ijk += num_mo;
	}
      }
#if !USE_LIBCHKPT
      wwritw(itap30,(char *) lagrangian,sizeof(double)*nx,locvec,&locvec); 
      pointers[m+11] = locvec; 
#else
      if(m==0) chkpt_wt_alpha_lagr(lagr);
      else chkpt_wt_beta_lagr(lagr);
#endif
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
	  for(j=0; j <= i ; j++) {
	    lagrangian[ioff[i+ijk]+j+ijk] = scf_info[k].fock_pac[ioff[i]+j];
	    lagr[i][j] = lagr[j][i] = scf_info[k].fock_pac[ioff[i]+j];
	  }
	ijk += num_mo;
      }
    }
#if !USE_LIBCHKPT
    wwritw(itap30,(char *) lagrangian,sizeof(double)*nx,locvec,&locvec);
    pointers[11]=0;
    pointers[12]=locvec;
#else
    chkpt_wt_lagr(lagr);
#endif
  }
  free(lagrangian);
  free_block(lagr);

#if !USE_LIBCHKPT
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
#endif
      
  if(ci_calc && iopen && irot) {
    fprintf(outfile,
	    "\n ci_typ is %s so mo vector will be rotated\n",ci_type);
    if (mo_print) print_mos();
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

  /* Print just the orbital eigenvalues --- CDS 4/01 */
  print_mo_eigvals();

  /*if(print & 1){
    print_mos_new();
    }*/
      
      
  fprintf(outfile,"\n%8cSCF total energy   = %20.12f\n",' ',etot);
      
  if(ksdft){
    fprintf(outfile,"%8ccoulomb energy     = %20.12f\n"
	    ,' ',coulomb_energy);  
    fprintf(outfile,"%8cexchange energy    = %20.12f\n"
	    ,' ' ,exch_energy);
    fprintf(outfile,"%8ccorrelation energy = %20.12f\n"
	    ,' ',corr_energy);
  }
      
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
      
#if USE_LIBCHKPT
  chkpt_close();
#else
  file30_close();
#endif

  if(!direct_scf){
    psio_close(Pmat.unit, 0);
    psio_close(PKmat.unit, 0);
  }
      
  if(!converged)
    fprintf(outfile,"\n%8cCalculation has not converged!\n",' ');
      
  tstop(outfile);
  psi_stop();
      
  if(!converged) exit(PSI_RETURN_FAILURE);
  exit(PSI_RETURN_SUCCESS);

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


/*
** print_mo_eigvals()
**
** Print out the MO eigenvalues, both RHF and UHF case
** C. David Sherrill
** April 2001
*/
void print_mo_eigvals(void)
{
  int i,irrep,done,lowest_irrep,printctr=0;
  int num_closed, num_open, num_mo, num_virt, intocc;
  int *counter;
  int *sorted_counter, **sorted_irreps, **sorted_index;
  double **sorted_evals;
  double tval,lowest,occup;
  double OCCTOL;

  OCCTOL = 1.0E-6;

  num_closed = num_open = num_mo = num_virt = 0;
  counter = init_int_array(num_ir);

  fprintf(outfile, "\nOrbital energies (a.u.):\n");

  /* TWOCON case */
  if (twocon) {

    /* Ok, just go through and pick out the lowest one each time */
    done = 0;
    printctr = 1;
    while (!done) {
      lowest = 1E9; lowest_irrep = 0;
      done = 1;
      for (irrep=0; irrep < num_ir; irrep++) {
        if (counter[irrep] == scf_info[irrep].num_mo) continue;
        done = 0;
        if ((tval = scf_info[irrep].fock_evals[counter[irrep]]) < lowest) {
          lowest = tval;
          lowest_irrep = irrep;
          occup = scf_info[lowest_irrep].occ_num[counter[lowest_irrep]];
        }
      }
      if (!done) {
        fprintf(outfile, " %3d%3s %9.4lf (%5.3lf) ",
                ++(counter[lowest_irrep]), scf_info[lowest_irrep].irrep_label,
                lowest, occup);
        if ((printctr % 3) == 0) fprintf(outfile, "\n");
        printctr++;
      } 
    }
    fprintf(outfile, "\n");
  }

  /* RHF/ROHF case */
  else if (!uhf) {

    /* Need to count up how many there are of each type */
    /* I'm not sure what the nopen/nhalf distinction is ... */
    /* It appears that nopen is used for both high-spin and twocon */
    for (irrep=0; irrep < num_ir; irrep++) {
      num_closed += scf_info[irrep].nclosed;    
      num_open += scf_info[irrep].nopen + scf_info[irrep].nhalf;
      num_mo += scf_info[irrep].num_mo;
    }

    num_virt = num_mo - num_closed - num_open;

    sorted_counter = init_int_array(3);
    sorted_index = init_int_matrix(3,num_mo);
    sorted_irreps = init_int_matrix(3,num_mo);
    sorted_evals = init_matrix(3,num_mo);

    /* Ok, just go through and pick out the lowest one each time */
    done = 0;
    while (!done) {
      lowest = 1E9; lowest_irrep = 0;
      done = 1;
      for (irrep=0; irrep < num_ir; irrep++) {
        if (counter[irrep] == scf_info[irrep].num_mo) continue;
        done = 0;
        if ((tval = scf_info[irrep].fock_evals[counter[irrep]]) < lowest) {
          lowest = tval;
          lowest_irrep = irrep;
          /* now figure out where to store it */
          occup = scf_info[lowest_irrep].occ_num[counter[lowest_irrep]];
        }
      }
      if (!done) {
        if (fabs(occup - 2.0) < OCCTOL) intocc = 2;
        else if (fabs(occup - 1.0) < OCCTOL) intocc = 1;
        else if (fabs(occup - 0.0) < OCCTOL) intocc = 0;
        else {
          fprintf(outfile, 
            "(print_mo_eigvals): I found an orbital with %f electrons...\n",
             occup);
          return;
        }
        sorted_irreps[intocc][sorted_counter[intocc]] = lowest_irrep;
        sorted_evals[intocc][sorted_counter[intocc]] = lowest;
        sorted_index[intocc][sorted_counter[intocc]] = counter[lowest_irrep]; 
        counter[lowest_irrep]++;
        sorted_counter[intocc]++;
      }
    }

    fprintf(outfile, "\n  Doubly occupied orbitals\n");
     for (i=0,printctr=1; i<num_closed; i++,printctr++) {
       fprintf(outfile, " %3d%3s %12.6lf  ", 
             sorted_index[2][i]+1, scf_info[sorted_irreps[2][i]].irrep_label,
             sorted_evals[2][i]);
       if ((printctr % 3) == 0) fprintf(outfile, "\n");
     } 
     if ((printctr-1) % 3) fprintf(outfile, "\n");

     if (num_open > 0) {
       fprintf(outfile, "\n  Singly occupied orbitals\n");
       for (i=0,printctr=1; i<num_open; i++,printctr++) {
         fprintf(outfile, " %3d%3s %12.6lf  ", 
               sorted_index[1][i]+1, scf_info[sorted_irreps[1][i]].irrep_label,
               sorted_evals[1][i]);
         if ((printctr % 3) == 0) fprintf(outfile, "\n");
       } 
     }
     if ((printctr-1) % 3) fprintf(outfile, "\n");
   
     fprintf(outfile, "\n  Unoccupied orbitals\n");
     for (i=0,printctr=1; i<num_virt; i++,printctr++) {
       fprintf(outfile, " %3d%3s %12.6lf  ", 
             sorted_index[0][i]+1, scf_info[sorted_irreps[0][i]].irrep_label,
             sorted_evals[0][i]);
       if ((printctr % 3) == 0) fprintf(outfile, "\n");
     } 
     if ((printctr-1) % 3) fprintf(outfile, "\n");

    free(sorted_counter);
    free_int_matrix(sorted_index,3);
    free_int_matrix(sorted_irreps,3);
    free_matrix(sorted_evals,3);
  }

  /* UHF case */
  else {

    /* I'm not sure what the nopen/nhalf distinction is ... */
    for (irrep=0; irrep < num_ir; irrep++) {
      num_closed += scf_info[irrep].nclosed;    
      num_open += scf_info[irrep].nopen + scf_info[irrep].nhalf;
      num_mo += scf_info[irrep].num_mo;
    }

    num_virt = num_mo - num_closed - num_open;

    sorted_counter = init_int_array(2);
    sorted_index = init_int_matrix(2,num_mo);
    sorted_irreps = init_int_matrix(2,num_mo);
    sorted_evals = init_matrix(2,num_mo);

    /* do alpha */

    /* Ok, just go through and pick out the lowest one each time */
    done = 0;
    while (!done) {
      lowest = 1E9; lowest_irrep = 0;
      done = 1;
      for (irrep=0; irrep < num_ir; irrep++) {
        if (counter[irrep] == scf_info[irrep].num_mo) continue;
        done = 0;
        if ((tval = spin_info[0].scf_spin[irrep].fock_evals[counter[irrep]]) 
            < lowest) {
          lowest = tval;
          lowest_irrep = irrep;
          /* now figure out where to store it */
    occup = spin_info[0].scf_spin[lowest_irrep].occ_num[counter[lowest_irrep]];
        }
      }
      if (!done) {
        if (fabs(occup - 1.0) < OCCTOL) intocc = 1;
        else if (fabs(occup - 0.0) < OCCTOL) intocc = 0;
        else {
          fprintf(outfile, 
            "(print_mo_eigvals): I found an orbital with %f electrons...\n",
             occup);
          return;
        }
        sorted_irreps[intocc][sorted_counter[intocc]] = lowest_irrep;
        sorted_evals[intocc][sorted_counter[intocc]] = lowest;
        sorted_index[intocc][sorted_counter[intocc]] = counter[lowest_irrep]; 
        counter[lowest_irrep]++;
        sorted_counter[intocc]++;
      }
    }

    fprintf(outfile, "\n  Alpha Occupied orbitals\n");
     for (i=0,printctr=1; i<sorted_counter[1]; i++,printctr++) {
       fprintf(outfile, " %3d%3s %12.6lf  ", 
             sorted_index[1][i]+1, scf_info[sorted_irreps[1][i]].irrep_label,
             sorted_evals[1][i]);
       if ((printctr % 3) == 0) fprintf(outfile, "\n");
     } 
     if ((printctr-1) % 3) fprintf(outfile, "\n");

     fprintf(outfile, "\n  Alpha Unoccupied orbitals\n");
     for (i=0,printctr=1; i<sorted_counter[0]; i++,printctr++) {
       fprintf(outfile, " %3d%3s %12.6lf  ", 
             sorted_index[0][i]+1, scf_info[sorted_irreps[0][i]].irrep_label,
             sorted_evals[0][i]);
       if ((printctr % 3) == 0) fprintf(outfile, "\n");
     } 
     if ((printctr-1) % 3) fprintf(outfile, "\n");



    /* do beta */
    zero_int_array(counter,num_ir);
    zero_int_array(sorted_counter,2);
    zero_int_matrix(sorted_index,2,num_mo);
    zero_int_matrix(sorted_irreps,2,num_mo);

    /* Ok, just go through and pick out the lowest one each time */
    done = 0;
    while (!done) {
      lowest = 1E9; lowest_irrep = 0;
      done = 1;
      for (irrep=0; irrep < num_ir; irrep++) {
        if (counter[irrep] == scf_info[irrep].num_mo) continue;
        done = 0;
        if ((tval = spin_info[1].scf_spin[irrep].fock_evals[counter[irrep]]) 
            < lowest) {
          lowest = tval;
          lowest_irrep = irrep;
          /* now figure out where to store it */
    occup = spin_info[1].scf_spin[lowest_irrep].occ_num[counter[lowest_irrep]];
        }
      }
      if (!done) {
        if (fabs(occup - 1.0) < OCCTOL) intocc = 1;
        else if (fabs(occup - 0.0) < OCCTOL) intocc = 0;
        else {
          fprintf(outfile, 
            "(print_mo_eigvals): I found an orbital with %f electrons...\n",
             occup);
          return;
        }
        sorted_irreps[intocc][sorted_counter[intocc]] = lowest_irrep;
        sorted_evals[intocc][sorted_counter[intocc]] = lowest;
        sorted_index[intocc][sorted_counter[intocc]] = counter[lowest_irrep]; 
        counter[lowest_irrep]++;
        sorted_counter[intocc]++;
      }
    }

    fprintf(outfile, "\n  Beta Occupied orbitals\n");
     for (i=0,printctr=1; i<sorted_counter[1]; i++,printctr++) {
       fprintf(outfile, " %3d%3s %12.6lf  ", 
             sorted_index[1][i]+1, scf_info[sorted_irreps[1][i]].irrep_label,
             sorted_evals[1][i]);
       if ((printctr % 3) == 0) fprintf(outfile, "\n");
     } 
     if ((printctr-1) % 3) fprintf(outfile, "\n");

     fprintf(outfile, "\n  Beta Unoccupied orbitals\n");
     for (i=0,printctr=1; i<sorted_counter[0]; i++,printctr++) {
       fprintf(outfile, " %3d%3s %12.6lf  ", 
             sorted_index[0][i]+1, scf_info[sorted_irreps[0][i]].irrep_label,
             sorted_evals[0][i]);
       if ((printctr % 3) == 0) fprintf(outfile, "\n");
     } 
     if ((printctr-1) % 3) fprintf(outfile, "\n");




    free(sorted_counter);
    free_int_matrix(sorted_index,2);
    free_int_matrix(sorted_irreps,2);
    free_matrix(sorted_evals,2);
  }

  fprintf(outfile, "\n");

  free(counter);

}



/*
 * write_scf_matrices():
 *
 * Write the full C matrix (including zeroes).
 *
 * libchkpt version only!  libfile30 is dead.
 *
 * This is taken out of the cleanup() code above.  I need to call just
 * this part when reorthogonalizing the C matrix.
 *
 * C. David Sherrill
 * October 2002
 */

void write_scf_matrices(void)
{
  struct symm *s;
  double **scr1;
  int row, col;
  int i, j, k, m;

  scr1 = block_matrix(nbfso,nmo);
  if(uhf) {
    for(m=0; m < 2; m++) {
      zero_mat(scr1, nbfso, nmo);
      for(k=0,row=0,col=0; k < num_ir; k++) {
	s=&scf_info[k];
	for(i=0; i < s->num_so; i++) {
	  for(j=0; j < s->num_mo; j++) {
	    scr1[i+row][j+col] = spin_info[m].scf_spin[k].cmat[i][j];
	  }
	}
	row += s->num_so;
	col += s->num_mo;
      }
      #if USE_LIBCHKPT
      if(m==0) chkpt_wt_alpha_scf(scr1);
      else chkpt_wt_beta_scf(scr1);
      #else
      if(m==0) file30_wt_alpha_scf(scr1);
      else file30_wt_beta_scf(scr1);
      #endif
    }
  }
  else {
    for(k=0,row=0,col=0; k < num_ir; k++) {
      s=&scf_info[k];
      for(i=0; i < s->num_so; i++) {
	for(j=0; j < s->num_mo; j++) {
	  scr1[i+row][j+col] = s->cmat[i][j];
	}
      }
      row += s->num_so;
      col += s->num_mo;
    }
    #if USE_LIBCHKPT
    chkpt_wt_scf(scr1);
    #else
    file30_wt_scf(scr1);
    #endif
  }
  free_block(scr1);
}


/*----------------------------------------------------------------------
  Figure out number of frozen DOCC's in each irrep from the eigenvalues
 ----------------------------------------------------------------------*/
static int* compute_frzcpi(int nfzc)
{
  int mo, nmo, docc, irrep;
  double last_lowest = -1.0E100;
  double lowest_eval;
  int lowest_eval_irrep;
  double *evals, eval;
  struct symm *s;
  int *frzcpi;

  frzcpi = init_int_array(num_ir);
  if (!uhf) {
    for(docc=0; docc<nfzc; docc++) {
      lowest_eval = 1.0E100;
      for(irrep=0; irrep<num_ir; irrep++) {
	s = &scf_info[irrep];
	evals = s->fock_evals;
	nmo = s->num_mo;
	for(mo=0; mo<nmo; mo++) {
	  eval = evals[mo];
	  if (eval < lowest_eval && eval > last_lowest) {
	    lowest_eval = eval;
	    lowest_eval_irrep = irrep;
	  }
	}
      }
      last_lowest = lowest_eval;
      frzcpi[lowest_eval_irrep]++;
    }
  }

  return frzcpi;
}

/*----------------------------------------------------------------------
  Figure out number of frozen UOCC's in each irrep from the eigenvalues
 ----------------------------------------------------------------------*/
static int* compute_frzvpi(int nfzv)
{
  int mo, nmo, uocc, irrep;
  double last_highest = 1.0E100;
  double highest_eval;
  int highest_eval_irrep;
  double *evals, eval;
  struct symm *s;
  int *frzvpi;

  frzvpi = init_int_array(num_ir);
  if (!uhf) {
    for(uocc=0; uocc<nfzv; uocc++) {
      highest_eval = -1.0E100;
      for(irrep=0; irrep<num_ir; irrep++) {
	s = &scf_info[irrep];
	evals = s->fock_evals;
	nmo = s->num_mo;
	for(mo=0; mo<nmo; mo++) {
	  eval = evals[mo];
	  if (eval > highest_eval && eval < last_highest) {
	    highest_eval = eval;
	    highest_eval_irrep = irrep;
	  }
	}
      }
      last_highest = highest_eval;
      frzvpi[highest_eval_irrep]++;
    }
  }

  return frzvpi;
}
