#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#if !USE_LIBCHKPT
#include <file30_params.h>
#endif
#include <psifiles.h>
#include "input.h"
#include "global.h"
#include "defines.h"

#if USE_LIBCHKPT
void write_scf_calc();
#else
void write_scf_calc(PSI_FPTR *ptr, int *pointers);
#endif

/*-----------------------------------------------------------------------------------------------------------------

 -----------------------------------------------------------------------------------------------------------------*/

void write_scf_to_file30()
{

  PSI_FPTR ptr, junk;
  int *constants,*pointers,*calcs;
  double *arr_double;
  int *scf_pointers;
  
#if USE_LIBCHKPT
  psio_open(PSIF_CHKPT, PSIO_OPEN_OLD);
#else
  constants = init_int_array(MCONST);
  calcs = init_int_array(MPOINT);
  rfile(CHECKPOINTFILE);
#endif

  /*-----------------
    Update constants
   -----------------*/
#if !USE_LIBCHKPT
  wreadw(CHECKPOINTFILE,(char *) constants, MCONST*(sizeof(int)),100*sizeof(int),&junk);
  constants[40] = num_so_typs;
  constants[41] = mxcoef;
  constants[42] = iopen;
  constants[44] = 1;       /* number of calculations present in the file */
  constants[45] = num_mo;
  constants[46] = ref;
  constants[50] = 0;
  wwritw(CHECKPOINTFILE,(char *) constants, MCONST*sizeof(int),100*sizeof(int),&junk);
#else
  chkpt_wt_nsymhf(num_so_typs);
  chkpt_wt_iopen(iopen);
  chkpt_wt_nmo(num_mo);
  chkpt_wt_ref(ref);
#endif

#if !USE_LIBCHKPT
  /* That's where the end of the file is
     and where the calculation-specific data will go */
  ptr = (constants[0] - 1)*sizeof(int);
  scf_pointers = init_int_array(20);
  /* write the data out */
  write_scf_calc(&ptr,scf_pointers);
  /* That's where pointers to calculation-specific data are */
  wreadw(CHECKPOINTFILE,(char *) calcs, MCALCS*(sizeof(int)),500*sizeof(int),&ptr);
  ptr = (calcs[0] + 60 - 1)*sizeof(int);
  wwritw(CHECKPOINTFILE,(char *) scf_pointers, 20*sizeof(int),ptr,&ptr);
  free(scf_pointers);
#else
  /* write the data out */
  write_scf_calc(&ptr,scf_pointers);
#endif
  
  /* Update energies */
#if !USE_LIBCHKPT
  arr_double = init_array(5);
  ptr += 3*num_atoms*sizeof(double);
  wreadw(CHECKPOINTFILE,(char *) arr_double, 5*sizeof(double), ptr, &junk);
  arr_double[1] = escf;   /* SCF energy */
  arr_double[2] = escf;   /* Ref. energy (if different from SCF - put SCF anyway) */
  wwritw(CHECKPOINTFILE,(char *) arr_double, 5*sizeof(double),ptr,&ptr);
  free(arr_double);
#else
  psio_write_entry(PSIF_CHKPT, "::SCF energy", (char *) &escf, sizeof(double));
  psio_write_entry(PSIF_CHKPT, "::Reference energy", (char *) &escf, sizeof(double));
#endif

#if USE_LIBCHKPT
  psio_close(PSIF_CHKPT, 1);
#else
  rclose(CHECKPOINTFILE,3);
  free(constants);
  free(calcs);
#endif
  return;
}

#if USE_LIBCHKPT
void write_scf_calc()
{
  double *zero_array;
  zero_array = init_array(num_mo*(num_mo+1)/2);
      
  if (spinrestr_ref) {
    /* SCF eigenvector */
    chkpt_wt_scf(scf_evect_so);
    
    /* SCF eigenvalues */
    chkpt_wt_evals(zero_array);
  }
  else {
    /* SCF eigenvectors */
    chkpt_wt_alpha_scf(scf_evect_so_alpha);
    chkpt_wt_beta_scf(scf_evect_so_beta);
    
    /* SCF eigenvalues */
    chkpt_wt_alpha_evals(zero_array);
    chkpt_wt_beta_evals(zero_array);
  }
      
  /* irrep labels for non-empty blocks */
  chkpt_wt_irr_labs(irr_labels);

  /* MOs per block */
  chkpt_wt_orbspi(orbspi);

  /* doubly-occupied MOs per block */
  chkpt_wt_clsdpi(clsdpi);

  /* singly-occupied MOs per block */
  chkpt_wt_openpi(openpi);
  
  return;
}
#else
void write_scf_calc(PSI_FPTR *ptr, int *pointers)
{
  int irrep, count, irr_count;
  int so_offset, mo_offset, so, mo, nso, nmo;
  int *cpi, *opi, *tpi;
  char **labs;
  double *scfvec1, *scfvec2;
  double *zero_array;

  cpi = init_int_array(num_so_typs);
  opi = init_int_array(num_so_typs);
  tpi = init_int_array(num_so_typs);
  labs = (char **) malloc(num_so_typs * sizeof(char *));
  zero_array = init_array(num_mo*(num_mo+1)/2);
      
  /* Get rid of non-symmetric blocks of the eigenvector(s) */
  scfvec1 = init_array(mxcoef);
  if (!spinrestr_ref) scfvec2 = init_array(mxcoef);
  count = 0;
  irr_count = 0;
  so_offset = 0;
  mo_offset = 0;
  for(irrep=0;irrep<nirreps;irrep++)
      if (orbspi[irrep]) {
	  nso = num_so_per_irrep[irrep];
	  nmo = orbspi[irrep];
	  if (spinrestr_ref)
	      for(mo=0;mo<nmo;mo++)
		  for(so=0;so<nso;so++)
		      scfvec1[count++] = scf_evect_so[so+so_offset][mo+mo_offset];
	  else
	      for(mo=0;mo<nmo;mo++)
		  for(so=0;so<nso;so++) {
		      scfvec1[count] = scf_evect_so_alpha[so+so_offset][mo+mo_offset];
		      scfvec2[count++] = scf_evect_so_beta[so+so_offset][mo+mo_offset];
		  }
	  so_offset += nso; mo_offset += nmo;
	  
	  cpi[irr_count] = clsdpi[irrep];
	  opi[irr_count] = openpi[irrep];
	  tpi[irr_count] = orbspi[irrep];
	  labs[irr_count++] = strdup(irr_labels[irrep]);
      }
  if (count != mxcoef)
      punt("Number of written SCF coefficients is invalid");

  if (spinrestr_ref) {
      /* SCF eigenvector */
      pointers[0] = (*ptr)/sizeof(int) + 1;
      wwritw(CHECKPOINTFILE,(char *) scfvec1, mxcoef*sizeof(double),*ptr,ptr);
      free(scfvec1);
      
      /* SCF eigenvalues */
      pointers[2] = (*ptr)/sizeof(int) + 1;
      wwritw(CHECKPOINTFILE,(char *) zero_array, num_mo*sizeof(double),*ptr,ptr);
  }
  else {
      /* SCF eigenvectors */
      pointers[0] = (*ptr)/sizeof(int) + 1;
      wwritw(CHECKPOINTFILE,(char *) scfvec1, mxcoef*sizeof(double),*ptr,ptr);
      free(scfvec1);
      pointers[1] = (*ptr)/sizeof(int) + 1;
      wwritw(CHECKPOINTFILE,(char *) scfvec2, mxcoef*sizeof(double),*ptr,ptr);
      free(scfvec2);
      
      /* SCF eigenvalues */
      pointers[2] = (*ptr)/sizeof(int) + 1;
      wwritw(CHECKPOINTFILE,(char *) zero_array, num_mo*sizeof(double),*ptr,ptr);
      pointers[3] = (*ptr)/sizeof(int) + 1;
      wwritw(CHECKPOINTFILE,(char *) zero_array, num_mo*sizeof(double),*ptr,ptr);
  }
      
  /* irrep labels for non-empty blocks */
  pointers[4] = (*ptr)/sizeof(int) + 1;
  for(irrep=0;irrep<num_so_typs;irrep++)
      wwritw(CHECKPOINTFILE,(char *) labs[irrep], 4*sizeof(char),*ptr,ptr);

  /* MOs per block */
  pointers[5] = (*ptr)/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) tpi, num_so_typs*sizeof(int),*ptr,ptr);

  /* doubly-occupied MOs per block */
  pointers[6] = (*ptr)/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) cpi, num_so_typs*sizeof(int),*ptr,ptr);

  /* singly-occupied MOs per block */
  pointers[7] = (*ptr)/sizeof(int) + 1;
  wwritw(CHECKPOINTFILE,(char *) opi, num_so_typs*sizeof(int),*ptr,ptr);
  
  return;
}
#endif
