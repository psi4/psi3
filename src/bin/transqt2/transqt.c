/*! \file 
    \ingroup (TRANSQT2)
    \brief Enter brief description of file here 
*/
/*
** TRANSQT:
**
** A program to transform one- and two-electron integrals from the
** symmetry-orbital basis to the molecular-orbital basis.
**
** This code replaces the original transqt code developed initially
** in 1995 by TDC, CDS, and JTF.  This version is designed to take
** advantage of libdpd's ability to handle easily four-index
** quantities, including symmetry.  This version requires
** significantly less disk space (ca. 1/2) than the original code,
** and is often much faster because of its reduced I/O requirements.
**
** This version of the code can do RHF, ROHF, and UHF transformations
** that are compatible with all the coupled cluster codes, including
** frozen orbitals.  
**
** Remaining tasks to achieve full replacement of transqt v.1:
**   (1) Add reordering arrays needed for DETCI.
**   (2) Add partial transforms for MP2 and MP2-R12.
**   (3) Replace the backtransformation.  (I want to do this with
**       symmetry, though, so there's no hurry here.)
**
** TDC, 7/06
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "globals.h"

/* Function prototypes */
void init_io(int argc, char *argv[]);
void init_ioff(void);
void title(void);
void get_params(void);
void get_moinfo(void);
void cleanup(void);
void exit_io(void);

int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_rhf(int **);

int file_build_presort(dpdfile4 *, int, double, long int, int, 
		       int, double *, double *, double *, double *, int);
void transtwo_rhf(void);
void transone(int m, int n, double *input, double *output, double **C, int nc, int *order, int *ioff);
void semicanonical_fock(void);

main(int argc, char *argv[])
{
  int nso, nmo, ntri_so, ntri_mo, nirreps;
  int **cachelist, *cachefiles;
  dpdfile4 I;
  int h, pq, p, q, i;
  double *H, *D, *F, *oei;
  double *H_a, *H_b, *D_a, *D_b, *F_a, *F_b;
  double ***C, ***C_a, ***C_b;
  int stat;
  int *offset;
  double efzc;

  init_io(argc,argv);
  init_ioff();
  title();
  get_params();
  get_moinfo();
  if(params.semicanonical) semicanonical_fock();

  timer_init();

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  ntri_so = nso*(nso+1)/2;
  ntri_mo = nmo*(nmo+1)/2;
  nirreps = moinfo.nirreps;

  cachefiles = init_int_array(PSIO_MAXUNIT);
  cachelist = cacheprep_rhf(params.cachelev, cachefiles); /* really just a placeholder */

  dpd_init(0, nirreps, params.memory, 0, cachefiles, cachelist,
	   NULL, 2, moinfo.sopi, moinfo.sosym, moinfo.mopi, moinfo.mosym);

  /*** Starting one-electron transforms and presort ***/

  /* For the one-electron integral transform, the full MO list is best */
  if(params.ref == 0 || params.ref == 1) {
    C = (double ***) malloc(1 * sizeof(double **));
    chkpt_init(PSIO_OPEN_OLD);
    C[0] = chkpt_rd_scf();
    chkpt_close();
  }
  else {
    C_a = (double ***) malloc(1 * sizeof(double **));
    C_b = (double ***) malloc(1 * sizeof(double **));
    chkpt_init(PSIO_OPEN_OLD);
    C_a[0] = chkpt_rd_alpha_scf();
    C_b[0] = chkpt_rd_beta_scf();
    chkpt_close();
  }

  /* build the frozen-core density (RHF) */
  offset = init_int_array(nirreps);
  for(h=1; h < nirreps; h++)
    offset[h] = offset[h-1] + moinfo.sopi[h-1];

  if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
    D = init_array(ntri_so);
    for(h=0; h < nirreps; h++)
      for(p=offset[h]; p < offset[h]+moinfo.sopi[h]; p++)
	for(q=offset[h]; q <=p; q++) {
	  pq = INDEX(p,q);
	  for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]; i++)
	    D[pq] += C[0][p][i] * C[0][q][i];
	}
    if(params.print_lvl > 2) {
      fprintf(outfile, "\n\tFrozen-core density (SO):\n");
      print_array(D, nso, outfile);
    }
  }
  else { /* UHF */
    D_a = init_array(ntri_so);
    D_b = init_array(ntri_so);
    for(h=0; h < nirreps; h++)
      for(p=offset[h]; p < offset[h]+moinfo.sopi[h]; p++)
	for(q=offset[h]; q <=p; q++) {
	  pq = INDEX(p,q);
	  for(i=offset[h]; i < offset[h]+moinfo.frdocc[h]; i++) {
	    D_a[pq] += C_a[0][p][i] * C_a[0][q][i];
	    D_b[pq] += C_b[0][p][i] * C_b[0][q][i];
	  }
	}
    if(params.print_lvl > 2) {
      fprintf(outfile, "\n\tAlpha Frozen-core density (SO):\n");
      print_array(D_a, nso, outfile);
      fprintf(outfile, "\n\tBeta Frozen-core density (SO):\n");
      print_array(D_b, nso, outfile);
    }
  }

  free(offset);

  /* pre-sort the SO-basis two-electron integrals and generate the fzc operator(s) */
  if(params.ref == 0 || params.ref == 1)
    F = init_array(ntri_so);
  else {
    F_a = init_array(ntri_so);
    F_b = init_array(ntri_so);
  }

  timer_on("presort");
  if(params.print_lvl) {
    fprintf(outfile, "\n\tPresorting SO-basis two-electron integrals.\n");
    fflush(outfile);
  }
  psio_open(PSIF_SO_PRESORT, 0);
  dpd_file4_init(&I, PSIF_SO_PRESORT, 0, 3, 3, "SO Ints (pq,rs)");
  if(params.ref == 0 || params.ref == 1) 
    file_build_presort(&I, PSIF_SO_TEI, params.tolerance, params.memory, 
	!params.delete_tei, moinfo.nfzc, D, NULL, F, NULL, params.ref);
  else 
    file_build_presort(&I, PSIF_SO_TEI, params.tolerance, params.memory, 
	!params.delete_tei, moinfo.nfzc, D_a, D_b, F_a, F_b, params.ref);
  dpd_file4_close(&I);
  psio_close(PSIF_SO_PRESORT, 1);
  timer_off("presort");

  /* read the bare one-electron integrals */
  oei = init_array(ntri_so);
  H = init_array(ntri_so);
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_T, H, ntri_so, 0, 0, outfile);
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_V, oei, ntri_so, 0, 0, outfile);
  for(pq=0; pq < ntri_so; pq++)
    H[pq] += oei[pq];

  /* add the remaining one-electron terms to the fzc operator(s) */
  if(params.ref == 0 || params.ref == 1) {
    for(pq=0; pq < ntri_so; pq++)
      F[pq] += H[pq];
  }
  else {
    for(pq=0; pq < ntri_so; pq++) {
      F_a[pq] += H[pq];
      F_b[pq] += H[pq];
    }
  }

  /* compute the frozen-core energy and write it to the chkpt file*/
  efzc = 0.0;
  if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
    for(p=0; p < nso; p++) {
      pq = INDEX(p,p);
      efzc += D[pq] * (H[pq] + F[pq]);
      for(q=0; q < p; q++) {
	pq = INDEX(p,q);
	efzc += 2.0 * D[pq] * (H[pq] + F[pq]);
      }
    }
  }
  else { /* UHF */
    for(p=0; p < nso; p++) {
      pq = INDEX(p,p);
      efzc += 0.5 * D_a[pq] * (H[pq] + F_a[pq]);
      efzc += 0.5 * D_b[pq] * (H[pq] + F_b[pq]);
      for(q=0; q < p; q++) {
	pq = INDEX(p,q);
	efzc += D_a[pq] * (H[pq] + F_a[pq]);
	efzc += D_b[pq] * (H[pq] + F_b[pq]);
      }
    }
  }
  if(params.print_lvl) {
    fprintf(outfile, "\tFrozen-core energy = %20.15f\n", efzc);
    fflush(outfile);
  }
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_efzc(efzc);
  chkpt_close();

  /* transform the bare one-electron integrals */
  if(params.ref == 0 || params.ref == 1) {
    transone(nso, nmo, H, oei, C[0], nmo, moinfo.pitzer2qt, ioff);
    if(params.print_lvl > 2) {
      fprintf(outfile, "\n\tOne-electron integrals (MO basis):\n");
      print_array(oei, nmo, outfile);
    }
    iwl_wrtone(PSIF_OEI, PSIF_MO_OEI, ntri_mo, oei);
  }
  else { /* UHF */
    /* alpha */
    transone(nso, nmo, H, oei, C_a[0], nmo, moinfo.pitzer2qt_A, ioff);
    if(params.print_lvl > 2) {
      fprintf(outfile, "\n\tAlpha one-electron integrals (MO basis):\n");
      print_array(oei, nmo, outfile);
    }
    iwl_wrtone(PSIF_OEI, PSIF_MO_A_OEI, ntri_mo, oei);

    /* beta */
    transone(nso, nmo, H, oei, C_b[0], nmo, moinfo.pitzer2qt_B, ioff);
    if(params.print_lvl > 2) {
      fprintf(outfile, "\n\tBeta one-electron integrals (MO basis):\n");
      print_array(oei, nmo, outfile);
    }
    iwl_wrtone(PSIF_OEI, PSIF_MO_B_OEI, ntri_mo, oei);
  }

  /* transform the frozen-core operator */
  if(params.ref == 0 || params.ref == 1) { /* RHF/ROHF */
    transone(nso, nmo, F, oei, C[0], nmo, moinfo.pitzer2qt, ioff);
    if(params.print_lvl > 2) {
      fprintf(outfile, "\n\tFrozen-core operator (MO basis):\n");
      print_array(oei, nmo, outfile);
    }
    iwl_wrtone(PSIF_OEI, PSIF_MO_FZC, ntri_mo, oei);
  }
  else { /* UHF */

    /* alpha */
    transone(nso, nmo, F_a, oei, C_a[0], nmo, moinfo.pitzer2qt_A, ioff);
    if(params.print_lvl > 2) {
      fprintf(outfile, "\n\tAlpha frozen-core operator (MO basis):\n");
      print_array(oei, nmo, outfile);
    }
    iwl_wrtone(PSIF_OEI, PSIF_MO_A_FZC, ntri_mo, oei);

    /* beta */
    transone(nso, nmo, F_b, oei, C_b[0], nmo, moinfo.pitzer2qt_B, ioff);
    if(params.print_lvl > 2) {
      fprintf(outfile, "\n\tBeta frozen-core operator (MO basis):\n");
      print_array(oei, nmo, outfile);
    }
    iwl_wrtone(PSIF_OEI, PSIF_MO_B_FZC, ntri_mo, oei);
  }

  free(oei);
  free(H);
  if(params.ref == 0 || params.ref == 1) {
    free(F);
    free(D);
    free_block(C[0]);
    free(C);
  }
  else {
    free(F_a);
    free(F_b);
    free(D_a);
    free(D_b);
    free_block(C_a[0]);
    free(C_a);
    free_block(C_b[0]);
    free(C_b);
  }

  /*** One-electron transforms complete ***/

  /*** Starting two-electron transforms ***/

  if(params.ref == 0 || params.ref == 1) transtwo_rhf();
  else transtwo_uhf();

  /*** Two-electron transforms complete ***/

  dpd_close(0);

  cachedone_rhf(cachelist);
  free(cachefiles);

  timer_done();

  cleanup();
  free(ioff);
  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

char *gprgid(void) { char *prgid = "TRANSQT"; return (prgid); }

void init_io(int argc, char *argv[])
{
  int i;
  extern char *gprgid(void);
  char *progid;
  int num_extra_args = 0;
  char **extra_args;
  extra_args = (char **) malloc(argc*sizeof(char *));

  params.print_lvl = 1;
  for (i=1; i<argc; i++) {
    if (!strcmp(argv[i], "--quiet"))
      params.print_lvl = 0;
    else
      extra_args[num_extra_args++] = argv[i];
  }

  psi_start(num_extra_args, extra_args, 0);

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s", gprgid());
  ip_cwk_add(progid);
  free(progid);

  if(params.print_lvl) tstart(outfile);
  psio_init();

  psio_open(CC_INFO, PSIO_OPEN_NEW);
}

void title(void)
{
  if(params.print_lvl) {
    fprintf(outfile, "\n");
    fprintf(outfile,"\t**************************************************\n");
    fprintf(outfile,"\t* TRANSQT:  Program to transform integrals from  *\n");
    fprintf(outfile,"\t*           the SO basis to the MO basis.        *\n");
    fprintf(outfile,"\t*                                                *\n");
    fprintf(outfile,"\t*            Daniel, David, & Justin             *\n");
    fprintf(outfile,"\t**************************************************\n");
    fprintf(outfile, "\n");
  }
}

void exit_io(void)
{
  psio_close(CC_INFO,1);
  psio_done();
  if(params.print_lvl) tstop(outfile);
  psi_stop();
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  for(i=1; i < IOFF_MAX; i++)
    ioff[i] = ioff[i-1] + i;
}
