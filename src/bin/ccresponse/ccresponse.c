/*
**  CCRESPONSE: Program to compute CC linear response properties.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <physconst.h>
#include <psifiles.h>
#include "globals.h"

/* Max length of ioff array */
#define IOFF_MAX 32641

/* Function prototypes */
void init_io(int argc, char *argv[]);
void init_ioff(void);
void title(void);
void get_moinfo(void);
void get_params(void);
void cleanup(void);
void exit_io(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);
void transmu(void);
void sortmu(void);
void hbar_extra(void);
void sort_lamps(void);
void compute_X(char *cart, int irrep, double omega);
double LCX(char *cart_c, int irrep_c, char *cart_x, int irrep_x, double omega);
double HXY(char *cart_x, int irrep_x, double omega_x, char *cart_y, int irrep_y, double omega_y);
double LHX1Y1(char *cart_x, int irrep_x, double omega_x, char *cart_y, int irrep_y, double omega_y);
double LHX2Y2(char *cart_x, int irrep_x, double omega_x, char *cart_y, int irrep_y, double omega_y);
double LHX1Y2(char *cart_x, int irrep_x, double omega_x, char *cart_y, int irrep_y, double omega_y);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;
  double polar, polar_LCX, polar_HXY, polar_LHX1Y1, polar_LHX2Y2, polar_LHX1Y2;
  char **cartcomp;
  int alpha, beta;
  double **tensor;

  init_io(argc, argv);
  init_ioff();
  title();
  get_moinfo();
  get_params();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /*** UHF references ***/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, 
	     NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi, moinfo.avir_sym,
	     moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  } 
  else { /*** RHF/ROHF references ***/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
	     2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);
  }

  transmu();
  sortmu();
  mubar();
  hbar_extra();
  sort_lamps();

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  tensor = block_matrix(3,3);

  /* Compute the dipole-perturbed CC wave functions */
  for(alpha=0; alpha < 3; alpha++) {
    compute_X(cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega);
    if(params.omega != 0.0) compute_X(cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega);
  }

  for(alpha=0; alpha < 3; alpha++) {
    for(beta=0; beta < 3; beta++) {
      if(!(moinfo.mu_irreps[alpha]^moinfo.mu_irreps[beta])) {

	if(params.omega != 0.0) {
	  polar_LCX = LCX(cartcomp[alpha], moinfo.mu_irreps[alpha], cartcomp[beta], 
			  moinfo.mu_irreps[beta], params.omega);
	  polar_LCX += LCX(cartcomp[beta], moinfo.mu_irreps[beta], cartcomp[alpha],
			   moinfo.mu_irreps[alpha], -params.omega);
	  polar_HXY = HXY(cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
			  cartcomp[beta], moinfo.mu_irreps[beta], params.omega);
	  polar_LHX1Y1 = LHX1Y1(cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				cartcomp[beta], moinfo.mu_irreps[beta], params.omega);
	  polar_LHX2Y2 = LHX2Y2(cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				cartcomp[beta], moinfo.mu_irreps[beta], params.omega);
	  polar_LHX1Y2 = LHX1Y2(cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				cartcomp[beta], moinfo.mu_irreps[beta], params.omega);
	  polar_LHX1Y2 += LHX1Y2(cartcomp[beta], moinfo.mu_irreps[beta], params.omega,
				 cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega);
	}
	else {
	  polar_LCX = LCX(cartcomp[alpha], moinfo.mu_irreps[alpha], cartcomp[beta],
			  moinfo.mu_irreps[beta], 0.0);
	  polar_LCX += LCX(cartcomp[beta], moinfo.mu_irreps[beta], cartcomp[alpha],
			  moinfo.mu_irreps[alpha], 0.0);
	  polar_HXY = HXY(cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
			  cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	  polar_LHX1Y1 = LHX1Y1(cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	  polar_LHX2Y2 = LHX2Y2(cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	  polar_LHX1Y2 = LHX1Y2(cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	  polar_LHX1Y2 += LHX1Y2(cartcomp[beta], moinfo.mu_irreps[beta], 0.0,
				 cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0);
	}

	polar = polar_LCX + polar_HXY + polar_LHX1Y1 + polar_LHX2Y2 + polar_LHX1Y2;
	/*
	fprintf(outfile, "polar_LCX    = %20.12f\n", polar_LCX);
	fprintf(outfile, "polar_HXY    = %20.12f\n", polar_HXY);
	fprintf(outfile, "polar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);
	fprintf(outfile, "polar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);
	fprintf(outfile, "polar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);
	*/

	tensor[alpha][beta] = -polar;
      }
      /*      fprintf(outfile, "%1s%1s polar = %20.12f\n", cartcomp[alpha], cartcomp[beta], polar); */
    }
  }

  fprintf(outfile, "\n                 CCSD Dipole Polarizability [(e^2 a0^2)/E_h]:\n");
  fprintf(outfile, "  -------------------------------------------------------------------------\n");
  fprintf(outfile,   "     Evaluated at omega = %5.3f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega,
	  (_c*_h*1e9)/(_hartree2J*params.omega), _hartree2ev*params.omega,
	  _hartree2wavenumbers*params.omega);
  fprintf(outfile, "  -------------------------------------------------------------------------\n");
  mat_print(tensor, 3, 3, outfile);

  for(alpha=0; alpha < 3; alpha++) free(cartcomp[alpha]);
  free(cartcomp);

  free_block(tensor);

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();

  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

void init_io(int argc, char *argv[])
{
  int i;
  extern char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(argc-1,argv+1,0); /* this assumes no cmdline args except filenames */
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);
  psio_init();

  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i, 1);
}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*       CCRESPONSE       *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;

  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i, 1);

  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
   char *prgid = "CCRESPONSE";

   return(prgid);
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}
