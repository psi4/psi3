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
#include <masses.h>
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
void mubar(void);
void transL(void);
void sortL(void);
void Lbar(void);
void hbar_extra(void);
void sort_lamps(void);
void compute_X(char *pert, char *cart, int irrep, double omega);
double LCX(char *pert_c, char *cart_c, int irrep_c, 
	   char *pert_x, char *cart_x, int irrep_x, double omega);
double HXY(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
	   char *pert_y, char *cart_y, int irrep_y, double omega_y);
double LHX1Y1(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
	      char *pert_y, char *cart_y, int irrep_y, double omega_y);
double LHX2Y2(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
	      char *pert_y, char *cart_y, int irrep_y, double omega_y);
double LHX1Y2(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
	      char *pert_y, char *cart_y, int irrep_y, double omega_y);

void local_init(void);
void local_done(void);

int main(int argc, char *argv[])
{
  int **cachelist, *cachefiles;
  double polar, polar_LCX, polar_HXY, polar_LHX1Y1, polar_LHX2Y2, polar_LHX1Y2;
  char **cartcomp;
  int alpha, beta, i;
  double **tensor;
  double TrG, M, nu, rotation;

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

  if(params.local) local_init();

  hbar_extra();
  sort_lamps();

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  /* prepare electric dipole integrals */
  if(!strcmp(params.prop,"POLARIZABILITY") || !strcmp(params.prop,"ALL")) {
    transmu();
    sortmu();
    mubar();

    tensor = block_matrix(3,3);

    /* Compute the electric-dipole-perturbed CC wave functions */
    for(alpha=0; alpha < 3; alpha++) {
      compute_X("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega);
      if(params.omega != 0.0) compute_X("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega);
    }

    for(alpha=0; alpha < 3; alpha++) {
      for(beta=0; beta < 3; beta++) {

	tensor[alpha][beta] = 0.0;
	polar_LCX = 0.0;
	polar_HXY = 0.0;
	polar_LHX1Y1 = 0.0;
	polar_LHX2Y2 = 0.0;
	polar_LHX1Y2 = 0.0;

	if(!(moinfo.mu_irreps[alpha]^moinfo.mu_irreps[beta])) {

	  if(params.omega != 0.0) {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "Mu", cartcomp[beta], 
			    moinfo.mu_irreps[beta], params.omega);
	    polar_LCX += LCX("Mu", cartcomp[beta], moinfo.mu_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], -params.omega);
	    polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
			    "Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega);
	    polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				  "Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega);
	    polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				  "Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega);
	    polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				  "Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega);
	    polar_LHX1Y2 += LHX1Y2("Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega,
				   "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega);
	  }
	  else {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "Mu", cartcomp[beta],
			    moinfo.mu_irreps[beta], 0.0);
	    polar_LCX += LCX("Mu", cartcomp[beta], moinfo.mu_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], 0.0);
	    polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
			    "Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	    polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	    polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	    polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	    polar_LHX1Y2 += LHX1Y2("Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0,
				   "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0);
	  }

	  polar = polar_LCX + polar_HXY + polar_LHX1Y1 + polar_LHX2Y2 + polar_LHX1Y2;
	  fprintf(outfile, "polar_LCX    = %20.12f\n", polar_LCX);
	  fprintf(outfile, "polar_HXY    = %20.12f\n", polar_HXY);
	  fprintf(outfile, "polar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);
	  fprintf(outfile, "polar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);
	  fprintf(outfile, "polar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);

	  tensor[alpha][beta] = -polar;
	}
      }
    }

    fprintf(outfile, "\n                 CCSD Dipole Polarizability [(e^2 a0^2)/E_h]:\n");
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    fprintf(outfile,   "     Evaluated at omega = %5.3f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega,
	    (_c*_h*1e9)/(_hartree2J*params.omega), _hartree2ev*params.omega,
	    _hartree2wavenumbers*params.omega);
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    mat_print(tensor, 3, 3, outfile);
  }

  /*** Optical rotation -- I STILL NEED TO CODE THE COMPLEX CONJUGATE OF <<mu;m>>!!! ***/

  /* prepare magnetic dipole integrals */
  if(!strcmp(params.prop,"ROTATION") || !strcmp(params.prop,"ALL")) {
    transL();
    sortL();
    Lbar();

    /* Compute the electric-dipole-perturbed CC wave functions, if necessary */
    if(!strcmp(params.prop,"ROTATION")) {
      for(alpha=0; alpha < 3; alpha++) {
	compute_X("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega);
	if(params.omega != 0.0) compute_X("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega);
      }
    }

    /* Compute the magnetic-dipole-perturbed CC wave functions */
    /* NB: For mixed response functions, we need only the +omega perturbed wfns */
    for(alpha=0; alpha < 3; alpha++)
      compute_X("L", cartcomp[alpha], moinfo.l_irreps[alpha], params.omega);

    for(alpha=0; alpha < 3; alpha++) {
      for(beta=0; beta < 3; beta++) {

	tensor[alpha][beta] = 0.0;
	polar_LCX = 0.0;
	polar_HXY = 0.0;
	polar_LHX1Y1 = 0.0;
	polar_LHX2Y2 = 0.0;
	polar_LHX1Y2 = 0.0;

	if(!(moinfo.mu_irreps[alpha]^moinfo.l_irreps[beta])) {

	  if(params.omega != 0.0) {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "L", cartcomp[beta], 
			    moinfo.l_irreps[beta], params.omega);
	    polar_LCX += LCX("L", cartcomp[beta], moinfo.l_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], -params.omega);
	    polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
			    "L", cartcomp[beta], moinfo.l_irreps[beta], params.omega);
	    polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], params.omega);
	    polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], params.omega);
	    polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], params.omega);
	    polar_LHX1Y2 += LHX1Y2("L", cartcomp[beta], moinfo.l_irreps[beta], params.omega,
				   "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega);
	  }
	  else {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "L", cartcomp[beta],
			    moinfo.l_irreps[beta], 0.0);
	    polar_LCX += LCX("L", cartcomp[beta], moinfo.l_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], 0.0);
	    polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
			    "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX1Y2 += LHX1Y2("L", cartcomp[beta], moinfo.l_irreps[beta], 0.0,
				   "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0);
	  }

	  polar = polar_LCX + polar_HXY + polar_LHX1Y1 + polar_LHX2Y2 + polar_LHX1Y2;
	  fprintf(outfile, "polar_LCX    = %20.12f\n", polar_LCX);
	  fprintf(outfile, "polar_HXY    = %20.12f\n", polar_HXY);
	  fprintf(outfile, "polar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);
	  fprintf(outfile, "polar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);
	  fprintf(outfile, "polar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);

	  tensor[alpha][beta] = -polar;
	}
	/*      fprintf(outfile, "%1s%1s polar = %20.12f\n", cartcomp[alpha], cartcomp[beta], polar); */
      }
    }

    fprintf(outfile, "\n                 CCSD Optical Rotation Tensor [(e^2 a0^2)/E_h]:\n");
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    fprintf(outfile,   "     Evaluated at omega = %5.3f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega,
	    (_c*_h*1e9)/(_hartree2J*params.omega), _hartree2ev*params.omega,
	    _hartree2wavenumbers*params.omega);
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    mat_print(tensor, 3, 3, outfile);

    /* compute the specific rotation */
    for(i=0,M=0.0; i < moinfo.natom ;i++) M += an2masses[(int) moinfo.zvals[i]];
    M *= _amu2g*_na; /* g/mol */
    TrG = -(2.0/3.0) * (tensor[0][0] + tensor[1][1] + tensor[2][2]);
    nu = params.omega * _hartree2wavenumbers;
    rotation = 1.343e-4 * TrG * nu * nu / M;
    fprintf(outfile, "[alpha]_(%5.3f) = %20.12f\n", params.omega, rotation);
  }

  for(alpha=0; alpha < 3; alpha++) free(cartcomp[alpha]);
  free(cartcomp);

  free_block(tensor);

  if(params.local) local_done();

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
  ip_cwk_add(":INPUT");
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
