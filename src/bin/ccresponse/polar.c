#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <physconst.h>
#define EXTERN
#include "globals.h"

void transmu(void);
void sortmu(void);
void mubar(void);

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

void polar(void)
{
  double **tensor;
  double polar, polar_LCX, polar_HXY, polar_LHX1Y1, polar_LHX2Y2, polar_LHX1Y2;
  char **cartcomp;
  int alpha, beta;
  double omega_nm, omega_ev, omega_cm;

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  tensor = block_matrix(3,3);

  transmu();
  sortmu();
  mubar();

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
	/*
	  fprintf(outfile, "polar_LCX    = %20.12f\n", polar_LCX);
	  fprintf(outfile, "polar_HXY    = %20.12f\n", polar_HXY);
	  fprintf(outfile, "polar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);
	  fprintf(outfile, "polar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);
	  fprintf(outfile, "polar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);
	*/

	tensor[alpha][beta] = -polar;
      }
    }
  }

  fprintf(outfile, "\n                 CCSD Dipole Polarizability [(e^2 a0^2)/E_h]:\n");
  fprintf(outfile, "  -------------------------------------------------------------------------\n");
  if(params.omega != 0.0) 
    omega_nm = (_c*_h*1e9)/(_hartree2J*params.omega);
  omega_ev = _hartree2ev*params.omega;
  omega_cm = _hartree2wavenumbers*params.omega;
  if(params.omega != 0.0)
    fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", 
	    params.omega, omega_nm, omega_ev, omega_cm);
  else
    fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (Inf nm, %5.3f eV, %8.2f cm-1)\n", 
	    params.omega, omega_ev, omega_cm);
  fprintf(outfile, "  -------------------------------------------------------------------------\n");
  mat_print(tensor, 3, 3, outfile);

  free_block(tensor);

  free(cartcomp[0]);
  free(cartcomp[1]);
  free(cartcomp[2]);
  free(cartcomp);
}
