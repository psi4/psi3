#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <physconst.h>
#define EXTERN
#include "globals.h"

void status(char *, FILE *);
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
double cc2_LHX1Y1(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
		  char *pert_y, char *cart_y, int irrep_y, double omega_y);
double cc2_LHX1Y2(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
		  char *pert_y, char *cart_y, int irrep_y, double omega_y);

void polar(void)
{
  double ***tensor;
  double polar, polar_LCX, polar_HXY, polar_LHX1Y1, polar_LHX2Y2, polar_LHX1Y2;
  char **cartcomp;
  int alpha, beta, i;
  double omega_nm, omega_ev, omega_cm, *trace;

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  tensor = (double ***) malloc(params.nomega * sizeof(double **));
  for(i=0; i < params.nomega; i++) 
    tensor[i] = block_matrix(3,3);

  trace = init_array(params.nomega);

  transmu();
  sortmu();
  mubar();

  for(i=0; i < params.nomega; i++) {

    /* Compute the electric-dipole-perturbed CC wave functions */
    for(alpha=0; alpha < 3; alpha++) {
      compute_X("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega[i]);
      if(params.omega[i] != 0.0) compute_X("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i]);
    }

    for(alpha=0; alpha < 3; alpha++) {
      for(beta=0; beta < 3; beta++) {

	tensor[i][alpha][beta] = 0.0;
	polar_LCX = 0.0;
	polar_HXY = 0.0;
	polar_LHX1Y1 = 0.0;
	polar_LHX2Y2 = 0.0;
	polar_LHX1Y2 = 0.0;

	if(!(moinfo.mu_irreps[alpha]^moinfo.mu_irreps[beta])) {

	  if(params.omega[i] != 0.0) {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "Mu", cartcomp[beta], 
			    moinfo.mu_irreps[beta], params.omega[i]);
	    polar_LCX += LCX("Mu", cartcomp[beta], moinfo.mu_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], -params.omega[i]);
	    if (!strcmp(params.wfn,"CC2")) {
	      polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
			      "Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega[i]);
	      polar_LHX1Y1 = cc2_LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
					"Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega[i]);
	      polar_LHX1Y2 = cc2_LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
					"Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega[i]);
	      polar_LHX1Y2 += cc2_LHX1Y2("Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega[i],
					 "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i]);
	    }
	    else {
	      polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
				    "Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega[i]);
	      polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
				    "Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega[i]);
	      polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
				    "Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega[i]);
	      polar_LHX1Y2 += LHX1Y2("Mu", cartcomp[beta], moinfo.mu_irreps[beta], params.omega[i],
				     "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i]);
	    }
	  }
	  else {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "Mu", cartcomp[beta],
			    moinfo.mu_irreps[beta], 0.0);
	    polar_LCX += LCX("Mu", cartcomp[beta], moinfo.mu_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], 0.0);
	    if (!strcmp(params.wfn,"CC2")) {
	      polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
			      "Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	      polar_LHX1Y1 = cc2_LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
					"Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	      polar_LHX1Y2 = cc2_LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
					"Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	      polar_LHX1Y2 += cc2_LHX1Y2("Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0,
					 "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0);
	    }
	    else {
	      polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				    "Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	      polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				    "Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	      polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				    "Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0);
	      polar_LHX1Y2 += LHX1Y2("Mu", cartcomp[beta], moinfo.mu_irreps[beta], 0.0,
				     "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0);
	    }
	  }

	  polar = polar_LCX + polar_HXY + polar_LHX1Y1 + polar_LHX2Y2 + polar_LHX1Y2;

	  if(params.print & 2) {
	    fprintf(outfile, "\tTensor[%s][%s]\n", cartcomp[alpha], cartcomp[beta]);
	    fprintf(outfile, "\tpolar_LCX    = %20.15f\n", polar_LCX);
	    fprintf(outfile, "\tpolar_HXY    = %20.15f\n", polar_HXY);
	    fprintf(outfile, "\tpolar_LHX1Y1 = %20.15f\n", polar_LHX1Y1);
	    fprintf(outfile, "\tpolar_LHX1Y2 = %20.15f\n", polar_LHX1Y2);
	    fprintf(outfile, "\tpolar_LHX2Y2 = %20.15f\n\n", polar_LHX2Y2);
	    fflush(outfile);
	  }

	  tensor[i][alpha][beta] = -polar;
	}
      }
    }

    psio_close(CC_LR, 0);
    psio_open(CC_LR, 0);

    if (!strcmp(params.wfn,"CC2"))
      fprintf(outfile, "\n                 CC2 Dipole Polarizability [(e^2 a0^2)/E_h]:\n");
    else
      fprintf(outfile, "\n                 CCSD Dipole Polarizability [(e^2 a0^2)/E_h]:\n");
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    if(params.omega[i] != 0.0) 
      omega_nm = (_c*_h*1e9)/(_hartree2J*params.omega[i]);
    omega_ev = _hartree2ev*params.omega[i];
    omega_cm = _hartree2wavenumbers*params.omega[i];
    if(params.omega[i] != 0.0)
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", 
	      params.omega[i], omega_nm, omega_ev, omega_cm);
    else
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (Inf nm, %5.3f eV, %8.2f cm-1)\n", 
	      params.omega[i], omega_ev, omega_cm);
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    mat_print(tensor[i], 3, 3, outfile);
    trace[i] = tensor[i][0][0] + tensor[i][1][1] + tensor[i][2][2];
    fprintf(outfile, "\n\talpha_(%5.3f) = %20.12f a.u.\n", params.omega[i], trace[i]);
  }

  if(params.nomega > 1) {  /* print a summary table for multi-wavelength calcs */

    fprintf(outfile, "\n\t-------------------------------\n");
    if (!strcmp(params.wfn,"CC2"))
      fprintf(outfile,   "\t      CC2 Polarizability\n");
    else
      fprintf(outfile,   "\t      CCSD Polarizability\n");
    fprintf(outfile,   "\t-------------------------------\n");
    fprintf(outfile,   "\t    Omega          alpha\n");
    fprintf(outfile,   "\t E_h      nm        a.u.        \n");
    fprintf(outfile,   "\t-----   ------ ----------------\n");
    for(i=0; i < params.nomega; i++)
      fprintf(outfile, "\t%5.3f   %6.2f      %10.5f\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]), trace[i]);
  }

  for(i=0; i < params.nomega; i++)
    free_block(tensor[i]);
  free(tensor);

  free(trace);

  free(cartcomp[0]);
  free(cartcomp[1]);
  free(cartcomp[2]);
  free(cartcomp);
}
