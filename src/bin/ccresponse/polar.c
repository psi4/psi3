#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <physconst.h>
#define EXTERN
#include "globals.h"

void transpert(char *pert);
void sort_pert(char *pert, double **pertX, double **pertY, double **pertZ,
	       int irrep_x, int irrep_y, int irrep_z);
void pertbar(char *pert, int irrep_x, int irrep_y, int irrep_z, int anti);
void compute_X(char *pert, char *cart, int irrep, double omega);
void linresp(double **tensor, double A, double B,
	     char *pert_x, int *x_irreps, double omega_x,
	     char *pert_y, int *y_irreps, double omega_y);

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

  for(i=0; i < params.nomega; i++) {

    transpert("Mu");
    sort_pert("Mu", moinfo.MUX, moinfo.MUY, moinfo.MUZ, moinfo.irrep_x,
	      moinfo.irrep_y, moinfo.irrep_z);
    pertbar("Mu", moinfo.irrep_x, moinfo.irrep_y, moinfo.irrep_z, 0);

    /* Compute the electric-dipole-perturbed CC wave functions */
    for(alpha=0; alpha < 3; alpha++) {
      compute_X("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega[i]);
      if(params.omega[i] != 0.0) compute_X("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i]);
    }

    linresp(tensor[i], -1.0, 0.0, "Mu", moinfo.mu_irreps, -params.omega[i],
	    "Mu", moinfo.mu_irreps, params.omega[i]);

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
