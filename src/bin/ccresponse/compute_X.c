#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void init_X(char *pert, char *cart, int irrep, double omega);
void sort_X(char *pert, char *cart, int irrep, double omega);
void X1_build(char *pert, char *cart, int irrep, double omega);
void X2_build(char *pert, char *cart, int irrep, double omega);
void cc2_X1_build(char *pert, char *cart, int irrep, double omega);
void cc2_X2_build(char *pert, char *cart, int irrep, double omega);
double converged(char *pert, char *cart, int irrep, double omega);
void save_X(char *pert, char *cart, int irrep, double omega);
void print_X(char *pert, char *cart, int irrep, double omega);
void update_X(char *pert, char *cart, int irrep, double omega);
void diis(int iter, char *pert, char *cart, int irrep, double omega);
double pseudopolar(char *pert, char *cart, int irrep, double omega);
void cleanup(void);
void exit_io(void);

void compute_X(char *pert, char *cart, int irrep, double omega)
{
  int iter=0, done=0;
  double rms, polar, X2_norm;
  char lbl[32];
  dpdbuf4 X2;

  fprintf(outfile, "\n\tComputing %s-%1s-Perturbed Wave Function (%5.3f E_h).\n", pert, cart, omega);
  fprintf(outfile, "\tIter   Pseudopolarizability       RMS \n");
  fprintf(outfile, "\t----   --------------------   -----------\n");
  fflush(outfile);
  init_X(pert, cart, irrep, omega);

  sort_X(pert, cart, irrep, omega);
  polar = -2.0*pseudopolar(pert, cart, irrep, omega);
  fprintf(outfile, "\t%4d   %20.12f\n", iter, polar);
  fflush(outfile);

  for(iter=1; iter <= params.maxiter; iter++) {

    sort_X(pert, cart, irrep, omega);
    if (!strcmp(params.wfn,"CC2")) {
      cc2_X1_build(pert, cart, irrep, omega);
      cc2_X2_build(pert, cart, irrep, omega);
    }
    else {
      X1_build(pert, cart, irrep, omega);
      X2_build(pert, cart, irrep, omega);
    }
    update_X(pert, cart, irrep, omega);
    rms = converged(pert, cart, irrep, omega);
    if(rms <= params.convergence) {
      done = 1;
      save_X(pert, cart, irrep, omega);
      sort_X(pert, cart, irrep, omega);
      fprintf(outfile, "\t-----------------------------------------\n");
      fprintf(outfile, "\tConverged %s-%1s-Perturbed Wfn to %4.3e\n", pert, cart, rms);
      if(params.print == 2) {
	sprintf(lbl, "X_%s_%1s_IjAb (%5.3f)", pert, cart, omega);
	dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
	X2_norm = dpd_buf4_dot_self(&X2);
	dpd_buf4_close(&X2);
	X2_norm = sqrt(X2_norm);
	fprintf(outfile, "\tNorm of the converged X2 amplitudes %20.15f\n", X2_norm);
      }
      fflush(outfile);
      break;
    }
    if(params.diis) diis(iter, pert, cart, irrep, omega);
    save_X(pert, cart, irrep, omega);
    sort_X(pert, cart, irrep, omega);

    polar = -2.0*pseudopolar(pert, cart, irrep, omega);
    fprintf(outfile, "\t%4d   %20.12f    %4.3e\n", iter, polar, rms);
    fflush(outfile);

  }
  if(!done) {
    fprintf(outfile, "\t *** Failed to Converge Perturbed Wave Function to %2.1e.\n", params.convergence);
    fflush(outfile);
    dpd_close(0);
    cleanup();
    exit_io();
    exit(PSI_RETURN_FAILURE);
  }

  /*  print_X(pert, cart, irrep, omega); */
}
