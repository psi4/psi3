#include <stdio.h>
#include <string.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void init_X(char *cart, int irrep, double omega);
void sort_X(char *cart, int irrep, double omega);
void X1_build(char *cart, int irrep, double omega, int iter);
void X2_build(char *cart, int irrep, double omega);
double converged(char *cart, int irrep, double omega);
void save_X(char *cart, int irrep, double omega);
void print_X(char *cart, int irrep, double omega);
void update_X(char *cart, int irrep, double omega);
void diis(int iter, char *cart, int irrep, double omega);
double pseudopolar(char *cart, int irrep, double omega);
void cleanup(void);
void exit_io(void);

void compute_X(char *cart, int irrep, double omega)
{
  int iter=0, done=0;
  double rms, polar;
  dpdfile2 X1;
  dpdbuf4 X2;
  char lbl[32];

  fprintf(outfile, "\n\tComputing Perturbed %1s Wave Function at Frequency %5.3f E_h.\n", cart, omega);
  fprintf(outfile, "\tIter   Pseudopolarizability       RMS \n");
  fprintf(outfile, "\t----   --------------------   -----------\n");
  fflush(outfile);
  init_X(cart, irrep, omega);

  sort_X(cart, irrep, omega);
  polar = -2.0*pseudopolar(cart, irrep, omega);
  fprintf(outfile, "\t%4d   %20.12f\n", iter, polar);
  fflush(outfile);

  for(iter=1; iter <= params.maxiter; iter++) {

    sort_X(cart, irrep, omega);
    X1_build(cart, irrep, omega, iter);
    X2_build(cart, irrep, omega);
    update_X(cart, irrep, omega);
    rms = converged(cart, irrep, omega);
    if(rms <= params.convergence) {
      done = 1;
      save_X(cart, irrep, omega);
      sort_X(cart, irrep, omega);
      fprintf(outfile, "\t-----------------------------------------\n");
      fprintf(outfile, "\tConverged Perturbed %1s Wave Function X(%5.3f) to %4.3e\n", cart, omega, rms);
      fflush(outfile);
      break;
    }
    if(params.diis) diis(iter, cart, irrep, omega);
    save_X(cart, irrep, omega);
    sort_X(cart, irrep, omega);

    polar = -2.0*pseudopolar(cart, irrep, omega);
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

  /*  print_X(cart, irrep, omega); */
}
