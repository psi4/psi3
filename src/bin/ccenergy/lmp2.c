#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void lmp2(void)
{
  int i, j, a, b, ij, ab, iter, conv, row, col, nocc, nvir;
  double energy, rms, weak_pair_energy;
  dpdbuf4 T2, newT2, D, D2;
  dpdfile2 fij, fab;

  nocc = moinfo.occpi[0];
  nvir = moinfo.virtpi[0];

  /* Clean out diagonal element of occ-occ Fock matrix */
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_copy(&fij, CC_OEI, "fIJ (non-diagonal)");
  dpd_file2_close(&fij);

  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fIJ (non-diagonal)");
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_rd(&fij);
  for(i=0; i < nocc; i++) fij.matrix[0][i][i] = 0.0;
  dpd_file2_mat_wrt(&fij);
  dpd_file2_close(&fij);

  /* Build initial LMP2 amplitudes */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_TAMPS, "LMP2 tIjAb");
  dpd_buf4_close(&D);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
  dpd_buf4_init(&D2, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
  dpd_buf4_dirprd(&D2, &T2);
  dpd_buf4_close(&D2);
  dpd_buf4_close(&T2);

  /* Compute the LMP2 energy */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
  energy = dpd_buf4_dot(&D, &T2);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&D);

  fprintf(outfile, "\n\tComputing LMP2 amplitudes:\n");
  fprintf(outfile,   "\t--------------------------\n");
  fprintf(outfile, "\titer = %d    LMP2 Energy = %20.14f\n", 0, energy);

  conv = 0;
  for(iter=1; iter < params.maxiter; iter++) {

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_copy(&D, CC_TAMPS, "New LMP2 tIjAb");
    dpd_buf4_close(&D);

    dpd_buf4_init(&newT2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
    dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fIJ (non-diagonal)");
    dpd_contract424(&T2, &fij, &newT2, 1, 0, 1, -1, 1);
    dpd_contract244(&fij, &T2, &newT2, 0, 0, 0, -1, 1);
    dpd_file2_close(&fij);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&D2, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    dpd_buf4_dirprd(&D2, &newT2);
    dpd_buf4_close(&D2);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    energy = dpd_buf4_dot(&D, &newT2);
    dpd_buf4_close(&D);

    dpd_buf4_close(&newT2);

    /* Check for convergence */
    dpd_buf4_init(&newT2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
    dpd_buf4_mat_irrep_init(&newT2, 0);
    dpd_buf4_mat_irrep_rd(&newT2, 0);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
    dpd_buf4_mat_irrep_init(&T2, 0);
    dpd_buf4_mat_irrep_rd(&T2, 0);

    rms = 0.0;
    for(row=0; row < T2.params->rowtot[0]; row++)
      for(col=0; col < T2.params->coltot[0]; col++)
	rms += (newT2.matrix[0][row][col] - T2.matrix[0][row][col]) *
	  (newT2.matrix[0][row][col] - T2.matrix[0][row][col]);

    dpd_buf4_mat_irrep_close(&T2, 0);
    dpd_buf4_mat_irrep_close(&newT2, 0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&newT2);

    rms = sqrt(rms);

    fprintf(outfile, "\titer = %d    LMP2 Energy = %20.14f   RMS = %4.3e\n", iter, energy, rms);

    if(rms < params.convergence) {
      conv = 1;
      fprintf(outfile, "\n\tLMP2 Iterations converged.\n");
      break;
    }
    else {
      dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
      dpd_buf4_copy(&T2, CC_TAMPS, "LMP2 tIjAb");
      dpd_buf4_close(&T2);
    }
  }

  if(!conv) {
    fprintf(outfile, "\n\tLMP2 Iterative procedure failed.\n");
    exit(2);
  }

  /* Compute the MP2 weak-pair energy */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  dpd_buf4_mat_irrep_init(&D, 0);
  dpd_buf4_mat_irrep_rd(&D, 0);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New LMP2 tIjAb");
  dpd_buf4_mat_irrep_init(&T2, 0);
  dpd_buf4_mat_irrep_rd(&T2, 0);

  weak_pair_energy = 0.0;
  for(ij=0; ij < nocc*nocc; ij++) {
    if(local.weak_pairs[ij]) {
      for(ab=0; ab < nvir*nvir; ab++) 
	weak_pair_energy += D.matrix[0][ij][ab] * T2.matrix[0][ij][ab];
    }
  }

  dpd_buf4_mat_irrep_close(&T2, 0);
  dpd_buf4_close(&T2);
  dpd_buf4_mat_irrep_close(&D, 0);
  dpd_buf4_close(&D);

  fprintf(outfile, "\n\tMP2 Weak Pair Energy = %20.14f\n\n", weak_pair_energy);

}
