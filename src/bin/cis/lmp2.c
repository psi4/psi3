#include <stdio.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"


/* mp2(): Computes the MP2 energy and associated first-order wfn.
** 
** Given a set of non-canonical occupied orbitals and canonical
** virtual orbitals, this routine computes the local-MP2 energy as
** described originally by Saebo and Pulay, J. Chem. Phys. 86, 914
** (1987).
**
** TDC, June 2002
*/

void mp2(void)
{
  int i, j, k, a, b, ij, ab, iter, conv, row, col;
  int nocc, nvir, natom, **domain, weak;
  double energy, rms, weak_pair_energy;
  dpdbuf4 T2, newT2, D, D2;
  dpdfile2 fIJ, fij, fAB, fab;

  /* Clean out diagonal element of occ-occ Fock matrix */
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_copy(&fij, CC_OEI, "fIJ (non-diagonal)");
  dpd_file2_close(&fij);

  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fIJ (non-diagonal)");
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_rd(&fij);
  for(h=0; h < nirreps; h++)
    for(i=0; i < fij.params->rowtot[h]; i++)
      fij.matrix[0][i][i] = 0.0;

  dpd_file2_mat_wrt(&fij);
  dpd_file2_close(&fij);

  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij (non-diagonal)");
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_rd(&fij);
  for(h=0; h < nirreps; h++)
    for(i=0; i < fij.params->rowtot[h]; i++)
      fij.matrix[0][i][i] = 0.0;

  dpd_file2_mat_wrt(&fij);
  dpd_file2_close(&fij);

  /* Build initial LMP2 amplitudes */
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_copy(&D, CC_TAMPS, "LMP2 tIjAb");
  dpd_buf4_close(&D);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "LMP2 tIjAb");
  local_filter_T2(&T2);
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

    local_filter_T2(&newT2);

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
    exit(PSI_RETURN_FAILURE);
  }

  /* Turn off weak pairs again for the LCCSD */
  for(i=0,ij=0; i < nocc; i++)
    for(j=0; j < nocc; j++,ij++) {
      weak = 1;
      for(k=0; k < local.natom; k++)
        if(local.domain[i][k] && local.domain[j][k]) weak = 0;

      if(weak) local.weak_pairs[ij] = 1;
      else local.weak_pairs[ij] = 0; 
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

  fprintf(outfile, "\n\tLMP2 Weak Pair Energy   = %20.14f\n\n", weak_pair_energy);
  fprintf(outfile, "\n\tLMP2 Correlation Energy = %20.14f\n\n", energy);
  fprintf(outfile, "\n\tLMP2 Total Energy       = %20.14f\n\n", energy+moinfo.eref);

  local.weak_pair_energy = weak_pair_energy;

}

