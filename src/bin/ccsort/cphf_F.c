#include <stdio.h>
#include <stdlib.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void cphf_F(void)
{
  int irrep, row, a, i, asym, isym, num_ai, info, *ipiv;
  double *vector, *vector2, polar;
  dpdbuf4 A;
  dpdfile2 mu;

  psio_open(PSIF_MO_HESS, 1);
  dpd_buf4_init(&A, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");

  irrep = moinfo.irrep_x;

  /* sort Mu elements into a single vector for lineq solver */
  dpd_file2_init(&mu, CC_OEI, irrep, 1, 0, "Mu_X_AI");
  dpd_file2_mat_init(&mu);
  dpd_file2_mat_rd(&mu);
  num_ai = A.params->rowtot[irrep];
  vector = init_array(num_ai);
  for(row=0; row < num_ai; row++) {
    a = A.params->roworb[irrep][row][0];
    i = A.params->roworb[irrep][row][1];
    asym = A.params->psym[a];
    isym = A.params->qsym[i];
    vector[row] = mu.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]];
  }
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  /* grab current irrep of MO Hessian */
  dpd_buf4_mat_irrep_init(&A, irrep);
  dpd_buf4_mat_irrep_rd(&A, irrep);


  /* solve CPHF equations */
  ipiv = init_int_array(num_ai);
  info = C_DGESV(num_ai, 1, &(A.matrix[irrep][0][0]), num_ai, ipiv, vector, num_ai);
  if(info) {
    fprintf(outfile, "CCSORT: cphf_F: Error in C_DGESV.  Info = %d.  Exiting.\n", info);
    exit(PSI_RETURN_FAILURE);
  }
  free(ipiv);

  dpd_buf4_mat_irrep_close(&A, irrep);

  /* sort CPHF solution to DPD format */
  dpd_file2_init(&mu, CC_OEI, irrep, 1, 0, "CPHF Uf_X_AI");
  dpd_file2_mat_init(&mu);
  for(row=0; row < num_ai; row++) {
    a = A.params->roworb[irrep][row][0];
    i = A.params->roworb[irrep][row][1];
    asym = A.params->psym[a];
    isym = A.params->qsym[i];
    mu.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]] = vector[row];
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_close(&mu);

  irrep = moinfo.irrep_y;

  /* sort Mu elements into a single vector for lineq solver */
  dpd_file2_init(&mu, CC_OEI, irrep, 1, 0, "Mu_Y_AI");
  dpd_file2_mat_init(&mu);
  dpd_file2_mat_rd(&mu);
  num_ai = A.params->rowtot[irrep];
  vector = init_array(num_ai);
  for(row=0; row < num_ai; row++) {
    a = A.params->roworb[irrep][row][0];
    i = A.params->roworb[irrep][row][1];
    asym = A.params->psym[a];
    isym = A.params->qsym[i];
    vector[row] = mu.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]];
  }
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  /* grab current irrep of MO Hessian */
  dpd_buf4_mat_irrep_init(&A, irrep);
  dpd_buf4_mat_irrep_rd(&A, irrep);


  /* solve CPHF equations */
  ipiv = init_int_array(num_ai);
  info = C_DGESV(num_ai, 1, &(A.matrix[irrep][0][0]), num_ai, ipiv, vector, num_ai);
  if(info) {
    fprintf(outfile, "CCSORT: cphf_F: Error in C_DGESV.  Info = %d.  Exiting.\n", info);
    exit(PSI_RETURN_FAILURE);
  }
  free(ipiv);

  dpd_buf4_mat_irrep_close(&A, irrep);

  /* sort CPHF solution to DPD format */
  dpd_file2_init(&mu, CC_OEI, irrep, 1, 0, "CPHF Uf_Y_AI");
  dpd_file2_mat_init(&mu);
  for(row=0; row < num_ai; row++) {
    a = A.params->roworb[irrep][row][0];
    i = A.params->roworb[irrep][row][1];
    asym = A.params->psym[a];
    isym = A.params->qsym[i];
    mu.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]] = vector[row];
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_close(&mu);

  irrep = moinfo.irrep_z;

  /* sort Mu elements into a single vector for lineq solver */
  dpd_file2_init(&mu, CC_OEI, irrep, 1, 0, "Mu_Z_AI");
  dpd_file2_mat_init(&mu);
  dpd_file2_mat_rd(&mu);
  num_ai = A.params->rowtot[irrep];
  vector = init_array(num_ai);
  for(row=0; row < num_ai; row++) {
    a = A.params->roworb[irrep][row][0];
    i = A.params->roworb[irrep][row][1];
    asym = A.params->psym[a];
    isym = A.params->qsym[i];
    vector[row] = mu.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]];
  }
  dpd_file2_mat_close(&mu);
  dpd_file2_close(&mu);

  /*
  vector2 = init_array(num_ai);
  for(row=0; row < num_ai; row++) vector2[row] = vector[row];
  */

  /* grab current irrep of MO Hessian */
  dpd_buf4_mat_irrep_init(&A, irrep);
  dpd_buf4_mat_irrep_rd(&A, irrep);

  /* solve CPHF equations */
  ipiv = init_int_array(num_ai);
  info = C_DGESV(num_ai, 1, &(A.matrix[irrep][0][0]), num_ai, ipiv, vector, num_ai);
  if(info) {
    fprintf(outfile, "CCSORT: cphf_F: Error in C_DGESV.  Info = %d.  Exiting.\n", info);
    exit(PSI_RETURN_FAILURE);
  }
  free(ipiv);

  dpd_buf4_mat_irrep_close(&A, irrep);

  /*
  polar = 0.0;
  for(row=0; row < num_ai; row++) polar += vector2[row] * vector[row];
  polar *= 4.0;
  fprintf(outfile, "\talpha_zz = %20.12f\n", polar);
  free(vector2);
  */

  /* sort CPHF solution to DPD format */
  dpd_file2_init(&mu, CC_OEI, irrep, 1, 0, "CPHF Uf_Z_AI");
  dpd_file2_mat_init(&mu);
  for(row=0; row < num_ai; row++) {
    a = A.params->roworb[irrep][row][0];
    i = A.params->roworb[irrep][row][1];
    asym = A.params->psym[a];
    isym = A.params->qsym[i];
    mu.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]] = vector[row];
  }
  dpd_file2_mat_wrt(&mu);
  dpd_file2_close(&mu);

  dpd_buf4_close(&A);
  psio_close(PSIF_MO_HESS, 1);
}
