#include <stdio.h>
#include <math.h>
#include <ip_libv1.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

int converged(void)
{
  int row,col,h,nirreps;
  double rms=0.0;
  struct oe_dpdfile L1, L1old;
  struct dpdbuf L2, L2old;

  nirreps = moinfo.nirreps;

  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "New LIA", 0, outfile);
  dpd_oe_file_mat_init(&L1);
  dpd_oe_file_mat_rd(&L1, 0, outfile);
  dpd_oe_file_init(&L1old, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_mat_init(&L1old);
  dpd_oe_file_mat_rd(&L1old, 0, outfile);
  for(h=0; h < nirreps; h++)
      for(row=0; row < L1.params->rowtot[h]; row++)
	  for(col=0; col < L1.params->coltot[h]; col++)
	      rms += (L1.matrix[h][row][col] - L1old.matrix[h][row][col]) *
		     (L1.matrix[h][row][col] - L1old.matrix[h][row][col]);

  dpd_oe_file_mat_close(&L1);
  dpd_oe_file_close(&L1);
  dpd_oe_file_mat_close(&L1old);
  dpd_oe_file_close(&L1old);

  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "New Lia", 0, outfile);
  dpd_oe_file_mat_init(&L1);
  dpd_oe_file_mat_rd(&L1, 0, outfile);
  dpd_oe_file_init(&L1old, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_oe_file_mat_init(&L1old);
  dpd_oe_file_mat_rd(&L1old, 0, outfile);
  for(h=0; h < nirreps; h++)
      for(row=0; row < L1.params->rowtot[h]; row++)
	  for(col=0; col < L1.params->coltot[h]; col++)
	      rms += (L1.matrix[h][row][col] - L1old.matrix[h][row][col]) *
		     (L1.matrix[h][row][col] - L1old.matrix[h][row][col]);

  dpd_oe_file_mat_close(&L1);
  dpd_oe_file_close(&L1);
  dpd_oe_file_mat_close(&L1old);
  dpd_oe_file_close(&L1old);

  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "New LIJAB", 0, outfile);
  dpd_buf_init(&L2old, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&L2, h);
      dpd_buf_mat_irrep_rd(&L2, h, 0, outfile);
      dpd_buf_mat_irrep_init(&L2old, h);
      dpd_buf_mat_irrep_rd(&L2old, h, 0, outfile);
      for(row=0; row < L2.params->rowtot[h]; row++)
	  for(col=0; col < L2.params->coltot[h]; col++)
	      rms += (L2.matrix[h][row][col] - L2old.matrix[h][row][col]) *
		     (L2.matrix[h][row][col] - L2old.matrix[h][row][col]);
      dpd_buf_mat_irrep_close(&L2, h);
      dpd_buf_mat_irrep_close(&L2old, h);
    }
  dpd_buf_close(&L2old);
  dpd_buf_close(&L2);

  dpd_buf_init(&L2, CC_LAMPS, 2, 7, 2, 7, 0, "New Lijab", 0, outfile);
  dpd_buf_init(&L2old, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&L2, h);
      dpd_buf_mat_irrep_rd(&L2, h, 0, outfile);
      dpd_buf_mat_irrep_init(&L2old, h);
      dpd_buf_mat_irrep_rd(&L2old, h, 0, outfile);
      for(row=0; row < L2.params->rowtot[h]; row++)
	  for(col=0; col < L2.params->coltot[h]; col++)
	      rms += (L2.matrix[h][row][col] - L2old.matrix[h][row][col]) *
		     (L2.matrix[h][row][col] - L2old.matrix[h][row][col]);
      dpd_buf_mat_irrep_close(&L2, h);
      dpd_buf_mat_irrep_close(&L2old, h);
    }
  dpd_buf_close(&L2old);
  dpd_buf_close(&L2);

  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb", 0, outfile);
  dpd_buf_init(&L2old, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&L2, h);
      dpd_buf_mat_irrep_rd(&L2, h, 0, outfile);
      dpd_buf_mat_irrep_init(&L2old, h);
      dpd_buf_mat_irrep_rd(&L2old, h, 0, outfile);
      for(row=0; row < L2.params->rowtot[h]; row++)
	  for(col=0; col < L2.params->coltot[h]; col++)
	      rms += (L2.matrix[h][row][col] - L2old.matrix[h][row][col]) *
		     (L2.matrix[h][row][col] - L2old.matrix[h][row][col]);
      dpd_buf_mat_irrep_close(&L2, h);
      dpd_buf_mat_irrep_close(&L2old, h);
    }
  dpd_buf_close(&L2old);
  dpd_buf_close(&L2);

  rms = sqrt(rms);
  moinfo.conv = rms;

  if(rms < params.convergence) return 1;
  else return 0;
}
