#include <stdio.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

int converged(void)
{
  int row,col,h,nirreps;
  double rms=0.0;
  dpdfile2 L1, L1old;
  dpdbuf4 L2, L2old;

  nirreps = moinfo.nirreps;

  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "New LIA");
  dpd_file2_mat_init(&L1);
  dpd_file2_mat_rd(&L1);
  dpd_file2_init(&L1old, CC_OEI, 0, 0, 1, "LIA");
  dpd_file2_mat_init(&L1old);
  dpd_file2_mat_rd(&L1old);
  for(h=0; h < nirreps; h++)
      for(row=0; row < L1.params->rowtot[h]; row++)
	  for(col=0; col < L1.params->coltot[h]; col++)
	      rms += (L1.matrix[h][row][col] - L1old.matrix[h][row][col]) *
		     (L1.matrix[h][row][col] - L1old.matrix[h][row][col]);

  dpd_file2_mat_close(&L1);
  dpd_file2_close(&L1);
  dpd_file2_mat_close(&L1old);
  dpd_file2_close(&L1old);

  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "New Lia");
  dpd_file2_mat_init(&L1);
  dpd_file2_mat_rd(&L1);
  dpd_file2_init(&L1old, CC_OEI, 0, 0, 1, "Lia");
  dpd_file2_mat_init(&L1old);
  dpd_file2_mat_rd(&L1old);
  for(h=0; h < nirreps; h++)
      for(row=0; row < L1.params->rowtot[h]; row++)
	  for(col=0; col < L1.params->coltot[h]; col++)
	      rms += (L1.matrix[h][row][col] - L1old.matrix[h][row][col]) *
		     (L1.matrix[h][row][col] - L1old.matrix[h][row][col]);

  dpd_file2_mat_close(&L1);
  dpd_file2_close(&L1);
  dpd_file2_mat_close(&L1old);
  dpd_file2_close(&L1old);

  dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 7, 2, 7, 0, "New LIJAB");
  dpd_buf4_init(&L2old, CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB");
  for(h=0; h < nirreps; h++) {
  dpd_buf4_mat_irrep_init(&L2, h); 0,
      dpd_buf4_mat_irrep_rd(&L2, h);
  dpd_buf4_mat_irrep_init(&L2old, h); 0,
      dpd_buf4_mat_irrep_rd(&L2old, h);
      for(row=0; row < L2.params->rowtot[h]; row++)
	  for(col=0; col < L2.params->coltot[h]; col++)
	      rms += (L2.matrix[h][row][col] - L2old.matrix[h][row][col]) *
		     (L2.matrix[h][row][col] - L2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&L2, h);
      dpd_buf4_mat_irrep_close(&L2old, h);
    }
  dpd_buf4_close(&L2old);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 7, 2, 7, 0, "New Lijab");
  dpd_buf4_init(&L2old, CC_LAMPS, 0, 2, 7, 2, 7, 0, "Lijab");
  for(h=0; h < nirreps; h++) {
  dpd_buf4_mat_irrep_init(&L2, h); 0,
      dpd_buf4_mat_irrep_rd(&L2, h);
  dpd_buf4_mat_irrep_init(&L2old, h); 0,
      dpd_buf4_mat_irrep_rd(&L2old, h);
      for(row=0; row < L2.params->rowtot[h]; row++)
	  for(col=0; col < L2.params->coltot[h]; col++)
	      rms += (L2.matrix[h][row][col] - L2old.matrix[h][row][col]) *
		     (L2.matrix[h][row][col] - L2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&L2, h);
      dpd_buf4_mat_irrep_close(&L2old, h);
    }
  dpd_buf4_close(&L2old);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");
  dpd_buf4_init(&L2old, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  for(h=0; h < nirreps; h++) {
  dpd_buf4_mat_irrep_init(&L2, h); 0,
      dpd_buf4_mat_irrep_rd(&L2, h);
  dpd_buf4_mat_irrep_init(&L2old, h); 0,
      dpd_buf4_mat_irrep_rd(&L2old, h);
      for(row=0; row < L2.params->rowtot[h]; row++)
	  for(col=0; col < L2.params->coltot[h]; col++)
	      rms += (L2.matrix[h][row][col] - L2old.matrix[h][row][col]) *
		     (L2.matrix[h][row][col] - L2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&L2, h);
      dpd_buf4_mat_irrep_close(&L2old, h);
    }
  dpd_buf4_close(&L2old);
  dpd_buf4_close(&L2);

  rms = sqrt(rms);
  moinfo.conv = rms;

  if(rms < params.convergence) return 1;
  else return 0;
}
