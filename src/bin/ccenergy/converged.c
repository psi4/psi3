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
  struct oe_dpdfile T1, T1old;
  struct dpdbuf T2, T2old;

  nirreps = moinfo.nirreps;

  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "New tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);
  dpd_oe_file_init(&T1old, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1old);
  dpd_oe_file_mat_rd(&T1old, 0, outfile);
  for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	  for(col=0; col < T1.params->coltot[h]; col++)
	      rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
		     (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&T1old);
  dpd_oe_file_close(&T1old);

  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "New tia", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);
  dpd_oe_file_init(&T1old, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1old);
  dpd_oe_file_mat_rd(&T1old, 0, outfile);
  for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	  for(col=0; col < T1.params->coltot[h]; col++)
	      rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
		     (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&T1old);
  dpd_oe_file_close(&T1old);

  dpd_buf_init(&T2, CC_TAMPS, 2, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&T2old, CC_TAMPS, 2, 7, 2, 7, 0, "tIJAB", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2, h);
      dpd_buf_mat_irrep_rd(&T2, h, 0, outfile);
      dpd_buf_mat_irrep_init(&T2old, h);
      dpd_buf_mat_irrep_rd(&T2old, h, 0, outfile);
      for(row=0; row < T2.params->rowtot[h]; row++)
	  for(col=0; col < T2.params->coltot[h]; col++)
	      rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
		     (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf_mat_irrep_close(&T2, h);
      dpd_buf_mat_irrep_close(&T2old, h);
    }
  dpd_buf_close(&T2old);
  dpd_buf_close(&T2);

  dpd_buf_init(&T2, CC_TAMPS, 2, 7, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&T2old, CC_TAMPS, 2, 7, 2, 7, 0, "tijab", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2, h);
      dpd_buf_mat_irrep_rd(&T2, h, 0, outfile);
      dpd_buf_mat_irrep_init(&T2old, h);
      dpd_buf_mat_irrep_rd(&T2old, h, 0, outfile);
      for(row=0; row < T2.params->rowtot[h]; row++)
	  for(col=0; col < T2.params->coltot[h]; col++)
	      rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
		     (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf_mat_irrep_close(&T2, h);
      dpd_buf_mat_irrep_close(&T2old, h);
    }
  dpd_buf_close(&T2old);
  dpd_buf_close(&T2);

  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);
  dpd_buf_init(&T2old, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&T2, h);
      dpd_buf_mat_irrep_rd(&T2, h, 0, outfile);
      dpd_buf_mat_irrep_init(&T2old, h);
      dpd_buf_mat_irrep_rd(&T2old, h, 0, outfile);
      for(row=0; row < T2.params->rowtot[h]; row++)
	  for(col=0; col < T2.params->coltot[h]; col++)
	      rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
		     (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_buf_mat_irrep_close(&T2, h);
      dpd_buf_mat_irrep_close(&T2old, h);
    }
  dpd_buf_close(&T2old);
  dpd_buf_close(&T2);

  rms = sqrt(rms);
  moinfo.conv = rms;

  if(rms < params.convergence) return 1;
  else return 0;
}
