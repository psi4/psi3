#include <stdio.h>
#include "dpd.h"
#include <qt.h>

int dpd_file2_scm(dpdfile2 *InFile, double alpha)
{
  int h, nirreps;
  int row, col, length;
  double *X;

  nirreps = InFile->params->nirreps;

  dpd_file2_mat_init(InFile);
  dpd_file2_mat_rd(InFile);

  for(h=0; h < nirreps; h++) {

      length = InFile->params->rowtot[h] * InFile->params->coltot[h];
      if(length) { 
         X = &(InFile->matrix[h][0][0]);
         C_DSCAL(length, alpha, X, 1);
       }

/*
      for(row=0; row < InFile->params->rowtot[h]; row++) {
	  for(col=0; col < InFile->params->coltot[h]; col++) {
	      InFile->matrix[h][row][col] *= alpha;
	    }
	}
*/
    }

  dpd_file2_mat_wrt(InFile);
  dpd_file2_mat_close(InFile);

  return 0;
}
