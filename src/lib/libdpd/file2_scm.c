#include <stdio.h>
#include "dpd.h"
#include <qt.h>

int dpd_file2_scm(dpdfile2 *InFile, double alpha)
{
  int h, nirreps, new_file2;
  int row, col, length;
  double *X;

  nirreps = InFile->params->nirreps;

  dpd_file2_mat_init(InFile);

  /* Look first for the TOC entry on disk */
  if(psio_tocscan(InFile->filenum, InFile->label) == NULL)
     new_file2 = 1;
  else new_file2 = 0;

  if(!new_file2) dpd_file2_mat_rd(InFile);

  for(h=0; h < nirreps; h++) {

      length = InFile->params->rowtot[h] * InFile->params->coltot[h];
      if(length) { 
         X = &(InFile->matrix[h][0][0]);
         C_DSCAL(length, alpha, X, 1);
       }
  }

  dpd_file2_mat_wrt(InFile);
  dpd_file2_mat_close(InFile);

  return 0;
}
