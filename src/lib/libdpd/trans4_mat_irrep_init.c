#include <stdio.h>
#include <libciomr/libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

int dpd_trans4_mat_irrep_init(dpdtrans4 *Trans, int irrep)
{
  int rowtot, coltot;

  rowtot = Trans->buf.params->coltot[irrep];
  coltot = Trans->buf.params->rowtot[irrep];

  if(rowtot * coltot) Trans->matrix[irrep] = dpd_block_matrix(rowtot,coltot);
  

  return 0;
}
