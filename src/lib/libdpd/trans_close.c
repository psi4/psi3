#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_trans_close(struct dpdtrans *Trans)
{
  int nirreps;

  nirreps = Trans->buf.params->nirreps;

  free(Trans->matrix);
  
  free_int_matrix(Trans->shift.rowtot, nirreps);
  free_int_matrix(Trans->shift.coltot, nirreps);
  free(Trans->shift.matrix);

  return 0;

}
