#include "dpd.h"

int dpd_oe_file_mat_close(struct oe_dpdfile *File)
{
  int h;

  for(h=0; h < File->params->nirreps; h++)
      dpd_oe_file_mat_irrep_close(File, h);

  return 0;
}
