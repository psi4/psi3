#include "dpd.h"

int dpd_file_mat_init(struct dpdfile *File)
{
  int h;

  for(h=0; h < File->params->nirreps; h++)
      dpd_file_mat_irrep_init(File, h);

  return 0;
}
