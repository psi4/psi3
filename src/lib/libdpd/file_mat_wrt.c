#include <stdio.h>
#include "dpd.h"

int dpd_file_mat_wrt(struct dpdfile *File, int print_flag, FILE *outfile)
{
  int h;

  for(h=0; h < File->params->nirreps; h++)
      dpd_file_mat_irrep_wrt(File, h, print_flag, outfile);

  return 0;
}
