#include <stdio.h>
#include "dpd.h"

int dpd_file_print(struct dpdfile *File, FILE *outfile)
{
  int h;

  fprintf(outfile, "\n\tDPD File: %s\n", File->label);
  dpd_params_print(File->params, outfile);

  for(h=0; h < File->params->nirreps; h++) {
      dpd_file_mat_irrep_init(File, h);
      dpd_file_mat_irrep_rd(File, h, 1, outfile);
      dpd_file_mat_irrep_close(File, h);
    }

  return 0;

}
