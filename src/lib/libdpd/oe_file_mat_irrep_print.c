#include <stdio.h>
#include <stdlib.h>
#include "dpd.h"

int dpd_oe_file_mat_irrep_print(struct oe_dpdfile *File, int irrep, FILE *outfile)
{
  fprintf(outfile, "\n\tFile %3d OE DPD: %s\n", File->filenum,
	  File->label);
  fprintf(outfile,   "\tMatrix for Irrep %1d\n", irrep);
  fprintf(outfile,   "\t----------------------------------------\n");

  dpd_oe_params_mat_irrep_print(File->matrix[irrep], File->params,
				irrep, outfile);

  return 0;

}
