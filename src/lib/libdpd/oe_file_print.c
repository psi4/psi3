#include <stdio.h>
#include "dpd.h"

/* dpd_oe_file_print(): Prints out data for all irreps of a dpd
** one-electron file.
**
** Arguments:
**   struct oe_dpdfile *File: A pointer to the dpdfile to be printed.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_oe_file_print(struct oe_dpdfile *File, FILE *outfile)
{
  int h;

  fprintf(outfile, "\n\tOE DPD File: %s\n", File->label);
  dpd_oe_params_print(File->params, outfile);

  for(h=0; h < File->params->nirreps; h++) {
      dpd_oe_file_mat_irrep_init(File, h);
      dpd_oe_file_mat_irrep_rd(File, h, 1, outfile);
      dpd_oe_file_mat_irrep_close(File, h);
    }

  return 0;

}
