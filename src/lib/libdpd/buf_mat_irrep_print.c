#include <stdio.h>
#include "dpd.h"

int dpd_buf_mat_irrep_print(struct dpdbuf *Buf, int irrep, FILE *outfile)
{
  fprintf(outfile, "\n\tDPD Buffer for file: %s\n", Buf->file.label);
  fprintf(outfile,   "\tMatrix for Irrep %1d\n", irrep);
  fprintf(outfile,   "\t-----------------------------------------\n");

  dpd_mat_irrep_print(Buf->matrix[irrep], Buf->params, irrep, outfile);

  return 0;

}
