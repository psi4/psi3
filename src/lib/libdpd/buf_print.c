#include <stdio.h>
#include "dpd.h"

/* dpd_buf_print(): Prints out data for all irreps of a dpd
** two-electron buffer.
**
** Arguments:
**   struct dpdbuf *Buf: A pointer to the dpdbuf to be printed.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_buf_print(struct dpdbuf *Buf, FILE *outfile)
{
  int h;

  fprintf(outfile, "\n\tDPD Buffer for file: %s\n", Buf->file.label);
  dpd_params_print(Buf->params, outfile);

  for(h=0; h < Buf->params->nirreps; h++) {
      dpd_buf_mat_irrep_init(Buf, h);
      dpd_buf_mat_irrep_rd(Buf, h, 1, outfile);
      dpd_buf_mat_irrep_close(Buf, h);
    }

  return 0;

}
