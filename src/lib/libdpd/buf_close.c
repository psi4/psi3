#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include "dpd.h"

/* dpd_buf_close(): Closes a dpd two-electron buffer.
**
** Arguments:
**   struct dpdbuf *Buf: A pointer to the dpdbuf to be closed.
*/

int dpd_buf_close(struct dpdbuf *Buf)
{
  int nirreps;

  nirreps = Buf->params->nirreps;
  
  dpd_file_close(&(Buf->file));

  free(Buf->matrix);

  free_int_matrix(Buf->shift.rowtot,nirreps);
  free_int_matrix(Buf->shift.coltot,nirreps);
  free(Buf->shift.matrix);

  return 0;
}
