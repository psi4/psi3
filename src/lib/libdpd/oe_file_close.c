#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_oe_file_close(): Closes a dpd one-electron file.
**
** Arguments:
**   struct oe_dpdfile *File: A pointer to the file to be closed.
*/

int dpd_oe_file_close(struct oe_dpdfile *File)
{
  free(File->lfiles);
  free(File->matrix);

  return 0;
}
