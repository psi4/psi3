/*! \file file2_close.c
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_file2_close(): Closes a two-index dpd file.
**
** Arguments:
**   dpdfile2 *File: A pointer to the file to be closed.
*/

int dpd_file2_close(dpdfile2 *File)
{
  free(File->lfiles);

  if(!File->incore) free(File->matrix);
  else File->matrix = NULL;

  return 0;
}
