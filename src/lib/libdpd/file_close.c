#include <stdio.h>
#include <stdlib.h>
#include "dpd.h"

int dpd_file_close(struct dpdfile *File)
{
  free(File->lfiles);
  free(File->matrix);
  return 0;
}
