#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libfile30/file30.h>
#define EXTERN
#include "globals.h"

void probable(void)
{
  double *exps;
  
  file30_init();
  file30_rd_exps();
  file30_close();

}
