#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#define EXTERN
#include "globals.h"

void probable(void)
{
  double *exps;
  
  chkpt_init();
  chkpt_rd_exps();
  chkpt_close();

}
