
#include "includes.h"
#include "param.h"
#include "types.h"

ram_t *
ram_ioopen(param,unit)
char *param;
int unit;
{
  fprintf(stderr,"no ram io yet\n");
  io_q_abort();
  }

void
ram_ioclos(ud,status)
ram_t *ud;
int status;
{
  fprintf(stderr,"no ram io yet\n");
  io_q_abort();
  }

void
ram_iordr(ud)
ram_t *ud;
{
  fprintf(stderr,"no ram io yet\n");
  io_q_abort();
  }

void
ram_iowrr(ud)
ram_t *ud;
{
  fprintf(stderr,"no ram io yet\n");
  io_q_abort();
  }
