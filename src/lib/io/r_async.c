
#include "includes.h"
#include "param.h"
#include "types.h"

r_async_t *
r_async_ioopen(param,unit)
char *param;
int unit;
{
  fprintf(stderr,"no r_async io yet\n");
  io_q_abort();
  return 0;
  }

void
r_async_ioclos(ud,status)
r_async_t *ud;
int status;
{
  fprintf(stderr,"no r_async io yet\n");
  io_q_abort();
  }

void
r_async_iordr(ud)
r_async_t *ud;
{
  fprintf(stderr,"no r_async io yet\n");
  io_q_abort();
  }

void
r_async_iowrr(ud)
r_async_t *ud;
{
  fprintf(stderr,"no r_async io yet\n");
  io_q_abort();
  }
