
#include "includes.h"
#include "param.h"
#include "types.h"

s_async_t *
s_async_ioopen(param,unit)
char *param;
int unit;
{
  fprintf(stderr,"no s_async io yet\n");
  io_q_abort();
  }

void
s_async_ioclos(ud,status)
s_async_t *ud;
int status;
{
  fprintf(stderr,"no s_async io yet\n");
  io_q_abort();
  }

void
s_async_iordr(ud)
s_async_t *ud;
{
  fprintf(stderr,"no s_async io yet\n");
  io_q_abort();
  }

void
s_async_iowrr(ud)
s_async_t *ud;
{
  fprintf(stderr,"no s_async io yet\n");
  io_q_abort();
  }
