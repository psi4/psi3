#include "iomrparam.h"
#include "includes.h"

extern int io_locate(FILE *, char[]);

int oldstyleinput(void)
{
  FILE *input;
  int ierr;

  input = fopen("input.dat","r");
  ierr = io_locate(input,"# FILES ##");
  fclose(input);
  if (ierr == 0) return(1);
  return 0;
  }
