#define EXTERN
#include "globals.h"

void sortone(void)
{
  if(params.ref == 0 || params.ref == 1) sortone_ROHF();
  else if(params.ref == 2) sortone_UHF();
}
