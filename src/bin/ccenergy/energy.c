#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

double rhf_energy(void);
double rohf_energy(void);

double energy(void)
{
  if(params.ref == 0) return(rhf_energy());
  else if(params.ref == 1) return(rohf_energy());
}
