#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void x_Gijab(void) {

  if ( (params.ref == 0) || (params.ref == 1) ) 
    x_Gijab_ROHF();
    /*
  else if (params.ref == 2)
    x_Gijab_UHF();
    */
    return;
}
