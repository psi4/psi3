#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/** Wabei intermediates are stored here as (ei,ab) **/

void Wabei_build(void)
{
  if(params.ref == 0) Wabei_RHF();
  else if(params.ref == 1) Wabei_ROHF();
  else if(params.ref == 2) {
    WABEI_UHF();
    Wabei_UHF();
    WAbEi_UHF();
    WaBeI_UHF();
  }
}
