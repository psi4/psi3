#include "iomrparam.h"
#include "includes.h"
#include "pointers.h"

extern PSI_FPTR sec2i(int);

/* sets file pointer for unit to sector address */

void rsetsa(unit,address)
   int unit;
   int address;
   {
      PSI_FPTR ipos;

      ipos = (PSI_FPTR) sec2i(--address);
      ptr.wptr[unit]=ipos;
    }
