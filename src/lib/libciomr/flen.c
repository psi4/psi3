#include "includes.h"
#include "pointers.h"
#include "types.h"

extern PSI_FPTR iosize_(int *);

/* get the number of bytes in a unit */

PSI_FPTR flen(itape)
   int itape;

   {
      PSI_FPTR length;
      length = iosize_(&itape);
      return(length);
   }
