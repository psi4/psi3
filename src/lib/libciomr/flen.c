/*!
  \file flen.c
*/
 
#include "includes.h"
#include "pointers.h"
#include "types.h"

extern PSI_FPTR iosize_(int *);

/*!
** flen(): Get the number of bytes for file number 'itape'.
*/
PSI_FPTR flen(int itape)
   {
      PSI_FPTR length;
      length = iosize_(&itape);
      return(length);
   }
