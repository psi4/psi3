/*!
   \file write_entry.c
   \ingroup (PSIO)
*/

#include "psio.h"

/*!
** PSIO_WRITE_ENTRY()
**
** \ingroup (PSIO)
*/
int psio_write_entry(unsigned int unit, char *key, char *buffer, ULI size)
{
  psio_address end;
  return psio_write(unit, key, buffer, size, PSIO_ZERO, &end);
}

