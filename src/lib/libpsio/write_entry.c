/*!
** \file write_entry.c
*/

#include "psio.h"

/*!
** PSIO_WRITE_ENTRY()
*/
int psio_write_entry(ULI unit, char *key, char *buffer, ULI size)
{
  psio_address end;
  return psio_write(unit, key, buffer, size, PSIO_ZERO, &end);
}

