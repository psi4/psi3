#include "psio.h"

int psio_write_entry(ULI unit, char *key, char *buffer, ULI size)
{
  psio_address end;
  return psio_write(unit, key, buffer, size, PSIO_ZERO, &end);
}

