#include "psio.h"

/*
** PSIO_OPEN_CHECK(): Check to see if a given PSI direct access file
** is already open.
**
** Argument: unit = the PSI unit number.
*/

int psio_open_check(ULI unit)
{
  psio_ud *this_unit;

  this_unit = &(psio_unit[unit]);

  if(this_unit->vol[0].stream != -1) return 1;
  else return 0;
}
