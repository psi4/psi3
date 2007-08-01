/*!
   \file open_check.c
   \ingroup (PSIO)
*/

#include <libpsio/psio.h>

int __psio_open_check(psio_lib* Lib, unsigned int unit)
{
  psio_ud *this_unit;

  this_unit = &(Lib->psio_unit[unit]);

  if(this_unit->vol[0].stream != -1) return 1;
  else return 0;
}

/*!
** PSIO_OPEN_CHECK(): Check to see if a given PSI direct access file
** is already open.
**
** \param unit = the PSI unit number.
**
** \ingroup (PSIO)
*/

int psio_open_check(unsigned int unit)
{
  return __psio_open_check(_default_psio_lib_,unit);
}
