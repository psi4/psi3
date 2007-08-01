/*!
   \file toclast.c
   \ingroup (PSIO)
*/

#include <libpsio/psio.h>

psio_tocentry* __psio_toclast(psio_lib* Lib, unsigned int unit)
{
  psio_tocentry *this_entry = Lib->psio_unit[unit].toc;

  while(this_entry->next != NULL) this_entry = this_entry->next;

  return(this_entry);
}

/*!
** PSIO_TOCLAST(): Returns the last TOC entry.
**
** \ingroup (PSIO)
*/

psio_tocentry* psio_toclast(unsigned int unit)
{
  return __psio_toclast(_default_psio_lib_,unit);
}
