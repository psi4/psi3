/*!
   \file toclen.c
   \ingroup (PSIO)
*/

#include "psio.h"

/*!
** PSIO_TOCLEN(): Compute the length of the TOC for a given unit.
**
** \ingroup (PSIO)
*/

unsigned int psio_toclen(unsigned int unit)
{
  unsigned int toclen=0;
  psio_tocentry *this_entry;

  this_entry = psio_unit[unit].toc;

  while(this_entry != NULL) {
      ++toclen;
      this_entry = this_entry->next;
    }

  return(toclen);
}
