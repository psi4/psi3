/*!
** \file toclen.c
*/

#include "psio.h"

/*!
** PSIO_TOCLEN(): Compute the length of the TOC for a given unit.
*/

ULI psio_toclen(ULI unit)
{
  ULI toclen=0;
  psio_tocentry *this_entry;

  this_entry = psio_unit[unit].toc;

  while(this_entry != NULL) {
      ++toclen;
      this_entry = this_entry->next;
    }

  return(toclen);
}
