#include "psio.h"

/*
** PSIO_TOCLAST(): Returns the last TOC entry.
*/

psio_tocentry *psio_toclast(ULI unit)
{
  psio_tocentry *this_entry;

  this_entry = psio_unit[unit].toc;

  while(this_entry->next != NULL) this_entry = this_entry->next;

  return(this_entry);
}

