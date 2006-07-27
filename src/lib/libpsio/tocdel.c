/*!
   \file tocdel.c
   \ingroup (PSIO)
*/

#include <stdlib.h>
#include <string.h>
#include "psio.h"

/*!
** PSIO_TOCDEL(): Delete an entry with identifier key from TOC.
** 
** \param unit = file number
** \param key = entry to delete from TOC
**
** \ingroup (PSIO)
*/

int psio_tocdel(unsigned int unit, char *key)
{
  psio_tocentry *this_entry, *last_entry, *next_entry;

  this_entry = psio_tocscan(unit, key);
  if(this_entry == NULL) {
    fprintf(stderr, "PSIO_ERROR: Can't find TOC Entry %s\n", key);
    psio_error(unit,PSIO_ERROR_NOTOCENT);
  }

  last_entry = this_entry->last;
  if(last_entry != NULL) last_entry->next = this_entry->next;

  next_entry = this_entry->next;
  if(next_entry != NULL) next_entry->last = this_entry->last;

  free(this_entry);

  /* Dump the new TOC to the file */
  psio_tocwrite(unit);

  return(1);
}
