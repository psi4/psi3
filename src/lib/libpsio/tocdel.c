/*!
   \file tocdel.c
   \ingroup (PSIO)
*/

#include <stdlib.h>
#include "psio.h"

/*!
** PSIO_TOCDEL(): Delete an entry with identifier key from TOC.
** 
** \param unit = file number
** \param key = entry to delete from TOC
**
** \ingroup (PSIO)
*/

int psio_tocdel(ULI unit, char *key)
{
  psio_tocentry *this_entry, *last_entry;

  /* Check the key length first */
  if((strlen(key)+1) > PSIO_KEYLEN) psio_error(unit,PSIO_ERROR_KEYLEN);

  this_entry = psio_tocscan(unit, key);
  if(this_entry == NULL) psio_error(unit,PSIO_ERROR_NOTOCENT);

  last_entry = this_entry->last;
  last_entry->next = this_entry->next;

  free(this_entry);

  return(0);
}
