/*!
   \file tocrename.c
   \ingroup (PSIO)
*/

#include <stdlib.h>
#include <string.h>
#include "psio.h"

/*!
** PSIO_TOCRENAME(): Rename an entry with identifier key from TOC.
** 
** \param unit = file number
** \param key = entry to rename
** \param newkey = new name of entry
**
** \ingroup (PSIO)
*/

int psio_tocrename(ULI unit, char *key, char *newkey)
{
  psio_tocentry *this_entry;

  /* Check the key length first */
  if((strlen(key)+1) > PSIO_KEYLEN) psio_error(unit,PSIO_ERROR_KEYLEN);
  if((strlen(newkey)+1) > PSIO_KEYLEN) psio_error(unit,PSIO_ERROR_KEYLEN);

  this_entry = psio_tocscan(unit, key);
  if(this_entry == NULL) {
      fprintf(stderr, "PSIO_ERROR: Can't find TOC Entry %s\n", key);
      psio_error(unit,PSIO_ERROR_NOTOCENT);
  }

  /* replace the old key with the new one */
  strcpy(this_entry->key,newkey);

  return(0);
}
