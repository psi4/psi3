/*!
** \file tocclean.c
*/

#include "psio.h"

/*!
** PSIO_TOCCLEAN(): Delete all TOC entries after the given key.
** If a blank key is given, the entire TOC will be wiped.
*/

int psio_tocclean(ULI unit, char *key)
{
  psio_tocentry *this_entry, *last_entry, *prev_entry;

  /* Check the key length first */
  if((strlen(key)+1) > PSIO_KEYLEN) psio_error(unit,PSIO_ERROR_KEYLEN);

  this_entry = psio_tocscan(unit, key);
  if(this_entry == NULL) {
      if(!strcmp(key,"")) this_entry = psio_unit[unit].toc;
      else psio_error(unit,PSIO_ERROR_NOTOCENT);
    }
  else this_entry = this_entry->next;

  /* Get the end of the TOC and work backwards */
  last_entry = psio_toclast(unit);

  while((last_entry != this_entry) && (last_entry != NULL)) { 
      /* Now free all the remaining members */
      prev_entry = last_entry->last;
      free(last_entry);
      last_entry = prev_entry;
    }
  return(0);
}
