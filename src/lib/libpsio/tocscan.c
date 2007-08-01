/*!
   \file tocscan.c
   \ingroup (PSIO)
*/

#include <string.h>
#include <libpsio/psio.h>

psio_tocentry* __psio_tocscan(psio_lib* Lib, unsigned int unit, char *key)
{
  psio_tocentry *this_entry;

  if(key == NULL) return(NULL);

  if((strlen(key)+1) > PSIO_KEYLEN) psio_error(unit,PSIO_ERROR_KEYLEN);

  this_entry = Lib->psio_unit[unit].toc;

  while(this_entry != NULL) {
      if(!strcmp(this_entry->key,key)) return(this_entry);
      this_entry = this_entry->next;
    }

  return(NULL);
}

/*!
** PSIO_TOCSCAN(): Scans the TOC for a particular keyword and returns either
** a pointer to the entry or NULL to the caller.
**
** \ingroup (PSIO)
*/

psio_tocentry *psio_tocscan(unsigned int unit, char *key)
{
  return __psio_tocscan(_default_psio_lib_,unit,key);
}
