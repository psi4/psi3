#include "psio.h"

/*
** PSIO_GET_ENTRY(): Returns the ending relative address of the tocentry
** for the given keyword.
*/

psio_address psio_get_entry_end(ULI unit, char *key)
{
  psio_tocentry *this_entry;

  this_entry = psio_tocscan(unit,key);

  if(this_entry ==  NULL) return(PSIO_ZERO);
  else return(this_entry->eadd);
}
