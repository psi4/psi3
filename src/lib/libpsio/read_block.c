#include <unistd.h>
#include <string.h>
#include "psio.h"

int psio_read_block(ULI unit, char *key, char *buffer, ULI blksiz,
	            ULI start_blk, ULI end_blk)
{
  ULI size, shift;
  psio_ud *this_unit;
  psio_address sadd, eadd;
  psio_tocentry *this_entry;

  this_unit = &(psio_unit[unit]);
  
  /* Find the entry in the TOC */
  this_entry = psio_tocscan(unit, key);

  if(this_entry == NULL) psio_error(unit,PSIO_ERROR_NOTOCENT);
  else {
      size = (end_blk - start_blk + 1) * blksiz; /* The total buffer size */
      shift = start_blk * blksiz; /* Number of bytes to shift from start */

      /* Compute the starting page and offset for the block */
      sadd = psio_get_address(this_entry->sadd, shift);

      /* Make sure the block starts and ends within the entry */
      if((sadd.page > this_entry->eadd.page))
	  psio_error(unit,PSIO_ERROR_BLKSTART);
      else if((sadd.page == this_entry->eadd.page) &&
	      (sadd.offset > this_entry->eadd.offset))
	  psio_error(unit,PSIO_ERROR_BLKSTART);

      eadd = psio_get_address(sadd, size);
      if((eadd.page > this_entry->eadd.page))
	  psio_error(unit,PSIO_ERROR_BLKEND);
      else if((eadd.page == this_entry->eadd.page) &&
	      (eadd.offset > this_entry->eadd.offset))
	  psio_error(unit,PSIO_ERROR_BLKEND);
    }

  /* Now read the actual data from the unit */
  psio_rw(unit, buffer, sadd, size, 0);

  return(0);
}
