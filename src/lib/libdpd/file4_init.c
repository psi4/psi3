#include <stdio.h>
#include <stdlib.h>
#include <libpsio/psio.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_file4_init(): Prepares a dpd four-index file on disk for
** reading/writing.
**
** Arguments:
**   dpdfile4 *File: A pointer to the dpdfile4 to be initialized.
**   int filenum: The PSI unit number for this file.
**   int irrep: The irrep of this quantity.
**   int pqnum: The index combination for the bra indices for the
**              data as it will be stored on disk.
**   int rsnum: The index combination for the ket indices for the
**              data as it will be stored on disk.
**   char *label: A string labelling for this buffer.
*/

int dpd_file4_init(dpdfile4 *File, int filenum, int irrep, int pqnum,
		   int rsnum,  char *label)
{
  int i;
  unsigned int priority;
  struct dpd_file4_cache_entry *this_entry;
  
  File->dpdnum = dpd_default;
  File->params = &(dpd_list[dpd_default].params4[pqnum][rsnum]);

  strcpy(File->label,label);
  File->filenum = filenum;
  File->my_irrep = irrep;

  this_entry = dpd_file4_cache_scan(filenum, irrep, pqnum, rsnum, label, dpd_default);
  if(this_entry != NULL) {
      File->incore = 1;
      File->matrix = this_entry->matrix;
    }
  else {
   File->incore = 0;
   File->matrix = (double ***) malloc(File->params->nirreps*sizeof(double **));
    }

  /* Construct logical subfile pointers */
  File->lfiles = (psio_address *) malloc(File->params->nirreps *
					 sizeof(psio_address));
  File->lfiles[0] = PSIO_ZERO;
  for(i=1; i < File->params->nirreps; i++) {
    File->lfiles[i] = psio_get_address(File->lfiles[i-1],
				       (File->params->rowtot[i-1] *
					File->params->coltot[(i-1)^irrep] *
					sizeof(double)));
  }

  /* Put this file4 into cache if requested */
  if(dpd_main.cachefiles[filenum] && dpd_main.cachelist[pqnum][rsnum]) 
      {
      /* Get the file4's cache priority */
      if(dpd_main.cachetype == 1)
          priority = dpd_file4_cache_get_priority(File);
      else priority = 0;

      dpd_file4_cache_add(File, priority); 

      /* Make sure this cache entry can't be deleted until we're done */
      dpd_file4_cache_lock(File);
    }

  return 0;
}
