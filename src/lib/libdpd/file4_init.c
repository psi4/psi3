#include <stdio.h>
#include <stdlib.h>
#include <psio.h>
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
  struct dpd_file4_cache_entry *this_entry;
  
  strcpy(File->label,label);
  File->filenum = filenum;
  File->params = &(dpd_default->params4[pqnum][rsnum]);
  File->my_irrep = irrep;

  this_entry = dpd_file4_cache_scan(filenum, irrep, pqnum, rsnum, label);
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
  for(i=1; i < File->params->nirreps; i++)
    File->lfiles[i] = psio_get_address(File->lfiles[i-1],
				       (File->params->rowtot[i-1] *
					File->params->coltot[(i-1)^irrep] *
					sizeof(double)));

  /* Put this file4 into cache if requested */
  if(dpd_default->cachefiles[filenum] && dpd_default->cachelist[pqnum][rsnum]) 
      dpd_file4_cache_add(File); 
/*
  dpd_file4_cache_add(File); 
*/

  return 0;
}
