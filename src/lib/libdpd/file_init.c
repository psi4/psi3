#include <stdio.h>
#include <stdlib.h>
#include <psio.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_file_init(): Prepares a dpd two-electron file on disk for
** reading/writing.
**
** Arguments:
**   struct dpdfile *File: A pointer to the dpdfile to be initialized.
**   int filenum: The PSI unit number for this file.
**   int pqnum: The index combination for the bra indices for the
**              data as it will be stored on disk.
**   int rsnum: The index combination for the ket indices for the
**              data as it will be stored on disk.
**   char *label: A string labelling for this buffer.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_file_init(struct dpdfile *File, int filenum, int pqnum,
		  int rsnum,  char *label, int print_flag,
		  FILE *outfile)
{
  int i;
  
  strcpy(File->label,label);
  File->filenum = filenum;
  File->params = &(dpd_params[pqnum][rsnum]);

  File->matrix = (double ***) malloc(File->params->nirreps*sizeof(double **));

  /* Construct logical subfile pointers */
  File->lfiles = (psio_address *) malloc(File->params->nirreps *
					 sizeof(psio_address));
  File->lfiles[0] = PSIO_ZERO;
  for(i=1; i < File->params->nirreps; i++)
    File->lfiles[i] = psio_get_address(File->lfiles[i-1],
				       (File->params->rowtot[i-1] *
					File->params->coltot[i-1] *
					sizeof(double)));

  if(print_flag) {
      fprintf(outfile, "\n\tDPD File %d: %s\n", filenum, label);
      dpd_params_print(File->params, outfile);
    }
  return 0;
}
		  
