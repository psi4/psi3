#include <stdio.h>
#include <stdlib.h>
#include <psio.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_oe_file_init(): Initializes a dpd one-electron file for reading
** or writing data.
**
** Arguments:
**   struct oe_dpdfile *File: A pointer to the one-electron dpdfile.
**   int filenum: The PSI unit number for this file.
**   int pnum: The orbital subspace number for the left index [see
**             dpd_init()].
**   int qnum: The orbital subspace number for the right index [see
**             dpd_init()].
**   char *label: A string labelling for this buffer.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_oe_file_init(struct oe_dpdfile *File, int filenum, int pnum, int qnum,
		     char *label, int print_flag, FILE *outfile)
{
  int i;

  strcpy(File->label,label);
  File->filenum = filenum;
  File->params = &(dpd_oe_params[pnum][qnum]);

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
      fprintf(outfile, "\n\tDPD OE File: %s\n", label);
      dpd_oe_params_print(File->params, outfile);
    }

  return 0;
}
