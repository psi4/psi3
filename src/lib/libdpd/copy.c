#include <stdio.h>
#include "dpd.h"

/* dpd_copy(): Copies an existing two-electron dpdbuf into another file.
**
** Arguments:
**   struct *dpdbuf InBuf: A pointer to the given dpd buffer.
**   int outfilenum: The PSI unit number for the new buffer.
**   char *label: A string labelling for this buffer.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_copy(struct dpdbuf *InBuf, int outfilenum, 
	     char *label, int print_flag, FILE *outfile)
{
  int h;
  double ***matrix;
  struct dpdfile OutFile;

  dpd_file_init(&OutFile, outfilenum, InBuf->params->pqnum,
		InBuf->params->rsnum, label, print_flag, outfile);

  /* Save the OutFile's matrix pointer */
  matrix = OutFile.matrix;

  /* Temporarily assign the OutFile's matrix pointer to InBuf */
  OutFile.matrix = InBuf->matrix;

  for(h=0; h < InBuf->params->nirreps; h++) {

      dpd_buf_mat_irrep_init(InBuf, h);
      dpd_buf_mat_irrep_rd(InBuf, h, print_flag, outfile);

      /* Since OutFile and InBuf point to the same memory, we can dump
	 the buffer's contents to OutFile */
      dpd_file_mat_irrep_wrt(&OutFile, h, print_flag, outfile);

      dpd_buf_mat_irrep_close(InBuf, h);
    }

  /* Clean up the OutFile's matrix pointer */
  OutFile.matrix = matrix;

  dpd_file_close(&OutFile);

  return 0;
}
