#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_buf_init(): Initializes a dpd two-electron buffer for reading or writing
** symmetry-packed data.
**
** Arguments:
**   struct dpdbuf *Buf: A pointer to the dpdbuf to be initialized.
**   int inputfile: The PSI unit number for the data on disk.
**   int pqnum: The index combination for the bra indices for the
**              data as it will be used in memory.
**   int rsnum: The index combination for the ket indices for the
**              data as it will be used in memory.
**   int file_pqnum: The index combination for the bra indices for the
**                   data as it will be stored on disk.
**   int file_rsnum: The index combination for the ket indices for the
**                   data as it will be stored on disk.
**   int anti: Boolean flag which indicates whether the data needs to
**             be antisymmetrized as it is read from disk.
**   char *label: A string labelling for this buffer.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_buf_init(struct dpdbuf *Buf, int inputfile, int pqnum, int rsnum,
		 int file_pqnum, int file_rsnum, int anti, char *label,
		 int print_flag, FILE *outfile)
{
  int nirreps;
  
  Buf->anti = anti;
  Buf->params = &(dpd_params[pqnum][rsnum]);

  dpd_file_init(&(Buf->file), inputfile, file_pqnum, file_rsnum,
		label, print_flag, outfile);

  if(print_flag) {
      fprintf(outfile, "\n\tDPD Buffer for file: %s\n", label);
      dpd_params_print(Buf->params, outfile);
    }

  Buf->matrix = (double ***) malloc(Buf->params->nirreps*sizeof(double **));

  /* Set up shifted matrix info */
  nirreps = Buf->params->nirreps;
  Buf->shift.shift_type = 0;
  Buf->shift.rowtot = init_int_matrix(nirreps, nirreps);
  Buf->shift.coltot = init_int_matrix(nirreps, nirreps);
  Buf->shift.matrix = (double ****) malloc(nirreps * sizeof(double ***));

  return 0;
}
