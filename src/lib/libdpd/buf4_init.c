#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

/* dpd_buf4_init(): Initializes a dpd four-index buffer for reading or writing
**   data.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the dpdbuf to be initialized.
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
**   char *label: The string labelling the PSIO TOC entry on disk.
*/

int dpd_buf4_init(dpdbuf4 *Buf, int inputfile, int irrep, int pqnum, int rsnum,
		 int file_pqnum, int file_rsnum, int anti, char *label)
{
  int nirreps;

  Buf->anti = anti;
  Buf->params = &(dpd_default->params4[pqnum][rsnum]);

  dpd_file4_init(&(Buf->file), inputfile, irrep, file_pqnum, file_rsnum, label);

  Buf->matrix = (double ***) malloc(Buf->params->nirreps*sizeof(double **));

  /* Set up shifted matrix info */
  nirreps = Buf->params->nirreps;
  Buf->shift.shift_type = 0;
  Buf->shift.rowtot = init_int_matrix(nirreps, nirreps);
  Buf->shift.coltot = init_int_matrix(nirreps, nirreps);
  Buf->shift.matrix = (double ****) malloc(nirreps * sizeof(double ***));

  return 0;
}
