#include <stdio.h>
#include <qt.h>
#include "dpd.h"

/* dpd_dirprd(): Computes the direct product between two dpd two-electron
** buffers.
**
** Arguments:
**   struct dpdbuf *BufA: A pointer to one of the dpd two-electron buffers.
**   struct dpdbuf *BufB: A pointer to the other dpd two-electron buffer.
**   int print_flag: A value for the print routines.
**   FILE *outfile: A pointer to the output file stream.
*/

int dpd_dirprd(struct dpdbuf *BufA, struct dpdbuf *BufB, 
	       int print_flag, FILE *outfile)
{
  int h, nirreps;
  int row, col;

  nirreps = BufA->params->nirreps;

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(BufA, h);
      dpd_buf_mat_irrep_init(BufB, h);
      dpd_buf_mat_irrep_rd(BufA, h, print_flag, outfile);
      dpd_buf_mat_irrep_rd(BufB, h, print_flag, outfile);

      dirprd_block(BufA->matrix[h], BufB->matrix[h],
		   BufA->params->rowtot[h], BufA->params->coltot[h]);

      dpd_buf_mat_irrep_wrt(BufB, h, print_flag, outfile);
      dpd_buf_mat_irrep_close(BufA, h);
      dpd_buf_mat_irrep_close(BufB, h);
    }

  return 0;
}
      
