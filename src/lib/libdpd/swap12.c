#include <stdio.h>
#include <libciomr.h>
#include <qt.h>
#include "dpd.h"

/* dpd_swap12(): Swaps indices 1 and 2 (p and q) in a dpd two-electron
** buffer and writes the data to a dpd two-electron file.
**
** Arguments:
**   struct dpdbuf *InBuf: A pointer to the already-initialzed input
**                         buffer.
**   int outfilenum: The PSI unit number for the target data.
**   int pqnum: The index combination for the bra indices for the new
**              dpd file.  N.B. this is NOT error checked for consistency
**              with the input buffer.
**   int rsnum: The index combination for the ket indices for the new
**              dpd file.  N.B. this is NOT error checked for consistency
**              with the input buffer.
**   char *label: A string labelling for this buffer.
**   int print_flag: A boolean for the print routines.
**   FILE *outfile: The formatted output file stream.
*/

int dpd_swap12(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile)
{
  int h,nirreps;
  int p, q, r, s, pq, rs, inbufqp, inbufrs;
  struct dpdfile OutFile;

  nirreps = InBuf->params->nirreps;

  dpd_file_init(&OutFile, outfilenum, pqnum, rsnum, label, print_flag, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_file_mat_irrep_init(&OutFile, h);
      dpd_buf_mat_irrep_init(InBuf, h);
      dpd_buf_mat_irrep_rd(InBuf, h, print_flag, outfile);

      for(pq=0; pq < OutFile.params->rowtot[h]; pq++) {
	  p = OutFile.params->roworb[h][pq][0];
	  q = OutFile.params->roworb[h][pq][1];
	  inbufqp = InBuf->params->rowidx[q][p];

	  for(rs=0; rs < OutFile.params->coltot[h]; rs++) {
	      r = OutFile.params->colorb[h][rs][0];
	      s = OutFile.params->colorb[h][rs][1];
	      inbufrs = InBuf->params->colidx[r][s];
	      
	      OutFile.matrix[h][pq][rs] = InBuf->matrix[h][inbufqp][inbufrs];
	    }
	}
      dpd_file_mat_irrep_wrt(&OutFile, h, print_flag, outfile);
      dpd_file_mat_irrep_close(&OutFile, h);
      dpd_buf_mat_irrep_close(InBuf, h);
    }
  dpd_file_close(&OutFile);
}
