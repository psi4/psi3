#include <stdio.h>
#include <libciomr.h>
#include <qt.h>
#include "dpd.h"

int dpd_swap34(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile)
{
  int h,nirreps;
  int p, q, r, s, pq, rs, inbufpq, inbufsr;
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
	  inbufpq = InBuf->params->rowidx[p][q];

	  for(rs=0; rs < OutFile.params->coltot[h]; rs++) {
	      r = OutFile.params->colorb[h][rs][0];
	      s = OutFile.params->colorb[h][rs][1];
	      inbufsr = InBuf->params->colidx[s][r];

	      OutFile.matrix[h][pq][rs] = InBuf->matrix[h][inbufpq][inbufsr];
	    }
	}
      dpd_file_mat_irrep_wrt(&OutFile, h, print_flag, outfile);
      dpd_file_mat_irrep_close(&OutFile, h);
      dpd_buf_mat_irrep_close(InBuf, h);
    }
  dpd_file_close(&OutFile);
}
