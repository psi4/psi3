#include <stdio.h>
#include <libciomr.h>
#include <qt.h>
#include "dpd.h"

/* dpd_swap24(): Swaps indices 2 and 4 (q and s) in a dpd two-electron
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

int dpd_swap24(struct dpdbuf *InBuf, int outfilenum, int pqnum, int rsnum,
	       char *label, int print_flag, FILE *outfile)
{
  int h,nirreps;
  int Gp, Gq, Gr, Gs, Gpq, Grs, Gps, Grq;
  int p, q, r, s, P, Q, R, S, pq, rs, ps, rq;
  struct dpdfile OutFile;

  nirreps = InBuf->params->nirreps;

  dpd_file_init(&OutFile, outfilenum, pqnum, rsnum, label, print_flag, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_file_mat_irrep_init(&OutFile, h);

      for(Gp=0; Gp < nirreps; Gp++) {
	  Gq = Gp^h;
	  for(Gr=0; Gr < nirreps; Gr++) {
	      Gs = Gr^h;

	      Grq = Gps = Gp^Gs;

	      dpd_buf_mat_irrep_init(InBuf, Gps);
	      dpd_buf_mat_irrep_rd(InBuf, Gps, print_flag, outfile);

	      for(p=0; p < OutFile.params->ppi[Gp]; p++) {
		  P = OutFile.params->poff[Gp] + p;
		  for(q=0; q < OutFile.params->qpi[Gq]; q++) {
		      Q = OutFile.params->qoff[Gq] + q;
		      pq = OutFile.params->rowidx[P][Q];

		      for(r=0; r < OutFile.params->rpi[Gr]; r++) {
			  R = OutFile.params->roff[Gr] + r;
			  rq = InBuf->params->colidx[R][Q];
			  
			  for(s=0; s < OutFile.params->spi[Gs]; s++) {
			      S = OutFile.params->soff[Gs] + s;
			      rs = OutFile.params->colidx[R][S];
			      ps = InBuf->params->rowidx[P][S];

		      OutFile.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][rq];
			      
			    }
			}
		    }
		}

	      dpd_buf_mat_irrep_close(InBuf, Gps);
	    }

	}

      dpd_file_mat_irrep_wrt(&OutFile, h, print_flag, outfile);
      dpd_file_mat_irrep_close(&OutFile, h);

    }

  dpd_file_close(&OutFile);

}
