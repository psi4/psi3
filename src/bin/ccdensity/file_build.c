#include <stdio.h>
#include <iwl.h>
#include "dpd.h"

void idx_permute(struct dpdfile *File, int p, int q, int r, int s,
		 int perm_pr, int perm_qs, int perm_prqs,
		 double value, FILE *outfile);

int file_build(struct dpdfile *File, int inputfile, double tolerance,
	       int perm_pr, int perm_qs, int perm_prqs,
	       int keep, int print_flag, FILE *outfile)
{
  struct iwlbuf InBuf;
  int lastbuf;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s;
  double value;

  /* This has got to change!!! */
  dpd_file_mat_init(File);

  iwl_buf_init(&InBuf, inputfile, tolerance, 1, 1);

  lblptr = InBuf.labels;
  valptr = InBuf.values;
  lastbuf = InBuf.lastbuf;

  for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      idx_permute(File,p,q,r,s,perm_pr,perm_qs,perm_prqs,value,outfile);

          } /* end loop through current buffer */

  /* Now run through the rest of the buffers in the file */
  while (!lastbuf) {
    iwl_buf_fetch(&InBuf);
    lastbuf = InBuf.lastbuf;

    for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      idx_permute(File,p,q,r,s,perm_pr,perm_qs,perm_prqs,value,outfile);
      
      } /* end loop through current buffer */
    } /* end loop over reading buffers */

  /* Dump the matrices out to disk */
  dpd_file_mat_wrt(File, print_flag, outfile);
  dpd_file_mat_close(File);

  iwl_buf_close(&InBuf, keep);

  return 0;
}
