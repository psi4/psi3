#include <stdio.h>
#include <libciomr.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void init_amps(void)
{
  struct oe_dpdfile tIA, tia, fIA, fia, dIA, dia;
  struct dpdbuf tIJAB, tijab, tIjAb, D, dIJAB, dijab, dIjAb;

  /* Restart from previous amplitudes if we can/should */
  /*  Still need to shift this to new I/O
  if(params.restart && psio_flen(CC_tIA) && psio_flen(CC_tia)
     && psio_flen(CC_tIJAB) && psio_flen(CC_tijab) && psio_flen(CC_tIjAb))
      return;
      */

  dpd_oe_file_init(&fIA, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_copy(&fIA, CC_OEI, "tIA", 0, outfile);
  dpd_oe_file_close(&fIA);
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
/*  dpd_oe_scm(&tIA, 0, 0, outfile);  */
  dpd_oe_file_close(&tIA);
  
  dpd_oe_file_init(&fia, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_copy(&fia, CC_OEI, "tia", 0, outfile);
  dpd_oe_file_close(&fia);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
/*  dpd_oe_scm(&tia, 0, 0, outfile); */
  dpd_oe_file_close(&tia);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&dIA, CC_OEI, 0, 1, "dIA", 0, outfile);
  dpd_oe_dirprd(&dIA, &tIA, 0, outfile);
  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&dIA);

  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_init(&dia, CC_OEI, 0, 1, "dia", 0, outfile);
  dpd_oe_dirprd(&dia, &tia, 0, outfile);
  dpd_oe_file_close(&tia);
  dpd_oe_file_close(&dia);

  dpd_buf_init(&D, CC_DINTS, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)",
	       0, outfile);
  dpd_copy(&D, CC_TAMPS, "tIJAB", 0, outfile);
  dpd_copy(&D, CC_TAMPS, "tijab", 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_init(&dIJAB, CC_DENOM, 1, 6, 1, 6, 0, "dIJAB", 0, outfile);
  dpd_buf_init(&tIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_dirprd(&dIJAB, &tIJAB, 0, outfile);
  dpd_buf_close(&tIJAB);
  dpd_buf_close(&dIJAB);

  dpd_buf_init(&dijab, CC_DENOM, 1, 6, 1, 6, 0, "dijab", 0, outfile);
  dpd_buf_init(&tijab, CC_TAMPS, 2, 7, 2, 7, 0, "tijab", 0, outfile);
  dpd_dirprd(&dijab, &tijab, 0, outfile);
  dpd_buf_close(&tijab);
  dpd_buf_close(&dijab);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_copy(&D, CC_TAMPS, "tIjAb", 0, outfile);
  dpd_buf_close(&D);
  
  dpd_buf_init(&dIjAb, CC_DENOM, 0, 5, 0, 5, 0, "dIjAb", 0, outfile);
  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_dirprd(&dIjAb, &tIjAb, 0, outfile);
  dpd_buf_close(&tIjAb);
  dpd_buf_close(&dIjAb);
}
