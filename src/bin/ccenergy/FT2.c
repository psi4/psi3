#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void FT2(void)
{
  struct oe_dpdfile tIA, tia;
  struct dpdbuf newtIJAB, newtijab, newtIjAb, t2;
  struct dpdbuf F_anti, F;

/*  timer_on("FT2", outfile); */

  dpd_buf_init(&newtIJAB, CC_TAMPS, 0, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&newtijab, CC_TAMPS, 0, 7, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&newtIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);

  dpd_buf_init(&F_anti, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_contract221(&F_anti, &tIA, &newtIJAB, 1, 1, 1, 1, 1, 0, outfile);
  dpd_contract221(&F_anti, &tia, &newtijab, 1, 1, 1, 1, 1, 0, outfile);
  dpd_contract221(&F, &tia, &newtIjAb, 1, 1, 1, 1, 1, 0, outfile);

  dpd_buf_init(&t2, CC_TMP0, 0, 7, 0, 7, 0, "tIJAB", 0, outfile);
  dpd_contract221(&F_anti, &tIA, &t2, 1, 1, 1, -1, 0, 0, outfile);
  dpd_swap12(&t2, CC_TMP1, 0, 7, "tJIAB", 0, outfile);
  dpd_buf_close(&t2);
  dpd_buf_init(&t2, CC_TMP1, 0, 7, 0, 7, 0, "tJIAB", 0, outfile);
  dpd_axpy(&t2, &newtIJAB, 1, 0, outfile);
  dpd_buf_close(&t2);


  dpd_buf_init(&t2, CC_TMP0, 0, 7, 0, 7, 0, "tijab", 0, outfile);
  dpd_contract221(&F_anti, &tia, &t2, 1, 1, 1, -1, 0, 0, outfile);
  dpd_swap12(&t2, CC_TMP1, 0, 7, "tjiab", 0, outfile);
  dpd_buf_close(&t2);
  dpd_buf_init(&t2, CC_TMP1, 0, 7, 0, 7, 0, "tjiab", 0, outfile);
  dpd_axpy(&t2, &newtijab, 1, 0, outfile);
  dpd_buf_close(&t2);
  

  dpd_buf_init(&t2, CC_TMP0, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_contract221(&F, &tIA, &t2, 1, 1, 1, 1, 0, 0, outfile);
  dpd_swap12(&t2, CC_TMP1, 0, 5, "tJiaB", 0, outfile);
  dpd_buf_close(&t2);
  dpd_buf_init(&t2, CC_TMP1, 0, 5, 0, 5, 0, "tJiaB", 0, outfile);
  dpd_swap34(&t2, CC_TMP0, 0, 5, "tJiBa", 0, outfile);
  dpd_buf_close(&t2);
  dpd_buf_init(&t2, CC_TMP0, 0, 5, 0, 5, 0, "tJiBa", 0, outfile);
  dpd_axpy(&t2, &newtIjAb, 1, 0, outfile);
  dpd_buf_close(&t2);

  dpd_oe_file_close(&tIA); dpd_oe_file_close(&tia);

  dpd_buf_close(&F_anti); dpd_buf_close(&F);
  
  dpd_buf_close(&newtIJAB); dpd_buf_close(&newtijab); dpd_buf_close(&newtIjAb);

/*  timer_off("FT2", outfile); */
}
