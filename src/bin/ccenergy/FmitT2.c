#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void FmitT2(void)
{
  struct oe_dpdfile FMIt, Fmit;
  struct dpdbuf newtIJAB, newtijab, newtIjAb;
  struct dpdbuf tIJAB, tijab, tIjAb;

  dpd_buf_init(&newtIJAB, CC_TAMPS, 0, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&newtijab, CC_TAMPS, 0, 7, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&newtIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);

  dpd_buf_init(&tIJAB, CC_TAMPS, 0, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&tijab, CC_TAMPS, 0, 7, 2, 7, 0, "tijab", 0, outfile);
  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);

  dpd_oe_file_init(&FMIt, CC_OEI, 0, 0, "FMIt", 0, outfile);
  dpd_oe_file_init(&Fmit, CC_OEI, 0, 0, "Fmit", 0, outfile);

  dpd_contract221(&tIJAB, &FMIt, &newtIJAB, 1, 0, 1, -1, 1, 0, outfile);
  dpd_contract212(&FMIt, &tIJAB, &newtIJAB, 0, 0, 0, -1, 1, 0, outfile);

  dpd_contract221(&tijab, &Fmit, &newtijab, 1, 0, 1, -1, 1, 0, outfile);
  dpd_contract212(&Fmit, &tijab, &newtijab, 0, 0, 0, -1, 1, 0, outfile);

  dpd_contract221(&tIjAb, &Fmit, &newtIjAb, 1, 0, 1, -1, 1, 0, outfile);
  dpd_contract212(&FMIt, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1, 0, outfile);

  dpd_oe_file_close(&FMIt); dpd_oe_file_close(&Fmit);

  dpd_buf_close(&tIJAB); dpd_buf_close(&tijab); dpd_buf_close(&tIjAb);

  dpd_buf_close(&newtIJAB); dpd_buf_close(&newtijab); dpd_buf_close(&newtIjAb);
}
