#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void FaetT2(void)
{
  struct oe_dpdfile FAEt, Faet;
  struct dpdbuf newtIJAB, newtijab, newtIjAb;
  struct dpdbuf tIJAB, tijab, tIjAb;

  dpd_buf_init(&newtIJAB, CC_TAMPS, 2, 5, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&newtijab, CC_TAMPS, 2, 5, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&newtIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);

  dpd_buf_init(&tIJAB, CC_TAMPS, 2, 5, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&tijab, CC_TAMPS, 2, 5, 2, 7, 0, "tijab", 0, outfile);
  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);

  dpd_oe_file_init(&FAEt, CC_OEI, 1, 1, "FAEt", 0, outfile);
  dpd_oe_file_init(&Faet, CC_OEI, 1, 1, "Faet", 0, outfile);

  dpd_contract221(&tIJAB, &FAEt, &newtIJAB, 3, 1, 0, 1, 1, 0, outfile);
  dpd_contract212(&FAEt, &tIJAB, &newtIJAB, 1, 2, 1, 1, 1, 0, outfile);

  dpd_contract221(&tijab, &Faet, &newtijab, 3, 1, 0, 1, 1, 0, outfile);
  dpd_contract212(&Faet, &tijab, &newtijab, 1, 2, 1, 1, 1, 0, outfile);

  dpd_contract221(&tIjAb, &Faet, &newtIjAb, 3, 1, 0, 1, 1, 0, outfile);
  dpd_contract212(&FAEt, &tIjAb, &newtIjAb, 1, 2, 1, 1, 1, 0, outfile);

  dpd_oe_file_close(&FAEt);  dpd_oe_file_close(&Faet);

  dpd_buf_close(&tIJAB);  dpd_buf_close(&tijab);  dpd_buf_close(&tIjAb);
  dpd_buf_close(&newtIJAB); dpd_buf_close(&newtijab);  dpd_buf_close(&newtIjAb);

}
