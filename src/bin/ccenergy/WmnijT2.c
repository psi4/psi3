#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void WmnijT2(void)
{
  struct dpdbuf newtIJAB, newtijab, newtIjAb;
  struct dpdbuf WMNIJ, Wmnij, WMnIj;
  struct dpdbuf tauIJAB, tauijab, tauIjAb;

  dpd_buf_init(&newtIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&newtijab, CC_TAMPS, 2, 7, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&newtIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);

  dpd_buf_init(&WMNIJ, CC_HBAR, 2, 2, 2, 2, 0, "WMNIJ", 0, outfile);
  dpd_buf_init(&Wmnij, CC_HBAR, 2, 2, 2, 2, 0, "Wmnij", 0, outfile);
  dpd_buf_init(&WMnIj, CC_HBAR, 0, 0, 0, 0, 0, "WMnIj", 0, outfile);

  dpd_buf_init(&tauIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&tauijab, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&tauIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);

  dpd_contract222(&WMNIJ, &tauIJAB, &newtIJAB, 1, 1, 1, 1, 0, outfile);
  dpd_contract222(&Wmnij, &tauijab, &newtijab, 1, 1, 1, 1, 0, outfile);
  dpd_contract222(&WMnIj, &tauIjAb, &newtIjAb, 1, 1, 1, 1, 0, outfile);

  dpd_buf_close(&tauIJAB); dpd_buf_close(&tauijab); dpd_buf_close(&tauIjAb);

  dpd_buf_close(&WMNIJ); dpd_buf_close(&Wmnij); dpd_buf_close(&WMnIj);

  dpd_buf_close(&newtIJAB); dpd_buf_close(&newtijab); dpd_buf_close(&newtIjAb);

}
