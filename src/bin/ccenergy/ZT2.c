#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void ZT2(void)
{
  struct dpdbuf ZIJMA, ZIJAM, Zijma, Zijam, ZIjMa, ZIjAm;
  struct dpdbuf newtIJAB, newtijab, newtIjAb;
  struct oe_dpdfile tIA, tia;

  dpd_buf_init(&ZIJMA, CC_MISC, 2, 10, 2, 10, 0, "ZIJMA", 0, outfile);
  dpd_buf_init(&ZIJAM, CC_MISC, 2, 11, 2, 11, 0, "ZIJAM", 0, outfile);
  dpd_buf_init(&Zijma, CC_MISC, 2, 10, 2, 10, 0, "Zijma", 0, outfile);
  dpd_buf_init(&Zijam, CC_MISC, 2, 11, 2, 11, 0, "Zijam", 0, outfile);
  dpd_buf_init(&ZIjMa, CC_MISC, 0, 10, 0, 10, 0,  "ZIjMa", 0, outfile);
  dpd_buf_init(&ZIjAm, CC_MISC, 0, 11, 0, 11, 0,  "ZIjAm", 0, outfile);

  dpd_buf_init(&newtIJAB, CC_TAMPS, 2, 5, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&newtijab, CC_TAMPS, 2, 5, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&newtIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_contract221(&ZIJAM, &tIA, &newtIJAB, 3, 0, 0, 1, 1, 0, outfile);
  dpd_contract212(&tIA, &ZIJMA, &newtIJAB, 0, 2, 1, -1, 1, 0, outfile);

  dpd_contract221(&Zijam, &tia, &newtijab, 3, 0, 0, 1, 1, 0, outfile);
  dpd_contract212(&tia, &Zijma, &newtijab, 0, 2, 1, -1, 1, 0, outfile);

  dpd_contract221(&ZIjAm, &tia, &newtIjAb, 3, 0, 0, -1, 1, 0, outfile);
  dpd_contract212(&tIA, &ZIjMa, &newtIjAb, 0, 2, 1, -1, 1, 0, outfile);

  dpd_oe_file_close(&tIA); dpd_oe_file_close(&tia);

  dpd_buf_close(&newtIJAB); dpd_buf_close(&newtijab); dpd_buf_close(&newtIjAb);

  dpd_buf_close(&ZIJMA); dpd_buf_close(&ZIJAM); dpd_buf_close(&Zijma);
  dpd_buf_close(&Zijam); dpd_buf_close(&ZIjMa); dpd_buf_close(&ZIjAm);
}
