#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void FaeL2(void)
{
  struct dpdbuf Lijab, LIJAB, LIjAb;
  struct dpdbuf newLIJAB, newLijab, newLIjAb;
  struct oe_dpdfile LFaet2, LFAEt2;

 /* RHS += P(ab)*Lijae*Feb */

  dpd_oe_file_init(&LFAEt2, CC_OEI, 1, 1, "FAEt", 0, outfile);
  dpd_oe_file_init(&LFaet2, CC_OEI, 1, 1, "Faet", 0, outfile);

  dpd_buf_init(&LIJAB, CC_LAMPS, 2, 5, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_buf_init(&newLIJAB, CC_LAMPS, 2, 5, 2, 7, 0, "New LIJAB",
               0, outfile);
  dpd_contract221(&LIJAB, &LFAEt2, &newLIJAB, 3, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_contract212(&LFAEt2, &LIJAB, &newLIJAB, 0, 2, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&LIJAB);
  dpd_buf_close(&newLIJAB);

  dpd_buf_init(&Lijab, CC_LAMPS, 2, 5, 2, 7, 0, "Lijab", 0, outfile);
  dpd_buf_init(&newLijab, CC_LAMPS, 2, 5, 2, 7, 0, "New Lijab",
              0, outfile);
  dpd_contract221(&Lijab, &LFaet2, &newLijab, 3, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_contract212(&LFaet2, &Lijab, &newLijab, 0, 2, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Lijab);
  dpd_buf_close(&newLijab);

  dpd_buf_init(&LIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_buf_init(&newLIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb",
              0, outfile);
  dpd_contract221(&LIjAb, &LFaet2, &newLIjAb, 3, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_contract212(&LFAEt2, &LIjAb, &newLIjAb, 0, 2, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&LIjAb);
  dpd_buf_close(&newLIjAb);

  dpd_oe_file_close(&LFaet2);
  dpd_oe_file_close(&LFAEt2);

}
