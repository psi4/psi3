#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void FaeL2(void)
{
  dpdbuf4 Lijab, LIJAB, LIjAb;
  dpdbuf4 newLIJAB, newLijab, newLIjAb;
  dpdfile2 LFaet2, LFAEt2;

 /* RHS += P(ab)*Lijae*Feb */

  dpd_file2_init(&LFAEt2, CC_OEI, 0, 1, 1, "FAEt");
  dpd_file2_init(&LFaet2, CC_OEI, 0, 1, 1, "Faet");

  dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&newLIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New LIJAB");
  dpd_contract424(&LIJAB, &LFAEt2, &newLIJAB, 3, 0, 0, 1.0, 1.0);
  dpd_contract244(&LFAEt2, &LIJAB, &newLIJAB, 0, 2, 1, 1.0, 1.0);
  dpd_buf4_close(&LIJAB);
  dpd_buf4_close(&newLIJAB);

  dpd_buf4_init(&Lijab, CC_LAMPS, 0, 2, 5, 2, 7, 0, "Lijab");
  dpd_buf4_init(&newLijab, CC_LAMPS, 0, 2, 5, 2, 7, 0, "New Lijab");
  dpd_contract424(&Lijab, &LFaet2, &newLijab, 3, 0, 0, 1.0, 1.0);
  dpd_contract244(&LFaet2, &Lijab, &newLijab, 0, 2, 1, 1.0, 1.0);
  dpd_buf4_close(&Lijab);
  dpd_buf4_close(&newLijab);

  dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_init(&newLIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "New LIjAb");
  dpd_contract424(&LIjAb, &LFaet2, &newLIjAb, 3, 0, 0, 1.0, 1.0);
  dpd_contract244(&LFAEt2, &LIjAb, &newLIjAb, 0, 2, 1, 1.0, 1.0);
  dpd_buf4_close(&LIjAb);
  dpd_buf4_close(&newLIjAb);

  dpd_file2_close(&LFaet2);
  dpd_file2_close(&LFAEt2);

}
