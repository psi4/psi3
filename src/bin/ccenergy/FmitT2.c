#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void FmitT2(void)
{
  dpdfile2 FMIt, Fmit;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 tIJAB, tijab, tIjAb;
  dpdbuf4 t2;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");

    dpd_contract424(&tIjAb, &FMIt, &newtIjAb, 1, 0, 1, -1, 1);
    dpd_contract244(&FMIt, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);

    dpd_file2_close(&FMIt);

    dpd_buf4_close(&tIjAb);

    dpd_buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_init(&Fmit, CC_OEI, 0, 0, 0, "Fmit");

    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract424(&tIJAB, &FMIt, &t2, 1, 0, 1, -1, 0);
    dpd_contract244(&FMIt, &tIJAB, &t2, 0, 0, 0, -1, 1);
    dpd_buf4_axpy(&t2, &newtIJAB, 1);
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract424(&tijab, &Fmit, &t2, 1, 0, 1, -1, 0);
    dpd_contract244(&Fmit, &tijab, &t2, 0, 0, 0, -1, 1);
    dpd_buf4_axpy(&t2, &newtijab, 1);
    dpd_buf4_close(&t2);

    dpd_contract424(&tIjAb, &Fmit, &newtIjAb, 1, 0, 1, -1, 1);
    dpd_contract244(&FMIt, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);

    dpd_file2_close(&FMIt); dpd_file2_close(&Fmit);

    dpd_buf4_close(&tIJAB);
    dpd_buf4_close(&tijab);
    dpd_buf4_close(&tIjAb);

    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&newtIjAb);
  }
}
