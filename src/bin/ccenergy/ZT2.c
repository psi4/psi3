#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void ZT2(void)
{
  dpdbuf4 ZIJMA, ZIJAM, Zijma, Zijam, ZIjMa, ZIjAm;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdfile2 tIA, tia;
  dpdbuf4 t2;

  dpd_buf4_init(&ZIJMA, CC_MISC, 0, 2, 10, 2, 10, 0, "ZIJMA");
  dpd_buf4_init(&ZIJAM, CC_MISC, 0, 2, 11, 2, 11, 0, "ZIJAM");
  dpd_buf4_init(&Zijma, CC_MISC, 0, 2, 10, 2, 10, 0, "Zijma");
  dpd_buf4_init(&Zijam, CC_MISC, 0, 2, 11, 2, 11, 0, "Zijam");
  dpd_buf4_init(&ZIjMa, CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjMa");
  dpd_buf4_init(&ZIjAm, CC_MISC, 0, 0, 11, 0, 11, 0, "ZIjAm");

  dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
  dpd_buf4_init(&newtijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
  dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
  dpd_contract424(&ZIJAM, &tIA, &t2, 3, 0, 0, 1, 0);
  dpd_contract244(&tIA, &ZIJMA, &t2, 0, 2, 1, -1, 1);
  dpd_buf4_axpy(&t2, &newtIJAB, 1);
  dpd_buf4_close(&t2);

  dpd_buf4_init(&t2, CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
  dpd_contract424(&Zijam, &tia, &t2, 3, 0, 0, 1, 0);
  dpd_contract244(&tia, &Zijma, &t2, 0, 2, 1, -1, 1);
  dpd_buf4_axpy(&t2, &newtijab, 1);
  dpd_buf4_close(&t2);

  dpd_contract424(&ZIjAm, &tia, &newtIjAb, 3, 0, 0, -1, 1);
  dpd_contract244(&tIA, &ZIjMa, &newtIjAb, 0, 2, 1, -1, 1);

  dpd_file2_close(&tIA); 
  dpd_file2_close(&tia);

  dpd_buf4_close(&newtIJAB); 
  dpd_buf4_close(&newtijab); 
  dpd_buf4_close(&newtIjAb); 

  dpd_buf4_close(&ZIJMA); 
  dpd_buf4_close(&ZIJAM); 
  dpd_buf4_close(&Zijma);
  dpd_buf4_close(&Zijam); 
  dpd_buf4_close(&ZIjMa); 
  dpd_buf4_close(&ZIjAm);
}
