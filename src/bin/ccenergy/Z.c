#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void Z_build(void)
{
  dpdbuf4 ZIJMA, Zijma, ZIjMa, ZIjmA;
  dpdbuf4 tauIJAB, tauijab, tauIjAb, tauIjbA, F_anti, F;

  timer_on("Z");

  dpd_buf4_init(&ZIJMA, CC_MISC, 0, 2, 10, 2, 10, 0, "ZIJMA");
  dpd_buf4_init(&Zijma, CC_MISC, 0, 2, 10, 2, 10, 0, "Zijma");
  dpd_buf4_init(&ZIjMa, CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjMa");
  dpd_buf4_init(&ZIjmA, CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjmA");

  dpd_buf4_init(&tauIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_buf4_init(&tauijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  dpd_buf4_init(&tauIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_buf4_init(&tauIjbA, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjbA");

  dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 7, 10, 7, 0, "F <ia||bc> (ia,b>c)");
  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");

  dpd_contract444(&tauIJAB, &F_anti, &ZIJMA, 0, 0, 1, 0);
  dpd_contract444(&tauijab, &F_anti, &Zijma, 0, 0, 1, 0);
  dpd_contract444(&tauIjAb, &F, &ZIjMa, 0, 0, 1, 0);
  dpd_contract444(&tauIjbA, &F, &ZIjmA, 0, 0, 1, 0);

  dpd_buf4_close(&tauIJAB); 
  dpd_buf4_close(&tauijab); 
  dpd_buf4_close(&tauIjAb);
  dpd_buf4_close(&tauIjbA);

  dpd_buf4_close(&F_anti); 
  dpd_buf4_close(&F);

  dpd_buf4_sort(&ZIJMA, CC_MISC, pqsr, 2, 11, "ZIJAM");
  dpd_buf4_sort(&Zijma, CC_MISC, pqsr, 2, 11, "Zijam");
  dpd_buf4_sort(&ZIjmA, CC_MISC, pqsr, 0, 11, "ZIjAm");

  dpd_buf4_close(&ZIJMA);  
  dpd_buf4_close(&Zijma);
  dpd_buf4_close(&ZIjMa);  
  dpd_buf4_close(&ZIjmA);

  timer_off("Z");
}

