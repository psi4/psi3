#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void Z_build(void)
{
  struct dpdbuf ZIJMA, Zijma, ZIjMa, ZIjmA;
  struct dpdbuf tauIJAB, tauijab, tauIjAb, tauIjbA, F_anti, F;

/*  timer_on("Z_build", outfile); */

  dpd_buf_init(&ZIJMA, CC_MISC, 2, 10, 2, 10, 0, "ZIJMA", 0, outfile);
  dpd_buf_init(&Zijma, CC_MISC, 2, 10, 2, 10, 0, "Zijma", 0, outfile);
  dpd_buf_init(&ZIjMa, CC_MISC, 0, 10, 0, 10, 0, "ZIjMa", 0, outfile);
  dpd_buf_init(&ZIjmA, CC_MISC, 0, 10, 0, 10, 0, "ZIjmA", 0, outfile);

  dpd_buf_init(&tauIJAB, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&tauijab, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&tauIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_buf_init(&tauIjbA, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjbA", 0, outfile);

  dpd_buf_init(&F_anti, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);

  dpd_contract222(&tauIJAB, &F_anti, &ZIJMA, 0, 0, 1, 0, 0, outfile);
  dpd_contract222(&tauijab, &F_anti, &Zijma, 0, 0, 1, 0, 0, outfile);
  dpd_contract222(&tauIjAb, &F, &ZIjMa, 0, 0, 1, 0, 0, outfile);
  dpd_contract222(&tauIjbA, &F, &ZIjmA, 0, 0, 1, 0, 0, outfile);

  dpd_buf_close(&tauIJAB);  dpd_buf_close(&tauijab);  dpd_buf_close(&tauIjAb);
  dpd_buf_close(&F_anti);  dpd_buf_close(&F);

  dpd_swap34(&ZIJMA, CC_MISC, 2, 11, "ZIJAM", 0, outfile);
  dpd_swap34(&Zijma, CC_MISC, 2, 11, "Zijam", 0, outfile);
  dpd_swap34(&ZIjmA, CC_MISC, 0, 11, "ZIjAm", 0, outfile);

  dpd_buf_close(&ZIJMA);  dpd_buf_close(&Zijma);
  dpd_buf_close(&ZIjMa);  dpd_buf_close(&ZIjmA);

/*  timer_off("Z_build", outfile); */
}

