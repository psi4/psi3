#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void ET2(void)
{
  struct oe_dpdfile tIA, tia;
  struct dpdbuf newtIJAB, newtijab, newtIjAb;
  struct dpdbuf E_anti, E, t2;

/*  timer_on("ET2", outfile); */

  dpd_buf_init(&newtIJAB, CC_TAMPS, 2, 5, 2, 7, 0, "New tIJAB", 0, outfile);
  dpd_buf_init(&newtijab, CC_TAMPS, 2, 5, 2, 7, 0, "New tijab", 0, outfile);
  dpd_buf_init(&newtIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "New tIjAb", 0, outfile);

  dpd_buf_init(&E_anti, CC_EINTS, 11, 2, 11, 0, 1, "E <ai|jk>", 0, outfile);
  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);
  
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_contract221(&E_anti, &tIA, &newtIJAB, 1, 0, 0, -1, 1, 0, outfile);
  dpd_contract221(&E_anti, &tia, &newtijab, 1, 0, 0, -1, 1, 0, outfile);
  dpd_contract221(&E, &tia, &newtIjAb, 1, 0, 0, -1, 1, 0, outfile);

  dpd_buf_init(&t2, CC_TMP0, 2, 5, 2, 5, 0, "Temp T (IJ,AB)", 0, outfile);
  dpd_contract221(&E_anti, &tIA, &t2, 1, 0, 0, 1, 0, 0, outfile);
  dpd_swap34(&t2, CC_TMP1, 2, 5, "Temp T (IJ,BA)", 0, outfile);
  dpd_buf_close(&t2);
  dpd_buf_init(&t2, CC_TMP1, 2, 5, 2, 5, 0, "Temp T (IJ,BA)", 0, outfile);
  dpd_axpy(&t2, &newtIJAB, 1, 0, outfile);
  dpd_buf_close(&t2);


  dpd_buf_init(&t2, CC_TMP0, 2, 5, 2, 5, 0, "Temp T (ij,ab)", 0, outfile);
  dpd_contract221(&E_anti, &tia, &t2, 1, 0, 0, 1, 0, 0, outfile);
  dpd_swap34(&t2, CC_TMP1, 2, 5, "Temp T (ij,ba)", 0, outfile);
  dpd_buf_close(&t2);
  dpd_buf_init(&t2, CC_TMP1, 2, 5, 2, 5, 0, "Temp T (ij,ba)", 0, outfile);
  dpd_axpy(&t2, &newtijab, 1, 0, outfile);
  dpd_buf_close(&t2);


  dpd_buf_init(&t2, CC_TMP0, 0, 5, 0, 5, 0, "Temp T (iJ,aB)", 0, outfile);
  dpd_contract221(&E, &tIA, &t2, 1, 0, 0, -1, 0, 0, outfile);
  dpd_swap34(&t2, CC_TMP1, 0, 5, "Temp T (iJ,Ba)", 0, outfile);
  dpd_buf_close(&t2);
  dpd_buf_init(&t2, CC_TMP1, 0, 5, 0, 5, 0, "Temp T (iJ,Ba)", 0, outfile);
  dpd_swap12(&t2, CC_TMP0, 0, 5, "Temp T (Ji,Ba)", 0, outfile);
  dpd_buf_close(&t2);
  dpd_buf_init(&t2, CC_TMP0, 0, 5, 0, 5, 0, "Temp T (Ji,Ba)", 0, outfile);
  dpd_axpy(&t2, &newtIjAb, 1, 0, outfile);
  dpd_buf_close(&t2);

  dpd_oe_file_close(&tIA); dpd_oe_file_close(&tia);

  dpd_buf_close(&E_anti); dpd_buf_close(&E);
  
  dpd_buf_close(&newtIJAB); dpd_buf_close(&newtijab); dpd_buf_close(&newtIjAb);

/*  timer_off("ET2", outfile); */
}
