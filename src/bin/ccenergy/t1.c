#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void t1_build(void)
{
  struct oe_dpdfile newtIA, newtia, tIA, tia, fIA, fia;
  struct oe_dpdfile FAE, Fae, FMI, Fmi, FME, Fme;
  struct oe_dpdfile dIA, dia;
  struct dpdbuf tIJAB, tijab, tIjAb, tiJaB;
  struct dpdbuf C_anti, D, F_anti, F, E_anti, E;

/*  timer_on("t1_build", outfile); */

  dpd_oe_file_init(&fIA, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_copy(&fIA, CC_OEI, "New tIA", 0, outfile);
  dpd_oe_file_close(&fIA);

  dpd_oe_file_init(&fia, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_copy(&fia, CC_OEI, "New tia", 0, outfile);
  dpd_oe_file_close(&fia);

  dpd_oe_file_init(&newtIA, CC_OEI, 0, 1, "New tIA", 0, outfile);
  dpd_oe_file_init(&newtia, CC_OEI, 0, 1, "New tia", 0, outfile);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_oe_file_init(&FAE, CC_OEI, 1, 1, "FAE", 0, outfile);
  dpd_oe_file_init(&Fae, CC_OEI, 1, 1, "Fae", 0, outfile);

  dpd_contract111(&tIA, &FAE, &newtIA, 0, 0, 1, 1, 0, outfile);
  dpd_contract111(&tia, &Fae, &newtia, 0, 0, 1, 1, 0, outfile);

  dpd_oe_file_close(&FAE);  dpd_oe_file_close(&Fae);

  dpd_oe_file_init(&FMI, CC_OEI, 0, 0, "FMI", 0, outfile);
  dpd_oe_file_init(&Fmi, CC_OEI, 0, 0, "Fmi", 0, outfile);

  dpd_contract111(&FMI, &tIA, &newtIA, 1, 1, -1, 1, 0, outfile);
  dpd_contract111(&Fmi, &tia, &newtia, 1, 1, -1, 1, 0, outfile);

  dpd_oe_file_close(&FMI);  dpd_oe_file_close(&Fmi);
  dpd_oe_file_close(&tIA);  dpd_oe_file_close(&tia);

  dpd_oe_file_init(&FME, CC_OEI, 0, 1, "FME", 0, outfile);
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "Fme", 0, outfile);

  dpd_buf_init(&tIJAB, CC_TAMPS, 0, 5, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&tijab, CC_TAMPS, 0, 5, 2, 7, 0, "tijab", 0, outfile);

  dpd_dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1, 0, outfile);
  dpd_dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1, 0, outfile);
  
  dpd_buf_close(&tIJAB);  dpd_buf_close(&tijab);

  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);

  dpd_dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1, 0, outfile);
  dpd_dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1, 0, outfile);
  
  dpd_buf_close(&tIjAb);
  
  dpd_oe_file_close(&FME);  dpd_oe_file_close(&Fme);

  dpd_buf_init(&C_anti, CC_CINTS, 10, 10, 10, 10, 0,
	       "C <ia||jb>", 0, outfile);
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1, 0, outfile);
  dpd_dot13(&tia, &D, &newtIA, 0, 0, 1, 1, 0, outfile);

  dpd_dot14(&tia, &C_anti, &newtia, 0, 1, -1, 1, 0, outfile);
  dpd_dot13(&tIA, &D, &newtia, 0, 0, 1, 1, 0, outfile);

  dpd_oe_file_close(&tIA);  dpd_oe_file_close(&tia);
  dpd_buf_close(&C_anti);  dpd_buf_close(&D);

  dpd_buf_init(&F_anti, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&tIJAB, CC_TAMPS, 0, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&tijab, CC_TAMPS, 0, 7, 2, 7, 0, "tijab", 0, outfile);

  dpd_contract122(&tIJAB, &F_anti, &newtIA, 1, 1, 1, 1, 0, outfile);
  dpd_contract122(&tijab, &F_anti, &newtia, 1, 1, 1, 1, 0, outfile);

  dpd_buf_close(&tIJAB);
  dpd_buf_close(&tijab);
  dpd_buf_close(&F_anti);

  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_buf_init(&tiJaB, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);

  dpd_contract122(&tiJaB, &F, &newtIA, 1, 1, 1, 1, 0, outfile);
  dpd_contract122(&tIjAb, &F, &newtia, 1, 1, 1, 1, 0, outfile);
  
  dpd_buf_close(&F);  dpd_buf_close(&tIjAb);  dpd_buf_close(&tiJaB);

  dpd_buf_init(&E_anti, CC_EINTS, 11, 2, 11, 0, 1, "E <ai|jk>", 0, outfile);
  dpd_buf_init(&tIJAB, CC_TAMPS, 2, 5, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&tijab, CC_TAMPS, 2, 5, 2, 7, 0, "tijab", 0, outfile);

  dpd_contract122(&E_anti, &tIJAB, &newtIA, 1, 3, -1, 1, 0, outfile);
  dpd_contract122(&E_anti, &tijab, &newtia, 1, 3, -1, 1, 0, outfile);

  dpd_buf_close(&E_anti);  dpd_buf_close(&tIJAB);  dpd_buf_close(&tijab);
  
  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);
  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_buf_init(&tiJaB, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);

  dpd_contract122(&E, &tiJaB, &newtIA, 1, 3, -1, 1, 0, outfile);
  dpd_contract122(&E, &tIjAb, &newtia, 1, 3, -1, 1, 0, outfile);

  dpd_buf_close(&E);  dpd_buf_close(&tIjAb);  dpd_buf_close(&tiJaB);

  dpd_oe_file_init(&dIA, CC_OEI, 0, 1, "dIA", 0, outfile);
  dpd_oe_dirprd(&dIA, &newtIA, 0, outfile);
  dpd_oe_file_close(&dIA);

  dpd_oe_file_init(&dia, CC_OEI, 0, 1, "dia", 0, outfile);
  dpd_oe_dirprd(&dia, &newtia, 0, outfile);
  dpd_oe_file_close(&dia);

  dpd_oe_file_close(&newtIA);
  dpd_oe_file_close(&newtia);

/*  timer_off("t1_build", outfile); */
}
