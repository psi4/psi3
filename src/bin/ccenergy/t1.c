#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void t1_build(void)
{
  dpdfile2 newtIA, newtia, tIA, tia, fIA, fia;
  dpdfile2 FAE, Fae, FMI, Fmi, FME, Fme;
  dpdfile2 dIA, dia;
  dpdbuf4 tIJAB, tijab, tIjAb, tiJaB;
  dpdbuf4 C_anti, D, F_anti, F, E_anti, E;

/*  timer_on("t1_build", outfile); */

  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_copy(&fIA, CC_OEI, "New tIA");
  dpd_file2_close(&fIA);

  dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
  dpd_file2_copy(&fia, CC_OEI, "New tia");
  dpd_file2_close(&fia);

  dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_init(&newtia, CC_OEI, 0, 0, 1, "New tia");

  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
  dpd_file2_init(&Fae, CC_OEI, 0, 1, 1, "Fae");

  dpd_contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
  dpd_contract222(&tia, &Fae, &newtia, 0, 0, 1, 1);

  dpd_file2_close(&FAE);  dpd_file2_close(&Fae);

  dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
  dpd_file2_init(&Fmi, CC_OEI, 0, 0, 0, "Fmi");

  dpd_contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
  dpd_contract222(&Fmi, &tia, &newtia, 1, 1, -1, 1);

  dpd_file2_close(&FMI);  dpd_file2_close(&Fmi);
  dpd_file2_close(&tIA);  dpd_file2_close(&tia);

  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");

  dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");

  dpd_dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1);
  dpd_dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1);
  
  dpd_buf4_close(&tIJAB);  
  dpd_buf4_close(&tijab);

  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

  dpd_dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1);
  dpd_dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1);
  
  dpd_buf4_close(&tIjAb);
  
  dpd_file2_close(&FME);  
  dpd_file2_close(&Fme);

  dpd_buf4_init(&C_anti, CC_CINTS, 0, 10, 10, 10, 10, 0,
	       "C <ia||jb>");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  dpd_dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1);
  dpd_dot13(&tia, &D, &newtIA, 0, 0, 1, 1);

  dpd_dot14(&tia, &C_anti, &newtia, 0, 1, -1, 1);
  dpd_dot13(&tIA, &D, &newtia, 0, 0, 1, 1);

  dpd_file2_close(&tIA);  
  dpd_file2_close(&tia);

  dpd_buf4_close(&C_anti);  
  dpd_buf4_close(&D);

  dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 7, 10, 7, 0, "F <ia||bc> (ia,b>c)");
  dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");

  dpd_contract442(&tIJAB, &F_anti, &newtIA, 1, 1, 1, 1);
  dpd_contract442(&tijab, &F_anti, &newtia, 1, 1, 1, 1);

  dpd_buf4_close(&tIJAB);
  dpd_buf4_close(&tijab);
  dpd_buf4_close(&F_anti);

  dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&tiJaB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

  dpd_contract442(&tiJaB, &F, &newtIA, 1, 1, 1, 1);
  dpd_contract442(&tIjAb, &F, &newtia, 1, 1, 1, 1);
  
  dpd_buf4_close(&F);  
  dpd_buf4_close(&tIjAb);  
  dpd_buf4_close(&tiJaB);

  dpd_buf4_init(&E_anti, CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
  dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");

  dpd_contract442(&E_anti, &tIJAB, &newtIA, 1, 3, -1, 1);
  dpd_contract442(&E_anti, &tijab, &newtia, 1, 3, -1, 1);

  dpd_buf4_close(&E_anti);  
  dpd_buf4_close(&tIJAB);  
  dpd_buf4_close(&tijab);
  
  dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&tiJaB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

  dpd_contract442(&E, &tiJaB, &newtIA, 1, 3, -1, 1);
  dpd_contract442(&E, &tIjAb, &newtia, 1, 3, -1, 1);

  dpd_buf4_close(&E);  
  dpd_buf4_close(&tIjAb);  
  dpd_buf4_close(&tiJaB);

  dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&dIA, &newtIA);
  dpd_file2_close(&dIA);

  dpd_file2_init(&dia, CC_OEI, 0, 0, 1, "dia");
  dpd_file2_dirprd(&dia, &newtia);
  dpd_file2_close(&dia);

  dpd_file2_close(&newtIA);
  dpd_file2_close(&newtia);

/*  timer_off("t1_build", outfile); */
}
