#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Fme_build(void)
{
  dpdfile2 FME, Fme, fIA, fia, tIA, tia;
  dpdbuf4 D_anti, D;

  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_copy(&fIA, CC_OEI, "FME");
  dpd_file2_close(&fIA);

  dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
  dpd_file2_copy(&fia, CC_OEI, "Fme");
  dpd_file2_close(&fia);
  
  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
  
  dpd_buf4_init(&D_anti, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  dpd_dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
  dpd_dot13(&tia, &D, &FME, 0, 0, 1.0, 1.0);

  dpd_dot13(&tia, &D_anti, &Fme, 0, 0, 1.0, 1.0);
  dpd_dot13(&tIA, &D, &Fme, 0, 0, 1.0, 1.0);

  dpd_file2_close(&tIA);
  dpd_file2_close(&tia);
  dpd_buf4_close(&D_anti);
  dpd_buf4_close(&D);

  dpd_file2_close(&FME);
  dpd_file2_close(&Fme);
}
