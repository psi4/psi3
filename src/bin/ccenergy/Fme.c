#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Fme_build(void)
{
  struct oe_dpdfile FME, Fme, fIA, fia, tIA, tia;
  struct dpdbuf D_anti, D;

  dpd_oe_file_init(&fIA, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_copy(&fIA, CC_OEI, "FME", 0, outfile);
  dpd_oe_file_close(&fIA);

  dpd_oe_file_init(&fia, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_copy(&fia, CC_OEI, "Fme", 0, outfile);
  dpd_oe_file_close(&fia);
  
  dpd_oe_file_init(&FME, CC_OEI, 0, 1, "FME", 0, outfile);
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "Fme", 0, outfile);
  
  dpd_buf_init(&D_anti, CC_DINTS, 0, 5, 0, 5, 0, "D <ij||ab>",
	       0, outfile);
  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot13(&tia, &D, &FME, 0, 0, 1.0, 1.0, 0, outfile);

  dpd_dot13(&tia, &D_anti, &Fme, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot13(&tIA, &D, &Fme, 0, 0, 1.0, 1.0, 0, outfile);

  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&tia);
  dpd_buf_close(&D_anti);
  dpd_buf_close(&D);

  dpd_oe_file_close(&FME);
  dpd_oe_file_close(&Fme);
}
