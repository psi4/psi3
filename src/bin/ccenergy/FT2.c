#include <stdio.h>
#include <dpd.h>
#include <psio.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void FT2(void)
{
  dpdfile2 tIA, tia;
  dpdbuf4 newtIJAB, newtijab, newtIjAb, t2, t2a, t2b;
  dpdbuf4 F_anti, F;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &tIA, &newtIjAb, 1, 1, 1, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_contract244(&tIA, &F, &newtIjAb, 1, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA); 

    dpd_buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /*** AA ***/

    dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 7, 10, 7, 0, "F <ia||bc> (ia,b>c)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract424(&F_anti, &tIA, &t2, 1, 1, 1, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_axpy(&t2a, &newtIJAB, 1);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&F_anti);

    /*** BB ***/

    dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 7, 10, 7, 0, "F <ia||bc> (ia,b>c)");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract424(&F_anti, &tia, &t2, 1, 1, 1, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_axpy(&t2a, &newtijab, 1);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&F_anti);

    /*** AB ***/

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &tia, &newtIjAb, 1, 1, 1, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_contract244(&tIA, &F, &newtIjAb, 1, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);

    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&newtIjAb);
  }

}
