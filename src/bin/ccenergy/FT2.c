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
  dpdbuf4 Z;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    /*
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &tIA, &newtIjAb, 1, 1, 1, 1, 1);
    dpd_buf4_close(&F);
    */

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z1(ij,ab)");
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&tIA, &F, &Z, 1, 0, 0, 1, 0);
    dpd_file2_close(&tIA); 
    dpd_buf4_close(&F);

    dpd_buf4_sort(&Z, CC_TMP0, qpsr, 0, 5, "Z2(ji,ba)");
    dpd_buf4_axpy(&Z, &newtIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2(ji,ba)");
    dpd_buf4_axpy(&Z, &newtIjAb, 1.0);
    dpd_buf4_close(&Z);

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
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&newtIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, CC_TAMPS, 0, 10, 17, 12, 17, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /*** AA ***/

    dpd_buf4_init(&F, CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
    dpd_buf4_init(&t2, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_contract424(&F, &tIA, &t2, 1, 1, 1, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_axpy(&t2a, &newtIJAB, 1);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&F);

    /*** BB ***/

    dpd_buf4_init(&F, CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
    dpd_buf4_init(&t2, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    dpd_contract424(&F, &tia, &t2, 1, 1, 1, 1, 0);
    dpd_buf4_sort(&t2, CC_TMP0, qprs, 10, 17, "T (ji,a>b)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    dpd_buf4_init(&t2b, CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ji,a>b)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_axpy(&t2a, &newtijab, 1);
    dpd_buf4_close(&t2b);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&F);

    /*** AB ***/

    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_contract424(&F, &tia, &newtIjAb, 1, 1, 1, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_init(&F, CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
    dpd_contract244(&tIA, &F, &newtIjAb, 1, 2, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);

    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&newtIjAb);

  }

}
