#include <stdio.h>
#include <math.h>
#include <string.h>
#define EXTERN
#include "globals.h"

void rzero(int C_irr) {
  double ra, rb, r2aa, r2bb, r2ab, rzero=0.0, energy, norm;
  dpdfile2 RIA, Ria, RIA2, Ria2, FIA, Fia;
  dpdbuf4 RIJAB, Rijab, RIjAb, D, R2;
  dpdbuf4 fRIJAB, fRijab, fRIjAb;

  if (C_irr == H_IRR) {
    dpd_file2_init(&FIA, CC_OEI, H_IRR, 0, 1, "FME");
    dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
    ra = dpd_file2_dot(&FIA, &RIA);
    dpd_file2_close(&RIA);
    dpd_file2_close(&FIA);

    dpd_file2_init(&Fia, CC_OEI, H_IRR, 0, 1, "Fme");
    dpd_file2_init(&Ria, CC_RAMPS, C_irr, 0, 1, "Ria");
    rb = dpd_file2_dot(&Fia, &Ria);
    dpd_file2_close(&Ria);
    dpd_file2_close(&Fia);

    dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "RIJAB");
    r2aa = dpd_buf4_dot(&D, &RIJAB);
    dpd_buf4_close(&RIJAB);
    dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "Rijab");
    r2bb = dpd_buf4_dot(&D, &Rijab);
    dpd_buf4_close(&Rijab);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
    r2ab = dpd_buf4_dot(&D, &RIjAb);
    dpd_buf4_close(&RIjAb);
    dpd_buf4_close(&D);

    psio_read_entry(CC_INFO, "CCEOM Energy", (char *) &energy, sizeof(double));
    rzero = (ra + rb + r2aa + r2bb + r2ab)/energy;
  }

  dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
  dpd_file2_init(&Ria, CC_RAMPS, C_irr, 0, 1, "Ria");
  dpd_buf4_init(&fRIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "RIJAB");
  dpd_buf4_init(&fRijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "Rijab");
  dpd_buf4_init(&fRIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");

  norm = norm_C(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb);
  norm *= norm;
  norm += rzero * rzero;
  norm = sqrt(norm);
  rzero = rzero / norm;
  scm_C(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb, 1.0/norm);

  dpd_file2_close(&RIA);
  dpd_file2_close(&Ria);
  dpd_buf4_close(&fRIJAB);
  dpd_buf4_close(&fRijab);
  dpd_buf4_close(&fRIjAb);

  psio_write_entry(CC_INFO, "EOM R0", (char *) &rzero, sizeof(double));

  /* Now make useful copies in RAMPS */
  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_sort(&R2, CC_RAMPS, qprs, 0, 5, "RjIAb");
  dpd_buf4_close(&R2);
  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RjIAb");
  dpd_buf4_sort(&R2, CC_RAMPS, pqsr, 0, 5, "RiJaB");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 2, 7, 0, "RIJAB");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RIAJB");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 2, 7, 0, "Rijab");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "Riajb");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RIAjb");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RiaJB");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 10, 10, 10, 10, 0, "RIAjb");
  dpd_buf4_sort(&R2, CC_RAMPS, psrq, 10, 10, "RIbjA");
  dpd_buf4_sort(&R2, CC_RAMPS, rqps, 10, 10, "RjAIb");
  dpd_buf4_close(&R2);

  return;
}

void rzero_rhf(int C_irr) {
  double r1, r2, rzero=0.0, energy, norm;
  dpdfile2 RIA, FIA;
  dpdbuf4 RIjAb, RIjbA, RIjAb1, RIjbA1, D, R2;

  dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_sort(&RIjAb, CC_RAMPS, pqsr, 0, 5, "RIjbA");
  dpd_buf4_close(&RIjAb);

  if (C_irr == H_IRR) {
    dpd_file2_init(&FIA, CC_OEI, H_IRR, 0, 1, "FME");
    dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
    r1 = 2.0 * dpd_file2_dot(&FIA, &RIA);
    dpd_file2_close(&RIA);
    dpd_file2_close(&FIA);

    dpd_buf4_init(&RIjAb1, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_scm(&RIjAb1, 2.0);
    dpd_buf4_init(&RIjbA, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjbA");
    dpd_buf4_axpy(&RIjbA, &RIjAb1, -1.0);
    dpd_buf4_close(&RIjbA);

    dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    r2 = dpd_buf4_dot(&D, &RIjAb1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&RIjAb1);

    psio_read_entry(CC_INFO, "CCEOM Energy", (char *) &energy, sizeof(double));
    rzero = (r1 + r2)/energy;
  }

  dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
  dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_init(&RIjbA, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjbA");

  norm = norm_C_rhf(&RIA, &RIjAb, &RIjbA);
  /*    norm *= norm;
        norm += rzero * rzero;
        norm = sqrt(norm);
        rzero = rzero / norm;
      
        dpd_file2_scm(&RIA, 1.0/norm);
        dpd_buf4_scm(&RIjAb, 1.0/norm);
  */
  dpd_file2_close(&RIA);
  dpd_buf4_close(&RIjAb);

  /*  psio_write_entry(CC_INFO, "EOM R0", (char *) &rzero, sizeof(double)); */

  /* Now make useful copies in RAMPS */
  /*
  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_sort(&R2, CC_RAMPS, qprs, 0, 5, "RjIAb");
  dpd_buf4_close(&R2);
  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RjIAb");
  dpd_buf4_sort(&R2, CC_RAMPS, pqsr, 0, 5, "RiJaB");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RIAjb");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RiaJB");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, C_irr, 10, 10, 10, 10, 0, "RIAjb");
  dpd_buf4_sort(&R2, CC_RAMPS, psrq, 10, 10, "RIbjA");
  dpd_buf4_sort(&R2, CC_RAMPS, rqps, 10, 10, "RjAIb");
  dpd_buf4_close(&R2);
  */

  return;
}
