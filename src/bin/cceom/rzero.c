#include <stdio.h>
#include <math.h>
#include <string.h>
#define EXTERN
#include "globals.h"

void rzero(void) {
  int irrep;
  double ra, rb, r2aa, r2bb, r2ab, rzero, energy, norm;
  dpdfile2 RIA, Ria, RIA2, Ria2, FIA, Fia;
  dpdbuf4 RIJAB, Rijab, RIjAb, D, R2;
  dpdbuf4 fRIJAB, fRijab, fRIjAb;

  psio_read_entry(CC_INFO, "CCEOM Irrep", (char *) &(irrep), sizeof(int));

  dpd_file2_init(&FIA, CC_OEI, irrep, 0, 1, "FME");
  dpd_file2_init(&RIA, CC_OEI, irrep, 0, 1, "RIA");
  ra = dpd_file2_dot(&FIA, &RIA);
  dpd_file2_close(&RIA);
  dpd_file2_close(&FIA);

  dpd_file2_init(&Fia, CC_OEI, irrep, 0, 1, "Fme");
  dpd_file2_init(&Ria, CC_OEI, irrep, 0, 1, "Ria");
  rb = dpd_file2_dot(&Fia, &Ria);
  dpd_file2_close(&Ria);
  dpd_file2_close(&Fia);

  dpd_buf4_init(&D, CC_DINTS, irrep, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_init(&RIJAB, CC_RAMPS, irrep, 2, 7, 2, 7, 0, "RIJAB");
  r2aa = dpd_buf4_dot(&D, &RIJAB);
  dpd_buf4_close(&RIJAB);
  dpd_buf4_init(&Rijab, CC_RAMPS, irrep, 2, 7, 2, 7, 0, "Rijab");
  r2bb = dpd_buf4_dot(&D, &Rijab);
  dpd_buf4_close(&Rijab);
  dpd_buf4_close(&D);

  dpd_buf4_init(&D, CC_DINTS, irrep, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&RIjAb, CC_RAMPS, irrep, 0, 5, 0, 5, 0, "RIjAb");
  r2ab = dpd_buf4_dot(&D, &RIjAb);
  dpd_buf4_close(&RIjAb);
  dpd_buf4_close(&D);

  psio_read_entry(CC_INFO, "CCEOM Energy", (char *) &energy, sizeof(double));
  rzero = (ra + rb + r2aa + r2bb + r2ab)/energy;

  /*
    fprintf(outfile, "One A R0 = %20.14f\n", ra);
    fprintf(outfile, "One B R0 = %20.14f\n", rb);
    fprintf(outfile, "Two AA R0 = %20.14f\n", r2aa);
    fprintf(outfile, "Two BB R0 = %20.14f\n", r2bb);
    fprintf(outfile, "Two AB R0 = %20.14f\n", r2ab);
    fprintf(outfile, "CCEOM Energy Read = %20.14f\n", energy);
    fprintf(outfile, "Value of R0 before norm = %20.14f\n", rzero);
  */

  dpd_file2_init(&RIA, CC_OEI, irrep, 0, 1, "RIA");
  dpd_file2_init(&Ria, CC_OEI, irrep, 0, 1, "Ria");
  dpd_buf4_init(&fRIJAB, CC_RAMPS, irrep, 2, 7, 2, 7, 0, "RIJAB");
  dpd_buf4_init(&fRijab, CC_RAMPS, irrep, 2, 7, 2, 7, 0, "Rijab");
  dpd_buf4_init(&fRIjAb, CC_RAMPS, irrep, 0, 5, 0, 5, 0, "RIjAb");

  /*
    norm = norm_v(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb);
    fprintf(outfile,"norm of R before clean: %18.13lf\n",norm);
    fflush(outfile);
    c_clean(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb);
  */

  norm = norm_C(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb);
  norm *= norm;
  norm += rzero * rzero;
  norm = sqrt(norm);
  rzero = rzero / norm;
  scm_C(&RIA, &Ria, &fRIJAB, &fRijab, &fRIjAb, 1.0/norm);

  fprintf(outfile, "Value of R0 in normalized R = %20.14f\n", rzero);

  dpd_file2_close(&RIA);
  dpd_file2_close(&Ria);
  dpd_buf4_close(&fRIJAB);
  dpd_buf4_close(&fRijab);
  dpd_buf4_close(&fRIjAb);

  psio_write_entry(CC_INFO, "EOM R0", (char *) &rzero, sizeof(double));

  /* Now make useful copies in RAMPS */
  dpd_buf4_init(&R2, CC_RAMPS, irrep, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_sort(&R2, CC_TMP0, qprs, 0, 5, "RjIAb");
  dpd_buf4_close(&R2);
  dpd_buf4_init(&R2, CC_TMP0, irrep, 0, 5, 0, 5, 0, "RjIAb");
  dpd_buf4_sort(&R2, CC_RAMPS, pqsr, 0, 5, "RiJaB");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, irrep, 0, 5, 2, 7, 0, "RIJAB");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RIAJB");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, irrep, 0, 5, 2, 7, 0, "Rijab");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "Riajb");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, irrep, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RIAjb");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, irrep, 0, 5, 0, 5, 0, "RiJaB");
  dpd_buf4_sort(&R2, CC_RAMPS, prqs, 10, 10, "RiaJB");
  dpd_buf4_close(&R2);

  dpd_buf4_init(&R2, CC_RAMPS, irrep, 10, 10, 10, 10, 0, "RIAjb");
  dpd_buf4_sort(&R2, CC_RAMPS, psrq, 10, 10, "RIbjA");
  dpd_buf4_sort(&R2, CC_RAMPS, rqps, 10, 10, "RjAIb");
  dpd_buf4_close(&R2);

  return;
}
