#include <stdio.h>
#include <math.h>
#include <string.h>
#define EXTERN
#include "globals.h"

/* this function determines R0, properly normalizes R, and checks orthogonality
 * with the ground state left eigenvector (1+lambda) */

void rzero(int C_irr) {
  double ra, rb, r2aa, r2bb, r2ab, rzero=0.0, energy, norm, dotval;
  dpdfile2 RIA, Ria, RIA2, Ria2, FIA, Fia, LIA, Lia;
  dpdbuf4 RIJAB, Rijab, RIjAb, D, R2, LIJAB, Lijab, LIjAb;
  dpdbuf4 fRIJAB, fRijab, fRIjAb;

  /* Calculate R0 consistent with R1 and R2 */
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
  else {
    rzero = 0.0;
  }

  /* normalize full R */
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
  fprintf(outfile,"R0 of normalized R = %15.10lf\n", rzero);

  /* testing */
  /*
     dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
     dpd_file2_init(&Ria, CC_RAMPS, C_irr, 0, 1, "Ria");
     dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "RIJAB");
     dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "Rijab");
     dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");

     dpd_file2_copy(&RIA, CC_OEI, "LIA");
     dpd_file2_copy(&Ria, CC_OEI, "Lia");
     dpd_buf4_copy(&RIJAB, CC_LAMPS, "LIJAB");
     dpd_buf4_copy(&Rijab, CC_LAMPS, "Lijab");
     dpd_buf4_copy(&RIjAb, CC_LAMPS, "LIjAb");

     dpd_file2_init(&LIA, CC_OEI, C_irr, 0, 1, "LIA");
     dpd_file2_init(&Lia, CC_OEI, C_irr, 0, 1, "Lia");
     dpd_buf4_init(&LIJAB, CC_LAMPS, C_irr, 2, 7, 2, 7, 0, "LIJAB");
     dpd_buf4_init(&Lijab, CC_LAMPS, C_irr, 2, 7, 2, 7, 0, "Lijab");
     dpd_buf4_init(&LIjAb, CC_LAMPS, C_irr, 0, 5, 0, 5, 0, "LIjAb");

     ra = dpd_file2_dot(&RIA, &LIA);
     rb = dpd_file2_dot(&Ria, &Lia);
     r2aa = dpd_buf4_dot(&RIJAB, &LIJAB);
     r2bb = dpd_buf4_dot(&Rijab, &Lijab);
     r2ab = dpd_buf4_dot(&RIjAb, &LIjAb);
     dotval = ra + rb + r2aa + r2bb + r2ab;
     fprintf(outfile,"R dot L (Rcopy) %15.10lf\n",dotval);

     dpd_file2_close(&LIA);
     dpd_file2_close(&Lia);
     dpd_buf4_close(&LIJAB);
     dpd_buf4_close(&Lijab);
     dpd_buf4_close(&LIjAb);

     dpd_file2_close(&RIA);
     dpd_file2_close(&Ria);
     dpd_buf4_close(&RIJAB);
     dpd_buf4_close(&Rijab);
     dpd_buf4_close(&RIjAb);
   */
  /* end testing junk */

  if (eom_params.dot_with_Lg) {
    /* evaluate check <Rx|Lg> == 0 */
    if (C_irr == H_IRR ) {
      dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
      dpd_file2_init(&Ria, CC_RAMPS, C_irr, 0, 1, "Ria");
      dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "RIJAB");
      dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "Rijab");
      dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");

      dpd_file2_init(&LIA, CC_OEI, H_IRR, 0, 1, "LIA");
      dpd_file2_init(&Lia, CC_OEI, H_IRR, 0, 1, "Lia");
      dpd_buf4_init(&LIJAB, CC_LAMPS, H_IRR, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_init(&Lijab, CC_LAMPS, H_IRR, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_init(&LIjAb, CC_LAMPS, H_IRR, 0, 5, 0, 5, 0, "LIjAb");

      ra = dpd_file2_dot(&RIA, &LIA);
      rb = dpd_file2_dot(&Ria, &Lia);
      r2aa = dpd_buf4_dot(&RIJAB, &LIJAB);
      r2bb = dpd_buf4_dot(&Rijab, &Lijab);
      r2ab = dpd_buf4_dot(&RIjAb, &LIjAb);
      dotval = rzero + ra + rb + r2aa + r2bb + r2ab;

      fprintf(outfile,"ra = %15.10lf\n", ra);
      fprintf(outfile,"rb = %15.10lf\n", rb);
      fprintf(outfile,"r2aa = %15.10lf\n", r2aa);
      fprintf(outfile,"r2bb = %15.10lf\n", r2bb);
      fprintf(outfile,"r2ab = %15.10lf\n", r2ab);

      dpd_file2_close(&LIA);
      dpd_file2_close(&Lia);
      dpd_buf4_close(&LIJAB);
      dpd_buf4_close(&Lijab);
      dpd_buf4_close(&LIjAb);
      dpd_file2_close(&RIA);
      dpd_file2_close(&Ria);
      dpd_buf4_close(&RIJAB);
      dpd_buf4_close(&Rijab);
      dpd_buf4_close(&RIjAb);
    }
    else {
      dotval = 0.0;
    }

    fprintf(outfile,"\n Check of orthogonality: <Rx|Lg> = %15.10lf\n", dotval);
  }

  return;
}

void rzero_rhf(int C_irr) {
  double r1, r2, rzero=0.0, energy, norm, dotval;
  double dot_IA, dot_ia, dot_IJAB, dot_ijab, dot_IjAb;
  dpdfile2 RIA, FIA, LIA, Lia, Ria;
  dpdbuf4 RIjAb, RIjbA, RIjAb1, RIjbA1, D, R2, LIjAb, RIJAB, Rijab;
  dpdbuf4 LIJAB, Lijab;

  /* produce RIjbA and 2RIjAb-RIjbA */
  dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_sort(&RIjAb, CC_RAMPS, pqsr, 0, 5, "RIjbA");
  dpd_buf4_copy(&RIjAb, CC_RAMPS, "2RIjAb - RIjbA");
  dpd_buf4_close(&RIjAb);

  dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
  dpd_buf4_scm(&RIjAb, 2.0);
  dpd_buf4_init(&RIjbA, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjbA");
  dpd_buf4_axpy(&RIjbA, &RIjAb, -1.0);
  dpd_buf4_close(&RIjbA);
  dpd_buf4_close(&RIjAb);

  /* calculate R0 consistent with R1 and R2 */
  if (C_irr == H_IRR) {
    dpd_file2_init(&FIA, CC_OEI, H_IRR, 0, 1, "FME");
    dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
    r1 = 2.0 * dpd_file2_dot(&FIA, &RIA);
    dpd_file2_close(&RIA);
    dpd_file2_close(&FIA);

    dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
    dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    r2 = dpd_buf4_dot(&D, &RIjAb);
    dpd_buf4_close(&D);
    dpd_buf4_close(&RIjAb);

    psio_read_entry(CC_INFO, "CCEOM Energy", (char *) &energy, sizeof(double));
    rzero = (r1 + r2)/energy;
  }
  else {
    rzero = 0.0;
  }

  /* normalize full R */
  dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
  dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_init(&RIjbA, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjbA");

  norm = norm_C_rhf(&RIA, &RIjAb, &RIjbA);
  norm *= norm;
  norm += rzero * rzero;
  norm = sqrt(norm);
  rzero = rzero / norm;
  dpd_file2_scm(&RIA, 1.0/norm);
  dpd_buf4_scm(&RIjAb, 1.0/norm);
  dpd_buf4_scm(&RIjbA, 1.0/norm);

  dpd_file2_close(&RIA);
  dpd_buf4_close(&RIjAb);
  dpd_buf4_close(&RIjbA);

  psio_write_entry(CC_INFO, "EOM R0", (char *) &rzero, sizeof(double));
  fprintf(outfile,"R0 of normalized R   = %15.10lf\n", rzero);

  /* produce Ria, RIJAB, Rijab, and 2RIjAb-RIjbA */
  dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
  dpd_file2_copy(&RIA, CC_RAMPS, "Ria");
  dpd_file2_close(&RIA);
  dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_buf4_copy(&RIjAb, CC_RAMPS, "2RIjAb - RIjbA");
  dpd_buf4_close(&RIjAb);
  dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 2, 7, 0, 5, 1, "RIjAb");
  dpd_buf4_copy(&RIjAb, CC_RAMPS, "RIJAB");
  dpd_buf4_copy(&RIjAb, CC_RAMPS, "Rijab");
  dpd_buf4_close(&RIjAb);

  dpd_buf4_init(&RIjbA, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjbA");
  dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
  dpd_buf4_scm(&RIjAb, 2.0);
  dpd_buf4_axpy(&RIjbA, &RIjAb, -1.0);
  dpd_buf4_close(&RIjAb);
  dpd_buf4_close(&RIjbA);

  if (eom_params.dot_with_Lg) {
    /* check orthogonality with ground state (1+lambda) */
    if (C_irr == H_IRR) {
      dpd_file2_init(&LIA, CC_OEI, H_IRR, 0, 1, "LIA");
      dpd_file2_init(&RIA, CC_RAMPS, H_IRR, 0, 1, "RIA");
      r1 = 2.0 * dpd_file2_dot(&LIA, &RIA);
      dpd_file2_close(&RIA);
      dpd_file2_close(&LIA);

      dpd_buf4_init(&LIjAb, CC_LAMPS, H_IRR, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
      r2 = dpd_buf4_dot(&LIjAb, &RIjAb);
      dpd_buf4_close(&RIjAb);
      dpd_buf4_close(&LIjAb);
    }
    else {
      r1 = r2 = 0.0;
    }
    dotval = r1 + r2 + rzero;
    fprintf(outfile,"<L0|R0>              = %15.10lf\n", 1.0*rzero);
    fprintf(outfile,"2*<LIAg|RIAx>        = %15.10lf\n",r1);
    fprintf(outfile,"<LIjAb|2RIjAb-RIjbA> = %15.10lf\n",r2);
    fprintf(outfile,"<Lg|Rx>              = %15.10lf\n", dotval);

     /* double check orthogonality rohf-like */
    fprintf(outfile,"\nROHF-like orthogonality test\n");
      dpd_file2_init(&LIA, CC_OEI, H_IRR, 0, 1, "LIA");
      dpd_file2_init(&RIA, CC_RAMPS, H_IRR, 0, 1, "RIA");
      dot_IA = dpd_file2_dot(&LIA, &RIA);
    fprintf(outfile,"<LIAg|RIAx>          = %15.10lf\n", dot_IA);
      dpd_file2_close(&RIA);
      dpd_file2_close(&LIA);
      dpd_file2_init(&Lia, CC_OEI, H_IRR, 0, 1, "Lia");
      dpd_file2_init(&Ria, CC_RAMPS, H_IRR, 0, 1, "Ria");
      dot_ia = dpd_file2_dot(&Lia, &Ria);
    fprintf(outfile,"<Liag|Riax>          = %15.10lf\n", dot_ia);
      dpd_file2_close(&Ria);
      dpd_file2_close(&Lia);

      dpd_buf4_init(&LIJAB, CC_LAMPS, H_IRR, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "RIJAB");
      dot_IJAB = dpd_buf4_dot(&LIJAB, &RIJAB);
    fprintf(outfile,"<LIJABg|RIJABx>      = %15.10lf\n", dot_IJAB);
      dpd_buf4_close(&RIJAB);
      dpd_buf4_close(&LIJAB);

      dpd_buf4_init(&Lijab, CC_LAMPS, H_IRR, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "Rijab");
      dot_ijab = dpd_buf4_dot(&Lijab, &Rijab);
    fprintf(outfile,"<Lijabg|Rijabx>      = %15.10lf\n", dot_ijab);
      dpd_buf4_close(&Rijab);
      dpd_buf4_close(&Lijab);

      dpd_buf4_init(&LIjAb, CC_LAMPS, H_IRR, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
      dot_IjAb = dpd_buf4_dot(&LIjAb, &RIjAb);
    fprintf(outfile,"<LIjAbg|RIjAbx>      = %15.10lf\n", dot_IjAb);
      dpd_buf4_close(&RIjAb);
      dpd_buf4_close(&LIjAb);
      fprintf(outfile,"rohf-like <Lg|Rx>    = %15.10lf\n", rzero + dot_IA +
          dot_ia + dot_IJAB + dot_ijab + dot_IjAb);
    }

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
