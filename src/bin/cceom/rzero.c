#include <stdio.h>
#include <math.h>
#include <string.h>
#define EXTERN
#include "globals.h"

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
            dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);
extern void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
            dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);

/* this function determines R0, properly normalizes R, and checks orthogonality
 * with the ground state left eigenvector (1+lambda) */

/* for ROHF and UHF */
void rzero(int C_irr) {
  double rzero=0.0, energy, norm, dotval;
  double dot_IA, dot_ia, dot_IJAB, dot_ijab, dot_IjAb;
  dpdfile2 RIA, Ria, RIA2, Ria2, FIA, Fia, LIA, Lia;
  dpdbuf4 RIJAB, Rijab, RIjAb, D, R2, LIJAB, Lijab, LIjAb;
  dpdbuf4 fRIJAB, fRijab, fRIjAb;
  int L_irr;
  int A_OCC, B_OCC, A_VIR, B_VIR;
  int AA_OCC, AA_VIR, BB_OCC, BB_VIR, AB_OCC, AB_VIR;

  A_OCC = 0; A_VIR = 1;
  AA_OCC = 2; AA_VIR = 7;
  if (params.eom_ref <= 1) {
    B_OCC = 0; B_VIR = 1;
    BB_OCC = 2; BB_VIR = 7;
    AB_OCC = 0; AB_VIR = 5;
  }
  else if (params.eom_ref == 2) {
    B_OCC = 2; B_VIR = 3;
    BB_OCC = 12; BB_VIR = 17;
    AB_OCC = 22; AB_VIR = 28;
  }

  L_irr = eom_params.L_irr;
  if ( psio_tocscan(CC_INFO, "CCEOM Energy") == NULL) {
    fprintf(outfile,"No CCEOM energy found in CC_INFO.  Not normalizing R.\n");
    return;
  };
  psio_read_entry(CC_INFO, "CCEOM Energy", (char *) &energy, sizeof(double));

  /* Calculate R0 consistent with R1 and R2 */
  if (C_irr == H_IRR) {
    dpd_file2_init(&FIA, CC_OEI, H_IRR, A_OCC, A_VIR, "FME");
    dpd_file2_init(&RIA, CC_RAMPS, C_irr, A_OCC, A_VIR, "RIA");
    dot_IA = dpd_file2_dot(&FIA, &RIA);
    dpd_file2_close(&RIA);
    dpd_file2_close(&FIA);

    dpd_file2_init(&Fia, CC_OEI, H_IRR, B_OCC, B_VIR, "Fme");
    dpd_file2_init(&Ria, CC_RAMPS, C_irr, B_OCC, B_VIR, "Ria");
    dot_ia = dpd_file2_dot(&Fia, &Ria);
    dpd_file2_close(&Ria);
    dpd_file2_close(&Fia);

    if (params.eom_ref == 1) {
      dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
      dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "RIJAB");
      dot_IJAB = dpd_buf4_dot(&D, &RIJAB);
      dpd_buf4_close(&RIJAB);
      dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "Rijab");
      dot_ijab = dpd_buf4_dot(&D, &Rijab);
      dpd_buf4_close(&Rijab);
      dpd_buf4_close(&D);
  
      dpd_buf4_init(&D, CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
      dot_IjAb = dpd_buf4_dot(&D, &RIjAb);
      dpd_buf4_close(&RIjAb);
      dpd_buf4_close(&D);
    }
    else if (params.eom_ref == 2) {
      dpd_buf4_init(&D, CC_DINTS, H_IRR, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
      dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "RIJAB");
      dot_IJAB = dpd_buf4_dot(&D, &RIJAB);
      dpd_buf4_close(&RIJAB);
      dpd_buf4_close(&D);

      dpd_buf4_init(&D, CC_DINTS, H_IRR, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
      dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 12, 17, 12, 17, 0, "Rijab");
      dot_ijab = dpd_buf4_dot(&D, &Rijab);
      dpd_buf4_close(&Rijab);
      dpd_buf4_close(&D);
                         
      dpd_buf4_init(&D, CC_DINTS, H_IRR, 22, 28, 22, 28, 0, "D <Ij|Ab>");
      dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 22, 28, 22, 28, 0, "RIjAb");
      dot_IjAb = dpd_buf4_dot(&D, &RIjAb);
      dpd_buf4_close(&RIjAb);
      dpd_buf4_close(&D);
    }

    rzero = (dot_IA + dot_ia + dot_IJAB + dot_ijab + dot_IjAb)/energy;
  }
  else {
    rzero = 0.0;
  }
  fprintf(outfile,"R0 of unnormalized R = %15.10lf\n", rzero);

  /* normalize full R */
  dpd_file2_init(&RIA, CC_RAMPS, C_irr, A_OCC, A_VIR, "RIA");
  dpd_file2_init(&Ria, CC_RAMPS, C_irr, B_OCC, B_VIR, "Ria");
  dpd_buf4_init(&fRIJAB, CC_RAMPS, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, "RIJAB");
  dpd_buf4_init(&fRijab, CC_RAMPS, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, "Rijab");
  dpd_buf4_init(&fRIjAb, CC_RAMPS, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, "RIjAb");

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

  fprintf(outfile,"Writing to CC_INFO...\n");
  psio_write_entry(CC_INFO, "EOM R Irrep", (char *) &C_irr, sizeof(int));
  fprintf(outfile,"Symmetry of R vector = %5s\n", moinfo.labels[C_irr]);
  psio_write_entry(CC_INFO, "EOM R0", (char *) &rzero, sizeof(double));
  /* energy actually written in diag */
  fprintf(outfile,"CCEOM Energy         = %15.10lf\n", energy);
  fprintf(outfile,"R0 of normalized R   = %15.10lf\n", rzero);

  if (eom_params.dot_with_L) {
    /* evaluate check <R|L> == 0 */
    if (C_irr == L_irr ) {
      dpd_file2_init(&RIA, CC_RAMPS, C_irr, A_OCC, A_VIR, "RIA");
      dpd_file2_init(&Ria, CC_RAMPS, C_irr, B_OCC, B_VIR, "Ria");
      dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, "RIJAB");
      dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, "Rijab");
      dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, "RIjAb");

      dpd_file2_init(&LIA, CC_OEI, L_irr, A_OCC, A_VIR, "LIA");
      dpd_file2_init(&Lia, CC_OEI, L_irr, B_OCC, B_VIR, "Lia");
      dpd_buf4_init(&LIJAB, CC_LAMPS, L_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, "LIJAB");
      dpd_buf4_init(&Lijab, CC_LAMPS, L_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, "Lijab");
      dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, "LIjAb");

      fprintf(outfile,"\nROHF orthogonality test\n");
      fprintf(outfile,"<L0|R0>            = %15.10lf\n", eom_params.L0 * rzero);
      dot_IA = dpd_file2_dot(&LIA, &RIA);
      fprintf(outfile,"<LIA|RIA>          = %15.10lf\n", dot_IA);
      dot_ia = dpd_file2_dot(&Lia, &Ria);
      fprintf(outfile,"<Lia|Ria>          = %15.10lf\n", dot_ia);
      dot_IJAB = dpd_buf4_dot(&LIJAB, &RIJAB);
      fprintf(outfile,"<LIJAB|RIJAB>      = %15.10lf\n", dot_IJAB);
      dot_ijab = dpd_buf4_dot(&Lijab, &Rijab);
      fprintf(outfile,"<Lijab|Rijab>      = %15.10lf\n", dot_ijab);
      dot_IjAb = dpd_buf4_dot(&LIjAb, &RIjAb);
      fprintf(outfile,"<LIjAb|RIjAb>      = %15.10lf\n", dot_IjAb);
      fprintf(outfile,"<L|R>              = %15.10lf\n", (eom_params.L0 * rzero)
          + dot_IA + dot_ia + dot_IJAB + dot_ijab + dot_IjAb);

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
      fprintf(outfile,"\nOverlap <R|L> zero by symmetry\n");
    }
  }

  return;
}

void rzero_rhf(int C_irr) {
  double r1, r2, rzero=0.0, energy, norm, dotval;
  double dot_IA, dot_ia, dot_IJAB, dot_ijab, dot_IjAb;
  dpdfile2 RIA, FIA, LIA, Lia, Ria;
  dpdbuf4 RIjAb, RIjbA, RIjAb1, RIjbA1, D, R2, LIjAb, RIJAB, Rijab;
  dpdbuf4 LIJAB, Lijab;
  int L_irr;

  L_irr = eom_params.L_irr;
  psio_read_entry(CC_INFO, "CCEOM Energy", (char *) &energy, sizeof(double));

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

    rzero = (r1 + r2)/energy;
  }
  else {
    rzero = 0.0;
  }

  fprintf(outfile,"R0 of unnormalized R = %15.10lf\n", rzero);
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

  fprintf(outfile,"Writing to CC_INFO...\n");
  psio_write_entry(CC_INFO, "EOM R Irrep", (char *) &C_irr, sizeof(int));
  fprintf(outfile,"Irrep of R vector    = %s\n", moinfo.labels[C_irr]);
  psio_write_entry(CC_INFO, "EOM R0", (char *) &rzero, sizeof(double));
  fprintf(outfile,"EOM energy           = %15.10lf\n", energy);
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

  dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
  norm = dpd_file2_dot_self(&RIA);
  dpd_file2_close(&RIA);
  dpd_file2_init(&Ria, CC_RAMPS, C_irr, 0, 1, "Ria");
  norm += dpd_file2_dot_self(&Ria);
  dpd_file2_close(&Ria);
  dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "RIJAB");
  norm += dpd_buf4_dot_self(&RIJAB);
  dpd_buf4_close(&RIJAB);
  dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "Rijab");
  norm += dpd_buf4_dot_self(&Rijab);
  dpd_buf4_close(&Rijab);
  dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");
  norm += dpd_buf4_dot_self(&RIjAb);
  dpd_buf4_close(&RIjAb);
  norm += rzero * rzero;
  /* if (fabs(norm - 1.0) > 1E-6) */
  fprintf(outfile,"checked normalization <R|R> = %15.10lf\n", norm);

  if (eom_params.dot_with_L) {
    /* check orthogonality with ground state (1+lambda) */
    if (C_irr == L_irr) {
      dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
      dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
      r1 = 2.0 * dpd_file2_dot(&LIA, &RIA);
      dpd_file2_close(&RIA);
      dpd_file2_close(&LIA);

      dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "2RIjAb - RIjbA");
      r2 = dpd_buf4_dot(&LIjAb, &RIjAb);
      dpd_buf4_close(&RIjAb);
      dpd_buf4_close(&LIjAb);

      dotval = r1 + r2 + (eom_params.L0 * rzero);
      fprintf(outfile,"Performing RHF orthogonality test\n");
      fprintf(outfile,"<L0|R0>              = %15.10lf\n", eom_params.L0 * rzero);
      fprintf(outfile,"2*<LIA|RIA>          = %15.10lf\n",r1);
      fprintf(outfile,"<LIjAb|2RIjAb-RIjbA> = %15.10lf\n",r2);
      fprintf(outfile,"<L|R>                = %15.10lf\n", dotval);
    }
    else {
      fprintf(outfile,"<L|R> zero by symmetry\n");
    }

    if (C_irr == L_irr ) {
     /* double check orthogonality rohf-like */
      dpd_file2_init(&RIA, CC_RAMPS, C_irr, 0, 1, "RIA");
      dpd_file2_init(&Ria, CC_RAMPS, C_irr, 0, 1, "Ria");
      dpd_buf4_init(&RIJAB, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "RIJAB");
      dpd_buf4_init(&Rijab, CC_RAMPS, C_irr, 2, 7, 2, 7, 0, "Rijab");
      dpd_buf4_init(&RIjAb, CC_RAMPS, C_irr, 0, 5, 0, 5, 0, "RIjAb");

      dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
      dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "Lia");
      dpd_buf4_init(&LIJAB, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_init(&Lijab, CC_LAMPS, L_irr, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");

      dot_IA = dpd_file2_dot(&LIA, &RIA);
      dot_ia = dpd_file2_dot(&Lia, &Ria);
      dot_IJAB = dpd_buf4_dot(&LIJAB, &RIJAB);
      dot_ijab = dpd_buf4_dot(&Lijab, &Rijab);
      dot_IjAb = dpd_buf4_dot(&LIjAb, &RIjAb);

      dpd_file2_close(&RIA);
      dpd_file2_close(&Ria);
      dpd_buf4_close(&RIJAB);
      dpd_buf4_close(&Rijab);
      dpd_buf4_close(&RIjAb);

      dpd_file2_close(&LIA);
      dpd_file2_close(&Lia);
      dpd_buf4_close(&LIJAB);
      dpd_buf4_close(&Lijab);
      dpd_buf4_close(&LIjAb);

      fprintf(outfile,"\nROHF-like orthogonality test\n");
      fprintf(outfile,"<L0|R0>              = %15.10lf\n", eom_params.L0 * rzero);
      fprintf(outfile,"<LIA|RIA>            = %15.10lf\n", dot_IA);
      fprintf(outfile,"<Lia|Ria>            = %15.10lf\n", dot_ia);
      fprintf(outfile,"<LIJAB|RIJAB>        = %15.10lf\n", dot_IJAB);
      fprintf(outfile,"<Lijab|Rijab>        = %15.10lf\n", dot_ijab);
      fprintf(outfile,"<LIjAb|RIjAb>        = %15.10lf\n", dot_IjAb);
      fprintf(outfile,"<L|R>                = %15.10lf\n", eom_params.L0 * rzero +
          dot_IA + dot_ia + dot_IJAB + dot_ijab + dot_IjAb);
    }
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
