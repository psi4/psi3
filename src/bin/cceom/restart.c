/* restart() collapses L vectors down to num vectors */

#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void restart(double **alpha, int L, int num, int C_irr, int ortho) {
  int i,I,j,h;
  char lbl[20];
  dpdfile2 C1, CME, Cme, SIA, Sia;
  dpdbuf4 C2, CMNEF, Cmnef, CMnEf, SIJAB, Sijab, SIjAb;
  double dotval, norm;

  /* Orthonormalize alpha[1] through alpha[num] */
  if (ortho) {
    for (I=1;I<num;++I) {
      for (i=0; i<I; i++) {
        dotval = 0.0;
        for (j=0;j<L;++j) {
          dotval += alpha[j][i] * alpha[j][I];
        }
        for (j=0; j<L; j++) alpha[j][I] -= dotval * alpha[j][i];
      }
      dotval = 0.0;
      for (j=0;j<L;++j) dotval += alpha[j][I] * alpha[j][I];
      norm = sqrt(dotval);
      for (j=0;j<L;++j) alpha[j][I] = alpha[j][I]/norm;
    }
  }

  /* Form restart vectors Ci = Sum_j(alpha[j][i]*Cj) */
  for (i=0; i<num; ++i) {
    sprintf(lbl, "%s %d", "CME", L+i);
    dpd_file2_init(&C1, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_scm(&C1, 0.0);
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "CME", j);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      dpd_file2_axpy(&CME, &C1, alpha[j][i], 0);
      dpd_file2_close(&CME);
    }
    dpd_file2_close(&C1);

    sprintf(lbl, "%s %d", "CMnEf", L+i);
    dpd_buf4_init(&C2, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&C2, 0.0);
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "CMnEf", j);
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_axpy(&CMnEf, &C2, alpha[j][i]);
      dpd_buf4_close(&CMnEf);
    }
    dpd_buf4_close(&C2);

    if (params.eom_ref > 0) {
      sprintf(lbl, "%s %d", "Cme", L+i);
      dpd_file2_init(&C1, EOM_Cme, C_irr, 0, 1, lbl);
      dpd_file2_scm(&C1, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "Cme", j);
        dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
        dpd_file2_axpy(&Cme, &C1, alpha[j][i], 0);
        dpd_file2_close(&Cme);
      }
      dpd_file2_close(&C1);

      sprintf(lbl, "%s %d", "CMNEF", L+i);
      dpd_buf4_init(&C2, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
      dpd_buf4_scm(&C2, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "CMNEF", j);
        dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
        dpd_buf4_axpy(&CMNEF, &C2, alpha[j][i]);
        dpd_buf4_close(&CMNEF);
      }
      dpd_buf4_close(&C2);

      sprintf(lbl, "%s %d", "Cmnef", L+i);
      dpd_buf4_init(&C2, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
      dpd_buf4_scm(&C2, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "Cmnef", j);
        dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
        dpd_buf4_axpy(&Cmnef, &C2, alpha[j][i]);
        dpd_buf4_close(&Cmnef);
      }
      dpd_buf4_close(&C2);
    }

    sprintf(lbl, "%s %d", "SIA", L+i);
    dpd_file2_init(&C1, EOM_SIA, C_irr, 0, 1, lbl);
    dpd_file2_scm(&C1, 0.0);
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "SIA", j);
      dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
      dpd_file2_axpy(&SIA, &C1, alpha[j][i], 0);
      dpd_file2_close(&SIA);
    }
    dpd_file2_close(&C1);

    sprintf(lbl, "%s %d", "SIjAb", L+i);
    dpd_buf4_init(&C2, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&C2, 0.0);
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "SIjAb", j);
      dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_axpy(&SIjAb, &C2, alpha[j][i]);
      dpd_buf4_close(&SIjAb);
    }
    dpd_buf4_close(&C2);

    if (params.eom_ref > 0) {
      sprintf(lbl, "%s %d", "Sia", L+i);
      dpd_file2_init(&C1, EOM_Sia, C_irr, 0, 1, lbl);
      dpd_file2_scm(&C1, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "Sia", j);
        dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
        dpd_file2_axpy(&Sia, &C1, alpha[j][i], 0);
        dpd_file2_close(&Sia);
      }
      dpd_file2_close(&C1);

      sprintf(lbl, "%s %d", "SIJAB", L+i);
      dpd_buf4_init(&C2, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
      dpd_buf4_scm(&C2, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "SIJAB", j);
        dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
        dpd_buf4_axpy(&SIJAB, &C2, alpha[j][i]);
        dpd_buf4_close(&SIJAB);
      }
      dpd_buf4_close(&C2);

      sprintf(lbl, "%s %d", "Sijab", L+i);
      dpd_buf4_init(&C2, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
      dpd_buf4_scm(&C2, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "Sijab", j);
        dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
        dpd_buf4_axpy(&Sijab, &C2, alpha[j][i]);
        dpd_buf4_close(&Sijab);
      }
      dpd_buf4_close(&C2);
    }

  }

  /* Copy restart vectors to beginning of file */
  for (i=0; i<num; ++i) {
    sprintf(lbl, "%s %d", "CME", L+i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_copy(&CME, EOM_CME, lbl);
    dpd_file2_close(&CME);
    sprintf(lbl, "%s %d", "CMnEf", L+i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_copy(&CMnEf, EOM_CMnEf, lbl);
    dpd_buf4_close(&CMnEf);

    if (params.eom_ref > 0) {
      sprintf(lbl, "%s %d", "Cme", L+i);
      dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "Cme", i);
      dpd_file2_copy(&Cme, EOM_Cme, lbl);
      dpd_file2_close(&Cme);
      sprintf(lbl, "%s %d", "CMNEF", L+i);
      dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "CMNEF", i);
      dpd_buf4_copy(&CMNEF, EOM_CMNEF, lbl);
      dpd_buf4_close(&CMNEF);
      sprintf(lbl, "%s %d", "Cmnef", L+i);
      dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "Cmnef", i);
      dpd_buf4_copy(&Cmnef, EOM_Cmnef, lbl);
      dpd_buf4_close(&Cmnef);
    }

    sprintf(lbl, "%s %d", "SIA", L+i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_copy(&SIA, EOM_SIA, lbl);
    dpd_file2_close(&SIA);
    sprintf(lbl, "%s %d", "SIjAb", L+i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_copy(&SIjAb, EOM_SIjAb, lbl);
    dpd_buf4_close(&SIjAb);
    
    if (params.eom_ref > 0) {
      sprintf(lbl, "%s %d", "Sia", L+i);
      dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "Sia", i);
      dpd_file2_copy(&Sia, EOM_Sia, lbl);
      dpd_file2_close(&Sia);
      sprintf(lbl, "%s %d", "SIJAB", L+i);
      dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "SIJAB", i);
      dpd_buf4_copy(&SIJAB, EOM_SIJAB, lbl);
      dpd_buf4_close(&SIJAB);
      sprintf(lbl, "%s %d", "Sijab", L+i);
      dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "Sijab", i);
      dpd_buf4_copy(&Sijab, EOM_Sijab, lbl); 
      dpd_buf4_close(&Sijab);
    }
  }

  return;
}

