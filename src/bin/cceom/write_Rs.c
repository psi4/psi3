/*
 *   write_Rs writes out all of the converged R's to RAMPS for the current irrep
 */

#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void write_Rs(int C_irr, double *evals, int *converged) {
  int i;
  dpdfile2 CME, Cme;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  static int R_index = 0;
  char C_lbl[32], R_lbl[32], E_lbl[32];

  for (i=0; i<eom_params.cs_per_irrep[C_irr]; ++i) {
    if (!converged[i]) continue; /* this root did not converge */
    ++R_index;

    sprintf(E_lbl, "EOM CCSD Energy for root %d", R_index);
    psio_write_entry(CC_INFO, E_lbl, (char *) &(evals[i]), sizeof(double));

    sprintf(C_lbl, "%s %d", "CME", i);
    sprintf(R_lbl, "%s %d", "RIA", R_index);

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, C_lbl);
    dpd_file2_copy(&CME, CC_RAMPS, R_lbl);
    dpd_file2_close(&CME);

    sprintf(C_lbl, "%s %d", "CMnEf", i);
    sprintf(R_lbl, "RIjAb %d", R_index);
    if (params.eom_ref <= 1)
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, C_lbl);
    else if (params.eom_ref == 2)
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, C_lbl);
    dpd_buf4_copy(&CMnEf, CC_RAMPS, R_lbl);
    dpd_buf4_close(&CMnEf);
          
    if(params.eom_ref > 0) {
      sprintf(C_lbl, "%s %d", "Cme", i);
      sprintf(R_lbl, "%s %d", "Ria", R_index);
      if (params.eom_ref == 1)
         dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, C_lbl);
      else if (params.eom_ref == 2)
         dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, C_lbl);
      dpd_file2_copy(&Cme, CC_RAMPS, R_lbl);
      dpd_file2_close(&Cme);

      sprintf(C_lbl, "%s %d", "CMNEF", i);
      sprintf(R_lbl, "%s %d", "RIJAB", R_index);

      dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, C_lbl);
      dpd_buf4_copy(&CMNEF, CC_RAMPS, R_lbl);
      dpd_buf4_close(&CMNEF);

      sprintf(C_lbl, "%s %d", "Cmnef", i);
      sprintf(R_lbl, "%s %d", "Rijab", R_index);

      if (params.eom_ref == 1)
        dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, C_lbl);
      else if (params.eom_ref ==2)
        dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, C_lbl);
        dpd_buf4_copy(&Cmnef, CC_RAMPS, R_lbl);
        dpd_buf4_close(&Cmnef);
    }
  }
}
