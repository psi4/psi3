#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function orthogonalizes the r residual vector with the set of
numCs C vectors and adds the new vector to the C list if its norm is greater
than params.residual_tol */

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
  dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

extern void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
  dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);

void schmidt_add(dpdfile2 *RIA, dpdfile2 *Ria,
  dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int *numCs, int irrep)
{
   double dotval;
   double norm;
   int i, I;
   dpdfile2 Cme, CME, Cme2, CME2;
   dpdbuf4 CMNEF, Cmnef, CMnEf, CMNEF2, Cmnef2, CMnEf2;
   dpdbuf4 CMnEf_buf;
   char CME_lbl[32], Cme_lbl[32], CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];

   for (i=0; i<*numCs; i++) {
      sprintf(CME_lbl, "%s %d", "CME", i);
      sprintf(Cme_lbl, "%s %d", "Cme", i);
      sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
      sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
      sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);

      dpd_file2_init(&CME, EOM_CME, irrep, 0, 1, CME_lbl);
      dpd_file2_init(&Cme, EOM_Cme, irrep, 0, 1, Cme_lbl);
      dpd_buf4_init(&CMNEF, EOM_CMNEF, irrep, 2, 7, 2, 7, 0, CMNEF_lbl);
      dpd_buf4_init(&Cmnef, EOM_Cmnef, irrep, 2, 7, 2, 7, 0, Cmnef_lbl);
      dpd_buf4_init(&CMnEf, EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);

      dotval  = dpd_file2_dot(RIA, &CME);
      dotval += dpd_file2_dot(Ria, &Cme);
      dotval += dpd_buf4_dot(RIJAB, &CMNEF);
      dotval += dpd_buf4_dot(Rijab, &Cmnef);
      dotval += dpd_buf4_dot(RIjAb, &CMnEf);

#ifdef EOM_DEBUG
 fprintf(outfile,"dotval in schmidt_add: %15.10lf\n",dotval);
#endif

      dpd_file2_axpy(&CME, RIA, -1.0*dotval, 0);
      dpd_file2_axpy(&Cme, Ria, -1.0*dotval, 0);
      dpd_buf4_axpy(&CMNEF, RIJAB, -1.0*dotval);
      dpd_buf4_axpy(&Cmnef, Rijab, -1.0*dotval);
      dpd_buf4_axpy(&CMnEf, RIjAb, -1.0*dotval);

      dpd_file2_close(&CME);
      dpd_file2_close(&Cme);
      dpd_buf4_close(&CMNEF);
      dpd_buf4_close(&Cmnef);
      dpd_buf4_close(&CMnEf);
   }

   norm = norm_C(RIA, Ria, RIJAB, Rijab, RIjAb);

   if (norm < eom_params.schmidt_add_residual_tol) {
      return;
   }
   else {
      scm_C(RIA, Ria, RIJAB, Rijab, RIjAb, 1.0/norm);
      sprintf(CME_lbl, "%s %d", "CME", *numCs);
      sprintf(Cme_lbl, "%s %d", "Cme", *numCs);
      sprintf(CMNEF_lbl, "%s %d", "CMNEF", *numCs);
      sprintf(Cmnef_lbl, "%s %d", "Cmnef", *numCs);
      sprintf(CMnEf_lbl, "%s %d", "CMnEf", *numCs);

      dpd_file2_copy(RIA, EOM_CME, CME_lbl);
      dpd_file2_copy(Ria, EOM_Cme, Cme_lbl);
      dpd_buf4_copy(RIJAB, EOM_CMNEF, CMNEF_lbl);
      dpd_buf4_copy(Rijab, EOM_Cmnef, Cmnef_lbl);
      dpd_buf4_copy(RIjAb, EOM_CMnEf, CMnEf_lbl);

      ++(*numCs);
   }
   return;
}

