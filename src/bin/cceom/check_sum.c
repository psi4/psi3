#include <stdio.h>
#include <math.h>
#include <string.h>
#define EXTERN
#include "globals.h"

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);

void check_sum(char *term_lbl, int index, int irrep) {
  int save_params_ref;
  dpdfile2 Sia, SIA;
  dpdbuf4 SIJAB, Sijab, SIjAb, SIjbA;
  static double old_norm=0, old_dotval=0;
  double norm,dotval;
  char lbl[80];

  if (!strcmp(term_lbl,"reset"))  {
    fprintf(outfile,"resetting norm\n");
  	old_norm = 0;
  	return;
  }

  //save_params_ref = params.ref;
  //params.ref = 0;

  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "SIA", index);
    dpd_file2_init(&SIA, EOM_SIA, irrep, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIjAb", index);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, irrep, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_sort(&SIjAb, EOM_SIjAb, pqsr, 0, 5, "SIjbA"); 
    dpd_buf4_init(&SIjbA, EOM_SIjAb, irrep, 0, 5, 0, 5, 0, "SIjbA");

    norm = norm_C_rhf(&SIA, &SIjAb, &SIjbA);

    dpd_file2_close(&SIA);
    dpd_buf4_close(&SIjAb);
    dpd_buf4_close(&SIjbA);
  }
  else {
    sprintf(lbl, "%s %d", "SIA", index);
    dpd_file2_init(&SIA, EOM_SIA, irrep, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", index);
    dpd_file2_init(&Sia, EOM_Sia, irrep, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIJAB", index);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, irrep, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "Sijab", index);
    dpd_buf4_init(&Sijab, EOM_Sijab, irrep, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "SIjAb", index);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, irrep, 0, 5, 0, 5, 0, lbl);

    norm = norm_C(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb); 

    dpd_file2_close(&SIA);
    dpd_file2_close(&Sia);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&SIjAb);
  }

  fprintf(outfile,"%7s, D(norm sigma)=%15.10lf\n", term_lbl, norm - old_norm);
  fflush(outfile);
  old_norm = norm;

  //params.ref = save_params_ref;

  return;
}
