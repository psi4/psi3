#include <stdio.h>
#include <math.h>
#include <string.h>
#define EXTERN
#include "globals.h"

double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
  dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

void check_sum(char *term_lbl, int index, int irrep) {
  dpdfile2 Sia, SIA;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  static double old_norm=0;
  double norm;
  char lbl[80];

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

/*
  c_clean(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb);
*/

  norm = norm_C(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb); 

/* #ifdef EOM_DEBUG */
  fprintf(outfile,"Change in norm from %s term: %15.13lf\n",
    term_lbl, norm - old_norm);
  fflush(outfile);
/* #endif */

  old_norm = norm;

  dpd_file2_close(&SIA);
  dpd_file2_close(&Sia);
  dpd_buf4_close(&SIJAB);
  dpd_buf4_close(&Sijab);
  dpd_buf4_close(&SIjAb);
  return;
}
