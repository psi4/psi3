#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar singles-doubles block contribution
to a Sigma vector stored at Sigma plus 'i' */

void FSD(int i, int irrep) {
  dpdfile2 Fme, FME;
  dpdfile2 SIA, Sia;
  dpdbuf4 CMNEF, Cmnef, CMnEf, CmNeF;
  char lbl[32];

  sprintf(lbl, "%s %d", "SIA", i);
  dpd_file2_init(&SIA, EOM_SIA, irrep, 0, 1, lbl);
  sprintf(lbl, "%s %d", "Sia", i);
  dpd_file2_init(&Sia, EOM_Sia, irrep, 0, 1, lbl);

 /* SIA += FME*CIMAE + Fme*CImAe */
  dpd_file2_init(&FME, CC_OEI, irrep, 0, 1, "FME");
  sprintf(lbl, "%s %d", "CMNEF", i);
  dpd_buf4_init(&CMNEF, EOM_CMNEF, irrep, 0, 5, 2, 7, 0, lbl);
  dpd_dot24(&FME,&CMNEF,&SIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&CMNEF);
  dpd_file2_close(&FME);

  dpd_file2_init(&Fme, CC_OEI, irrep, 0, 1, "Fme");
  sprintf(lbl, "%s %d", "CMnEf", i);
  dpd_buf4_init(&CMnEf, EOM_CMnEf, irrep, 0, 5, 0, 5, 0, lbl);
  dpd_dot24(&Fme,&CMnEf,&SIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&CMnEf);
  dpd_file2_close(&Fme);

 /* Sia += Fme*Cimae + FME*CiMaE */
  dpd_file2_init(&Fme, CC_OEI, irrep, 0, 1, "Fme");
  sprintf(lbl, "%s %d", "Cmnef", i);
  dpd_buf4_init(&Cmnef, EOM_Cmnef, irrep, 0, 5, 2, 7, 0, lbl);
  dpd_dot24(&Fme,&Cmnef,&Sia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Cmnef);
  dpd_file2_close(&Fme);

  dpd_file2_init(&FME, CC_OEI, irrep, 0, 1, "FME");
  dpd_buf4_init(&CmNeF, EOM_TMP, irrep, 0, 5, 0, 5, 0, "CmNeF");
  dpd_dot24(&FME,&CmNeF,&Sia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&CmNeF);
  dpd_file2_close(&Fme);

  dpd_file2_close(&SIA);
  dpd_file2_close(&Sia);

#ifdef EOM_DEBUG
  check_sum("FSD",i,irrep);
#endif

  return;
}
