#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

/* This function computes the H-bar singles-doubles block contribution
of Wamef to a Sigma vector stored at Sigma plus 'i' */

void WamefSD(int i, int C_irr) {
  dpdfile2 SIA, Sia;
  dpdbuf4 WAMEF, Wamef, WAmEf, WaMeF;
  dpdbuf4 CMNEF, Cmnef, CMnEf, CmNeF;
  char lbl[32];

  sprintf(lbl, "%s %d", "SIA", i);
  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
  sprintf(lbl, "%s %d", "Sia", i);
  dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);

 /* SIA += 0.5 CIMEF* WAMEF + CImEf * WAmEf */
  dpd_buf4_init(&WAMEF, CC_HBAR, H_IRR, 11, 7, 11, 7, 0, "WAMEF (AM,E>F)");
  sprintf(lbl, "%s %d", "CMNEF", i);
  dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 0, 7, 2, 7, 0, lbl);
  dpd_contract442(&CMNEF, &WAMEF, &SIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&CMNEF);
  dpd_buf4_close(&WAMEF);

  dpd_buf4_init(&WAmEf, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf (Am,Ef)");
  sprintf(lbl, "%s %d", "CMnEf", i);
  dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
  dpd_contract442(&CMnEf, &WAmEf, &SIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&CMnEf);
  dpd_buf4_close(&WAmEf);

 /* Sia += 0.5 Cimef * Wamef + CiMeF * WaMeF */
  dpd_buf4_init(&Wamef, CC_HBAR, H_IRR, 11, 7, 11, 7, 0, "Wamef (am,e>f)");
  sprintf(lbl, "%s %d", "Cmnef", i);
  dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 0, 7, 2, 7, 0, lbl);
  dpd_contract442(&Cmnef, &Wamef, &Sia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Cmnef);
  dpd_buf4_close(&Wamef);

  dpd_buf4_init(&WaMeF, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WaMeF (aM,eF)");
  dpd_buf4_init(&CmNeF, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CmNeF");
  dpd_contract442(&CmNeF, &WaMeF, &Sia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&CmNeF);
  dpd_buf4_close(&WaMeF);

  dpd_file2_close(&SIA);
  dpd_file2_close(&Sia);

#ifdef EOM_DEBUG
  check_sum("WamefSD",i,C_irr);
#endif

  return;
}
