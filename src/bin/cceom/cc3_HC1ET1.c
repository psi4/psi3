#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* compute [ [H,eT1] , C1 ]
 *
 * [H,eT^1] matrix elements are in CC3_HET1
 *
 * C1 has symmetry C_irr and is read from EOM_CME with name "CME i"
 *
 * illustrative code that computes [H,C1] is in cc3_HC1.c
 *
 * output matrix elements are written to CC_HC1ET1
 *
*/


void cc3_HC1ET1 (int i, int C_irr) {

  /* only need new Wmbij and Wabei for CC3 EOM energies */

  HC1ET1_Wmbij(i, C_irr);
  HC1ET1_Wabei(i, C_irr);

  return;
}

void HC1ET1_Wmbij(int i, int C_irr)
{
  dpdbuf4 C, D, E, F, W, W1, Z;
  dpdfile2 CME, Cme;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_close(&CME);

  }
  else if (params.ref == 1) { /** ROHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }
  else if (params.ref == 2) { /** UHF **/

    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }

  return;
}

void HC1ET1_Wabei(int i, int C_irr)
{
  dpdfile2 CME, Cme;
  dpdbuf4 Z, Z1, Z2, Z3, B, C, D, E, F, W;
  char CME_lbl[32], Cme_lbl[32];
  sprintf(CME_lbl, "%s %d", "CME", i);
  sprintf(Cme_lbl, "%s %d", "Cme", i);

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_close(&CME);
  }
  else if (params.ref == 1) { /** ROHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }
  else if (params.ref == 2) { /** UHF **/
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);

    dpd_file2_close(&CME);
    dpd_file2_close(&Cme);
  }

  return;
}
