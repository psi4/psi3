#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void DL2(int L_irr, int root_L_irr)
{
  dpdbuf4 D, Dold;

  if (params.ground) {
    /* RHS = <ij||ab> */
    if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
      dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
      dpd_buf4_copy(&D, CC_LAMBDA, "New LIJAB");
      dpd_buf4_copy(&D, CC_LAMBDA, "New Lijab");
      dpd_buf4_close(&D);
  
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_buf4_copy(&D, CC_LAMBDA, "New LIjAb");
      dpd_buf4_close(&D);
    }
    else if(params.ref == 2) { /** UHF **/
      dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
      dpd_buf4_copy(&D, CC_LAMBDA, "New LIJAB");
      dpd_buf4_close(&D);
  
      dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
      dpd_buf4_copy(&D, CC_LAMBDA, "New Lijab");
      dpd_buf4_close(&D);
  
      dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
      dpd_buf4_copy(&D, CC_LAMBDA, "New LIjAb");
      dpd_buf4_close(&D);
    }
  }
  else { /* excited state - no homogeneous term, first term is E*L */
    if (params.ref == 0 || params.ref == 1 ) { /** RHF/ROHF **/
      dpd_buf4_init(&D, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
      dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_axpy(&Dold, &D, -1.0 * params.cceom_energy[L_irr][root_L_irr]);
      dpd_buf4_close(&Dold);
      dpd_buf4_close(&D);
      dpd_buf4_init(&D, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
      dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
      dpd_buf4_axpy(&Dold, &D, -1.0 * params.cceom_energy[L_irr][root_L_irr]);
      dpd_buf4_close(&Dold);
      dpd_buf4_close(&D);
      dpd_buf4_init(&D, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
      dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      dpd_buf4_axpy(&Dold, &D, -1.0 * params.cceom_energy[L_irr][root_L_irr]);
      dpd_buf4_close(&Dold);
      dpd_buf4_close(&D);
    }
    else { /** UHF **/
      dpd_buf4_init(&D, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
      dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
      dpd_buf4_axpy(&Dold, &D, -1.0 * params.cceom_energy[L_irr][root_L_irr]);
      dpd_buf4_close(&Dold);
      dpd_buf4_close(&D);
      dpd_buf4_init(&D, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
      dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
      dpd_buf4_axpy(&Dold, &D, -1.0 * params.cceom_energy[L_irr][root_L_irr]);
      dpd_buf4_close(&Dold);
      dpd_buf4_close(&D);
      dpd_buf4_init(&D, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
      dpd_buf4_init(&Dold, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
      dpd_buf4_axpy(&Dold, &D, -1.0 * params.cceom_energy[L_irr][root_L_irr]);
      dpd_buf4_close(&Dold);
      dpd_buf4_close(&D);
    }
  }
}
