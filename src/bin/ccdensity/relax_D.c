#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* relax_D(): Add the orbital-response contributions to the Dia block
** of the one-electron density matrix:
**
** D(A,I) = D(amp)(A,I) + D(orb)(A,I)
**
** D(I,A) = D(amp)(I,A) + D(orb)(A,I)
**
** */

void relax_D(struct RHO_Params rho_params)
{
  dpdfile2 D1, D2;

  if(params.ref == 0 || params.ref == 1) {   
    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);
  }
  else if(params.ref == 2) {

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
    dpd_file2_init(&D2, CC_OEI, 0, 3, 2, "D(orb)(a,i)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);

    dpd_file2_init(&D1, CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
    dpd_file2_init(&D2, CC_OEI, 0, 3, 2, "D(orb)(a,i)");
    dpd_file2_axpy(&D2, &D1, 1.0, 1);
    dpd_file2_close(&D2);
    dpd_file2_close(&D1);


  }
}
