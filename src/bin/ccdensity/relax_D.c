#include <dpd.h>
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

void relax_D(void)
{
  dpdfile2 D1, D2;
  
  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "DAI");
  dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_file2_axpy(&D2, &D1, 1.0, 1);
  dpd_file2_close(&D2);
  dpd_file2_close(&D1);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "DIA");
  dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  dpd_file2_axpy(&D2, &D1, 1.0, 1);
  dpd_file2_close(&D2);
  dpd_file2_close(&D1);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "Dai");
  dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  dpd_file2_axpy(&D2, &D1, 1.0, 1);
  dpd_file2_close(&D2);
  dpd_file2_close(&D1);

  dpd_file2_init(&D1, CC_OEI, 0, 0, 1, "Dia");
  dpd_file2_init(&D2, CC_OEI, 0, 1, 0, "D(orb)(a,i)");
  dpd_file2_axpy(&D2, &D1, 1.0, 1);
  dpd_file2_close(&D2);
  dpd_file2_close(&D1);
}
