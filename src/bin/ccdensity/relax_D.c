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
  struct oe_dpdfile D1, D2;
  
  dpd_oe_file_init(&D1, CC_OEI, 0, 1, "DAI", 0, outfile);
  dpd_oe_file_init(&D2, CC_OEI, 1, 0, "D(orb)(A,I)", 0 , outfile);
  dpd_oe_axpy(&D2, &D1, 1.0, 1, 0, outfile);
  dpd_oe_file_close(&D2);
  dpd_oe_file_close(&D1);

  dpd_oe_file_init(&D1, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_oe_file_init(&D2, CC_OEI, 1, 0, "D(orb)(A,I)", 0 , outfile);
  dpd_oe_axpy(&D2, &D1, 1.0, 1, 0, outfile);
  dpd_oe_file_close(&D2);
  dpd_oe_file_close(&D1);

  dpd_oe_file_init(&D1, CC_OEI, 0, 1, "Dai", 0, outfile);
  dpd_oe_file_init(&D2, CC_OEI, 1, 0, "D(orb)(a,i)", 0 , outfile);
  dpd_oe_axpy(&D2, &D1, 1.0, 1, 0, outfile);
  dpd_oe_file_close(&D2);
  dpd_oe_file_close(&D1);

  dpd_oe_file_init(&D1, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_oe_file_init(&D2, CC_OEI, 1, 0, "D(orb)(a,i)", 0 , outfile);
  dpd_oe_axpy(&D2, &D1, 1.0, 1, 0, outfile);
  dpd_oe_file_close(&D2);
  dpd_oe_file_close(&D1);
}
