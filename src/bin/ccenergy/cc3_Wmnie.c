#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* cc3_Wmnie(): Compute the Wmnie matrix from CC3 theory, which is
** given in spin-orbitals as:
** 
** Wmnie = <mn||ie> + t_i^f <mn||fe>
**
** TDC, Feb 2004
*/

void cc3_Wmnie(void)
{
  dpdbuf4 E, D, W;
  dpdfile2 t1;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_buf4_copy(&E, CC_MISC, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_close(&E);

    dpd_buf4_init(&W, CC_MISC, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&t1, &D, &W, 1, 2, 1, 1, 1);
    dpd_file2_close(&t1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);
  }
}
