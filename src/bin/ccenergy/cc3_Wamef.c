#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* cc3_Wamef(): Compute the Wamef matrix from CC3 theory, which is
** given in spin-orbitals as:
** 
** Wamef = <am||ef> - t_n^a <nm||ef>
**
** TDC, Feb 2004
*/

void cc3_Wmnie(void)
{
  dpdbuf4 F, D, W;
  dpdfile2 t1;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&F, CC_EINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_sort(&F, CC_MISC, qpsr, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_close(&F);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&t1, &D, &W, 1, 1, 0, -1, 1);
    dpd_file2_close(&t1);
    dpd_buf4_close(&D);
  }
}
