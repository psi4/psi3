#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* cc3_Wmnij(): Compute the Wmnij matrix from CC3 theory, which is
** given in spin-orbitals as:
**
** Wmnij = <mn||ij> + P(ij) t_j^e <mn||ie> + t_i^e t_j^f <mn||ef>
**
** TDC, Feb 2004
*/

void cc3_Wmnij(void)
{
  dpdbuf4 A, E, D, Z, W;
  dpdfile2 t1;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_copy(&A, CC_MISC, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_close(&A);

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 0, 0, 0, 0, "CC3 ZMnIj (Mn,Ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_contract424(&E, &t1, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);
    dpd_buf4_init(&W, CC_MISC, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC_MISC, qpsr, 0, 0, "CC3 WMnIj (Mn,Ij)", 1);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 10, 0, 10, 0, "CC3 ZMnIf (Mn,If)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&t1, &D, &Z, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC_MISC, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract424(&Z, &t1, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    dpd_file2_close(&t1);
  }
}
