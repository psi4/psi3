#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* cc3_Wmbij(): Compute the Wmbij matrix from CC3 theory, which is
** given in spin-orbitals as:
**
** Wmbij = <mb||ij> + P(ij) t_i^e <mb||ej> - t_n^b Wmnij + t_i^e t_j^f <mb||ef>
**
** where the Wmnij intermediate is described in cc3_Wmnij.c.
**
** TDC, Feb 2004
*/

void cc3_Wmbij(void)
{
  dpdbuf4 C, D, E, F, W, W1, Z;
  dpdfile2 t1;

  if(params.ref == 0) { /** RHF **/

    /* W(Mb,Ij) <-- <Mb|Ij> */
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0,"E <ij|ka>");
    dpd_buf4_sort(&E, CC_MISC, rspq, 10, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_close(&E);

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    /* W(Mb,Ij) <-- - t(n,b) W(Mn,Ij) */
    dpd_buf4_init(&W, CC_MISC, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&W1, CC_MISC, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract424(&W1, &t1, &W, 1, 0, 1, -1, 1);
    dpd_buf4_close(&W1);
    dpd_buf4_close(&W);


    /* W(Mb,Ij) <-- + t(j,e) <Mb|Ie) */
    dpd_buf4_init(&W, CC_MISC, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_contract424(&C, &t1, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&C);
    dpd_buf4_close(&W);

    /* Z(Mb,Ej) = <Mb|Ej> + t(j,f) * <Mb|Ef> */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ib|aj> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "CC3 ZMbEj (Mb,Ej)");
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC3 ZMbEj (Mb,Ej)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &t1, &Z, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    /* W(Mb,Ij) <-- t(I,E) * Z(Mb,Ej) */
    dpd_buf4_init(&W, CC_MISC, 0, 10, 0, 10, 0, 0, "CC3 WMbIj (Mb,Ij)");
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC3 ZMbEj (Mb,Ej)");
    dpd_contract244(&t1, &Z, &W, 1, 2, 1, 1, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_sort(&W, CC_MISC, rspq, 0, 10, "CC3 WMbIj (Ij,Mb)");
    dpd_buf4_close(&W);

    dpd_file2_close(&t1);

  }
}
