#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* cc3_Wabei(): Compute the Wabei matrix from CC3 theory, which is
** given in spin orbitals as:
**
** Wabei = <ab||ei> - P(ab) t_m^a <mb||ei> + t_i^f <ab||ef> 
**         - P(ab) t_i^f t_m^b <am||ef> + t_m^a t_n^b <mn||ei>
**         + t_m^a t_i^f t_n^b <mn||ef>
**
** The basic strategy for this code is to generate two intermediate
** quantities, Z1(Ab,EI) and Z2(Ei,Ab), which are summed in the final
** step to give the complete W(Ei,Ab) intermediate.  This is sorted
** to W(iE,Ab) storage for use in the triples equations.
**
** TDC, Feb 2004
*/

void cc3_Wabei(void)
{
  dpdfile2 t1;
  dpdbuf4 Z, Z1, Z2, Z3;
  dpdbuf4 B, C, D, E, F, W;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_buf4_copy(&F, CC_TMP0, "CC3 Z(Ei,Ab)");
    dpd_buf4_close(&F);

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "CC3 Z(Ab,Ei)");
    dpd_buf4_init(&Z2, CC_TMP0, 0, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");

    /* Z1(Ab,Ei) <-- <Ab|Ef> * t(i,f) */
    dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    dpd_contract424(&B, &t1, &Z1, 3, 1, 0, 1, 0);
    dpd_buf4_close(&B);

    /* Z(Mb,Ei) <-- <Mb|Ei> + <Mb|Ef> t(i,f) */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ib|aj> (ib,aj)");
    dpd_buf4_copy(&D, CC_TMP0, "CC3 Z(Mb,Ei)");
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, CC_TMP0, 0, 10, 11, 10, 11, 0, "CC3 Z(Mb,Ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &t1, &Z, 3, 1, 0, 1, 1);
    dpd_buf4_close(&F);
    /* Z1(Ab,Ei) <-- - t(M,A) * Z(Mb,Ei) */
    dpd_contract244(&t1, &Z, &Z1, 0, 0, 0, -1, 1);
    dpd_buf4_close(&Z);

    /* Z(Ei,Am) <-- <Ei|Am> + <Am|Ef> t(i,f) */
    dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_buf4_sort(&C, CC_TMP0, qpsr, 11, 11, "CC3 Z(Ei,Am)");
    dpd_buf4_close(&C);
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "CC3 Z(Am,Ei)");
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_contract424(&F, &t1, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&F);
    dpd_buf4_sort_axpy(&Z, CC_TMP0, rspq, 11, 11, "CC3 Z(Ei,Am)", 1);
    dpd_buf4_close(&Z);
    /* Z2(Ei,Ab) <-- - Z(Ei,Am) t(m,b) */
    dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "CC3 Z(Ei,Am)");
    dpd_contract424(&Z, &t1, &Z2, 3, 0, 0, -1, 1);
    dpd_buf4_close(&Z);

    /* Z(Mn,Ei) = <Mn|Ei> + <Mn|Ef> t_i^f */
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    dpd_buf4_copy(&E, CC_TMP0, "CC3 Z(Mn,Ei)");
    dpd_buf4_close(&E);
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "CC3 Z(Mn,Ei)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &t1, &Z, 3, 1, 0, 1, 1);
    dpd_buf4_close(&D);
    /* Z'(An,Ei) = t_M^A Z(Mn,Ei) */
    dpd_buf4_init(&Z3, CC_TMP0, 0, 11, 11, 11, 11, 0, "CC3 Z(An,Ei)");
    dpd_contract244(&t1, &Z, &Z3, 0, 0, 0, 1, 0);
    dpd_buf4_close(&Z);
    /* Z2(Ei,Ab) <-- Z'(An,Ei) t_n^b */
    dpd_contract424(&Z3, &t1, &Z2, 1, 0, 0, 1, 1);
    dpd_buf4_close(&Z3);

    dpd_buf4_close(&Z2);
    dpd_buf4_close(&Z1);

    /* W(Ab,Ei) = Z1(Ab,Ei) + Z2(Ei,Ab) */
    dpd_buf4_init(&Z1, CC_TMP0, 0, 5, 11, 5, 11, 0, "CC3 Z(Ab,Ei)");
    dpd_buf4_sort_axpy(&Z1, CC_TMP0, rspq, 11, 5, "CC3 Z(Ei,Ab)", 1);
    dpd_buf4_close(&Z1);
    dpd_buf4_init(&Z2, CC_TMP0, 0, 11, 5, 11, 5, 0, "CC3 Z(Ei,Ab)");
    dpd_buf4_sort(&Z2, CC_MISC, qpsr, 10, 5, "CC3 WAbEi (Ie,Ab)");
    dpd_buf4_close(&Z2);

    dpd_file2_close(&t1);
  }
}
