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

void cc3_Wamef(void)
{
  dpdbuf4 F, D, W;
  dpdfile2 t1,tia,tIA;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_buf4_sort(&F, CC3_HET1, qpsr, 11, 5, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_close(&F);

    dpd_buf4_init(&W, CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract244(&t1, &D, &W, 0, 0, 0, -1, 1);
    dpd_file2_close(&t1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);
  }

  else if (params.ref == 1) { /** ROHF **/

    /** W(AM,E>F) <--- <AM||EF> **/
    /** W(am,e>f) <--- <am||ef> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 7, 11, 5, 1, "F <ai|bc>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WAMEF (AM,E>F)");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 Wamef (am,e>f)");
    dpd_buf4_close(&F);

    /** W(Am,Ef) <--- <Am|Ef> **/
    /** W(aM,eF) <--- <aM|eF> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WaMeF (aM,eF)");
    dpd_buf4_close(&F);

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /* t(N,A) <NM||EF> --> W(AM,E>F) */
    dpd_buf4_init(&W, CC3_HET1, 0, 11, 7, 11, 7, 0, "CC3 WAMEF (AM,E>F)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* t(n,a) <nm||ef> --> W(am,e>f) */
    dpd_buf4_init(&W, CC3_HET1, 0, 11, 7, 11, 7, 0, "CC3 Wamef (am,e>f)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* t(N,A) <Nm|Ef> --> W(Am,Ef) */
    dpd_buf4_init(&W, CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    /* t(n,a) <nM|eF> --> W(aM,eF) */
    dpd_buf4_init(&W, CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WaMeF (aM,eF)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1.0);
    dpd_buf4_close(&D);
    dpd_buf4_close(&W);

    dpd_file2_close(&tia);
    dpd_file2_close(&tIA);
  }
  
  else if (params.ref == 2) {

    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");

    /** W(AM,E>F) <--- <AM||EF> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WAMEF (AM,E>F)");
    dpd_buf4_close(&F);

    /** W(am,e>f) <--- <am||ef> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 Wamef (am,e>f)");
    dpd_buf4_close(&F);

    /** W(Am,Ef) <--- <Am|Ef> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_close(&F);

    /** W(aM,eF) <--- <aM|eF> **/
    dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    dpd_buf4_copy(&F, CC3_HET1, "CC3 WaMeF (aM,eF)");
    dpd_buf4_close(&F);

    /** W(AM,E>F) <--- tNA * <NM||EF> **/
    dpd_buf4_init(&W, CC3_HET1, 0, 21, 7, 21, 7, 0, "CC3 WAMEF (AM,E>F)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&D);

    /** W(am,e>f) <--- tna * <nm||ef> **/
    dpd_buf4_init(&W, CC3_HET1, 0, 31, 17, 31, 17, 0, "CC3 Wamef (am,e>f)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&D);

    /** W(Am,Ef) <--- tNA * <Nm|Ef> **/
    dpd_buf4_init(&W, CC3_HET1, 0, 26, 28, 26, 28, 0, "CC3 WAmEf (Am,Ef)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract244(&tIA, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&D);

    /** W(aM,eF) <--- tna * <nM|eF> **/
    dpd_buf4_init(&W, CC3_HET1, 0, 25, 29, 25, 29, 0, "CC3 WaMeF (aM,eF)");
    dpd_buf4_init(&D, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    dpd_contract244(&tia, &D, &W, 0, 0, 0, -1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&D);

    dpd_file2_close(&tia);
    dpd_file2_close(&tIA);
  }
}
