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
  dpdbuf4 A, E, D, Z, W, Z1, X;
  dpdfile2 t1, tIA, tia;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_copy(&A, CC3_HET1, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_close(&A);

    dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 0, 0, 0, 0, "CC3 ZMnIj (Mn,Ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_contract424(&E, &t1, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);
    dpd_buf4_init(&W, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&W);
    dpd_buf4_sort_axpy(&Z, CC3_HET1, qpsr, 0, 0, "CC3 WMnIj (Mn,Ij)", 1);
    dpd_buf4_close(&Z);

    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 10, 0, 10, 0, "CC3 ZMnIf (Mn,If)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract244(&t1, &D, &Z, 1, 2, 1, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&W, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract424(&Z, &t1, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    dpd_file2_close(&t1);
  }

  else if (params.ref == 1) {
    /** W(M>N,I>J) <--- <MN||IJ> **/
    /** W(m>n,i>j) <--- <mn||ij> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 2, 2, 0, 0, 1, "A <ij|kl>");
    dpd_buf4_copy(&A, CC3_HET1, "CC3 WMNIJ (M>N,I>J)");
    dpd_buf4_copy(&A, CC3_HET1, "CC3 Wmnij (m>n,i>j)");
    dpd_buf4_close(&A);

    /** W(Mn,Ij) <--- <Mn|Ij> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_copy(&A, CC3_HET1, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_close(&A);

    /** term 2 **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    /**** W(M>N,I>J) <-- ZMNIJ <-- P(I/J)( <MN||IE> * t1[J][E] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,IJ)");
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_contract424(&E, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 2, 0, "Z (M>N,JI)");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,JI)");
    dpd_buf4_axpy(&Z1, &Z, -1);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(m>n,i>j) <-- Zmnij <-- P(i/j)( <mn||ie> * t1[j][e] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (m>n,ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_contract424(&E, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 2, 0, "Z (m>n,ji)");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (m>n,ji)");
    dpd_buf4_axpy(&Z1, &Z, -1);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 Wmnij (m>n,i>j)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ie> * t1[j][e] ****/
    dpd_buf4_init(&W, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_contract424(&E, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ej> * t1[I][E] ****/
    dpd_buf4_init(&W, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 0, 11, 0, 11, 0, "E <ij|ak>");
    dpd_contract244(&tIA, &E, &W, 1, 2, 1, 1.0, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    /** term 3 **/

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 11, 2, 11, 0, "Z (M>N,EJ)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1.0, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 11, 2, 11, 0, "Z (m>n,ej)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1.0, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 Wmnij (m>n,i>j)");
    dpd_contract244(&tia, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    /**** W(Mn,Ij) <-- tIE tjf <Mn|Ef> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 0, 11, 0, 11, 0, "Z (Mn,Ej)");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1.0, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  }

  else if (params.ref == 2) {

    /** W(M>N,I>J) <--- <MN||IJ> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 2, 2, 0, 0, 1, "A <IJ|KL>");
    dpd_buf4_copy(&A, CC3_HET1, "CC3 WMNIJ (M>N,I>J)");
    dpd_buf4_close(&A);

    /** W(m>n,i>j) <--- <mn||ij> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 12, 12, 10, 10, 1, "A <ij|kl>");
    dpd_buf4_copy(&A, CC3_HET1, "CC3 Wmnij (m>n,i>j)");
    dpd_buf4_close(&A);

    /** W(Mn,Ij) <--- <Mn|Ij> **/
    dpd_buf4_init(&A, CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    dpd_buf4_copy(&A, CC3_HET1, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_close(&A);

    /** term 2 **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    /**** W(M>N,I>J) <-- ZMNIJ <-- P(I/J)( <MN||IE> * t1[J][E] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,IJ)");
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_contract424(&E, &tIA, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 2, 0, "Z (M>N,JI)");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 2, 0, 2, 0, 0, "Z (M>N,JI)");
    dpd_buf4_axpy(&Z1, &Z, -1);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(m>n,i>j) <-- Zmnij <-- P(i/j)( <mn||ie> * t1[j][e] ) ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 12, 10, 12, 10, 0, "Z (m>n,ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_contract424(&E, &tia, &Z, 3, 1, 0, 1, 0);
    dpd_buf4_close(&E);

    dpd_buf4_sort(&Z, CC_TMP0, pqsr, 12, 10, "Z (m>n,ji)");
    dpd_buf4_init(&Z1, CC_TMP0, 0, 12, 10, 12, 10, 0, "Z (m>n,ji)");
    dpd_buf4_axpy(&Z1, &Z, -1);
    dpd_buf4_close(&Z1);

    dpd_buf4_init(&W, CC3_HET1, 0, 12, 10, 12, 12, 0, "CC3 Wmnij (m>n,i>j)");
    dpd_buf4_axpy(&Z, &W, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ie> * t1[j][e] ****/
    dpd_buf4_init(&W, CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_contract424(&E, &tia, &W, 3, 1, 0, 1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    /**** W(Mn,Ij) <-- <Mn|Ej> * t1[I][E] ****/
    dpd_buf4_init(&W, CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_contract244(&tIA, &E, &W, 1, 2, 1, 1.0, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&W);

    /** term 3 **/

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 2, 21, 2, 21, 0, "Z (M>N,EJ)");
    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_contract424(&D, &tIA, &Z, 3, 1, 0, 1.0, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC3_HET1, 0, 2, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    /**** W(M>N,I>J) <-- tIE tJF <MN||EF> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 12, 31, 12, 31, 0, "Z (m>n,ej)");
    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1.0, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC3_HET1, 0, 12, 10, 12, 12, 0, "CC3 Wmnij (m>n,i>j)");
    dpd_contract244(&tia, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    /**** W(Mn,Ij) <-- tIE tjf <Mn|Ef> ****/
    dpd_buf4_init(&Z, CC_TMP0, 0, 22, 26, 22, 26, 0, "Z (Mn,Ej)");
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_contract424(&D, &tia, &Z, 3, 1, 0, 1.0, 0);
    dpd_buf4_close(&D);

    dpd_buf4_init(&W, CC3_HET1, 0, 22, 22, 22, 22, 0, "CC3 WMnIj (Mn,Ij)");
    dpd_contract244(&tIA, &Z, &W, 1, 2, 1, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
  }
}
