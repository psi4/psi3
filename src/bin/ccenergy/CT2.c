#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

/* CT2(): Contributions of C-class integrals to T2.
**
** t(ij,ab) <--- P(ij) P(ab) t(m,a) t(i,e) <mb||je>
**
** This term is evaluated in two N^5 steps:
**  (1)  Y(mb,ji) = t(i,e) <mb||je>  (o^3 v^2)
**  (2)  t(ij,ab) <--- P(ij) P(ab) t(m,a) Y(mb,ji) (4 * o^3 v^2)
**
** Spin cases for UHF or ROHF orbitals:
** ------------------------------------
**                 *** AA ***
** + t(M,A) t(I,E) <MB||JE> - t(M,B) t(I,E) <MA||JE>
** - t(M,A) t(J,E) <MB||IE> + t(M,B) t(J,E) <MA||IE>
**
**                 *** BB ***
** + t(m,a) t(i,e) <mb||je> - t(m,b) t(i,e) <ma||je>
** - t(m,a) t(j,e) <me||ie> + t(m,b) t(j,e) <ma||ie>
**
**                 *** AB ***
** - t(M,A) t(I,E) <Mj|Eb> - t(m,b) t(I,E) <mA|jE>
** - t(M,A) t(j,e) <Mb|Ie> - t(m,b) t(j,e) <mI|eA>
**
** For the AA and BB spin cases, only the first term needs to be evaluated,
** while for the AB case, all four terms are different.
**
** This code was rewritten to eliminate all buf4_sort calls involving 
** mixing indices between bra and ket (i.e., a "complex" sort).  
** The current version requires six pairs of contractions (one each for AA 
** and BB and four for AB) and twelve simple sorts (four each for AA, BB, 
** and AB).
** TDC
** May 2000
*/

void CT2(void)
{
  dpdfile2 tIA, tia;
  dpdbuf4 Y, C, D, T2new, T2;

  dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

  /*** AA ***/

  /* C(MB||JE) * T(I,E) --> Y(MB,JI) */
  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (MB,JI)");
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
  dpd_buf4_close(&C);

  /* T(M,A) * Y(MB,JI) --> T(AB,JI) */
  dpd_buf4_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, 0, "T2 (AB,JI)");
  dpd_contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
  dpd_buf4_close(&Y);

  /* T(AB,JI) --> T(IJ,AB) */
  dpd_buf4_sort(&T2new, CC_TMP0, srpq, 0, 5, "T2 (IJ,AB)");

  /* P(IJ) P(AB) T2(IJ,AB) */
  dpd_buf4_init(&T2new, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (IJ,AB)");
  dpd_buf4_sort(&T2new, CC_TMP0, qprs, 0, 5, "T2 (JI,AB)");
  dpd_buf4_sort(&T2new, CC_TMP0, pqsr, 0, 5, "T2 (IJ,BA)");
  dpd_buf4_sort(&T2new, CC_TMP0, qpsr, 0, 5, "T2 (JI,BA)");
  
  /* T2(IJ,AB) - T2(JI,AB) - T2(IJ,BA) - T2(JI,BA) */
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (JI,AB)");
  dpd_buf4_axpy(&T2, &T2new, -1);
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (IJ,BA)");
  dpd_buf4_axpy(&T2, &T2new, -1);
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (JI,BA)");
  dpd_buf4_axpy(&T2, &T2new, 1);
  dpd_buf4_close(&T2);

  /* T2(IJ,AB) --> T2new (IJ,AB) */
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
  dpd_buf4_axpy(&T2new, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);
  

  /*** BB ***/

  /* C(mb||je) * T(i,e) --> Y(mb,ji) */
  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (MB,JI)");
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
  dpd_buf4_close(&C);

  /* T(m,a) * Y(mb,ji) --> T(ab,ji) */
  dpd_buf4_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, 0, "T2 (AB,JI)");
  dpd_contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
  dpd_buf4_close(&Y);

  /* T(ab,ji) --> T(ij,ab) */
  dpd_buf4_sort(&T2new, CC_TMP0, srpq, 0, 5, "T2 (IJ,AB)");

  /* P(ij) P(ab) T2(ij,ab) */
  dpd_buf4_init(&T2new, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (IJ,AB)");
  dpd_buf4_sort(&T2new, CC_TMP0, qprs, 0, 5, "T2 (JI,AB)");
  dpd_buf4_sort(&T2new, CC_TMP0, pqsr, 0, 5, "T2 (IJ,BA)");
  dpd_buf4_sort(&T2new, CC_TMP0, qpsr, 0, 5, "T2 (JI,BA)");
 
  /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) - T2(ji,ba) */
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (JI,AB)");
  dpd_buf4_axpy(&T2, &T2new, -1);
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (IJ,BA)");
  dpd_buf4_axpy(&T2, &T2new, -1);
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (JI,BA)");
  dpd_buf4_axpy(&T2, &T2new, 1);
  dpd_buf4_close(&T2);

  /* T2(ij,ab) --> T2new (ij,ab) */
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tijab");
  dpd_buf4_axpy(&T2new, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);


  /*** AB ***/

  /* C(mA|jE) * T(I,E) --> Y(mA,jI) */
  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (mA,jI)");
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
  dpd_buf4_close(&C);

  /* T(m,b) * Y(mA,jI) --> T2(bA,jI) */
  dpd_buf4_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, 0, "T2 (bA,jI)");
  dpd_contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
  dpd_buf4_close(&Y);

  /* T(bA,jI) --> Tnew(Ij,Ab) */
  dpd_buf4_sort(&T2new, CC_TMP0, srqp, 0, 5, "T2 (Ij,Ab)");
  dpd_buf4_close(&T2new);
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab)");
  dpd_buf4_init(&T2new, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
  dpd_buf4_axpy(&T2, &T2new, -1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);


  /* C(Mb|Ie) * T(j,e) --> Y(Mb,Ij) */
  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (Mb,Ij)");
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
  dpd_buf4_close(&C);

  /* T(M,A) * Y(Mb,Ij) --> T2(Ab,Ij) */
  dpd_buf4_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, 0, "T2 (Ab,Ij)");
  dpd_contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
  dpd_buf4_close(&Y);

  /* T(Ab,Ij) --> Tnew(Ij,Ab) */
  dpd_buf4_sort(&T2new, CC_TMP0, rspq, 0, 5, "T2 (Ij,Ab)");
  dpd_buf4_close(&T2new);
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab)");
  dpd_buf4_init(&T2new, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
  dpd_buf4_axpy(&T2, &T2new, -1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);


  /* D(Mb,jE) * T(I,E) --> Y(Mb,jI) */
  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(Mb,jI)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
  dpd_buf4_close(&D);

  /* T(M,A) * Y(Mb,jI) --> T2(Ab,jI) */
  dpd_buf4_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, 0, "T2 (Ab,jI)");
  dpd_contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
  dpd_buf4_close(&Y);
  
  /* T2(Ab,jI) --> Tnew(Ij,Ab) */
  dpd_buf4_sort(&T2new, CC_TMP0, srpq, 0, 5, "T2 (Ij,Ab)");
  dpd_buf4_close(&T2new);
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab)");
  dpd_buf4_init(&T2new, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
  dpd_buf4_axpy(&T2, &T2new, -1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);


  /* D(mA,Ie) * T(j,e) --> Y(mA,Ij) */
  dpd_buf4_init(&Y, CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(mA,Ij)");
  dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_contract424(&D, &tia, &Y, 3, 1, 0, 1, 0);
  dpd_buf4_close(&D);
 
  /* T(m,b) * Y(mA,Ij) --> T2(bA,Ij) */
  dpd_buf4_init(&T2new, CC_TMP0, 0, 5, 0, 5, 0, 0, "T2 (bA,Ij)");
  dpd_contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
  dpd_buf4_close(&Y);
  
  /* T2(bA,Ij) --> Tnew(Ij,Ab) */
  dpd_buf4_sort(&T2new, CC_TMP0, rsqp, 0, 5, "T2 (Ij,Ab)");
  dpd_buf4_close(&T2new);
  dpd_buf4_init(&T2, CC_TMP0, 0, 0, 5, 0, 5, 0, "T2 (Ij,Ab)");
  dpd_buf4_init(&T2new, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
  dpd_buf4_axpy(&T2, &T2new, -1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);

  dpd_file2_close(&tIA); dpd_file2_close(&tia);
}
