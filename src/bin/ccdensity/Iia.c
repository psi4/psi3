#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* Iia(): Build the occupied-virtual block of the orbital Lagrangian
** using the expression given in lag.c.
** */

#define CC_FOCK CC_OEI
/*
#define CC_OEI CC_OEI_NEW
#define CC_AINTS CC_AINTS_NEW
#define CC_BINTS CC_BINTS_NEW
#define CC_CINTS CC_CINTS_NEW
#define CC_DINTS CC_DINTS_NEW
#define CC_EINTS CC_EINTS_NEW
#define CC_FINTS CC_FINTS_NEW
*/

void Iia(void)
{
  int h, nirreps, a, i;
  dpdfile2 F, D, I;
  dpdbuf4 G, Aints, Fints, Eints, Dints, Cints;

  nirreps = moinfo.nirreps;

  /* I'IA <-- sum_J fIJ (DAJ + DJA) + sum_B fIB (DAB + DBA) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

  dpd_file2_init(&F, CC_FOCK, 0, 0, 0, "fIJ");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 0.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_close(&F);

  dpd_file2_init(&F, CC_FOCK, 0, 0, 1, "fIA");
  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_close(&F);

  dpd_file2_close(&I);

  /* I'ia <-- sum_j fij (Daj + Dja) + sum_b fib (Dab + Dba) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

  dpd_file2_init(&F, CC_FOCK, 0, 0, 0, "fij");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 0.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_close(&F);

  dpd_file2_init(&F, CC_FOCK, 0, 0, 1, "fia");
  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_close(&F);

  dpd_file2_close(&I);

  /* I'IA <-- sum_JKL <LK||JI> G(LK,JA) + 2 sum_jKl <lK|jI> G(lK,jA) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

  dpd_buf4_init(&Aints, CC_AINTS, 0, 2, 0, 0, 0, 1, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");
  dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Aints);

  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Aints);

  dpd_file2_close(&I);

  /* I'ia <-- sum_jkl <lk||ji> G(lk,ja) + 2 sum_JkL <Lk|Ji> G(Lk,Ja) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

  dpd_buf4_init(&Aints, CC_AINTS, 0, 2, 0, 0, 0, 1, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");
  dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Aints);

  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_contract442(&Aints, &G, &I, 3, 3, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Aints);

  dpd_file2_close(&I);

  /* I'IA <-- sum_BCD <IB||CD> G(AB,CD) + 2 sum_bCd <Ib|Cd> G(Ab,Cd) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 5, 7, 7, 7, 0, "GABCD");
  dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
  dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_file2_close(&I);

  /* I'ia <-- sum_bcd <ib||cd> G(ab,cd) + 2 sum_BcD <Dc|Bi> G(Dc,Ba) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 5, 7, 7, 7, 0, "Gabcd");
  dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_sort(&Fints, CC_TMP0, srqp, 5, 11, "F <cb|ai>");
  dpd_buf4_close(&Fints);
  dpd_buf4_init(&Fints, CC_TMP0, 0, 5, 11, 5, 11, 0, "F <cb|ai>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
  dpd_contract442(&Fints, &G, &I, 3, 3, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_file2_close(&I);

  /* I'IA <-- 2 sum_JKB <JI||KB> G(JA,KB) + 2 sum_jKb <Ij|Kb> G(Aj,Kb) +
              2 sum_jkB <jI|kB>> G(jA,kB) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
  dpd_buf4_sort(&G, CC_TMP0, qprs, 11, 10, "GiBJa (Bi,Ja)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 11, 10, 11, 10, 0, "GiBJa (Bi,Ja)");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_contract442(&Eints, &G, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&G);

  dpd_file2_close(&I);

  /* I'ia <-- 2 sum_jkb <ji||kb> G(ja,kb) + 2 sum_JkB <iJ|kB> G(aJ,kB) +
              2 sum_JKb <Ji|kB> G(Ja,Kb) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
  dpd_buf4_sort(&G, CC_TMP0, qprs, 11, 10, "GIbjA (bI,jA)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 11, 10, 11, 10, 0, "GIbjA (bI,jA)");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_contract442(&Eints, &G, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&G);

  dpd_file2_close(&I);

  /* I'IA <-- sum_BJK <JK||IB> G(JK,AB) + 2 sum_bJk <Jk|Ib> G(Jk,Ab) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "GIJAB");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'ia <-- sum_bjk <jk||ib> G(jk,ab) + 2 sum_BjK <Kj|Bi> G(Kj,Ba) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "Gijab");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 0, 5, "GjIbA");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GjIbA");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'IA <-- sum_JBC <IJ||BC> G(AJ,BC) + 2 sum_jBc <Ij|Bc> G(Aj,Bc) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'ia <-- sum_jbc <ij||bc> G(aj,bc) + 2 sum_JbC <iJ|bC> G(aJ,bC) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'IA <-- 2 sum_BJC <JC||IB> G(JC,AB) + 2 sum_bJc <Jc|Ib> G(Jc,Ab) +
              2 sum_bjC <Cj|Ib> G(Cj,Ab) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 7, 0, "GCIAB");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GICBA");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GICBA");
  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GIcBa");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GIcBa");
  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&Dints, CC_TMP0, spqr, 11, 10, "D <ij|ab> (bi,ja)");
  dpd_buf4_close(&Dints);
  dpd_buf4_init(&Dints, CC_TMP0, 0, 11, 10, 11, 10, 0, "D <ij|ab> (bi,ja)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
  dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'ia <-- 2 sum_bjc <jc||ib> G(jc,ab) + 2 sum_BjC <jC|iB> G(jC,aB) +
              2 sum_BJc <cJ|iB> G(cJ,aB) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 7, 0, "Gciab");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "Gicba");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "Gicba");
  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GiCbA");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GiCbA");
  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_contract442(&Cints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_buf4_close(&G);

  /* This set of sorted D-integrals is generated in the previous code block */
  dpd_buf4_init(&Dints, CC_TMP0, 0, 11, 10, 11, 10, 0, "D <ij|ab> (bi,ja)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
  dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);
}
