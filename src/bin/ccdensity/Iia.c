#include <dpd.h>
#define EXTERN
#include "globals.h"

/* Iia(): Build the occupied-virtual block of the orbital Lagrangian
** using the expression given in lag.c.
** */

void Iia(void)
{
  int h, nirreps, a, i;
  struct oe_dpdfile F, D, I;
  struct dpdbuf G, Aints, Fints, Eints, Dints, Cints;

  nirreps = moinfo.nirreps;

  /* I'IA <-- sum_J fIJ (DAJ + DJA) + sum_B fIB (DAB + DBA) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);

  dpd_oe_file_init(&F, CC_OEI, 0, 0, "fIJ", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DAI", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 1, 1, "DAB", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_close(&I);

  /* I'ia <-- sum_j fij (Daj + Dja) + sum_b fib (Dab + Dba) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);

  dpd_oe_file_init(&F, CC_OEI, 0, 0, "fij", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dai", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 1, 1, "Dab", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_close(&I);

  /* I'IA <-- sum_JKL <LK||JI> G(LK,JA) + 2 sum_jKl <lK|jI> G(lK,jA) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);

  dpd_buf_init(&Aints, CC_AINTS, 2, 0, 0, 0, 1, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "GIJKA", 0, outfile);
  dpd_contract122(&Aints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Aints);

  dpd_buf_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  dpd_contract122(&Aints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Aints);

  dpd_oe_file_close(&I);

  /* I'ia <-- sum_jkl <lk||ji> G(lk,ja) + 2 sum_JkL <Lk|Ji> G(Lk,Ja) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);

  dpd_buf_init(&Aints, CC_AINTS, 2, 0, 0, 0, 1, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "Gijka", 0, outfile);
  dpd_contract122(&Aints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Aints);

  dpd_buf_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  dpd_contract122(&Aints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Aints);

  dpd_oe_file_close(&I);

  /* I'IA <-- sum_BCD <IB||CD> G(AB,CD) + 2 sum_bCd <Ib|Cd> G(Ab,Cd) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);

  dpd_buf_init(&Fints, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 5, 7, 7, 7, 0, "GABCD", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_buf_init(&Fints, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 5, 5, 5, 5, 0, "GAbCd", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_oe_file_close(&I);

  /* I'ia <-- sum_bcd <ib||cd> G(ab,cd) + 2 sum_BcD <Dc|Bi> G(Dc,Ba) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);

  dpd_buf_init(&Fints, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 5, 7, 7, 7, 0, "Gabcd", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_buf_init(&Fints, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_sort(&Fints, CC_TMP0, srqp, 5, 11, "F <cb|ai>", 0, outfile);
  dpd_buf_close(&Fints);
  dpd_buf_init(&Fints, CC_TMP0, 5, 11, 5, 11, 0, "F <cb|ai>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 5, 5, 5, 5, 0, "GAbCd", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_oe_file_close(&I);

  /* I'IA <-- 2 sum_JKB <JI||KB> G(JA,KB) + 2 sum_jKb <Ij|Kb> G(Aj,Kb) +
              2 sum_jkB <jI|kB>> G(jA,kB) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIBJA", 0, outfile);
  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_contract122(&Eints, &G, &I, 1, 1, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&Eints);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GiBJa", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qprs, 11, 10, "GiBJa (Bi,Ja)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 11, 10, 11, 10, 0, "GiBJa (Bi,Ja)", 0, outfile);
  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&Eints);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GiBjA", 0, outfile);
  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 1, 1, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&Eints);
  dpd_buf_close(&G);

  dpd_oe_file_close(&I);

  /* I'ia <-- 2 sum_jkb <ji||kb> G(ja,kb) + 2 sum_JkB <iJ|kB> G(aJ,kB) +
              2 sum_JKb <Ji|kB> G(Ja,Kb) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "Gibja", 0, outfile);
  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_contract122(&Eints, &G, &I, 1, 1, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&Eints);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIbjA", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qprs, 11, 10, "GIbjA (bI,jA)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 11, 10, 11, 10, 0, "GIbjA (bI,jA)", 0, outfile);
  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&Eints);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIbJa", 0, outfile);
  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 1, 1, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&Eints);
  dpd_buf_close(&G);

  dpd_oe_file_close(&I);

  /* I'IA <-- sum_BJK <JK||IB> G(JK,AB) + 2 sum_bJk <Jk|Ib> G(Jk,Ab) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 5, 2, 7, 0, "GIJAB", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_oe_file_close(&I);

  /* I'ia <-- sum_bjk <jk||ib> G(jk,ab) + 2 sum_BjK <Kj|Bi> G(Kj,Ba) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 5, 2, 7, 0, "Gijab", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qpsr, 0, 5, "GjIbA", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 0, 5, 0, 5, 0, "GjIbA", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_oe_file_close(&I);

  /* I'IA <-- sum_JBC <IJ||BC> G(AJ,BC) + 2 sum_jBc <Ij|Bc> G(Aj,Bc) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);

  dpd_buf_init(&Dints, CC_DINTS, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 7, 11, 7, 0, "GCIAB", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GCiAb", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'ia <-- sum_jbc <ij||bc> G(aj,bc) + 2 sum_JbC <iJ|bC> G(aJ,bC) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);

  dpd_buf_init(&Dints, CC_DINTS, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 7, 11, 7, 0, "Gciab", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GcIaB", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'IA <-- 2 sum_BJC <JC||IB> G(JC,AB) + 2 sum_bJc <Jc|Ib> G(Jc,Ab) +
              2 sum_bjC <Cj|Ib> G(Cj,Ab) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 7, 0, "GCIAB", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qpsr, 10, 5, "GICBA", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 10, 5, 10, 5, 0, "GICBA", 0, outfile);
  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&Cints);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GcIaB", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qpsr, 10, 5, "GIcBa", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 10, 5, 10, 5, 0, "GIcBa", 0, outfile);
  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&Cints);
  dpd_buf_close(&G);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_sort(&Dints, CC_TMP0, spqr, 11, 10, "D <ij|ab> (bi,ja)", 0, outfile);
  dpd_buf_close(&Dints);
  dpd_buf_init(&Dints, CC_TMP0, 11, 10, 11, 10, 0, "D <ij|ab> (bi,ja)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GCiAb", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'ia <-- 2 sum_bjc <jc||ib> G(jc,ab) + 2 sum_BjC <jC|iB> G(jC,aB) +
              2 sum_BJc <cJ|iB> G(cJ,aB) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 7, 0, "Gciab", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qpsr, 10, 5, "Gicba", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 10, 5, 10, 5, 0, "Gicba", 0, outfile);
  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&Cints);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GCiAb", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qpsr, 10, 5, "GiCbA", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 10, 5, 10, 5, 0, "GiCbA", 0, outfile);
  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&Cints);
  dpd_buf_close(&G);

  /* This set of sorted D-integrals is generated in the previous code block */
  dpd_buf_init(&Dints, CC_TMP0, 11, 10, 11, 10, 0, "D <ij|ab> (bi,ja)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GcIaB", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);
}
