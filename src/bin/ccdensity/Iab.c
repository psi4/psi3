#include <dpd.h>
#define EXTERN
#include "globals.h"

/* Iab(): Build the virtual-virtual block of the orbital Lagrangian
** using the expression given in lag.c.
** */

void Iab(void)
{
  int h, nirreps, a, b;
  struct oe_dpdfile F, D, I;
  struct dpdbuf G, Bints, Cints, Dints, Eints, Fints;

  nirreps = moinfo.nirreps;

  /* I'AB <-- sum_I fAI (DBI + DIB) + sum_C fAC (DBC + DCB) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'AB", 0, outfile);

  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DAI", 0, outfile);
  dpd_contract111(&F, &D, &I, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_contract111(&F, &D, &I, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_init(&F, CC_OEI, 1, 1, "fAB", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 1, 1, "DAB", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_close(&I);

  /* I'ab <-- sum_i fai (Dbi + Dib) + sum_c fac (Dbc + Dcb) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'ab", 0, outfile);

  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dai", 0, outfile);
  dpd_contract111(&F, &D, &I, 1, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_contract111(&F, &D, &I, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_init(&F, CC_OEI, 1, 1, "fab", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 1, 1, "Dab", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_close(&I);

  /* I'AB <-- sum_JKI <JK||IA> G(JK,IB) + 2 sum_jKi <jK|iA> G(jK,iB) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'AB", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "GIJKA", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_oe_file_close(&I);

  /* I'ab <-- sum_jki <jk||ia> G(jk,ib) + 2 sum_JkI <Jk|Ia> G(Jk,Ib) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'ab", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "Gijka", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_oe_file_close(&I);

  /* I'AB <-- sum_CDE <AC||DE> G(BC,DE) + 2 sum_cDe <Ac|De> G(Bc,De) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'AB", 0, outfile);

  dpd_buf_init(&Bints, CC_BINTS, 5, 7, 5, 5, 1, "B <ab|cd>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 5, 7, 7, 7, 0, "GABCD", 0, outfile);
  dpd_contract122(&Bints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Bints);

  dpd_buf_init(&Bints, CC_BINTS, 5, 5, 5, 5, 0, "B <ab|cd>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 5, 5, 5, 5, 0, "GAbCd", 0, outfile);
  dpd_contract122(&Bints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Bints);

  dpd_oe_file_close(&I);

  /* I'ab <-- sum_cde <ac||de> G(bc,de) + 2 sum_CdE <Ed|Ca> G(Ed,Cb) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'ab", 0, outfile);

  dpd_buf_init(&Bints, CC_BINTS, 5, 7, 5, 5, 1, "B <ab|cd>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 5, 7, 7, 7, 0, "Gabcd", 0, outfile);
  dpd_contract122(&Bints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Bints);

  dpd_buf_init(&Bints, CC_BINTS, 5, 5, 5, 5, 0, "B <ab|cd>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 5, 5, 5, 5, 0, "GAbCd", 0, outfile);
  dpd_contract122(&Bints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Bints);

  dpd_oe_file_close(&I);

  /* I'AB <-- sum_ICD <AI||CD> G(BI,CD) + 2 sum_iCd <Ai|Cd> G(Bi,Cd) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'AB", 0, outfile);

  dpd_buf_init(&Fints, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_sort(&Fints, CC_TMP0, qprs, 11, 7, "F(CI,AB)", 0, outfile);
  dpd_buf_close(&Fints);
  dpd_buf_init(&Fints, CC_TMP0, 11, 7, 11, 7, 0, "F(CI,AB)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 7, 11, 7, 0, "GCIAB", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_buf_init(&Fints, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_sort(&Fints, CC_TMP0, qpsr, 11, 5, "F(Ai,Cd)", 0, outfile);
  dpd_buf_close(&Fints);
  dpd_buf_init(&Fints, CC_TMP0, 11, 5, 11, 5, 0, "F(Ai,Cd)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GCiAb", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_oe_file_close(&I);

  /* I'ab <-- sum_icd <ai||cd> G(bi,cd) + 2 sum_IcD <aI|cD> G(bI,cD) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'ab", 0, outfile);

  dpd_buf_init(&Fints, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_sort(&Fints, CC_TMP0, qprs, 11, 7, "F(ci,ab)", 0, outfile);
  dpd_buf_close(&Fints);
  dpd_buf_init(&Fints, CC_TMP0, 11, 7, 11, 7, 0, "F(ci,ab)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 7, 11, 7, 0, "Gciab", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_buf_init(&Fints, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_sort(&Fints, CC_TMP0, qpsr, 11, 5, "F(aI,cD)", 0, outfile);
  dpd_buf_close(&Fints);
  dpd_buf_init(&Fints, CC_TMP0, 11, 5, 11, 5, 0, "F(aI,cD)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GcIaB", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_oe_file_close(&I);

  /* I'AB <-- 2 sum_CDI <DI||CA> G(DI,CB) + 2 sum_cDi <Di|Ac> G(Di,Bc)
             + 2 sum_cdI <dI|cA> G(dI,cB) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'AB", 0, outfile);

  dpd_buf_init(&Fints, CC_FINTS, 10, 5, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_sort(&Fints, CC_TMP0, qprs, 11, 5, "F (DI,CA)", 0, outfile);
  dpd_buf_close(&Fints);
  dpd_buf_init(&Fints, CC_TMP0, 11, 5, 11, 5, 0, "F (DI,CA)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 7, 0, "GCIAB", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 3, 3, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_buf_init(&Fints, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_sort(&Fints, CC_TMP0, qpsr, 11, 5, "F (Di,Ac)", 0, outfile);
  dpd_buf_close(&Fints);
  dpd_buf_init(&Fints, CC_TMP0, 11, 5, 11, 5, 0, "F (Di,Ac)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GCiAb", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GcIaB", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_oe_file_close(&I);

  /* I'ab <-- 2 sum_cdi <di||ca> G(di,cb) + 2 sum_CdI <dI|aC> G(dI,bC)
             + 2 sum_CDi <Di|Ca> G(Di,Cb) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'ab", 0, outfile);

  /* Both sorted F-blocks used here were generated above */
  dpd_buf_init(&Fints, CC_TMP0, 11, 5, 11, 5, 0, "F (DI,CA)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 7, 0, "Gciab", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 3, 3, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_buf_init(&Fints, CC_TMP0, 11, 5, 11, 5, 0, "F (Di,Ac)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GcIaB", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GCiAb", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_oe_file_close(&I);

  /* I'AB <-- 2 sum_IJC <JC||IA> G(JC,IB) + 2 sum_jCi <jC|iA> G(jC,iB)
            - 2 sum_Jci <Jc|Ai> G(Jc,Bi) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'AB", 0, outfile);

  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIBJA", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Cints);

  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GiBjA", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Cints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIbjA", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, rpqs, 0, 5, "GIbjA (jI,bA)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 0, 5, 0, 5, 0, "GIbjA (jI,bA)", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 3, 3, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'ab <-- 2 sum_jci <jc||ia> G(jc,ib) + 2 sum_JcI <Jc|Ia> G(Jc,Ib)
            - 2 sum_jCI <Ij|Ca> GjCbI (IC,jb) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'ab", 0, outfile);

  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "Gibja", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Cints);

  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIbJa", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Cints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GiBJa", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, rpqs, 0, 5, "GiBJa (Ji,Ba)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 0, 5, 0, 5, 0, "GiBJa (Ji,Ba)", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 3, 3, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'AB <-- sum_CIJ <IJ||CA> G(IJ,CB) + 2 sum_Ijc <Ij|Ac> G(Ij,Bc) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'AB", 0, outfile);

  dpd_buf_init(&Dints, CC_DINTS, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 5, 2, 7, 0, "GIJAB", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'ab <-- sum_cij <ij||ca> G(ij,cb) + 2 sum_IjC <Ij|Ca> G(Ij,Cb) */
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'ab", 0, outfile);

  dpd_buf_init(&Dints, CC_DINTS, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 5, 2, 7, 0, "Gijab", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 3, 3, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);
}
