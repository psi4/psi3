#include <dpd.h>
#define EXTERN
#include "globals.h"

/* Iij(): Build the occupied-occupied block of the orbital Lagrangian
** using the expression given in lag.c.  Note that we include an
** addition term here referred to as the reference contribution,
** 2*fij.  This comes from the general spin-orbital SCF gradient
** expression and is present for all reference types (though for
** canonical unperturbed orbitals, only the diagonal elements
** contribute, of course).  */

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

void Iij(void)
{
  int h, nirreps, i,j;
  dpdfile2 I, F, D;
  dpdbuf4 G, Aints, Fints, Dints, Cints, Eints;

  nirreps = moinfo.nirreps;
  /* I'IJ <-- sum_K fIK (DJK + DKJ) + sum_A fIA (DJA + DAJ) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_file2_init(&F, CC_FOCK, 0, 0, 0, "fIJ");
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);

  /* Add reference contribution: I'IJ <-- 2 fIJ */
  dpd_file2_axpy(&F, &I, 2.0, 0);
  dpd_file2_close(&F);

  dpd_file2_init(&F, CC_FOCK, 0, 0, 1, "fIA");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_close(&F);

  dpd_file2_close(&I);

  /* I'ij <-- sum_k fik (Djk + Dkj)  + sum_a fia (Dja + Daj) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_file2_init(&F, CC_FOCK, 0, 0, 0, "fij");
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 0.0);
  dpd_contract222(&F, &D, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);

  /* Add reference contribution: I'ij <-- 2 fij */
  dpd_file2_axpy(&F, &I, 2.0, 0);
  dpd_file2_close(&F);

  dpd_file2_init(&F, CC_FOCK, 0, 0, 1, "fia");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_close(&F);

  dpd_file2_close(&I);

  /* I'IJ <-- sum_KL <IK||JL> (D_KL + D_LK) + sum_kl <Ik|Jl> (D_kl + D_lk) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Aints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Aints);
  dpd_file2_close(&D);
  
  dpd_file2_close(&I);

  /* I'ij <-- sum_kl <ik||jl> (D_kl + D_lk) + sum_KL <iK|jL> (D_KL + D_LK) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Aints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Aints);
  dpd_file2_close(&D);
  
  dpd_file2_close(&I);

  /* I'IJ <-- sum_KA <IK||JA> (D_KA + D_AK) + sum_ka <Ik|Ja> (D_ka + D_ak) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);
  
  dpd_file2_close(&I);

  /* I'ij <-- sum_ka <ik||ja> (D_ka + D_ak) + sum_KA <iK|jA> (D_KA + D_AK) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);
  
  dpd_file2_close(&I);

  /* I'IJ <-- sum_AK <JK||IA> (D_AK + D_KA) + sum_ak <Jk|Ia> (D_ak + D_ka) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);
  
  dpd_file2_close(&I);

  /* I'ij <-- sum_ak <jk||ia> (D_ak + D_ka) + sum_AK <jK|iA> (D_AK + D_KA) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Eints);
  
  dpd_file2_close(&I);

  /* I'IJ <-- sum_AB <IA||JB> (D_AB + D_BA) + sum_ab <Ia|Jb> (D_ab + D_ba) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'ij <-- sum_ab <ia||jb> (D_ab + D_ba) + sum_AB <iA|jB> (D_AB + D_BA) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Cints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'IJ <-- sum_KLM <IK||LM> G(JK,LM) + 2 sum_kLm <Ik|Lm> G(Jk,Lm) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 2, 2, 2, 0, "GIJKL");
  dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Aints);

  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Aints);

  dpd_file2_close(&I);

  /* I'ij <-- sum_klm <ik||lm> G(jk,lm) + 2 sum_KlM <Ki|Ml> G(Kj,Ml) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 2, 0, 0, 1, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 2, 2, 2, 0, "Gijkl");
  dpd_contract442(&Aints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Aints);

  dpd_buf4_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  dpd_contract442(&Aints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Aints);

  dpd_file2_close(&I);

  /* I'IJ <-- sum_ABC <IA||BC> G(JA,BC) + 2 sum_AbC <aI|bC> G(aJ,bC) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_sort(&G, CC_TMP0, qprs, 10, 7, "GICAB");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 7, 10, 7, 0, "GICAB");
  dpd_contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GIcBa");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GIcBa");
  dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_file2_close(&I);

  /* I'ij <-- sum_abc <ia||bcC> G(ja,bc) + 2 sum_AbC <Ai|Bc> G(Aj,Bc) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_sort(&G, CC_TMP0, qprs, 10, 7, "Gicab");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 7, 10, 7, 0, "Gicab");
  dpd_contract442(&Fints, &G, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GiCbA");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GiCbA");
  dpd_contract442(&Fints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_file2_close(&I);

  /* I'IJ <-- sum_KAB <IK||AB> G(JK,AB) + 2 sum_kAb <Ik|Ab> G(Jk,Ab) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'ij <-- sum_kab <ik||ab> G(jk,ab) + 2 sum_KaB <iK|aB> G(jK,aB) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "Gijab");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_contract442(&Dints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'IJ <-- 2 sum_AKB <IA||KB> G(JA,KB) + 2 sum_aKb <Ia|Kb> G(Ja,Kb) -
	      2 sum_akB <Ik|Ba> GJakB(Jk,Ba) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Cints);

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Cints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
  dpd_buf4_sort(&G, CC_TMP0, prsq, 0, 5, "GIbjA (Ij,Ab)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GIbjA (Ij,Ab)");
  dpd_contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'ij <-- 2 sum_akb <ia||kb> G(ja,kb) + 2 sum_AkB <iA|kB> G(jA,kB) +
	      2 sum_AKb <iK|bA> GjAKb(jK,bA) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Cints);

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  dpd_contract442(&Cints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Cints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
  dpd_buf4_sort(&G, CC_TMP0, prsq, 0, 5, "GiBJa (iJ,aB)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GiBJa (iJ,aB)");
  dpd_contract442(&Dints, &G, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'IJ <-- 2 sum_KLA <IK||LA> G(JK,LA) + 2 sum_kLa <Ik|La> G(Jk,La)
              + 2 sum_kAl <kI|lA> G(kJ,lA) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
  dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'ij <-- 2 sum_kla <ik||la> G(jk,la) + 2 sum_KlA <iK|lA> G(jK,lA)
              + 2 sum_KaL <Ki|La> G(Kj,La) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
  dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_contract442(&Eints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'IJ <-- 2 sum_AKL <K>L||IA> G(K>L,JA) + 2 sum_aKl <Kl|Ia> G(Kl,Ja) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'IJ");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'ij <-- 2 sum_akl <k>l||ia> G(k>l,ja) + 2 sum_AkL <kL|iA> G(kL,jA) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 0, "I'ij");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);
}
