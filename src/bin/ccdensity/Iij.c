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

void Iij(void)
{
  int h, nirreps, i,j;
  struct oe_dpdfile I, F, D;
  struct dpdbuf G, Aints, Fints, Dints, Cints, Eints;

  nirreps = moinfo.nirreps;
  /* I'IJ <-- sum_K fIK (DJK + DKJ) + sum_A fIA (DJA + DAJ) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_oe_file_init(&F, CC_OEI, 0, 0, "fIJ", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 0, "DIJ", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);

  /* Add reference contribution: I'IJ <-- 2 fIJ */
  dpd_oe_axpy(&F, &I, 2.0, 0, 0, outfile);
  dpd_oe_file_close(&F);

  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DAI", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_close(&I);

  /* I'ij <-- sum_k fik (Djk + Dkj)  + sum_a fia (Dja + Daj) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_oe_file_init(&F, CC_OEI, 0, 0, "fij", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 0, "Dij", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);

  /* Add reference contribution: I'ij <-- 2 fij */
  dpd_oe_axpy(&F, &I, 2.0, 0, 0, outfile);
  dpd_oe_file_close(&F);

  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dai", 0, outfile);
  dpd_contract111(&F, &D, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&F);

  dpd_oe_file_close(&I);

  /* I'IJ <-- sum_KL <IK||JL> (D_KL + D_LK) + sum_kl <Ik|Jl> (D_kl + D_lk) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_oe_file_init(&D, CC_OEI, 0, 0, "DIJ", 0, outfile);
  dpd_buf_init(&Aints, CC_AINTS, 0, 0, 0, 0, 1, "A <ij|kl>", 0, outfile);
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Aints);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 0, 0, "Dij", 0, outfile);
  dpd_buf_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, "A <ij|kl>", 0, outfile);
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Aints);
  dpd_oe_file_close(&D);
  
  dpd_oe_file_close(&I);

  /* I'ij <-- sum_kl <ik||jl> (D_kl + D_lk) + sum_KL <iK|jL> (D_KL + D_LK) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_oe_file_init(&D, CC_OEI, 0, 0, "Dij", 0, outfile);
  dpd_buf_init(&Aints, CC_AINTS, 0, 0, 0, 0, 1, "A <ij|kl>", 0, outfile);
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Aints);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 0, 0, "DIJ", 0, outfile);
  dpd_buf_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, "A <ij|kl>", 0, outfile);
  dpd_dot24(&D, &Aints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot24(&D, &Aints, &I, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Aints);
  dpd_oe_file_close(&D);
  
  dpd_oe_file_close(&I);

  /* I'IJ <-- sum_KA <IK||JA> (D_KA + D_AK) + sum_ka <Ik|Ja> (D_ka + D_ak) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DAI", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dai", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_buf_close(&Eints);
  
  dpd_oe_file_close(&I);

  /* I'ij <-- sum_ka <ik||ja> (D_ka + D_ak) + sum_KA <iK|jA> (D_KA + D_AK) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dai", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DAI", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_buf_close(&Eints);
  
  dpd_oe_file_close(&I);

  /* I'IJ <-- sum_AK <JK||IA> (D_AK + D_KA) + sum_ak <Jk|Ia> (D_ak + D_ka) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DAI", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dai", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_buf_close(&Eints);
  
  dpd_oe_file_close(&I);

  /* I'ij <-- sum_ak <jk||ia> (D_ak + D_ka) + sum_AK <jK|iA> (D_AK + D_KA) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dai", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DAI", 0, outfile);
  dpd_dot24(&D, &Eints, &I, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&D);
  dpd_buf_close(&Eints);
  
  dpd_oe_file_close(&I);

  /* I'IJ <-- sum_AB <IA||JB> (D_AB + D_BA) + sum_ab <Ia|Jb> (D_ab + D_ba) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_oe_file_init(&D, CC_OEI, 1, 1, "DAB", 0, outfile);
  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Cints);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 1, 1, "Dab", 0, outfile);
  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Cints);
  dpd_oe_file_close(&D);

  dpd_oe_file_close(&I);

  /* I'ij <-- sum_ab <ia||jb> (D_ab + D_ba) + sum_AB <iA|jB> (D_AB + D_BA) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_oe_file_init(&D, CC_OEI, 1, 1, "Dab", 0, outfile);
  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Cints);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 1, 1, "DAB", 0, outfile);
  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_dot24(&D, &Cints, &I, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot24(&D, &Cints, &I, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Cints);
  dpd_oe_file_close(&D);

  dpd_oe_file_close(&I);

  /* I'IJ <-- sum_KLM <IK||LM> G(JK,LM) + 2 sum_kLm <Ik|Lm> G(Jk,Lm) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_buf_init(&Aints, CC_AINTS, 0, 2, 0, 0, 1, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 2, 2, 2, 0, "GIJKL", 0, outfile);
  dpd_contract122(&Aints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Aints);

  dpd_buf_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, "GIjKl", 0, outfile);
  dpd_contract122(&Aints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Aints);

  dpd_oe_file_close(&I);

  /* I'ij <-- sum_klm <ik||lm> G(jk,lm) + 2 sum_KlM <Ki|Ml> G(Kj,Ml) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_buf_init(&Aints, CC_AINTS, 0, 2, 0, 0, 1, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 2, 2, 2, 0, "Gijkl", 0, outfile);
  dpd_contract122(&Aints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Aints);

  dpd_buf_init(&Aints, CC_AINTS, 0, 0, 0, 0, 0, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, "GIjKl", 0, outfile);
  dpd_contract122(&Aints, &G, &I, 1, 1, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Aints);

  dpd_oe_file_close(&I);

  /* I'IJ <-- sum_ABC <IA||BC> G(JA,BC) + 2 sum_AbC <aI|bC> G(aJ,bC) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_buf_init(&Fints, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 7, 11, 7, 0, "GCIAB", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qprs, 10, 7, "GICAB", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 10, 7, 10, 7, 0, "GICAB", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_buf_init(&Fints, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GcIaB", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qpsr, 10, 5, "GIcBa", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 10, 5, 10, 5, 0, "GIcBa", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_oe_file_close(&I);

  /* I'ij <-- sum_abc <ia||bcC> G(ja,bc) + 2 sum_AbC <Ai|Bc> G(Aj,Bc) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_buf_init(&Fints, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 7, 11, 7, 0, "Gciab", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qprs, 10, 7, "Gicab", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 10, 7, 10, 7, 0, "Gicab", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_buf_init(&Fints, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GCiAb", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, qpsr, 10, 5, "GiCbA", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 10, 5, 10, 5, 0, "GiCbA", 0, outfile);
  dpd_contract122(&Fints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Fints);

  dpd_oe_file_close(&I);

  /* I'IJ <-- sum_KAB <IK||AB> G(JK,AB) + 2 sum_kAb <Ik|Ab> G(Jk,Ab) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_buf_init(&Dints, CC_DINTS, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 7, 2, 7, 0, "GIJAB", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'ij <-- sum_kab <ik||ab> G(jk,ab) + 2 sum_KaB <iK|aB> G(jK,aB) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_buf_init(&Dints, CC_DINTS, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 7, 2, 7, 0, "Gijab", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 1, 1, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'IJ <-- 2 sum_AKB <IA||KB> G(JA,KB) + 2 sum_aKb <Ia|Kb> G(Ja,Kb) -
	      2 sum_akB <Ik|Ba> GJakB(Jk,Ba) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIBJA", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Cints);

  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIbJa", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Cints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIbjA", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, prsq, 0, 5, "GIbjA (Ij,Ab)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 0, 5, 0, 5, 0, "GIbjA (Ij,Ab)", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'ij <-- 2 sum_akb <ia||kb> G(ja,kb) + 2 sum_AkB <iA|kB> G(jA,kB) +
	      2 sum_AKb <iK|bA> GjAKb(jK,bA) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "Gibja", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Cints);

  dpd_buf_init(&Cints, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GiBjA", 0, outfile);
  dpd_contract122(&Cints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Cints);

  dpd_buf_init(&Dints, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GiBJa", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, prsq, 0, 5, "GiBJa (iJ,aB)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 0, 5, 0, 5, 0, "GiBJa (iJ,aB)", 0, outfile);
  dpd_contract122(&Dints, &G, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Dints);

  dpd_oe_file_close(&I);

  /* I'IJ <-- 2 sum_KLA <IK||LA> G(JK,LA) + 2 sum_kLa <Ik|La> G(Jk,La)
              + 2 sum_kAl <kI|lA> G(kJ,lA) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 2, 10, 0, "GIJKA", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 1, 1, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_oe_file_close(&I);

  /* I'ij <-- 2 sum_kla <ik||la> G(jk,la) + 2 sum_KlA <iK|lA> G(jK,lA)
              + 2 sum_KaL <Ki|La> G(Kj,La) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 2, 10, 0, "Gijka", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 0, 0, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 1, 1, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_oe_file_close(&I);

  /* I'IJ <-- 2 sum_AKL <K>L||IA> G(K>L,JA) + 2 sum_aKl <Kl|Ia> G(Kl,Ja) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "GIJKA", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_oe_file_close(&I);

  /* I'ij <-- 2 sum_akl <k>l||ia> G(k>l,ja) + 2 sum_AkL <kL|iA> G(kL,jA) */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);

  dpd_buf_init(&Eints, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)",
	       0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "Gijka", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_buf_init(&Eints, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  dpd_contract122(&Eints, &G, &I, 2, 2, 2.0, 1.0, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&Eints);

  dpd_oe_file_close(&I);
}
