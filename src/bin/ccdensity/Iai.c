#include <dpd.h>
#define EXTERN
#include "globals.h"

/* Iai(): Build the virtual-occupied block of the orbital Lagrangian
** using the expression given in lag.c.  Note that we include an
** addition term here referred to as the reference contribution,
** 2*fai.  This comes from the general spin-orbital SCF gradient
** expression, but for unperturbed canonical Hartree-Fock orbitals
** (i.e., RHF and UHF only) this contribution is zero.  However, since
** the code to include the terms is trivial, we go ahead and do the
** work for all reference types. */

#define CC_FOCK CC_OEI
#define CC_OEI CC_OEI_NEW
#define CC_AINTS CC_AINTS_NEW
#define CC_BINTS CC_BINTS_NEW
#define CC_CINTS CC_CINTS_NEW
#define CC_DINTS CC_DINTS_NEW
#define CC_EINTS CC_EINTS_NEW
#define CC_FINTS CC_FINTS_NEW

void Iai(void)
{
  int h, nirreps, a, i;
  dpdfile2 F, D, I;
  dpdbuf4 G, Eints, Dints, Cints, Fints, Bints;

  nirreps = moinfo.nirreps;

  /* I'AI <-- sum_J fAJ (DIJ + DJI) + sum_B fAB (DIB + DBI) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_file2_init(&F, CC_FOCK, 0, 0, 1, "fIA");
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_contract222(&F, &D, &I, 1, 0, 1.0, 0.0);
  dpd_contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
  dpd_file2_close(&D);

  /* Add reference contribution: I'AI <-- 2 fAI */
  dpd_file2_axpy(&F, &I, 2.0, 1);
  dpd_file2_close(&F);

  dpd_file2_init(&F, CC_FOCK, 0, 1, 1, "fAB");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_close(&F);

  dpd_file2_close(&I);

  /* I'ai <-- sum_j faj (Dij + Dji) + sum_b fab (Dib + Dbi) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_file2_init(&F, CC_FOCK, 0, 0, 1, "fia");
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_contract222(&F, &D, &I, 1, 0, 1.0, 0.0);
  dpd_contract222(&F, &D, &I, 1, 1, 1.0, 1.0);
  dpd_file2_close(&D);

  /* Add reference contribution: I'ai <-- 2 fai */
  dpd_file2_axpy(&F, &I, 2.0, 1);
  dpd_file2_close(&F);

  dpd_file2_init(&F, CC_FOCK, 0, 1, 1, "fab");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_contract222(&F, &D, &I, 0, 0, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_close(&F);

  dpd_file2_close(&I);

  /* I'AI <-- sum_JK <AJ||IK> (D_JK + D_KJ) + sum_jk <Aj|Ik> (D_jk + D_kj) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'ai <-- sum_jk <aj||ik> (D_jk + D_kj) + sum_jk <aJ|iK> (D_JK + D_KJ) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'AI <-- - sum_JB <JA||IB> (D_JB + D_BJ) + sum_jb <Ij|Ab> (D_jb + D_bj) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Cints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'ai <-- - sum_jb <ja||ib> (D_jb + D_bj) + sum_JB <iJ|aB> (D_JB + D_BJ) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_dot14(&D, &Cints, &I, 0, 0, -1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Cints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'AI <-- sum_BJ <IJ||AB> (D_BJ + D_JB) + sum_bj <Ij|Ab> (D_bj + D_jb) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'ai <-- sum_bj <ij||ab> (D_bj + D_jb) + sum_BJ <iJ|aB> (D_BJ + D_JB) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_dot24(&D, &Dints, &I, 0, 1, 1.0, 1.0);
  dpd_file2_close(&D);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'AI <-- sum_BC <IC||AB> (D_BC + D_CB) + sum_bc <Ib|Ac>(D_bc + D_cb) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  dpd_dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
  dpd_dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
  dpd_dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'ai <-- sum_bc <ic||ab> (D_bc + D_cb) + sum_BC <iB|aC>(D_BC + D_CB) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  dpd_dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
  dpd_dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
  dpd_dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'AI <-- sum_JKL <AJ||KL> G(IJ,KL) + 2 sum_jKl <Aj|Kl> G(Ij,Kl) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 2, 2, 2, 0, "GIJKL");
  dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'ai <-- sum_jkl <aj||kl> G(ij,kl) + 2 sum_JkL <Lk|Ja> G(Lk,Ji) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 2, 2, 2, 0, "Gijkl");
  dpd_contract442(&Eints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  dpd_contract442(&Eints, &G, &I, 3, 3, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'AI <-- 2 sum_JKB <JA||KB> G(JI,KB) + 2 sum_jkB <jA|kB> G(jI,kB) +
              2 sum_jKb <Kj|Ab> (Aj,Kb) G(Ij,Kb) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "GIJKA");
  dpd_contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Cints);

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  dpd_contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Cints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_sort(&Dints, CC_TMP0, rqps, 11, 10, "D <ij|ab> (aj,ib)");
  dpd_buf4_close(&Dints);
  dpd_buf4_init(&Dints, CC_TMP0, 0, 11, 10, 11, 10, 0, "D <ij|ab> (aj,ib)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'ai <-- 2 sum_jkb <ja||kb> G(ji,kb) + 2 sum_JKb <jA|Kb> G(Ji,Kb) +
              2 sum_JkB <kJ|aB> (aJ,kB) G(iJ,kB) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 2, 10, 0, "Gijka");
  dpd_contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Cints);

  dpd_buf4_init(&Cints, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_contract442(&Cints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Cints);

  /* This sorted D-group is formed in the last code block */
  dpd_buf4_init(&Dints, CC_TMP0, 0, 11, 10, 11, 10, 0, "D <ij|ab> (aj,ib)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  dpd_contract442(&Dints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'AI <-- sum_BJK <JK||AB> G(JK,IB) + 2 sum_Jkb <Jk|Ab> G(Jk,Ib) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_buf4_init(&Dints, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");
  dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'ai <-- sum_bjk <jk||ab> G(jk,ib) + 2 sum_jKB <jK|aB> G(jK,iB) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_buf4_init(&Dints, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");
  dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  dpd_contract442(&Dints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Dints);

  dpd_file2_close(&I);

  /* I'AI <-- sum_JBC <JA||BC> G(JI,BC) + 2 sum_jbC <jA|bC> G(jI,bC) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
  dpd_contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 0, 5, "GiJaB");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GiJaB");
  dpd_contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_file2_close(&I);

  /* I'ai <-- sum_jbc <ja||bc> G(ji,bc) + 2 sum_JBc <Ja|Bc> G(Ji,Bc) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "Gijab");
  dpd_contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_file2_close(&I);

  /* I'AI <-- 2 sum_BJC <JC||AB> G(JC,IB) + 2 sum_bJc <Jc||Ab> G(Jc,Ib) +
              2 sum_bjC <jC|bA> GjCIb(jC,bI) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  dpd_contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
  dpd_buf4_sort(&G, CC_TMP0, rsqp, 10, 11, "GIbjA (jA,bI)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 11, 10, 11, 0, "GIbjA (jA,bI)");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_contract442(&Fints, &G, &I, 3, 3, -2.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&G);

  dpd_file2_close(&I);

  /* I'ai <-- 2 sum_bjc <jc||ab> G(jc,ib) + 2 sum_BjC <jC||aB> G(jC,iB) +
              2 sum_BJc <Jc|Ba> GJciB(Jc,Bi) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  dpd_contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_contract442(&Fints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
  dpd_buf4_sort(&G, CC_TMP0, rsqp, 10, 11, "GiBJa (Ja,Bi)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 11, 10, 11, 0, "GiBJa (Ja,Bi)");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_contract442(&Fints, &G, &I, 3, 3, -2.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&G);

  dpd_file2_close(&I);

  /* I'AI <-- sum_BCD <AB||CD> G(IB,CD) + 2 sum_bCd <Ab|Cd> G(Ib,Cd) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_sort(&G, CC_TMP0, qprs, 10, 7, "GICAB");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 7, 10, 7, 0, "GICAB");
  dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 7, 5, 5, 1, "B <ab|cd>");
  dpd_contract442(&Bints, &G, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&Bints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GIcAb");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GIcAb");
  dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  dpd_contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&Bints);
  dpd_buf4_close(&G);

  dpd_file2_close(&I);

  /* I'ai <-- sum_bcd <ab||cd> G(ib,cd) + 2 sum_BcD <aB|cD> G(iB,cD) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_sort(&G, CC_TMP0, qprs, 10, 7, "Gicab");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 7, 10, 7, 0, "Gicab");
  dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 7, 5, 5, 1, "B <ab|cd>");
  dpd_contract442(&Bints, &G, &I, 0, 0, -2.0, 1.0);
  dpd_buf4_close(&Bints);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 10, 5, "GiCaB");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 10, 5, 10, 5, 0, "GiCaB");
  dpd_buf4_init(&Bints, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  dpd_contract442(&Bints, &G, &I, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&Bints);
  dpd_buf4_close(&G);

  dpd_file2_close(&I);
  
}
