#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

/* ENERGY(): Compute the CC energy using the one- and two-particle
** density matrices.
**
** E = sum_pq Dpq fpq + 1/4 sum_pqrs Gpqrs <pq||rs>
**
** The individual two-electron components are:
**
** E(ijkl) = 1/4 sum_ijkl Gijkl <ij||kl>
**
** E(ijka) = 1/4 sum_ijka [Gijka <ij||ka> + Gijak <ij||ak> +
**                         Giajk <ia||jk> + Gaijk <ai||jk>]
**         = sum_ijka Gijka <ij||ka>
**
** E(ijab) = 1/4 sum_ijab [Gijab <ij||ab> + Gabij <ab||ij>]
**         = 1/2 sum_ijab Gijab <ij||ab>
**
** E(iajb) = 1/4 sum_iajb [Giajb <ia||jb> + Giabj <ia||bj> +
**                         Gaijb <ai||jb> + Gaibj <ai||bj>]
**         = sum_iajb Giajb <ia||jb>
**
** E(abci) = 1/4 sum_abci [Gabci <ab||ci> + Gabic <ab||ic> +
**                         Gciab <ic||ab> + Gicab <ic||ab>]
**         = sum_abci Gabci <ab||ci>
** E(abcd) = 1/4 sum_abcd Gabcd <ab||cd>
**
** Individual spin cases are handled below.
*/

void energy(void)
{
  dpdfile2 D, F;
  dpdbuf4 G, A, B, C, DInts, E, FInts;
  double one_energy=0.0, two_energy=0.0, total_two_energy = 0.0;

  fprintf(outfile, "\n\tEnergies re-computed from CC density:\n");
  fprintf(outfile,   "\t-------------------------------------\n");

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fIJ");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_file2_init(&F, CC_OEI, 0, 0, 0, "fij");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fAB");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_file2_init(&F, CC_OEI, 0, 1, 1, "fab");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fia");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fIA");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_file2_init(&F, CC_OEI, 0, 0, 1, "fia");
  one_energy += dpd_file2_dot(&D, &F);
  dpd_file2_close(&F);
  dpd_file2_close(&D);

  fprintf(outfile, "\tOne-electron energy        = %20.15f\n", one_energy);
  fflush(outfile);

  total_two_energy = 0.0;

  two_energy = 0.0;
  dpd_buf4_init(&A, CC_AINTS, 0, 2, 2, 0, 0, 1, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 2, 2, 2, 0, "GIJKL");
  two_energy += dpd_buf4_dot(&G, &A);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 2, 2, 2, 0, "Gijkl");
  two_energy += dpd_buf4_dot(&G, &A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&A);
  dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  two_energy += dpd_buf4_dot(&G, &A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&A);

  total_two_energy += two_energy;
  fprintf(outfile, "\tIJKL energy                = %20.15f\n", two_energy);
  fflush(outfile);

  two_energy = 0.0;
  dpd_buf4_init(&E, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "GIJKA");
  two_energy += dpd_buf4_dot(&G, &E);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 10, 2, 10, 0, "Gijka");
  two_energy += dpd_buf4_dot(&G, &E);
  dpd_buf4_close(&G);
  dpd_buf4_close(&E);
  dpd_buf4_init(&E, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  two_energy += dpd_buf4_dot(&G, &E);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GiJkA");
  two_energy += dpd_buf4_dot(&G, &E);
  dpd_buf4_close(&G);
  dpd_buf4_close(&E);

  two_energy *= 2;
  total_two_energy += two_energy;
  fprintf(outfile, "\tIJKA energy                = %20.15f\n", two_energy);
  fflush(outfile);

  two_energy = 0.0;
  dpd_buf4_init(&DInts, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 7, 2, 7, 0, "GIJAB");
  two_energy += dpd_buf4_dot(&G, &DInts);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 7, 2, 7, 0, "Gijab");
  two_energy += dpd_buf4_dot(&G, &DInts);
  dpd_buf4_close(&G);
  dpd_buf4_close(&DInts);
  dpd_buf4_init(&DInts, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  two_energy += dpd_buf4_dot(&G, &DInts);
  dpd_buf4_close(&G);
  dpd_buf4_close(&DInts);

  two_energy *= 2;
  total_two_energy += two_energy;
  fprintf(outfile, "\tIJAB energy                = %20.15f\n", two_energy);
  fflush(outfile);

  /*
  ** Compute the Gibja contribution to the two-electron energy.  By
  ** spin-case this contribution looks like:
  **
  **  E(AA) <-- sum_IBJA G(IB,JA) <JA||IB>
  **  E(BB) <-- sum_ibja G(ib,ja) <ja||ib>
  **  E(AB) <-- sum_IbJa ( G(Ib,Ja) <Ja|Ib> + G(iB,jA) <jA|iB> -
  **                         G(Ib,jA) <jA|bI> - G(iB,Ja) <Ja|Bi> )
  **
  **  See Gibja.c for the definition of G here.
  */
  two_energy = 0.0;
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  two_energy += dpd_buf4_dot(&G, &C);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "Gibja");
  two_energy += dpd_buf4_dot(&G, &C);
  dpd_buf4_close(&G);
  dpd_buf4_close(&C);
  dpd_buf4_init(&C, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
  two_energy += dpd_buf4_dot(&G, &C);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBjA");
  two_energy += dpd_buf4_dot(&G, &C);
  dpd_buf4_close(&G);
  dpd_buf4_close(&C);
  dpd_buf4_init(&DInts, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbjA");
  two_energy -= dpd_buf4_dot(&G, &DInts);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GiBJa");
  two_energy -= dpd_buf4_dot(&G, &DInts);
  dpd_buf4_close(&G);
  dpd_buf4_close(&DInts);

  total_two_energy += two_energy;
  fprintf(outfile, "\tIBJA energy                = %20.15f\n", two_energy);
  fflush(outfile);

  two_energy = 0.0;
  dpd_buf4_init(&FInts, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_sort(&FInts, CC_TMP0, qprs, 11, 7, "F(CI,AB)");
  dpd_buf4_close(&FInts);
  dpd_buf4_init(&FInts, CC_TMP0, 0, 11, 7, 11, 7, 0, "F(CI,AB)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "GCIAB");
  two_energy -= dpd_buf4_dot(&G, &FInts);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 7, 11, 7, 0, "Gciab");
  two_energy -= dpd_buf4_dot(&G, &FInts);
  dpd_buf4_close(&G);
  dpd_buf4_close(&FInts);
  dpd_buf4_init(&FInts, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_sort(&FInts, CC_TMP0, qprs, 11, 5, "F(cI,Ba)");
  dpd_buf4_close(&FInts);
  dpd_buf4_init(&FInts, CC_TMP0, 0, 11, 5, 11, 5, 0, "F(cI,Ba)");
  dpd_buf4_sort(&FInts, CC_TMP1, pqsr, 11, 5, "F(cI,aB)");
  dpd_buf4_close(&FInts);
  dpd_buf4_init(&FInts, CC_TMP1, 0, 11, 5, 11, 5, 0, "F(cI,aB)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GcIaB");
  two_energy += dpd_buf4_dot(&G, &FInts);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
  two_energy += dpd_buf4_dot(&G, &FInts); 
  dpd_buf4_close(&G);
  dpd_buf4_close(&FInts);

  two_energy *= 2;
  total_two_energy += two_energy;
  fprintf(outfile, "\tCIAB energy                = %20.15f\n", two_energy);
  fflush(outfile);

  two_energy = 0.0;
  dpd_buf4_init(&B, CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 7, 7, 7, 7, 0, "GABCD");
  two_energy += dpd_buf4_dot(&G, &B);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_GAMMA, 0, 7, 7, 7, 7, 0, "Gabcd");
  two_energy += dpd_buf4_dot(&G, &B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&B);
  dpd_buf4_init(&B, CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
  two_energy += dpd_buf4_dot(&G, &B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&B);

  total_two_energy += two_energy;
  fprintf(outfile, "\tABCD energy                = %20.15f\n", two_energy);
  fprintf(outfile, "\tTotal two-electron energy  = %20.15f\n", total_two_energy);
  fprintf(outfile, "\tCCSD correlation energy    = %20.15f\n",
	  one_energy + total_two_energy);
  fprintf(outfile, "\tTotal CCSD energy          = %20.15f\n",
	  one_energy + total_two_energy + moinfo.eref);
}
