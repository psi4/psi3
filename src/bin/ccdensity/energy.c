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
  struct oe_dpdfile D, F;
  struct dpdbuf G, A, B, C, DInts, E, FInts;
  double one_energy=0.0, two_energy=0.0, total_two_energy = 0.0;

  fprintf(outfile, "\n\tEnergies re-computed from CC density:\n");
  fprintf(outfile,   "\t-------------------------------------\n");

  dpd_oe_file_init(&D, CC_OEI, 0, 0, "DIJ", 0, outfile);
  dpd_oe_file_init(&F, CC_OEI, 0, 0, "fIJ", 0, outfile);
  one_energy += dpd_oe_dot(&D, &F, 0, outfile);
  dpd_oe_file_close(&F);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 0, 0, "Dij", 0, outfile);
  dpd_oe_file_init(&F, CC_OEI, 0, 0, "fij", 0, outfile);
  one_energy += dpd_oe_dot(&D, &F, 0, outfile);
  dpd_oe_file_close(&F);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 1, 1, "DAB", 0, outfile);
  dpd_oe_file_init(&F, CC_OEI, 1, 1, "fAB", 0, outfile);
  one_energy += dpd_oe_dot(&D, &F, 0, outfile);
  dpd_oe_file_close(&F);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 1, 1, "Dab", 0, outfile);
  dpd_oe_file_init(&F, CC_OEI, 1, 1, "fab", 0, outfile);
  one_energy += dpd_oe_dot(&D, &F, 0, outfile);
  dpd_oe_file_close(&F);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fIA", 0, outfile);
  one_energy += dpd_oe_dot(&D, &F, 0, outfile);
  dpd_oe_file_close(&F);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fia", 0, outfile);
  one_energy += dpd_oe_dot(&D, &F, 0, outfile);
  dpd_oe_file_close(&F);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DAI", 0, outfile);
  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fIA", 0, outfile);
  one_energy += dpd_oe_dot(&D, &F, 0, outfile);
  dpd_oe_file_close(&F);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dai", 0, outfile);
  dpd_oe_file_init(&F, CC_OEI, 0, 1, "fia", 0, outfile);
  one_energy += dpd_oe_dot(&D, &F, 0, outfile);
  dpd_oe_file_close(&F);
  dpd_oe_file_close(&D);

  fprintf(outfile, "\tOne-electron energy        = %20.15f\n", one_energy);
  fflush(outfile);

  total_two_energy = 0.0;

  two_energy = 0.0;
  dpd_buf_init(&A, CC_AINTS, 2, 2, 0, 0, 1, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 2, 2, 2, 0, "GIJKL", 0, outfile);
  two_energy += dpd_dot(&G, &A, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 2, 2, 2, 2, 0, "Gijkl", 0, outfile);
  two_energy += dpd_dot(&G, &A, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&A);
  dpd_buf_init(&A, CC_AINTS, 0, 0, 0, 0, 0, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, "GIjKl", 0, outfile);
  two_energy += dpd_dot(&G, &A, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&A);

  total_two_energy += two_energy;
  fprintf(outfile, "\tIJKL energy                = %20.15f\n", two_energy);
  fflush(outfile);

  two_energy = 0.0;
  dpd_buf_init(&E, CC_EINTS, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)", 0,
	       outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "GIJKA", 0, outfile);
  two_energy += dpd_dot(&G, &E, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "Gijka", 0, outfile);
  two_energy += dpd_dot(&G, &E, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&E);
  dpd_buf_init(&E, CC_EINTS, 0, 10, 0, 10, 0, "E <ij|ka>", 0,
		outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  two_energy += dpd_dot(&G, &E, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  two_energy += dpd_dot(&G, &E, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&E);

  two_energy *= 2;
  total_two_energy += two_energy;
  fprintf(outfile, "\tIJKA energy                = %20.15f\n", two_energy);
  fflush(outfile);

  two_energy = 0.0;
  dpd_buf_init(&DInts, CC_DINTS, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)", 0,
	       outfile);
  dpd_buf_init(&G, CC_GAMMA, 2, 7, 2, 7, 0, "GIJAB", 0, outfile);
  two_energy += dpd_dot(&G, &DInts, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 2, 7, 2, 7, 0, "Gijab", 0, outfile);
  two_energy += dpd_dot(&G, &DInts, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&DInts);
  dpd_buf_init(&DInts, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0,
	       outfile);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  two_energy += dpd_dot(&G, &DInts, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&DInts);

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
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia||jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIBJA", 0, outfile);
  two_energy += dpd_dot(&G, &C, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "Gibja", 0, outfile);
  two_energy += dpd_dot(&G, &C, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&C);
  dpd_buf_init(&C, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIbJa", 0, outfile);
  two_energy += dpd_dot(&G, &C, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GiBjA", 0, outfile);
  two_energy += dpd_dot(&G, &C, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&C);
  dpd_buf_init(&DInts, CC_DINTS, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)", 0,
	       outfile);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIbjA", 0, outfile);
  two_energy -= dpd_dot(&G, &DInts, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GiBJa", 0, outfile);
  two_energy -= dpd_dot(&G, &DInts, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&DInts);

  total_two_energy += two_energy;
  fprintf(outfile, "\tIBJA energy                = %20.15f\n", two_energy);
  fflush(outfile);

  two_energy = 0.0;
  dpd_buf_init(&FInts, CC_FINTS, 10, 7, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_swap12(&FInts, CC_TMP0, 11, 7, "F(CI,AB)", 0, outfile);
  dpd_buf_close(&FInts);
  dpd_buf_init(&FInts, CC_TMP0, 11, 7, 11, 7, 0, "F(CI,AB)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 7, 11, 7, 0, "GCIAB", 0, outfile);
  two_energy -= dpd_dot(&G, &FInts, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 11, 7, 11, 7, 0, "Gciab", 0, outfile);
  two_energy -= dpd_dot(&G, &FInts, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&FInts);
  dpd_buf_init(&FInts, CC_FINTS, 10, 5, 10, 5, 0, "F <ia|bc>", 0, outfile);
  dpd_swap12(&FInts, CC_TMP0, 11, 5, "F(cI,Ba)", 0, outfile);
  dpd_buf_close(&FInts);
  dpd_buf_init(&FInts, CC_TMP0, 11, 5, 11, 5, 0, "F(cI,Ba)", 0, outfile);
  dpd_swap34(&FInts, CC_TMP1, 11, 5, "F(cI,aB)", 0, outfile);
  dpd_buf_close(&FInts);
  dpd_buf_init(&FInts, CC_TMP1, 11, 5, 11, 5, 0, "F(cI,aB)", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GcIaB", 0, outfile);
  two_energy += dpd_dot(&G, &FInts, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GCiAb", 0, outfile);
  two_energy += dpd_dot(&G, &FInts, 0, outfile); 
  dpd_buf_close(&G);
  dpd_buf_close(&FInts);

  two_energy *= 2;
  total_two_energy += two_energy;
  fprintf(outfile, "\tCIAB energy                = %20.15f\n", two_energy);
  fflush(outfile);

  two_energy = 0.0;
  dpd_buf_init(&B, CC_BINTS, 7, 7, 5, 5, 1, "B <ab|cd>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 7, 7, 7, 7, 0, "GABCD", 0, outfile);
  two_energy += dpd_dot(&G, &B, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 7, 7, 7, 7, 0, "Gabcd", 0, outfile);
  two_energy += dpd_dot(&G, &B, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&B);
  dpd_buf_init(&B, CC_BINTS, 5, 5, 5, 5, 0, "B <ab|cd>", 0, outfile);
  dpd_buf_init(&G, CC_GAMMA, 5, 5, 5, 5, 0, "GAbCd", 0, outfile);
  two_energy += dpd_dot(&G, &B, 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_close(&B);

  total_two_energy += two_energy;
  fprintf(outfile, "\tABCD energy                = %20.15f\n", two_energy);
  fprintf(outfile, "\tTotal two-electron energy  = %20.15f\n", total_two_energy);
  fprintf(outfile, "\tCCSD correlation energy    = %20.15f\n",
	  one_energy + total_two_energy);
  fprintf(outfile, "\tTotal CCSD energy          = %20.15f\n",
	  one_energy + total_two_energy + moinfo.eref);
}
