#include <dpd.h>
#define EXTERN
#include "globals.h"

/* relax_I(): Add the orbital-response contributions from the
** one-electron density matrix to the I(I,J) and I(I,A) blocks of the
** Lagrangian.  These terms arise from the first-order CPHF
** equations.  I *think* the following code is general enough to deal
** with both RHF and ROHF cases.
** */

void relax_I(void)
{
  struct oe_dpdfile I, D, f;
  struct dpdbuf E;
  int h, nirreps, i, j, e, *occpi, *virtpi, *openpi;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;

  /*** occupied-virtual relaxation terms */

  /* I(I,A) = I'(I,A) - sum_M f(I,M) D(orb)(A,M) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);
  dpd_oe_copy(&I, CC_OEI, "I(I,A)", 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I(I,A)", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 1, 0, "D(orb)(A,I)", 0, outfile);
  dpd_oe_file_init(&f, CC_OEI, 0, 0, "fIJ", 0, outfile);
  dpd_contract111(&f, &D, &I, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&f);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&I);

  /* I(i,a) = I'(i,a) - sum_m f(i,m) D(orb)(a,m) */
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);
  dpd_oe_copy(&I, CC_OEI, "I(i,a)", 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I(i,a)", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 1, 0, "D(orb)(a,i)", 0, outfile);
  dpd_oe_file_init(&f, CC_OEI, 0, 0, "fij", 0, outfile);
  dpd_contract111(&f, &D, &I, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&f);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&I);

  /*** occupied-occupied relaxtion terms */

  /* I(I,J) <-- I'(I,J) - sum_E,M D(orb)(E,M) [<EI||MJ> + <EJ||MI>]
                      - 2 sum_e,m D(orb)(e,m) <eI|mJ> */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);
  dpd_oe_copy(&I, CC_OEI, "I(I,J)", 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I(I,J)", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 1, 0, "D(orb)(A,I)", 0, outfile);
  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 1, "E <ai|jk>", 0, outfile);
  dpd_dot13(&D, &E, &I, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_dot13(&D, &E, &I, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&E);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 1, 0, "D(orb)(a,i)", 0, outfile);
  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);
  dpd_dot13(&D, &E, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&E);
  dpd_oe_file_close(&D);

  /* I(I,J) <-- - 2 sum_E  f(E,I) D(orb)(E,J) (J,socc)

   Note that this same term is not needed in the I(i,j) block since J
   is required to be a singly occupied orbital */
  dpd_oe_file_mat_init(&I);
  dpd_oe_file_mat_rd(&I, 0, outfile);

  dpd_oe_file_init(&f, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_file_mat_init(&f);
  dpd_oe_file_mat_rd(&f, 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 1, 0, "D(orb)(A,I)", 0, outfile);
  dpd_oe_file_mat_init(&D);
  dpd_oe_file_mat_rd(&D, 0, outfile);

  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++)
	  for(j=(occpi[h] - openpi[h]); j < occpi[h]; j++)
	      for(e=0; e < virtpi[h]; e++)
		  I.matrix[h][i][j] -= 2 * f.matrix[h][i][e] * D.matrix[h][e][j];
    }

  dpd_oe_file_mat_close(&D);
  dpd_oe_file_close(&D);
  dpd_oe_file_mat_close(&f);
  dpd_oe_file_close(&f);

  dpd_oe_file_mat_wrt(&I, 0, outfile);
  dpd_oe_file_mat_close(&I);
  dpd_oe_file_close(&I);

  /* I(i,j) <-- I'(i,j) - sum_e,m D(orb)(e,m) [<ei||mj> + <ej||mi>]
                      - 2 sum_E,M D(orb)(E,M) <Ei|Mj> */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);
  dpd_oe_copy(&I, CC_OEI, "I(i,j)", 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I(i,j)", 0, outfile);
  dpd_oe_file_init(&D, CC_OEI, 1, 0, "D(orb)(a,i)", 0, outfile);
  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 1, "E <ai|jk>", 0, outfile);
  dpd_dot13(&D, &E, &I, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_dot13(&D, &E, &I, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&E);
  dpd_oe_file_close(&D);
  dpd_oe_file_init(&D, CC_OEI, 1, 0, "D(orb)(A,I)", 0, outfile);
  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);
  dpd_dot13(&D, &E, &I, 0, 0, -2.0, 1.0, 0, outfile);
  dpd_buf_close(&E);
  dpd_oe_file_close(&D);
  dpd_oe_file_close(&I);

  /* Clean the I(i,j) block yet again */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I(i,j)", 0, outfile);
  dpd_oe_file_mat_init(&I);
  dpd_oe_file_mat_rd(&I, 0, outfile);

  for(h=0; h < nirreps; h++)
      for(i=0; i < occpi[h]; i++)
	  for(j=0; j < occpi[h]; j++)
	      if((i >= (occpi[h] - openpi[h])) ||
		 (j >= (occpi[h] - openpi[h])) )
		  I.matrix[h][i][j] = 0.0;

  dpd_oe_file_mat_wrt(&I, 0, outfile);
  dpd_oe_file_mat_close(&I);
  dpd_oe_file_close(&I);
}
