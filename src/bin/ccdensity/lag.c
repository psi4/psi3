#include <dpd.h>
#define EXTERN
#include "globals.h"

/* lag(): Build the orbital Lagrangian, I'pq, defined in spin-orbitals
** as:
**
** I'pq = sum_r fpr (Dqr + Drq) + sum_rs <pr||qs> (Drs + Dsr) (q,occ)
**        + sum_rst <pr||st> Gqrst + 2 fpq (q,occ)
**
** The orbital-response component of the gradient (for non-correlated
** orbitals) is defined as
**
** dE/dx <--- sum_pq I'pq U(x)pq
**
** where U(x)pq is the usual CPHF coefficient.  Note, however, that
** the final expression we want involves not CPHF coefficients, but
** overlap derivatives.  For example, in the occupied-occupied and
** virtual-virtual blocks, a choice of non-canonical perturbed
** orbitals allows the assignments
**
** U(x)ij = -1/2 S(x)ij      and      U(x)ab = -1/2 S(x)ab
**
** to be made.  We also choose to incorporate the -1/2 prefactor
** into the Largrangian itself so the final orbital response
** expression will appear as
**
** dE/dx <--- sum_pq Ipq S(x)pq
**
** where Ipq is the "relaxed" Lagrangian (see relax_I.c).
**
** The final set of loops force the appropriate open-shell terms to
** zero. (See the description of the treatment of open-shells in this
** ROHF-CCSD code as discussed in CCSORT for an explanation of why this
** is necessary.)
** */

void Iij(void);
void Iab(void);
void Iai(void);
void Iia(void);

void lag(void)
{
  int h, nirreps, i, j, a, b;
  int *occpi, *virtpi, *openpi;
  struct oe_dpdfile I;
  
  Iij();
  Iab();
  Iai();
  Iia();

  /* Multiply all I'pq components by -1/2 for compatibility with the
     final gradient expression */
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);
  dpd_oe_scm(&I, -0.5, 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);
  dpd_oe_scm(&I, -0.5, 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'AB", 0, outfile);
  dpd_oe_scm(&I, -0.5, 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'ab", 0, outfile);
  dpd_oe_scm(&I, -0.5, 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);
  dpd_oe_scm(&I, -0.5, 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);
  dpd_oe_scm(&I, -0.5, 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 1, 0, "I'AI", 0, outfile);
  dpd_oe_scm(&I, -0.5, 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_init(&I, CC_OEI, 1, 0, "I'ai", 0, outfile);
  dpd_oe_scm(&I, -0.5, 0, outfile);
  dpd_oe_file_close(&I);

  /* Now go through all terms involving open-shell orbitals and force
     the appropriate spin cases to zero. */
  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi; openpi = moinfo.openpi;

  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'IJ", 0, outfile);
  dpd_oe_file_close(&I);

  dpd_oe_file_init(&I, CC_OEI, 0, 0, "I'ij", 0, outfile);
  dpd_oe_file_mat_init(&I);
  dpd_oe_file_mat_rd(&I, 0, outfile);
  for(h=0; h < nirreps; h++) {
      for(i=(occpi[h]-openpi[h]); i < occpi[h]; i++) {
	  for(j=(occpi[h]-openpi[h]); j < occpi[h]; j++) {
	      I.matrix[h][i][j] = 0.0;
	    }
	}
      for(i=(occpi[h]-openpi[h]); i < occpi[h]; i++) {
	  for(j=0; j < occpi[h]; j++) {
	      I.matrix[h][i][j] = 0.0;
	    }
	}
      for(i=0; i < occpi[h]; i++) {
	  for(j=(occpi[h]-openpi[h]); j < occpi[h]; j++) {
	      I.matrix[h][i][j] = 0.0;
	    }
	}
    }
  dpd_oe_file_mat_wrt(&I, 0, outfile);
  dpd_oe_file_mat_close(&I);
  dpd_oe_file_close(&I);

  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'AB", 0, outfile);
  dpd_oe_file_mat_init(&I);
  dpd_oe_file_mat_rd(&I, 0, outfile);
  for(h=0; h < nirreps; h++) {
      for(a=(virtpi[h]-openpi[h]); a < virtpi[h]; a++) {
	  for(b=(virtpi[h]-openpi[h]); b < virtpi[h]; b++) {
	      I.matrix[h][a][b] = 0.0;
	    }
	}
      for(a=(virtpi[h]-openpi[h]); a < virtpi[h]; a++) {
	  for(b=0; b < virtpi[h]; b++) {
	      I.matrix[h][a][b] = 0.0;
	    }
	}
      for(a=0; a < virtpi[h]; a++) {
	  for(b=(virtpi[h]-openpi[h]); b < virtpi[h]; b++) {
	      I.matrix[h][a][b] = 0.0;
	    }
	}
    }
  dpd_oe_file_mat_wrt(&I, 0, outfile);
  dpd_oe_file_mat_close(&I);
  dpd_oe_file_close(&I);

  dpd_oe_file_init(&I, CC_OEI, 1, 1, "I'ab", 0, outfile);
  dpd_oe_file_close(&I);

  dpd_oe_file_init(&I, CC_OEI, 1, 0, "I'AI", 0, outfile);
  dpd_oe_file_mat_init(&I);
  dpd_oe_file_mat_rd(&I, 0, outfile);
  for(h=0; h < nirreps; h++) {
      for(a=(virtpi[h]-openpi[h]); a < virtpi[h]; a++) {
	  for(i=0; i < occpi[h]; i++) {
	      I.matrix[h][a][i] = 0.0;
	    }
	}
    }
  dpd_oe_file_mat_wrt(&I, 0, outfile);
  dpd_oe_file_mat_close(&I);
  dpd_oe_file_close(&I);

  dpd_oe_file_init(&I, CC_OEI, 1, 0, "I'ai", 0, outfile);
  dpd_oe_file_mat_init(&I);
  dpd_oe_file_mat_rd(&I, 0, outfile);
  for(h=0; h < nirreps; h++) {
      for(a=0; a < virtpi[h]; a++) {
	  for(i=(occpi[h] - openpi[h]); i < occpi[h]; i++) {
	      I.matrix[h][a][i] = 0.0;
	    }
	}
    }
  dpd_oe_file_mat_wrt(&I, 0, outfile);
  dpd_oe_file_mat_close(&I);
  dpd_oe_file_close(&I);

  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);
  dpd_oe_file_mat_init(&I);
  dpd_oe_file_mat_rd(&I, 0, outfile);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
	  for(a=(virtpi[h] - openpi[h]); a < virtpi[h]; a++) {
	      I.matrix[h][i][a] = 0.0;
	    }
	}
    }
  dpd_oe_file_mat_wrt(&I, 0, outfile);
  dpd_oe_file_mat_close(&I);
  dpd_oe_file_close(&I);

  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);
  dpd_oe_file_mat_init(&I);
  dpd_oe_file_mat_rd(&I, 0, outfile);
  for(h=0; h < nirreps; h++) {
      for(i=(occpi[h] - openpi[h]); i < occpi[h]; i++) {
	  for(a=0; a < virtpi[h]; a++) {
	      I.matrix[h][i][a] = 0.0;
	    }
	}
    }
  dpd_oe_file_mat_wrt(&I, 0, outfile);
  dpd_oe_file_mat_close(&I);
  dpd_oe_file_close(&I);
}
