#include <stdio.h>
#include <psio.h>
#include <libciomr.h>
#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

/* fock(): Build the alpha and beta Fock matrices from the
** one-electron integrals/frozen-core operator and active two-electron
** integrals on disk.
**
** Daniel Crawford, 1996
**
** Since the full spin-orbital Fock matrices are also built using an
** out-of-core algorithm in fock_build() [part of the two-electron
** integral sort in sort_tei()], this code is actually redundant.
** However, I'm leaving it in place since it already produces the DPD
** blocks of Fock matrix elements needed by CCENERGY.
** */

void fock(void)
{
  int h, Gi, Gj, Ga, Gb, Gm;
  int i,j,a,b,m;
  int I, J, A, B, M;
  int IM, JM, MA, MB, AM;
  int nirreps;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *occ_sym, *vir_sym;
  int *openpi;
  struct oe_dpdfile fIJ, fij, fAB, fab, fIA, fia;
  struct dpdbuf AInts_anti, AInts, CInts, CInts_anti, EInts_anti, EInts;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  /* Prepare the alpha and beta occ-occ Fock matrix files */
  dpd_oe_file_init(&fIJ, CC_OEI, 0, 0, "fIJ", 0, outfile);
  dpd_oe_file_init(&fij, CC_OEI, 0, 0, "fij", 0, outfile);
  dpd_oe_file_mat_init(&fIJ);
  dpd_oe_file_mat_init(&fij);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(i=0; i < occpi[h]; i++) 
          for(j=0; j < occpi[h]; j++) 
              fIJ.matrix[h][i][j] = hoo[h][i][j];   

      for(i=0; i < (occpi[h]-openpi[h]); i++) 
          for(j=0; j < (occpi[h]-openpi[h]); j++) 
              fij.matrix[h][i][j] = hoo[h][i][j];   
    }

  /* Two-electron contributions */

  /* Prepare the A integral buffers */
  dpd_buf_init(&AInts_anti, CC_AINTS, 0, 0, 0, 0, 1, "A <ij|kl>", 0, outfile);
  dpd_buf_init(&AInts, CC_AINTS, 0, 0, 0, 0, 0, "A <ij|kl>", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&AInts_anti, h);
      dpd_buf_mat_irrep_rd(&AInts_anti, h, 0, outfile);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
          Gj=Gi; Gm=h^Gi;

          /* Loop over the orbitals of the target */
          for(i=0; i < occpi[Gi]; i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < occpi[Gj]; j++) {
                  J = occ_off[Gj] + j;
                  for(m=0; m < occpi[Gm]; m++) {
                      M = occ_off[Gm] + m;

                      IM = AInts_anti.params->rowidx[I][M];
                      JM = AInts_anti.params->colidx[J][M];

                      fIJ.matrix[Gi][i][j] += AInts_anti.matrix[h][IM][JM];

                    }
                }
            }

          /* Loop over the orbitals of the target */
          for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < (occpi[Gj] - openpi[Gj]); j++) {
                  J = occ_off[Gj] + j;
                  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
                      M = occ_off[Gm] + m;

                      IM = AInts_anti.params->rowidx[I][M];
                      JM = AInts_anti.params->colidx[J][M];

                      fij.matrix[Gi][i][j] += AInts_anti.matrix[h][IM][JM];

                    }
                }
            }

        }
      
      dpd_buf_mat_irrep_close(&AInts_anti, h);

      dpd_buf_mat_irrep_init(&AInts, h);
      dpd_buf_mat_irrep_rd(&AInts, h, 0, outfile);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
          Gj=Gi; Gm=h^Gi;

          /* Loop over the orbitals of the target */
          for(i=0; i < occpi[Gi]; i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < occpi[Gj]; j++) {
                  J = occ_off[Gj] + j;
                  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
                      M = occ_off[Gm] + m;

                      IM = AInts.params->rowidx[I][M];
                      JM = AInts.params->colidx[J][M];

                      fIJ.matrix[Gi][i][j] += AInts.matrix[h][IM][JM];

                    }
                }
            }

          /* Loop over the orbitals of the target */
          for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
              I = occ_off[Gi] + i;
              for(j=0; j < (occpi[Gj] - openpi[Gj]); j++) {
                  J = occ_off[Gj] + j;
                  for(m=0; m < occpi[Gm]; m++) {
                      M = occ_off[Gm] + m;

                      IM = AInts.params->rowidx[I][M];
                      JM = AInts.params->colidx[J][M];

                      fij.matrix[Gi][i][j] += AInts.matrix[h][IM][JM];

                    }
                }
            }

        }
      
      dpd_buf_mat_irrep_close(&AInts, h);

    }

  /* Close the A Integral buffers */
  dpd_buf_close(&AInts_anti);
  dpd_buf_close(&AInts);
  
  /* Close the alpha and beta occ-occ Fock matrix files */
  dpd_oe_file_mat_wrt(&fIJ, 0, outfile);
  dpd_oe_file_mat_wrt(&fij, 0, outfile);
  dpd_oe_file_mat_close(&fIJ);
  dpd_oe_file_mat_close(&fij);
  dpd_oe_file_close(&fIJ);
  dpd_oe_file_close(&fij);

  /* Prepare the alpha and beta vir-vir Fock matrix files */
  dpd_oe_file_init(&fAB, CC_OEI, 1, 1, "fAB", 0, outfile);
  dpd_oe_file_init(&fab, CC_OEI, 1, 1, "fab", 0, outfile);
  dpd_oe_file_mat_init(&fAB);
  dpd_oe_file_mat_init(&fab);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(a=0; a < (virtpi[h] - openpi[h]); a++) 
          for(b=0; b < (virtpi[h] - openpi[h]); b++) 
              fAB.matrix[h][a][b] = hvv[h][a][b];   

      for(a=0; a < virtpi[h]; a++) 
          for(b=0; b < virtpi[h]; b++) 
              fab.matrix[h][a][b] = hvv[h][a][b];   
    }

  /* Two-electron contributions */

  /* Prepare the C integral buffers */
  dpd_buf_init(&CInts_anti, CC_CINTS, 10, 10, 10, 10, 0, 
               "C <ia||jb>", 0, outfile);
  dpd_buf_init(&CInts, CC_CINTS, 10, 10, 10, 10, 0, "C <ia|jb>", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&CInts_anti, h);
      dpd_buf_mat_irrep_rd(&CInts_anti, h, 0, outfile);

      /* Loop over irreps of the target */
      for(Ga=0; Ga < nirreps; Ga++) {
	  Gb = Ga; Gm = h^Ga;

	  /* Loop over orbitals of the target */
	  for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
	      A = vir_off[Ga] + a;
	      for(b=0; b < (virtpi[Gb] - openpi[Gb]); b++) {
		  B = vir_off[Gb] + b;

		  for(m=0; m < occpi[Gm]; m++) {
		      M = occ_off[Gm] + m;

		      MA = CInts_anti.params->rowidx[M][A];
		      MB = CInts_anti.params->colidx[M][B];

		      fAB.matrix[Ga][a][b] += CInts_anti.matrix[h][MA][MB];

		    }
		}
	    }

	  /* Loop over orbitals of the target */
	  for(a=0; a < virtpi[Ga]; a++) {
	      A = vir_off[Ga] + a;
	      for(b=0; b < virtpi[Gb]; b++) {
		  B = vir_off[Gb] + b;

		  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
		      M = occ_off[Gm] + m;

		      MA = CInts_anti.params->rowidx[M][A];
		      MB = CInts_anti.params->colidx[M][B];

		      fab.matrix[Ga][a][b] += CInts_anti.matrix[h][MA][MB];
		    }
		}
	    }
	}

      dpd_buf_mat_irrep_close(&CInts_anti, h);

      dpd_buf_mat_irrep_init(&CInts, h);
      dpd_buf_mat_irrep_rd(&CInts, h, 0, outfile);

      /* Loop over irreps of the target */
      for(Ga=0; Ga < nirreps; Ga++) {
	  Gb = Ga; Gm = h^Ga;

	  /* Loop over orbitals of the target */
	  for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
	      A = vir_off[Ga] + a;
	      for(b=0; b < (virtpi[Gb] - openpi[Gb]); b++) {
		  B = vir_off[Gb] + b;

		  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
		      M = occ_off[Gm] + m;

		      MA = CInts.params->rowidx[M][A];
		      MB = CInts.params->colidx[M][B];

		      fAB.matrix[Ga][a][b] += CInts.matrix[h][MA][MB];

		    }
		}
	    }

	  /* Loop over orbitals of the target */
	  for(a=0; a < virtpi[Ga]; a++) {
	      A = vir_off[Ga] + a;
	      for(b=0; b < virtpi[Gb]; b++) {
		  B = vir_off[Gb] + b;

		  for(m=0; m < occpi[Gm]; m++) {
		      M = occ_off[Gm] + m;

		      MA = CInts.params->rowidx[M][A];
		      MB = CInts.params->colidx[M][B];

		      fab.matrix[Ga][a][b] += CInts.matrix[h][MA][MB];
		    }
		}
	    }
	}

      dpd_buf_mat_irrep_close(&CInts, h);
    }

  /* Close the C integral buffers */
  dpd_buf_close(&CInts_anti);
  dpd_buf_close(&CInts);

  /* Close the alpha and beta vir-vir Fock matrix files */
  dpd_oe_file_mat_wrt(&fAB, 0, outfile);
  dpd_oe_file_mat_wrt(&fab, 0, outfile);
  dpd_oe_file_mat_close(&fAB);
  dpd_oe_file_mat_close(&fab);
  dpd_oe_file_close(&fAB);
  dpd_oe_file_close(&fab);

  /* Prepare the alpha and beta occ-vir Fock matrix files */
  dpd_oe_file_init(&fIA, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_file_init(&fia, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_file_mat_init(&fIA);
  dpd_oe_file_mat_init(&fia);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(i=0; i < occpi[h]; i++) 
          for(a=0; a < (virtpi[h] - openpi[h]); a++) 
              fIA.matrix[h][i][a] = hov[h][i][a];   

      for(i=0; i < (occpi[h] - openpi[h]); i++) 
          for(a=0; a < virtpi[h]; a++) 
              fia.matrix[h][i][a] = hov[h][i][a];   
    }

  /* Two-electron contributions */

  /* Prepare the E integral buffers */
  dpd_buf_init(&EInts_anti, CC_EINTS, 11, 0, 11, 0, 1,
	       "E <ai|jk>", 0, outfile);
  dpd_buf_init(&EInts, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&EInts_anti, h);
      dpd_buf_mat_irrep_rd(&EInts_anti, h, 0, outfile);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
	  Ga = Gi; Gm = h^Gi;

	  /* Loop over orbitals of the target */
	  for(i=0; i < occpi[Gi]; i++) {
	      I = occ_off[Gi] + i;
	      for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
		  A = vir_off[Ga] + a;

		  for(m=0; m < occpi[Gm]; m++) {
		      M = occ_off[Gm] + m;

		      AM = EInts_anti.params->rowidx[A][M];
		      IM = EInts_anti.params->colidx[I][M];

		      fIA.matrix[Gi][i][a] += EInts_anti.matrix[h][AM][IM];

		    }
		}
	    }

	  /* Loop over orbitals of the target */
	  for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
	      I = occ_off[Gi] + i;
	      for(a=0; a < virtpi[Ga]; a++) {
		  A = vir_off[Ga] + a;

		  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
		      M = occ_off[Gm] + m;

		      AM = EInts_anti.params->rowidx[A][M];
		      IM = EInts_anti.params->colidx[I][M];

		      fia.matrix[Gi][i][a] += EInts_anti.matrix[h][AM][IM];

		    }
		}
	    }
	}

      dpd_buf_mat_irrep_close(&EInts_anti, h);

      dpd_buf_mat_irrep_init(&EInts, h);
      dpd_buf_mat_irrep_rd(&EInts, h, 0, outfile);

      /* Loop over irreps of the target */
      for(Gi=0; Gi < nirreps; Gi++) {
	  Ga = Gi; Gm = h^Gi;

	  /* Loop over orbitals of the target */
	  for(i=0; i < occpi[Gi]; i++) {
	      I = occ_off[Gi] + i;
	      for(a=0; a < (virtpi[Ga] - openpi[Ga]); a++) {
		  A = vir_off[Ga] + a;

		  for(m=0; m < (occpi[Gm] - openpi[Gm]); m++) {
		      M = occ_off[Gm] + m;

		      AM = EInts.params->rowidx[A][M];
		      IM = EInts.params->colidx[I][M];

		      fIA.matrix[Gi][i][a] += EInts.matrix[h][AM][IM];

		    }
		}
	    }

	  /* Loop over orbitals of the target */
	  for(i=0; i < (occpi[Gi] - openpi[Gi]); i++) {
	      I = occ_off[Gi] + i;
	      for(a=0; a < virtpi[Ga]; a++) {
		  A = vir_off[Ga] + a;

		  for(m=0; m < occpi[Gm]; m++) {
		      M = occ_off[Gm] + m;

		      AM = EInts.params->rowidx[A][M];
		      IM = EInts.params->colidx[I][M];

		      fia.matrix[Gi][i][a] += EInts.matrix[h][AM][IM];

		    }
		}
	    }
	}

      dpd_buf_mat_irrep_close(&EInts, h);

    }

  /* Close the E integral buffers */
  dpd_buf_close(&EInts_anti);
  dpd_buf_close(&EInts);

  /* Close the alpha and beta occ-vir Fock matrix files */
  dpd_oe_file_mat_wrt(&fIA, 0, outfile);
  dpd_oe_file_mat_wrt(&fia, 0, outfile);
  dpd_oe_file_mat_close(&fIA);
  dpd_oe_file_mat_close(&fia);
  dpd_oe_file_close(&fIA);
  dpd_oe_file_close(&fia);
}
