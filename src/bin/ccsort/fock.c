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
** An alternative but currently unused algorithm may be found in fock_build.c
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
  dpdfile2 fIJ, fij, fAB, fab, fIA, fia, Hoo, Hvv, Hov;
  dpdbuf4 AInts_anti, AInts, CInts, CInts_anti, EInts_anti, EInts;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  dpd_file2_init(&Hoo, CC_OEI, 0, 0, 0, "h(i,j)");
  dpd_file2_init(&Hvv, CC_OEI, 0, 1, 1, "h(a,b)");
  dpd_file2_init(&Hov, CC_OEI, 0, 0, 1, "h(i,a)");
  dpd_file2_mat_init(&Hoo);
  dpd_file2_mat_init(&Hvv);
  dpd_file2_mat_init(&Hov);
  dpd_file2_mat_rd(&Hoo);
  dpd_file2_mat_rd(&Hvv);
  dpd_file2_mat_rd(&Hov);
  
  /* Prepare the alpha and beta occ-occ Fock matrix files */
  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fij);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(i=0; i < occpi[h]; i++) 
          for(j=0; j < occpi[h]; j++) 
              fIJ.matrix[h][i][j] = Hoo.matrix[h][i][j];   

      for(i=0; i < (occpi[h]-openpi[h]); i++) 
          for(j=0; j < (occpi[h]-openpi[h]); j++) 
              fij.matrix[h][i][j] = Hoo.matrix[h][i][j];   
    }

  /* Two-electron contributions */

  /* Prepare the A integral buffers */
  dpd_buf4_init(&AInts_anti, CC_AINTS, 0, 0, 0, 0, 0, 1, "A <ij|kl>");
  dpd_buf4_init(&AInts, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&AInts_anti, h);
      dpd_buf4_mat_irrep_rd(&AInts_anti, h);

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
      
      dpd_buf4_mat_irrep_close(&AInts_anti, h);

      dpd_buf4_mat_irrep_init(&AInts, h);
      dpd_buf4_mat_irrep_rd(&AInts, h);

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
      
      dpd_buf4_mat_irrep_close(&AInts, h);

    }

  /* Close the A Integral buffers */
  dpd_buf4_close(&AInts_anti);
  dpd_buf4_close(&AInts);
  
  /* Close the alpha and beta occ-occ Fock matrix files */
  dpd_file2_mat_wrt(&fIJ);
  dpd_file2_mat_wrt(&fij);
  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fij);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fij);

  /* Prepare the alpha and beta vir-vir Fock matrix files */
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_init(&fab);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(a=0; a < (virtpi[h] - openpi[h]); a++) 
          for(b=0; b < (virtpi[h] - openpi[h]); b++) 
              fAB.matrix[h][a][b] = Hvv.matrix[h][a][b];   

      for(a=0; a < virtpi[h]; a++) 
          for(b=0; b < virtpi[h]; b++) 
              fab.matrix[h][a][b] = Hvv.matrix[h][a][b];   
    }

  /* Two-electron contributions */

  /* Prepare the C integral buffers */
  dpd_buf4_init(&CInts_anti, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
  dpd_buf4_init(&CInts, CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&CInts_anti, h);
      dpd_buf4_mat_irrep_rd(&CInts_anti, h);

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

      dpd_buf4_mat_irrep_close(&CInts_anti, h);

      dpd_buf4_mat_irrep_init(&CInts, h);
      dpd_buf4_mat_irrep_rd(&CInts, h);

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

      dpd_buf4_mat_irrep_close(&CInts, h);
    }

  /* Close the C integral buffers */
  dpd_buf4_close(&CInts_anti);
  dpd_buf4_close(&CInts);

  /* Close the alpha and beta vir-vir Fock matrix files */
  dpd_file2_mat_wrt(&fAB);
  dpd_file2_mat_wrt(&fab);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fab);

  /* Prepare the alpha and beta occ-vir Fock matrix files */
  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
  dpd_file2_mat_init(&fIA);
  dpd_file2_mat_init(&fia);

  /* One-electron (frozen-core) contributions */
  for(h=0; h < nirreps; h++) {

      for(i=0; i < occpi[h]; i++) 
          for(a=0; a < (virtpi[h] - openpi[h]); a++) 
              fIA.matrix[h][i][a] = Hov.matrix[h][i][a];   

      for(i=0; i < (occpi[h] - openpi[h]); i++) 
          for(a=0; a < virtpi[h]; a++) 
              fia.matrix[h][i][a] = Hov.matrix[h][i][a];   
    }

  /* Two-electron contributions */

  /* Prepare the E integral buffers */
  dpd_buf4_init(&EInts_anti, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  dpd_buf4_init(&EInts, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(&EInts_anti, h);
      dpd_buf4_mat_irrep_rd(&EInts_anti, h);

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

      dpd_buf4_mat_irrep_close(&EInts_anti, h);

      dpd_buf4_mat_irrep_init(&EInts, h);
      dpd_buf4_mat_irrep_rd(&EInts, h);

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

      dpd_buf4_mat_irrep_close(&EInts, h);

    }

  /* Close the E integral buffers */
  dpd_buf4_close(&EInts_anti);
  dpd_buf4_close(&EInts);

  /* Close the alpha and beta occ-vir Fock matrix files */
  dpd_file2_mat_wrt(&fIA);
  dpd_file2_mat_wrt(&fia);
  dpd_file2_mat_close(&fIA);
  dpd_file2_mat_close(&fia);
  dpd_file2_close(&fIA);
  dpd_file2_close(&fia);

  dpd_file2_mat_close(&Hoo);
  dpd_file2_mat_close(&Hvv);
  dpd_file2_mat_close(&Hov);
  dpd_file2_close(&Hoo);
  dpd_file2_close(&Hvv);
  dpd_file2_close(&Hov);
}
