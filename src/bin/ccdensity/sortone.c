#include <stdio.h>
#include <dpd.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <psio.h>
#include <iwl.h>
#define EXTERN
#include "globals.h"
#include <psifiles.h>

/*
** SORTONE(): Place all the components of the 1pdm into a large
** matrix, O (moinfo.opdm), which we also symmetrize by computing
** Opq = 1/2 (Opq + Oqp).  This matrix is later written to disk in
** dump() for subsequent backtransformation.  Note that the components
** of the 1pdm computed into the DIJ, Dij, DAB, Dab, DAI, Dai, DIA,
** and Dia matrices remain non-symmetric (e.g., DIJ neq DJI).
*/

void sortone(void)
{
  int h, nirreps, nmo, nfzv, nfzc, nclsd, nopen;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *occpi, *virtpi, *occ_off, *vir_off; 
  int *occ_sym, *vir_sym, *openpi;
  int *qt_occ, *qt_vir;
  double **O, chksum, value;
  dpdfile2 D;
  psio_address next;

  nmo = moinfo.nmo;
  nfzc = moinfo.nfzc;
  nfzv = moinfo.nfzv;
  nclsd = moinfo.nclsd;
  nopen = moinfo.nopen;
  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;
  qt_occ = moinfo.qt_occ; qt_vir = moinfo.qt_vir;

  O = init_matrix(nmo-nfzv,nmo-nfzv);

  /* Sort A components first */
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(j=0; j < occpi[h]; j++) {
              J = qt_occ[occ_off[h] + j];
              O[I][J] += D.matrix[h][i][j];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(a=0; a < (virtpi[h] - openpi[h]); a++) {
          A = qt_vir[vir_off[h] + a];
          for(b=0; b < (virtpi[h] - openpi[h]); b++) {
              B = qt_vir[vir_off[h] + b];

              O[A][B] += D.matrix[h][a][b];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < (virtpi[h] - openpi[h]); a++) {
              A = qt_vir[vir_off[h] + a];

              O[A][I] += D.matrix[h][i][a];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < (virtpi[h] - openpi[h]); a++) {
              A = qt_vir[vir_off[h] + a];

              O[I][A] += D.matrix[h][i][a];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Sort B components */
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_file2_mat_init(&D); 
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < (occpi[h] - openpi[h]); i++) { 
          I = qt_occ[occ_off[h] + i];
          for(j=0; j < (occpi[h] - openpi[h]); j++) {
              J = qt_occ[occ_off[h] + j];
              O[I][J] += D.matrix[h][i][j];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(a=0; a < virtpi[h]; a++) {
          A = qt_vir[vir_off[h] + a];
          for(b=0; b < virtpi[h]; b++) {
              B = qt_vir[vir_off[h] + b];

              O[A][B] += D.matrix[h][a][b];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < virtpi[h]; a++) {
              A = qt_vir[vir_off[h] + a];

              O[A][I] += D.matrix[h][i][a];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < virtpi[h]; a++) {
              A = qt_vir[vir_off[h] + a];

              O[I][A] += D.matrix[h][i][a];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Symmetrize the onepdm */
  for(p=0; p < (nmo-nfzv); p++) {
      for(q=0; q < p; q++) {
          value = 0.5 * (O[p][q] + O[q][p]);
          O[p][q] = O[q][p] = value;
        }
    }

  moinfo.opdm = O;

  /* Write the correlation-only opdm to disk */
/*
  psio_open(PSIF_OEI, PSIO_OPEN_OLD);
  next = PSIO_ZERO;
  for(p=0; p < (nmo - nfzv); p++)
    psio_write(PSIF_OEI, "CC OPDM", (char *) &(O[p][0]), (nmo-nfzv)*sizeof(double), next, &next);
  psio_close(PSIF_OEI, 1);
*/

}
