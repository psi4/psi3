/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libdpd/dpd.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {
#include <psifiles.h>

/*
** sortone_rhf(): Place all the components of the 1pdm into a large
** matrix, O (moinfo.opdm), which we also symmetrize by computing Opq
** = 1/2 (Opq + Oqp).  This matrix is later written to disk in dump()
** for subsequent backtransformation.  Note that the components of the
** 1pdm computed into the DIJ, Dij, DAB, Dab, DAI, Dai, DIA, and Dia
** matrices remain non-symmetric (e.g., DIJ neq DJI).
**
** This version doesn't work with frozen orbitals yet.
**
** TDC, 1/03
**
** NB: For now, I just multiply the components by two to account for
** spin adaptation.  The factors of two should later move into the I
** build itself.
**
** TDC, 2/2008
*/

void sortone_RHF(struct RHO_Params rho_params)
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

  /* O = block_matrix(nmo-nfzc,nmo-nfzc); */
  O = block_matrix(nmo-nfzv, nmo-nfzv);

  /* Sort A components first */
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(j=0; j < occpi[h]; j++) {
              J = qt_occ[occ_off[h] + j];
              O[I][J] += 2.0 * D.matrix[h][i][j];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(a=0; a < virtpi[h]; a++) {
          A = qt_vir[vir_off[h] + a];
          for(b=0; b < virtpi[h]; b++) {
              B = qt_vir[vir_off[h] + b];

              O[A][B] += 2.0 * D.matrix[h][a][b];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < virtpi[h]; a++) {
              A = qt_vir[vir_off[h] + a];

              O[A][I] += 2.0 * D.matrix[h][i][a];
            }
        }
    }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
      for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          for(a=0; a < virtpi[h]; a++) {
              A = qt_vir[vir_off[h] + a];

              O[I][A] += 2.0 * D.matrix[h][i][a];
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

}

}} // namespace psi::ccdensity
