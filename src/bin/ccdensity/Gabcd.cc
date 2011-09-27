/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include <libqt/qt.h>
#include "globals.h"

namespace psi { namespace ccdensity {

void Gabcd(void)
{
  dpdbuf4 G, L, T;
  int G_irr;
  G_irr = params.G_irr;

  int n, h, nirreps, nbuckets;
  long int size_L, size_T, memoryd, num_rows_G, rows_left_G;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    if(strcmp(params.wfn, "OOCCD") || !params.ooccd_ooc)
      dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    else {
//    replacing above contract444 call with following code to perform an out of
//    core contraction for large orbital optimized optical rotation calculations
      fprintf(outfile, "\tRunning out-of-core Gabcd <-- L*T\n");
      dpd_buf4_sort(&L, CC_GLG, rspq, 5, 0, "LAbIj");
      dpd_buf4_close(&L);
      dpd_buf4_init(&L, CC_GLG, G_irr, 5, 0, 5, 0, 0, "LAbIj");

      nirreps = L.params->nirreps;

      for(h=0; h < nirreps; h++) {

        size_L = ((long) L.params->rowtot[h]) * ((long) L.params->coltot[h]);
        size_T = ((long) T.params->rowtot[h]) * ((long) T.params->coltot[h]);

        memoryd = dpd_memfree() - (size_L + size_T);
        num_rows_G = memoryd/G.params->coltot[h];
        nbuckets = (int) G.params->rowtot[h] / num_rows_G;
        
        rows_left_G = G.params->rowtot[h] % num_rows_G;

        dpd_buf4_mat_irrep_init(&L, h);
        dpd_buf4_mat_irrep_rd(&L, h);
      
        dpd_buf4_mat_irrep_init(&T, h);
        dpd_buf4_mat_irrep_rd(&T, h);

        dpd_buf4_mat_irrep_init_block(&G, h, num_rows_G);

        for(n=0; n < nbuckets; n++) {
//          dpd_buf4_mat_irrep_rd_block(&G, h, n*num_rows_G, num_rows_G);
          C_DGEMM('n', 'n', num_rows_G, G.params->coltot[h], L.params->coltot[h], 1,
                  &(L.matrix[h][n*num_rows_G][0]), L.params->coltot[h],
                  &(T.matrix[h][0][0]), T.params->coltot[h], 0,
                  &(G.matrix[h][0][0]), G.params->coltot[h]);
          dpd_buf4_mat_irrep_wrt_block(&G, h, n*num_rows_G, num_rows_G);
        }

        if(rows_left_G) {
//          dpd_buf4_mat_irrep_rd_block(&G, h, n*num_rows_G, rows_left_G);
          C_DGEMM('n', 'n', rows_left_G, G.params->coltot[h], L.params->coltot[h], 1,
                  &(L.matrix[h][n*num_rows_G][0]), L.params->coltot[h],
                  &(T.matrix[h][0][0]), T.params->coltot[h], 0,
                  &(G.matrix[h][0][0]), G.params->coltot[h]); 
          dpd_buf4_mat_irrep_wrt_block(&G, h, n*num_rows_G, rows_left_G);
        }

        dpd_buf4_mat_irrep_close(&L, h);
        dpd_buf4_mat_irrep_close(&T, h);
        dpd_buf4_mat_irrep_close_block(&G, h, num_rows_G);

      }
    }

    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
//    dpd_buf4_print(&G,outfile,0);
//    dpd_buf4_print(&G,outfile,1);
    if(strcmp(params.wfn, "OOCCD")) {
      if (params.ground)
        dpd_buf4_symm(&G);
    }
    dpd_buf4_close(&G);
  }
  else if(params.ref == 1) { /** RHF/ROHF **/

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_buf4_init(&L, CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    if (params.ground)
      dpd_buf4_symm(&G);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "Gabcd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    if (params.ground)
      dpd_buf4_symm(&G);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 5, 5, 5, 5, 0, "GAbCd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    if (params.ground)
      dpd_buf4_symm(&G);
    dpd_buf4_close(&G);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&G, CC_GAMMA, G_irr, 7, 7, 7, 7, 0, "GABCD");
    dpd_buf4_init(&L, CC_GLG, G_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 17, 17, 17, 17, 0, "Gabcd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, G_irr, 28, 28, 28, 28, 0, "GAbCd");
    dpd_buf4_init(&L, CC_GLG, G_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_contract444(&L, &T, &G, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&T);
    dpd_buf4_close(&L);
    dpd_buf4_symm(&G);
    dpd_buf4_close(&G);
  }
}

}} // namespace psi::ccdensity
