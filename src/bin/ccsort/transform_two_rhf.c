/* transform_two_rhf(): Carry out a four-index transform of the
** SO-basis two-electron integrals using DPD constructs.
**
** NB: Two key DPD pairs for the SO/MO subspaces:
**  pair   spaces
**    0    All SO,SO
**    3    SO >= SO (+)
**    5    All SO,SO
**    8    MO >= MO (+)
**
** TDC, March 2006
*/

#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include <ccfiles.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

int file_build_presort(dpdfile4 *, int, double, int);

void transform_two_rhf(void)
{
  dpdfile4 I;
  dpdbuf4 J, K;
  int h, pq;
  double **TMP;
  int *sopi, *mopi;
  int nso, nmo, nrows, ncols, nlinks, J_rs, K_rs, Gs, Gr;
  double ***C;

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  sopi = moinfo.sopi;
  mopi = moinfo.mopi;
  C = (double ***) malloc(moinfo.nirreps * sizeof(double **));
  chkpt_init(PSIO_OPEN_OLD);
  for(h=0; h < moinfo.nirreps; h++)
    C[h] = chkpt_rd_scf_irrep(h);
  chkpt_close();

  dpd_file4_init(&I, PSIF_SO_PRESORT, 0, 3, 3, "SO Ints (pq,rs)");
  file_build_presort(&I, PSIF_SO_TEI, params.tolerance, 1);
  dpd_file4_close(&I);

  TMP = block_matrix(nso,nso);

  dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (pq,rs)");
  dpd_buf4_init(&K, CC_TEI_HALFT, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  for(h=0; h < moinfo.nirreps; h++) {
    dpd_buf4_mat_irrep_row_init(&J, h);
    dpd_buf4_mat_irrep_row_init(&K, h);

    for(pq=0; pq < J.params->rowtot[h]; pq++) {
      dpd_buf4_mat_irrep_row_rd(&J, h, pq);

      for(Gr=0; Gr < moinfo.nirreps; Gr++) {
	Gs = h^Gr;
 
	nrows = sopi[Gr];
	ncols = mopi[Gs];
	nlinks = sopi[Gs];
	J_rs = J.col_offset[h][Gr];
	if(nrows && ncols && nlinks)
	  C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][0][J_rs],nlinks,C[Gs][0],
		  nlinks,0.0,TMP[0],nso);

	nrows = mopi[Gr];
	ncols = mopi[Gs];
	nlinks = sopi[Gr];
	K_rs = K.col_offset[h][Gr];
	if(nrows && ncols && nlinks)
	  C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[Gr][0],nlinks,TMP[0],nso,
		  0.0,&K.matrix[h][0][K_rs],ncols);
      }

      dpd_buf4_mat_irrep_row_wrt(&K, h, pq);
    }

    dpd_buf4_mat_irrep_row_close(&J, h);
    dpd_buf4_mat_irrep_row_close(&K, h);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  dpd_buf4_init(&K, CC_TEI_HALFT, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  dpd_buf4_sort(&K, CC_TEI_HALFT, rspq, 8, 3, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_close(&K);

  dpd_buf4_init(&J, CC_TEI_HALFT, 0, 8, 0, 8, 3, 0, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_init(&K, CC_MISC, 0, 8, 5, 8, 8, 0, "MO Ints (ij,kl)");
  for(h=0; h < moinfo.nirreps; h++) {
    dpd_buf4_mat_irrep_row_init(&J, h);
    dpd_buf4_mat_irrep_row_init(&K, h);

    for(pq=0; pq < J.params->rowtot[h]; pq++) {
      dpd_buf4_mat_irrep_row_rd(&J, h, pq);

      for(Gr=0; Gr < moinfo.nirreps; Gr++) {
	Gs = h^Gr;
 
	nrows = sopi[Gr];
	ncols = mopi[Gs];
	nlinks = sopi[Gs];
	J_rs = J.col_offset[h][Gr];
	if(nrows && ncols && nlinks)
	  C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][0][J_rs],nlinks,C[Gs][0],
		  nlinks,0.0,TMP[0],nso);

	nrows = mopi[Gr];
	ncols = mopi[Gs];
	nlinks = sopi[Gr];
	K_rs = K.col_offset[h][Gr];
	if(nrows && ncols && nlinks)
	  C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[Gr][0],nlinks,TMP[0],nso,
		  0.0,&K.matrix[h][0][K_rs],ncols);
      }

      dpd_buf4_mat_irrep_row_wrt(&K, h, pq);
    }

    dpd_buf4_mat_irrep_row_close(&J, h);
    dpd_buf4_mat_irrep_row_close(&K, h);
  }
  dpd_buf4_print(&K, outfile, 1);
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  free_block(TMP);
  for(h=0; h < moinfo.nirreps; h++)
    free_block(C[h]);
  free(C);
}
