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

#define CC_TEI_HALFT0 91
#define CC_TEI_HALFT1 92

int file_build_presort(dpdfile4 *, int, double, int);

void transform_two_rhf(void)
{
  dpdfile4 I;
  dpdbuf4 J, K;
  int h, pq, p, q;
  double **TMP;
  int *sopi, *mopi;
  int nso, nmo, nrows, ncols, nlinks, J_rs, K_rs, Gs, Gr;
  double ***C;
  int *pitz2dpd_occ, *pitz2dpd_vir;

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  sopi = moinfo.sopi;
  mopi = moinfo.mopi;
  C = (double ***) malloc(moinfo.nirreps * sizeof(double **));
  chkpt_init(PSIO_OPEN_OLD);
  for(h=0; h < moinfo.nirreps; h++)
    C[h] = chkpt_rd_scf_irrep(h);
  chkpt_close();

  /* We'll need to map from Pitzer directly to CC/DPD ordering, so let's build the lookups */
  pitz2dpd_occ = init_int_array(moinfo.nactive);
  pitz2dpd_vir = init_int_array(moinfo.nactive);
  for(i=0; i < moinfo.nactive; i++) {
    pitz2dpd_occ[i] = moinfo.cc_occ[moinfo.pitz2qt[i]];
    pitz2dpd_vir[i] = moinfo.cc_vir[moinfo.pitz2qt[i]];
  }

  psio_open(PSIF_SO_PRESORT, 0); /* presorted integrals */

  dpd_file4_init(&I, PSIF_SO_PRESORT, 0, 3, 3, "SO Ints (pq,rs)");
  file_build_presort(&I, PSIF_SO_TEI, params.tolerance, 1);  /* still need to add frozen-core operator in here */
  dpd_file4_close(&I);

  TMP = block_matrix(nso,nso);

  psio_open(CC_TEI_HALFT0, 0); /* half-transformed integrals */

  dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (pq,rs)");
  dpd_buf4_init(&K, CC_TEI_HALFT0, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
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

  psio_close(PSIF_SO_PRESORT, 0); /* delete the presorted integrals */
  psio_open(CC_TEI_HALFT1, 0); /* transposed half-transformed integrals */

  dpd_buf4_init(&K, CC_TEI_HALFT0, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  dpd_buf4_sort(&K, CC_TEI_HALFT1, rspq, 8, 3, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_close(&K);

  psio_close(CC_TEI_HALFT0, 0); /* delete the original half-transformed integrals */

  dpd_buf4_init(&J, CC_TEI_HALFT1, 0, 8, 0, 8, 3, 0, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_init(&K, CC_MISC, 0, 8, 5, 8, 8, 0, "MO Ints (ij,kl)"); /* place-holder buffer */
  for(h=0; h < moinfo.nirreps; h++) {
    dpd_buf4_mat_irrep_row_init(&J, h);
    buffer = init_array(K.params->coltot[h]);

    for(pq=0; pq < J.params->rowtot[h]; pq++) {
      if(occ[] && occ[])
	p = J.params->roworb[h][pq][0];
      q = J.params->roworb[h][pq][1];

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
		  0.0,&buffer[K_rs],ncols);

	/* dump the result straight into the distribute function */
      }
      for(rs=0; rs < J.params->coltot[h]; rs++) {
	r = J.params->colworb[h][rs][0];
	s = J.params->colworb[h][rs][1];


      }
    }

    dpd_buf4_mat_irrep_row_close(&J, h);
    dpd_buf4_mat_irrep_row_close(&K, h);
  }
  dpd_buf4_print(&K, outfile, 1);
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  psio_close(CC_TEI_HALFT1, 0); /* delete the transposed half-transformed integrals */

  free_block(TMP);
  for(h=0; h < moinfo.nirreps; h++)
    free_block(C[h]);
  free(C);

  free(pitz2dpd_occ);
  free(pitz2dpd_vir);
}
