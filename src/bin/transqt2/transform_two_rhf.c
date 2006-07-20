#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include <ccfiles.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void transform_two_rhf(void)
{
  dpdfile4 I;
  dpdbuf4 J, K;
  int h, pq, p, q, rs, r, s, P;
  double **TMP;
  int *sopi, *mopi;
  int *occ, *vir, *cc_occ, *cc_vir;
  int nso, nmo, nrows, ncols, nlinks, J_rs, K_rs, Gs, Gr;
  double ***C;
  double *buffer, value;
  struct iwlbuf ABuf, BBuf, CBuf, DBuf, EBuf, FBuf;
  dpdfile4 A;
  dpdbuf4 Abuf;

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  sopi = moinfo.sopi;
  mopi = moinfo.mopi;
  C = (double ***) malloc(moinfo.nirreps * sizeof(double **));
  chkpt_init(PSIO_OPEN_OLD);
  for(h=0; h < moinfo.nirreps; h++)
    C[h] = chkpt_rd_scf_irrep(h);
  chkpt_close();

  /* Build some sorting arrays specific to this transformation. */
  /* We're sorting from Pitzer to DPD/CC ordering, so we can skip the QT stuff */
  occ = init_int_array(moinfo.nactive);
  vir = init_int_array(moinfo.nactive);
  cc_occ = init_int_array(moinfo.nactive);
  cc_vir = init_int_array(moinfo.nactive);
  for(p=0; p < moinfo.nactive; p++) {  /* Pitzer index */
    P = moinfo.pitzer2qt[p];
    occ[p] = moinfo.occ[P];
    vir[p] = moinfo.vir[P];
    cc_occ[p] = moinfo.cc_occ[P];
    cc_vir[p] = moinfo.cc_vir[P];
  }

  psio_open(PSIF_SO_PRESORT, 0);

  dpd_file4_init(&I, PSIF_SO_PRESORT, 0, 3, 3, "SO Ints (pq,rs)");
  file_build_presort(&I, PSIF_SO_TEI, params.tolerance, 1);
  dpd_file4_close(&I);

  TMP = block_matrix(nso,nso);

  psio_open(CC_TEI_HALFT0, 0);

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

  psio_close(PSIF_SO_PRESORT, 0);
  psio_open(CC_TEI_HALFT1, 0);

  dpd_buf4_init(&K, CC_TEI_HALFT0, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  dpd_buf4_sort(&K, CC_TEI_HALFT1, rspq, 8, 3, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_close(&K);

  psio_close(CC_TEI_HALFT0, 0); 

  /* prepare target IWL buffers */
  iwl_buf_init(&ABuf, NEWTMP, params.tolerance, 0, 0);
  iwl_buf_init(&BBuf, NEWTMP+1, params.tolerance, 0, 0);
  iwl_buf_init(&CBuf, NEWTMP+2, params.tolerance, 0, 0);
  iwl_buf_init(&DBuf, NEWTMP+3, params.tolerance, 0, 0);
  iwl_buf_init(&EBuf, NEWTMP+4, params.tolerance, 0, 0);
  iwl_buf_init(&FBuf, NEWTMP+5, params.tolerance, 0, 0);

  dpd_buf4_init(&J, CC_TEI_HALFT1, 0, 8, 0, 8, 3, 0, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_init(&K, CC_MISC, 0, 8, 5, 8, 8, 0, "MO Ints (ij,kl)"); /* place-holder buffer */
  for(h=0; h < moinfo.nirreps; h++) {
    dpd_buf4_mat_irrep_row_init(&J, h);
    dpd_buf4_mat_irrep_row_init(&K, h); 
    buffer = init_array(K.params->coltot[h]);

    for(pq=0; pq < J.params->rowtot[h]; pq++) {
      p = K.params->roworb[h][pq][0];
      q = K.params->roworb[h][pq][1];

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

      for(rs=0; rs < K.params->coltot[h]; rs++) {
	r = K.params->colorb[h][rs][0];
	s = K.params->colorb[h][rs][1];

	value = K.matrix[h][0][rs];

	if(occ[p] && occ[q] && occ[r] && occ[s])
	  iwl_buf_wrt_val(&ABuf, cc_occ[p], cc_occ[q], cc_occ[r], cc_occ[s], 
			  value, 1, outfile, 1);
	else if(vir[p] && vir[q] && vir[r] && vir[s])
	  iwl_buf_wrt_val(&BBuf, cc_vir[p], cc_vir[q], cc_vir[r], cc_vir[s], 
			  value, 0, outfile, 1);
	else if(occ[p] && occ[q] && vir[r] && vir[s])
	  iwl_buf_wrt_val(&CBuf, cc_occ[p], cc_occ[q], cc_vir[r], cc_vir[s],
			  value, 0, outfile, 1);
	/* 	else if(occ[p] && vir[q] && occ[r] && vir[s]) { */
	/* 	  iwl_buf_wrt_val(&CBuf, cc_occ[p], cc_vir[q], cc_occ[r], cc_vir[s]); */
	/* 	  iwl_buf_wrt_val(&CBuf, cc_occ[q], cc_vir[p], cc_occ[r], cc_vir[s]); */
	/* 	  iwl_buf_wrt_val(&CBuf, cc_occ[p], cc_vir[q], cc_occ[s], cc_vir[r]); */
	/* 	  iwl_buf_wrt_val(&CBuf, cc_occ[q], cc_vir[p], cc_occ[s], cc_vir[r]); */
	/* 	} */
	/* 	else if(occ[p] && occ[q] && occ[r] && vir[s]) { */
	/* 	} */
	/* 	else if(occ[p] && ) { */
	/* 	} */
      }
    }

    dpd_buf4_mat_irrep_row_close(&J, h);
    dpd_buf4_mat_irrep_row_close(&K, h);
    free(buffer);
  }
  dpd_buf4_print(&K, outfile, 1);
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  iwl_buf_flush(&ABuf, 1);
  iwl_buf_flush(&BBuf, 1);
  iwl_buf_flush(&CBuf, 1);
  iwl_buf_flush(&DBuf, 1);
  iwl_buf_flush(&EBuf, 1);
  iwl_buf_flush(&FBuf, 1);

  iwl_buf_close(&ABuf, 1);
  iwl_buf_close(&BBuf, 1);
  iwl_buf_close(&CBuf, 1);
  iwl_buf_close(&DBuf, 1);
  iwl_buf_close(&EBuf, 1);
  iwl_buf_close(&FBuf, 1);

  psio_close(CC_TEI_HALFT1, 0);

  free(cc_occ); free(cc_vir);
  free(occ); free(vir);

  free_block(TMP);
  for(h=0; h < moinfo.nirreps; h++)
    free_block(C[h]);
  free(C);

  dpd_file4_init_nocache(&Afile, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  file_build(&Afile, NEWTMP, params.tolerance, 1, 1, 1, 0);
  dpd_file4_close(&Afile);

  dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
  dpd_buf4_print(&A, outfile, 1);
  dpd_buf4_close(&A);
}
