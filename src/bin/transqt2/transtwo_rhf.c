#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void transtwo_rhf(void)
{
  int nirreps, nso, nmo;
  double ***C, **scratch, **TMP;
  int h, p, q, r, s, pq, rs, Gs, Gr, PQ, RS;
  int nrows, ncols, nlinks;
  unsigned long int memfree, rows_per_bucket, rows_left; 
  int nbuckets, n;
  dpdbuf4 J, K;
  struct iwlbuf MBuff;
  int stat;

  nirreps = moinfo.nirreps;
  nso = moinfo.nso;
  nmo = moinfo.nmo;

  /* grab MOs and remove frozen core/virt */
  C = (double ***) malloc(nirreps * sizeof(double **));
  chkpt_init(PSIO_OPEN_OLD);
  for(h=0; h < nirreps; h++) {
    scratch = chkpt_rd_scf_irrep(h);
    C[h] = block_matrix(moinfo.sopi[h],moinfo.mopi[h]);
    for(q=0; q < moinfo.mopi[h]; q++)
      for(p=0; p < moinfo.sopi[h]; p++)
	C[h][p][q] = scratch[p][q+moinfo.frdocc[h]];
    if(params.print_lvl > 2) {
      fprintf(outfile, "\n\tMOs for irrep %d:\n",h);
      mat_print(C[h], moinfo.sopi[h], moinfo.mopi[h], outfile);
    }
    free_block(scratch);
  }
  chkpt_close();

  TMP = block_matrix(nso,nso);

  if(params.print_lvl) {
    fprintf(outfile, "\tStarting first half-transformation.\n");
    fflush(outfile);
  }

  psio_open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
  psio_open(PSIF_HALFT0, PSIO_OPEN_NEW);

  dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (pq,rs)");
  dpd_buf4_init(&K, PSIF_HALFT0, 0, 3, 5, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  for(h=0; h < nirreps; h++) {

    memfree = (unsigned long int) (dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
    rows_per_bucket = memfree/(2 * J.params->coltot[h]);
    if(rows_per_bucket > J.params->rowtot[h]) rows_per_bucket = (unsigned long int) J.params->rowtot[h];
    nbuckets = (int) ceil(((double) J.params->rowtot[h])/((double) rows_per_bucket));
    rows_left = (unsigned long int) (J.params->rowtot[h] % rows_per_bucket);
    if(params.print_lvl > 1) {
      fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memfree);
      fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rows_per_bucket);
      fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rows_left);
      fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nbuckets);
      fflush(outfile);
    }

    dpd_buf4_mat_irrep_init_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_init_block(&K, h, rows_per_bucket);

    for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {
      dpd_buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, rows_per_bucket);
      for(pq=0; pq < rows_per_bucket; pq++) {
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = h^Gr;
	  nrows = moinfo.sopi[Gr];
	  ncols = moinfo.mopi[Gs];
	  nlinks = moinfo.sopi[Gs];
	  rs = J.col_offset[h][Gr];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][pq][rs],nlinks,
		    C[Gs][0],ncols,0.0,TMP[0],nso);

	  nrows = moinfo.mopi[Gr];
	  ncols = moinfo.mopi[Gs];
	  nlinks = moinfo.sopi[Gr];
	  rs = K.col_offset[h][Gr];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[Gr][0],nrows,TMP[0],nso,
		    0.0,&K.matrix[h][pq][rs],ncols);
	} /* Gr */
      } /* pq */
      dpd_buf4_mat_irrep_wrt_block(&K, h, n*rows_per_bucket, rows_per_bucket);
    }
    if(rows_left) {
      dpd_buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, rows_left);
      for(pq=0; pq < rows_left; pq++) {
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = h^Gr;

	  nrows = moinfo.sopi[Gr];
	  ncols = moinfo.mopi[Gs];
	  nlinks = moinfo.sopi[Gs];
	  rs = J.col_offset[h][Gr];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][pq][rs],nlinks,
		    C[Gs][0],ncols,0.0,TMP[0],nso);

	  nrows = moinfo.mopi[Gr];
	  ncols = moinfo.mopi[Gs];
	  nlinks = moinfo.sopi[Gr];
	  rs = K.col_offset[h][Gr];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[Gr][0],nrows,TMP[0],nso,
		    0.0,&K.matrix[h][pq][rs],ncols);
	} /* Gr */
      } /* pq */

      dpd_buf4_mat_irrep_wrt_block(&K, h, n*rows_per_bucket, rows_left);
    }

    dpd_buf4_mat_irrep_close_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_close_block(&K, h, rows_per_bucket);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  psio_close(PSIF_SO_PRESORT, 0);

  if(params.print_lvl) {
    fprintf(outfile, "\tSorting half-transformed integrals.\n");
    fflush(outfile);
  }

  psio_open(PSIF_HALFT1, PSIO_OPEN_NEW);

  dpd_buf4_init(&K, PSIF_HALFT0, 0, 3, 8, 3, 8, 0, "Half-Transformed Ints (pq,ij)");
  dpd_buf4_sort(&K, PSIF_HALFT1, rspq, 8, 3, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_close(&K);

  psio_close(PSIF_HALFT0, 0);

  if(params.print_lvl) {
    fprintf(outfile, "\tStarting second half-transformation.\n");
    fflush(outfile);
  }
  iwl_buf_init(&MBuff, PSIF_MO_TEI, params.tolerance, 0, 0);

  dpd_buf4_init(&J, PSIF_HALFT1, 0, 8, 0, 8, 3, 0, "Half-Transformed Ints (ij,pq)");
  dpd_buf4_init(&K, CC_MISC, 0, 8, 5, 8, 8, 0, "MO Ints (ij,kl)");
  for(h=0; h < nirreps; h++) {

    memfree = (unsigned long int) (dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
    rows_per_bucket = memfree/(2 * J.params->coltot[h]);
    if(rows_per_bucket > J.params->rowtot[h]) rows_per_bucket = (unsigned long int) J.params->rowtot[h];
    nbuckets = (int) ceil(((double) J.params->rowtot[h])/((double) rows_per_bucket));
    rows_left = (unsigned long int) (J.params->rowtot[h] % rows_per_bucket);

    if(params.print_lvl > 1) {
      fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memfree);
      fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rows_per_bucket);
      fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rows_left);
      fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nbuckets);
      fflush(outfile);
    }

    dpd_buf4_mat_irrep_init_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_init_block(&K, h, rows_per_bucket);

    for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {
      dpd_buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, rows_per_bucket);
      for(pq=0; pq < rows_per_bucket; pq++) {
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = h^Gr;
	  nrows = moinfo.sopi[Gr];
	  ncols = moinfo.mopi[Gs];
	  nlinks = moinfo.sopi[Gs];
	  rs = J.col_offset[h][Gr];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][pq][rs],nlinks,
		    C[Gs][0],ncols,0.0,TMP[0],nso);

	  nrows = moinfo.mopi[Gr];
	  ncols = moinfo.mopi[Gs];
	  nlinks = moinfo.sopi[Gr];
	  rs = K.col_offset[h][Gr];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[Gr][0],nrows,TMP[0],nso,
		    0.0,&K.matrix[h][pq][rs],ncols);
	} /* Gr */

	p = moinfo.act2qt[K.params->roworb[h][pq+n*rows_per_bucket][0]];
	q = moinfo.act2qt[K.params->roworb[h][pq+n*rows_per_bucket][1]];
	PQ = INDEX(p,q);
	for(rs=0; rs < K.params->coltot[h]; rs++) {
	  r = moinfo.act2qt[K.params->colorb[h][rs][0]];
	  s = moinfo.act2qt[K.params->colorb[h][rs][1]];
	  RS = INDEX(r,s);
	  if(r >= s && RS <= PQ)
	    iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][pq][rs], params.print_tei, outfile, 0);
	} /* rs */
      } /* pq */
    }
    if(rows_left) {
      dpd_buf4_mat_irrep_rd_block(&J, h, n*rows_per_bucket, rows_left);
      for(pq=0; pq < rows_left; pq++) {
	for(Gr=0; Gr < nirreps; Gr++) {
	  Gs = h^Gr;
	  nrows = moinfo.sopi[Gr];
	  ncols = moinfo.mopi[Gs];
	  nlinks = moinfo.sopi[Gs];
	  rs = J.col_offset[h][Gr];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('n','n',nrows,ncols,nlinks,1.0,&J.matrix[h][pq][rs],nlinks,
		    C[Gs][0],ncols,0.0,TMP[0],nso);

	  nrows = moinfo.mopi[Gr];
	  ncols = moinfo.mopi[Gs];
	  nlinks = moinfo.sopi[Gr];
	  rs = K.col_offset[h][Gr];
	  if(nrows && ncols && nlinks)
	    C_DGEMM('t','n',nrows,ncols,nlinks,1.0,C[Gr][0],nrows,TMP[0],nso,
		    0.0,&K.matrix[h][pq][rs],ncols);
	} /* Gr */

	p = moinfo.act2qt[K.params->roworb[h][pq+n*rows_per_bucket][0]];
	q = moinfo.act2qt[K.params->roworb[h][pq+n*rows_per_bucket][1]];
	PQ = INDEX(p,q);
	for(rs=0; rs < K.params->coltot[h]; rs++) {
	  r = moinfo.act2qt[K.params->colorb[h][rs][0]];
	  s = moinfo.act2qt[K.params->colorb[h][rs][1]];
	  RS = INDEX(r,s);
	  if(r >= s && RS <= PQ)
	    iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][pq][rs], params.print_tei, outfile, 0);
	} /* rs */
      } /* pq */
    }
    dpd_buf4_mat_irrep_close_block(&J, h, rows_per_bucket);
    dpd_buf4_mat_irrep_close_block(&K, h, rows_per_bucket);
  }
  dpd_buf4_close(&K);
  dpd_buf4_close(&J);

  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  free_block(TMP);

  for(h=0; h < nirreps; h++)
    free_block(C[h]);
  free(C);

  psio_close(PSIF_HALFT1, 0);

  if(params.print_lvl) {
    fprintf(outfile, "\tTwo-electron integral transformation complete.\n");
    fflush(outfile);
  }
}