#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include <psio.h>
#include <iwl.h>
#include <dpd.h>
#include <qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

void halftrans(dpdbuf4 *Buf1, int dpdnum1, dpdbuf4 *Buf2, int dpdnum2, double ***C, int nirreps, 
	       int **mo_row, int **so_row, int *mospi, int *sospi, int type, double alpha, double beta);

void AO_contribute(struct iwlbuf *InBuf, dpdbuf4 *tau1_AO, dpdbuf4 *tau2_AO);

void BT2_AO(void)
{
  int h, nirreps, i, Gc, Gd, Ga, Gb, ij;
  double ***C, **X;
  int *orbspi, *virtpi;
  int **T2_cd_row_start, **T2_pq_row_start, offset, cd, pq;
  dpdbuf4 tau, t2, tau1_AO, tau2_AO;
  psio_address next;
  struct iwlbuf InBuf;
  int lastbuf;
  double tolerance=1e-14;
  double **integrals;
  int **tau1_cols, **tau2_cols, *num_ints;

  nirreps = moinfo.nirreps;
  orbspi = moinfo.orbspi;
  virtpi = moinfo.virtpi;
  C = moinfo.C;

  T2_cd_row_start = init_int_matrix(nirreps,nirreps);
  for(h=0; h < nirreps; h++) {
    for(Gc=0,offset=0; Gc < nirreps; Gc++) {
      Gd = Gc ^ h;
      T2_cd_row_start[h][Gc] = offset;
      offset += virtpi[Gc] * virtpi[Gd];
    }
  }

  T2_pq_row_start = init_int_matrix(nirreps,nirreps);
  for(h=0; h < nirreps; h++) {
    for(Gc=0,offset=0; Gc < nirreps; Gc++) {
      Gd = Gc ^ h;
      T2_pq_row_start[h][Gc] = offset;
      offset += orbspi[Gc] * orbspi[Gd];
    }
  }

  if(params.ref == 0) { /** RHF **/
    dpd_set_default(1);
    dpd_buf4_init(&tau1_AO, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (1)");
    dpd_buf4_scm(&tau1_AO, 0.0);

    dpd_set_default(0);
    dpd_buf4_init(&tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

    halftrans(&tau, 0, &tau1_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start, 
	      virtpi, orbspi, 0, 1.0, 0.0);

    dpd_buf4_close(&tau);
    dpd_buf4_close(&tau1_AO);


    /* Transpose tau1_AO for better memory access patterns */
    dpd_set_default(1);
    dpd_buf4_init(&tau1_AO, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (1)");
    dpd_buf4_sort(&tau1_AO, CC_TMP0, rspq, 5, 0, "tauPqIj (1)");
    dpd_buf4_close(&tau1_AO);

    dpd_buf4_init(&tau1_AO, CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (1)");
    dpd_buf4_init(&tau2_AO, CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (2)");
    dpd_buf4_scm(&tau2_AO, 0.0);

    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&tau1_AO, h);
      dpd_buf4_mat_irrep_rd(&tau1_AO, h);
      dpd_buf4_mat_irrep_init(&tau2_AO, h);
    }

    iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

    lastbuf = InBuf.lastbuf;

    AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

    while(!lastbuf) {
      iwl_buf_fetch(&InBuf);
      lastbuf = InBuf.lastbuf;

      AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
    }

    iwl_buf_close(&InBuf, 1);

    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_wrt(&tau2_AO, h);
      dpd_buf4_mat_irrep_close(&tau2_AO, h);
      dpd_buf4_mat_irrep_close(&tau1_AO, h);
    }
    dpd_buf4_close(&tau1_AO);
    dpd_buf4_close(&tau2_AO);

    /* Transpose tau2_AO for the half-backtransformation */
    dpd_set_default(1);
    dpd_buf4_init(&tau2_AO, CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (2)");
    dpd_buf4_sort(&tau2_AO, CC_TAMPS, rspq, 0, 5, "tauIjPq (2)");
    dpd_buf4_close(&tau2_AO);

    dpd_buf4_init(&tau2_AO, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (2)");

    dpd_set_default(0);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    halftrans(&t2, 0, &tau2_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start, 
	      virtpi, orbspi, 1, 1.0, 1.0);

    dpd_buf4_close(&t2);
    dpd_buf4_close(&tau2_AO);

  }
  else if(params.ref == 1) { /** ROHF **/

    /************************************* AA *****************************************/

    dpd_set_default(1);
    dpd_buf4_init(&tau1_AO, CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauIJPQ (1)");
    dpd_buf4_scm(&tau1_AO, 0.0);

    dpd_set_default(0);
    dpd_buf4_init(&tau, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");

    halftrans(&tau, 0, &tau1_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start, 
	      virtpi, orbspi, 0, 1.0, 0.0);

    dpd_buf4_close(&tau);
    dpd_buf4_close(&tau1_AO);

    /* Transpose tau1_AO for better memory access patterns */
    dpd_set_default(1);
    dpd_buf4_init(&tau1_AO, CC_TAMPS, 0, 2, 5, 2, 5, 1, "tauIJPQ (1)");
    dpd_buf4_sort(&tau1_AO, CC_TMP0, rspq, 5, 2, "tauPQIJ (1)");
    dpd_buf4_close(&tau1_AO);

    dpd_buf4_init(&tau1_AO, CC_TMP0, 0, 5, 2, 5, 2, 0, "tauPQIJ (1)");
    dpd_buf4_init(&tau2_AO, CC_TMP0, 0, 5, 2, 5, 2, 0, "tauPQIJ (2)");
    dpd_buf4_scm(&tau2_AO, 0.0);

    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&tau1_AO, h);
      dpd_buf4_mat_irrep_rd(&tau1_AO, h);
      dpd_buf4_mat_irrep_init(&tau2_AO, h);
    }

    iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

    lastbuf = InBuf.lastbuf;

    AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

    while(!lastbuf) {
      iwl_buf_fetch(&InBuf);
      lastbuf = InBuf.lastbuf;

      AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
    }

    iwl_buf_close(&InBuf, 1);

    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_wrt(&tau2_AO, h);
      dpd_buf4_mat_irrep_close(&tau2_AO, h);
      dpd_buf4_mat_irrep_close(&tau1_AO, h);
    }
    dpd_buf4_close(&tau1_AO);
    dpd_buf4_close(&tau2_AO);


    /* Transpose tau2_AO for the half-backtransformation */
    dpd_set_default(1);
    dpd_buf4_init(&tau2_AO, CC_TMP0, 0, 5, 2, 5, 2, 0, "tauPQIJ (2)");
    dpd_buf4_sort(&tau2_AO, CC_TAMPS, rspq, 2, 5, "tauIJPQ (2)");
    dpd_buf4_close(&tau2_AO);

    dpd_buf4_init(&tau2_AO, CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauIJPQ (2)");

    dpd_set_default(0);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");

    halftrans(&t2, 0, &tau2_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start, 
	      virtpi, orbspi, 1, 0.5, 1.0);

    dpd_buf4_close(&t2);
    dpd_buf4_close(&tau2_AO);

    /************************************* BB *****************************************/

    dpd_set_default(1);
    dpd_buf4_init(&tau1_AO, CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauijpq (1)");
    dpd_buf4_scm(&tau1_AO, 0.0);

    dpd_set_default(0);
    dpd_buf4_init(&tau, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauijab");

    halftrans(&tau, 0, &tau1_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start, 
	      virtpi, orbspi, 0, 1.0, 0.0);

    dpd_buf4_close(&tau);
    dpd_buf4_close(&tau1_AO);

    /* Transpose tau1_AO for better memory access patterns */
    dpd_set_default(1);
    dpd_buf4_init(&tau1_AO, CC_TAMPS, 0, 2, 5, 2, 5, 1, "tauijpq (1)");
    dpd_buf4_sort(&tau1_AO, CC_TMP0, rspq, 5, 2, "taupqij (1)");
    dpd_buf4_close(&tau1_AO);

    dpd_buf4_init(&tau1_AO, CC_TMP0, 0, 5, 2, 5, 2, 0, "taupqij (1)");
    dpd_buf4_init(&tau2_AO, CC_TMP0, 0, 5, 2, 5, 2, 0, "taupqij (2)");
    dpd_buf4_scm(&tau2_AO, 0.0);

    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&tau1_AO, h);
      dpd_buf4_mat_irrep_rd(&tau1_AO, h);
      dpd_buf4_mat_irrep_init(&tau2_AO, h);
    }

    iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

    lastbuf = InBuf.lastbuf;

    AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

    while(!lastbuf) {
      iwl_buf_fetch(&InBuf);
      lastbuf = InBuf.lastbuf;

      AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
    }

    iwl_buf_close(&InBuf, 1);

    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_wrt(&tau2_AO, h);
      dpd_buf4_mat_irrep_close(&tau2_AO, h);
      dpd_buf4_mat_irrep_close(&tau1_AO, h);
    }
    dpd_buf4_close(&tau1_AO);
    dpd_buf4_close(&tau2_AO);


    /* Transpose tau2_AO for the half-backtransformation */
    dpd_set_default(1);
    dpd_buf4_init(&tau2_AO, CC_TMP0, 0, 5, 2, 5, 2, 0, "taupqij (2)");
    dpd_buf4_sort(&tau2_AO, CC_TAMPS, rspq, 2, 5, "tauijpq (2)");
    dpd_buf4_close(&tau2_AO);

    dpd_buf4_init(&tau2_AO, CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauijpq (2)");

    dpd_set_default(0);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");

    halftrans(&t2, 0, &tau2_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start, 
	      virtpi, orbspi, 1, 0.5, 1.0);

    dpd_buf4_close(&t2);
    dpd_buf4_close(&tau2_AO);

    /************************************* AB *****************************************/

    dpd_set_default(1);
    dpd_buf4_init(&tau1_AO, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (1)");
    dpd_buf4_scm(&tau1_AO, 0.0);

    dpd_set_default(0);
    dpd_buf4_init(&tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

    halftrans(&tau, 0, &tau1_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start, 
	      virtpi, orbspi, 0, 1.0, 0.0);

    dpd_buf4_close(&tau);
    dpd_buf4_close(&tau1_AO);


    /* Transpose tau1_AO for better memory access patterns */
    dpd_set_default(1);
    dpd_buf4_init(&tau1_AO, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (1)");
    dpd_buf4_sort(&tau1_AO, CC_TMP0, rspq, 5, 0, "tauPqIj (1)");
    dpd_buf4_close(&tau1_AO);

    dpd_buf4_init(&tau1_AO, CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (1)");
    dpd_buf4_init(&tau2_AO, CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (2)");
    dpd_buf4_scm(&tau2_AO, 0.0);

    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&tau1_AO, h);
      dpd_buf4_mat_irrep_rd(&tau1_AO, h);
      dpd_buf4_mat_irrep_init(&tau2_AO, h);
    }

    iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

    lastbuf = InBuf.lastbuf;

    AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

    while(!lastbuf) {
      iwl_buf_fetch(&InBuf);
      lastbuf = InBuf.lastbuf;

      AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
    }

    iwl_buf_close(&InBuf, 1);

    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_wrt(&tau2_AO, h);
      dpd_buf4_mat_irrep_close(&tau2_AO, h);
      dpd_buf4_mat_irrep_close(&tau1_AO, h);
    }
    dpd_buf4_close(&tau1_AO);
    dpd_buf4_close(&tau2_AO);

    /* Transpose tau2_AO for the half-backtransformation */
    dpd_set_default(1);
    dpd_buf4_init(&tau2_AO, CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (2)");
    dpd_buf4_sort(&tau2_AO, CC_TAMPS, rspq, 0, 5, "tauIjPq (2)");
    dpd_buf4_close(&tau2_AO);

    dpd_buf4_init(&tau2_AO, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (2)");

    dpd_set_default(0);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    halftrans(&t2, 0, &tau2_AO, 1, C, nirreps, T2_cd_row_start, T2_pq_row_start, 
	      virtpi, orbspi, 1, 1.0, 1.0);

    dpd_buf4_close(&t2);
    dpd_buf4_close(&tau2_AO);

  }

  free(T2_cd_row_start);
  free(T2_pq_row_start);

  /* Reset the default dpd back to 0 --- this stuff gets really ugly */
  dpd_set_default(0);
}
