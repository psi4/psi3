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

/* Max length of ioff array */
#define IOFF_MAX 32641

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void AO_contribute(int p, int q, int r, int s, double value, double ***tau1_AO, 
		   double ***tau2_AO, int *orbspi, int *orbsym, int *rowtot, int *ioff);

void BT2_AO(void)
{
  int h, nirreps, i, *ioff, Gc, Gd, Ga, Gb, ij;
  double ***C, ***tau1_AO, ***tau2_AO, **X;
  int *orbspi, *virtpi, *orbsym, *rowtot, *coltot;
  dpdbuf4 tau, t2;
  psio_address next;
  struct iwlbuf InBuf;
  int idx, p, q, r, s, filenum;
  int lastbuf;
  double value, tolerance=1e-14;
  Value *valptr;
  Label *lblptr;

  return;

  ioff = (int *) malloc(IOFF_MAX * sizeof(int));
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;

  nirreps = moinfo.nirreps;
  orbspi = moinfo.orbspi;
  virtpi = moinfo.virtpi;

  C = (double ***) malloc(nirreps * sizeof(double **));
  next = PSIO_ZERO;
  for(h=0; h < nirreps; h++) {
    if(virtpi[h]) {
      C[h] = block_matrix(orbspi[h],virtpi[h]);
      psio_read(CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *) C[h][0],
		orbspi[h]*virtpi[h]*sizeof(double), next, &next);
    }
  }

  orbsym = init_int_array(moinfo.nso);
  for(h=0,q=0; h < nirreps; h++)
    for(p=0; p < orbspi[h]; p++)
      orbsym[q++] = h;

  /*** AA ***/

  dpd_buf4_init(&tau, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tauIJAB");
  rowtot = tau.params->rowtot;

  coltot = init_int_array(nirreps);
  for(h=0; h < nirreps; h++) {
    for(Gc=0; Gc < nirreps; Gc++) {
      Gd = h ^ Gc;
      coltot[h] += orbspi[Gc] * orbspi[Gd];
    }
  }

  tau1_AO = (double ***) malloc(nirreps * sizeof(double **));
  tau2_AO = (double ***) malloc(nirreps * sizeof(double **));

  for(h=0; h < nirreps; h++) {

    dpd_buf4_mat_irrep_init(&tau, h);
    dpd_buf4_mat_irrep_rd(&tau, h);

    tau1_AO[h] = dpd_block_matrix(rowtot[h],coltot[h]);
    tau2_AO[h] = dpd_block_matrix(rowtot[h],coltot[h]);

    for(Gc=0; Gc < nirreps; Gc++) {
      Gd = h^Gc;

      if(virtpi[Gc] && virtpi[Gd] && orbspi[Gc] && orbspi[Gd]) {
	X = block_matrix(orbspi[Gc],virtpi[Gd]);

	for(ij=0; ij < tau.params->rowtot[h]; ij++) {

	  C_DGEMM('n','n', orbspi[Gc], virtpi[Gd], virtpi[Gc], 1.0,
		  &(C[Gc][0][0]), virtpi[Gc], &(tau.matrix[h][ij][0]), virtpi[Gd],
		  0.0, &(X[0][0]), virtpi[Gd]);

	  C_DGEMM('n','t', orbspi[Gc], orbspi[Gd], virtpi[Gd], 1.0, 
		  &(X[0][0]), virtpi[Gd], &(C[Gd][0][0]), virtpi[Gd], 
		  0.0, &(tau1_AO[h][ij][0]), orbspi[Gd]);

	}

	free_block(X);
      }
    }

    dpd_buf4_mat_irrep_close(&tau, h);
  }

  dpd_buf4_close(&tau);

  iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

  lblptr = InBuf.labels;
  valptr = InBuf.values;
  lastbuf = InBuf.lastbuf;

  for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
    p = abs((int) lblptr[idx++]);
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    value = (double) valptr[InBuf.idx];

    /*    fprintf(outfile, "<%d %d %d %d = %20.10lf\n", p, q, r, s, value); */

    AO_contribute(p, q, r, s, value, tau1_AO, tau2_AO, orbspi, orbsym, rowtot, ioff);

  }
  while(!lastbuf) {
    iwl_buf_fetch(&InBuf);
    lastbuf = InBuf.lastbuf;
    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      /*      fprintf(outfile, "<%d %d %d %d = %20.10lf\n", p, q, r, s, value); */

      AO_contribute(p, q, r, s, value, tau1_AO, tau2_AO, orbspi, orbsym, rowtot, ioff);

    }
  }

  iwl_buf_close(&InBuf, 1);

  dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");

  for(h=0; h < nirreps; h++) {

    dpd_buf4_mat_irrep_init(&t2, h);
    dpd_buf4_mat_irrep_rd(&t2, h);

    for(Ga=0; Ga < nirreps; Ga++) {
      Gb = h^Ga;

      if(virtpi[Ga] && virtpi[Gb] && orbspi[Ga] && orbspi[Gb]) {
	X = block_matrix(orbspi[Ga],virtpi[Gb]);

	for(ij=0; ij < t2.params->rowtot[h]; ij++) {

	  C_DGEMM('n','n', orbspi[Ga], virtpi[Gb], orbspi[Gb], 1.0,
		  &(tau2_AO[h][ij][0]), orbspi[Gb], &(C[Gb][0][0]), virtpi[Gb],
		  0.0, &(X[0][0]), virtpi[Gb]);

	  C_DGEMM('t','n', virtpi[Ga], virtpi[Gb], orbspi[Ga], 1.0, 
		  &(C[Ga][0][0]), virtpi[Ga], &(X[0][0]), virtpi[Gb],
		  1.0, &(t2.matrix[h][ij][0]), virtpi[Gb]);

	}
	free_block(X);
      }
    }

    dpd_free_block(tau1_AO[h], rowtot[h], coltot[h]);
    dpd_free_block(tau2_AO[h], rowtot[h], coltot[h]);

    dpd_buf4_mat_irrep_wrt(&t2, h);
    dpd_buf4_mat_irrep_close(&t2, h);
  }

  dpd_buf4_close(&t2);

  for(h=0; h < nirreps; h++) 
    if(virtpi[h]) free_block(C[h]);

  free(tau1_AO);
  free(tau2_AO);
  free(coltot);
  free(ioff);
}

void AO_contribute(int p, int q, int r, int s, double value, double ***tau1_AO, 
		   double ***tau2_AO, int *orbspi, int *orbsym, int *rowtot, int *ioff)
{
  int Gp, Gq, Gr, Gs, Gpr, Gps, Gqr, Gqs, Grp, Gsp, Grq, Gsq;
  int pr, ps, qr, qs, rp, rq, sp, sq, pq, rs;
  int row;

  Gp = orbsym[p]; Gq = orbsym[q]; Gr = orbsym[r]; Gs = orbsym[s];

  pq = INDEX(p,q);  rs = INDEX(r,s);

  if(p!=q && r!=s) {

    /* (pq|rs) */
    Gpr = Gp ^ Gr;
    pr = p * orbspi[Gr] + r;
    qs = q * orbspi[Gs] + s;
    sq = s * orbspi[Gq] + q;

    for(row=0; row < rowtot[Gpr]; row++) {
      tau2_AO[Gpr][row][pr] += value * tau1_AO[Gpr][row][qs];
      tau2_AO[Gpr][row][pr] -= value * tau1_AO[Gpr][row][sq];
    }

    /* (pq|sr) */
    Gps = Gp ^ Gs;
    ps = p * orbspi[Gs] + s;
    qr = q * orbspi[Gr] + r;
    rq = r * orbspi[Gq] + q;

    for(row=0; row < rowtot[Gps]; row++) {
      tau2_AO[Gps][row][ps] += value * tau1_AO[Gps][row][qr];
      tau2_AO[Gps][row][ps] -= value * tau1_AO[Gps][row][rq];
    }

    /* (qp|rs) */
    Gqr = Gq ^ Gr;
    qr = q * orbspi[Gr] + r;
    ps = p * orbspi[Gs] + s;
    sp = s * orbspi[Gp] + p;

    for(row=0; row < rowtot[Gqr]; row++) {
      tau2_AO[Gqr][row][qr] += value * tau1_AO[Gqr][row][ps];
      tau2_AO[Gqr][row][qr] -= value * tau1_AO[Gqr][row][sp];
    }

    /* (qp|sr) */
    Gqs = Gq ^ Gs;
    qs = q * orbspi[Gs] + s;
    pr = p * orbspi[Gr] + r;
    rp = r * orbspi[Gp] + p;

    for(row=0; row < rowtot[Gqs]; row++) {
      tau2_AO[Gqs][row][qs] += value * tau1_AO[Gqs][row][pr];
      tau2_AO[Gqs][row][qs] -= value * tau1_AO[Gqs][row][rp];
    }

    if(pq != rs) {
      /* (rs|pq) */
      Grp = Gp ^ Gr;
      rp = r * orbspi[Gp] + p;
      sq = s * orbspi[Gq] + q;
      qs = q * orbspi[Gs] + s;

      for(row=0; row < rowtot[Grp]; row++) {
	tau2_AO[Grp][row][rp] += value * tau1_AO[Grp][row][sq];
	tau2_AO[Grp][row][rp] -= value * tau1_AO[Grp][row][qs];
      }

      /* (sr|pq) */
      Gsp = Gp ^ Gs;
      sp = s * orbspi[Gp] + p;
      qr = q * orbspi[Gr] + r;
      rq = r * orbspi[Gq] + q;

      for(row=0; row < rowtot[Gsp]; row++) {
	tau2_AO[Gsp][row][sp] += value * tau1_AO[Gsp][row][rq];
	tau2_AO[Gsp][row][sp] -= value * tau1_AO[Gsp][row][qr];
      }

      /* (rs|qp) */
      Grq = Gq ^ Gr;
      rq = r * orbspi[Gq] + q;
      ps = p * orbspi[Gs] + s;
      sp = s * orbspi[Gp] + p;

      for(row=0; row < rowtot[Grq]; row++) {
	tau2_AO[Grq][row][rq] += value * tau1_AO[Grq][row][sp];
	tau2_AO[Grq][row][rq] -= value * tau1_AO[Grq][row][ps];
      }

      /* (sr|qp) */
      Gsq = Gq ^ Gs;
      sq = s * orbspi[Gq] + q;
      pr = p * orbspi[Gr] + r;
      rp = r * orbspi[Gp] + p;

      for(row=0; row < rowtot[Gsq]; row++) {
	tau2_AO[Gsq][row][sq] += value * tau1_AO[Gsq][row][rp];
	tau2_AO[Gsq][row][sq] -= value * tau1_AO[Gsq][row][pr];
      }
    }

  }
  else if(p!=q && r==s) {

    /* (pq|rs) */
    Gpr = Gp ^ Gr;
    pr = p * orbspi[Gr] + r;
    qs = q * orbspi[Gs] + s;
    sq = s * orbspi[Gq] + q;

    for(row=0; row < rowtot[Gpr]; row++) {
      tau2_AO[Gpr][row][pr] += value * tau1_AO[Gpr][row][qs];
      tau2_AO[Gpr][row][pr] -= value * tau1_AO[Gpr][row][sq];
    }

    /* (qp|rs) */
    Gqr = Gq ^ Gr;
    qr = q * orbspi[Gr] + r;
    ps = p * orbspi[Gs] + s;
    sp = s * orbspi[Gp] + p;

    for(row=0; row < rowtot[Gqr]; row++) {
      tau2_AO[Gqr][row][qr] += value * tau1_AO[Gqr][row][ps];
      tau2_AO[Gqr][row][qr] -= value * tau1_AO[Gqr][row][sp];
    }

    if(pq != rs) {
      /* (rs|pq) */
      Grp = Gp ^ Gr;
      rp = r * orbspi[Gp] + p;
      qs = q * orbspi[Gs] + s;
      sq = s * orbspi[Gq] + q;

      for(row=0; row < rowtot[Grp]; row++) {
	tau2_AO[Grp][row][rp] += value * tau1_AO[Grp][row][sq];
	tau2_AO[Grp][row][rp] -= value * tau1_AO[Grp][row][qs];
      }

      /* (rs|qp) */
      Grq = Gq ^ Gr;
      rq = r * orbspi[Gq] + q;
      ps = p * orbspi[Gs] + s;
      sp = s * orbspi[Gp] + p;

      for(row=0; row < rowtot[Grq]; row++) {
	tau2_AO[Grq][row][rq] += value * tau1_AO[Grq][row][sp];
	tau2_AO[Grq][row][rq] -= value * tau1_AO[Grq][row][ps];
      }
    }

  }

  else if(p==q && r!=s) {

    /* (pq|rs) */
    Gpr = Gp ^ Gr;
    pr = p * orbspi[Gr] + r;
    qs = q * orbspi[Gs] + s;
    sq = s * orbspi[Gq] + q;

    for(row=0; row < rowtot[Gpr]; row++) {
      tau2_AO[Gpr][row][pr] += value * tau1_AO[Gpr][row][qs];
      tau2_AO[Gpr][row][pr] -= value * tau1_AO[Gpr][row][sq];
    }

    /* (pq|sr) */
    Gps = Gp ^ Gs;
    ps = p * orbspi[Gs] + s;
    qr = q * orbspi[Gr] + r;
    rq = r * orbspi[Gq] + q;

    for(row=0; row < rowtot[Gps]; row++) {
      tau2_AO[Gps][row][ps] += value * tau1_AO[Gps][row][qr];
      tau2_AO[Gps][row][ps] -= value * tau1_AO[Gps][row][rq];
    }

    if(pq != rs) {
      /* (rs|pq) */
      Grp = Gp ^ Gr;
      rp = r * orbspi[Gp] + p;
      qs = q * orbspi[Gs] + s;
      sq = s * orbspi[Gq] + q;

      for(row=0; row < rowtot[Grp]; row++) {
	tau2_AO[Grp][row][rp] += value * tau1_AO[Grp][row][sq];
	tau2_AO[Grp][row][rp] -= value * tau1_AO[Grp][row][qs];
      }

      /* (sr|pq) */
      Gsp = Gp ^ Gs;
      sp = s * orbspi[Gp] + p;
      qr = q * orbspi[Gr] + r;
      rq = r * orbspi[Gq] + q;

      for(row=0; row < rowtot[Gsp]; row++) {
	tau2_AO[Gsp][row][sp] += value * tau1_AO[Gsp][row][rq];
	tau2_AO[Gsp][row][sp] -= value * tau1_AO[Gsp][row][qr];
      }
    }

  }
  return;
}



