/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void rhf_fock_build(double **fock, double  **D)
{
  int i, j, ij;
  int nso, ntri, stat;
  double *scratch;
  int lastbuf, idx, p, q, r, s, pq, rs;
  double value;
  Value *valptr;
  Label *lblptr;
  struct iwlbuf InBuf;

  nso = moinfo.nso;
  ntri = nso * (nso+1)/2;

  /* one-electron contributions */
  scratch = init_array(ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_T, scratch, ntri, 0, 0, outfile);
  for(i=0, ij=0; i < nso; i++)
    for(j=0; j <= i; j++, ij++)
      fock[i][j] = fock[j][i] = scratch[ij];
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_V, scratch, ntri, 0, 0, outfile);
  for(i=0, ij=0; i < nso; i++)
    for(j=0; j <= i; j++, ij++) {
      fock[i][j] += scratch[ij];
      if(i!=j) fock[j][i] += scratch[ij];
    }
  free(scratch);

  iwl_buf_init(&InBuf, PSIF_SO_TEI, 0.0, 1, 1);
  do {

    lastbuf = InBuf.lastbuf;
    lblptr = InBuf.labels;
    valptr = InBuf.values;

    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];
      value = (double) valptr[InBuf.idx];

      pq = INDEX(p,q);
      rs = INDEX(r,s);

      /*
      fprintf(outfile, "%d %d %d %d [%d] [%d] %20.15f\n", p, q, r, s, pq, rs, value);
      */

      /* (pq|rs) */
      fock[p][q] += 2.0 * D[r][s] * value;
      fock[p][r] -= D[q][s] * value;

      if(p!=q && r!=s && pq != rs) {

	/* (pq|sr) */
        fock[p][q] += 2.0 * D[s][r] * value;
	fock[p][s] -= D[q][r] * value;

	/* (qp|rs) */
	fock[q][p] += 2.0 * D[r][s] * value;
	fock[q][r] -= D[p][s] * value;

	/* (qp|sr) */
	fock[q][p] += 2.0 * D[s][r] * value;
        fock[q][s] -= D[p][r] * value;

	/* (rs|pq) */
	fock[r][s] += 2.0 * D[p][q] * value;
	fock[r][p] -= D[s][q] * value;

	/* (rs|qp) */
	fock[r][s] += 2.0 * D[q][p] * value;
	fock[r][q] -= D[s][p] * value;

	/* (sr|pq) */
	fock[s][r] += 2.0 * D[p][q] * value;
	fock[s][p] -= D[r][q] * value;

	/* (sr|qp) */
	fock[s][r] += 2.0 * D[q][p] * value;
	fock[s][q] -= D[r][p] * value;
      }
      else if(p!=q && r!=s && pq==rs) {

        /* (pq|sr) */
        fock[p][q] += 2.0 * D[s][r] * value;
        fock[p][s] -= D[q][r] * value;

        /* (qp|rs) */
        fock[q][p] += 2.0 * D[r][s] * value;
        fock[q][r] -= D[p][s] * value;

        /* (qp|sr) */
        fock[q][p] += 2.0 * D[s][r] * value;
        fock[q][s] -= D[p][r] * value;

      }
      else if(p!=q && r==s) {

        /* (qp|rs) */
        fock[q][p] += 2.0 * D[r][s] * value;
        fock[q][r] -= D[p][s] * value;

        /* (rs|pq) */
        fock[r][s] += 2.0 * D[p][q] * value;
        fock[r][p] -= D[s][q] * value;

        /* (rs|qp) */
        fock[r][s] += 2.0 * D[q][p] * value;
        fock[r][q] -= D[s][p] * value;

      }
      else if(p==q && r!=s) {

        /* (pq|sr) */
        fock[p][q] += 2.0 * D[s][r] * value;
        fock[p][s] -= D[q][r] * value;

        /* (rs|pq) */
        fock[r][s] += 2.0 * D[p][q] * value;
        fock[r][p] -= D[s][q] * value;

        /* (sr|pq) */
        fock[s][r] += 2.0 * D[p][q] * value;
        fock[s][p] -= D[r][q] * value;

      }
      else if(p==q && r==s && pq!=rs) {

        /* (rs|pq) */
        fock[r][s] += 2.0 * D[p][q] * value;
        fock[r][p] -= D[s][q] * value;

      }
    }

    if(!lastbuf) iwl_buf_fetch(&InBuf);

  } while (!lastbuf);
  iwl_buf_close(&InBuf, 1);
}

void uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b)
{
  int i, j, ij;
  int nso, ntri, ntei, stat;
  double *scratch;
  int lastbuf, idx, p, q, r, s, pq, rs;
  double value;
  Value *valptr;
  Label *lblptr;
  struct iwlbuf InBuf;
  double **Dt;

  nso = moinfo.nso;
  ntri = nso * (nso+1)/2;
  ntei = ntri * (ntri+1)/2;

  Dt = block_matrix(nso, nso);
  for(p=0; p < nso; p++)
    for(q=0; q < nso; q++)
      Dt[p][q] = D_a[p][q] + D_b[p][q];

  /* one-electron contributions */
  scratch = init_array(ntri);
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_T, scratch, ntri, 0, 0, outfile);
  for(i=0, ij=0; i < nso; i++)
    for(j=0; j <= i; j++, ij++) {
      fock_a[i][j] = fock_a[j][i] = scratch[ij];
      fock_b[i][j] = fock_b[j][i] = scratch[ij];
    }
  stat = iwl_rdone(PSIF_OEI, PSIF_SO_V, scratch, ntri, 0, 0, outfile);
  for(i=0, ij=0; i < nso; i++)
    for(j=0; j <= i; j++, ij++) {
      fock_a[i][j] += scratch[ij];
      if(i!=j) fock_a[j][i] += scratch[ij];
      fock_b[i][j] += scratch[ij];
      if(i!=j) fock_b[j][i] += scratch[ij];
    }
  free(scratch);

  iwl_buf_init(&InBuf, PSIF_SO_TEI, 0.0, 1, 1);
  do {

    lastbuf = InBuf.lastbuf;
    lblptr = InBuf.labels;
    valptr = InBuf.values;

    for(idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = abs((int) lblptr[idx++]);
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];
      value = (double) valptr[InBuf.idx];

      pq = INDEX(p,q);
      rs = INDEX(r,s);

      /*
      fprintf(outfile, "%d %d %d %d [%d] [%d] %20.15f\n", p, q, r, s, pq, rs, value);
      */

      /* (pq|rs) */
      fock_a[p][q] += Dt[r][s] * value;
      fock_a[p][r] -= D_a[q][s] * value;
      fock_b[p][q] += Dt[r][s] * value;
      fock_b[p][r] -= D_b[q][s] * value;

      if(p!=q && r!=s && pq != rs) {

	/* (pq|sr) */
        fock_a[p][q] += Dt[s][r] * value;
	fock_a[p][s] -= D_a[q][r] * value;
        fock_b[p][q] += Dt[s][r] * value;
	fock_b[p][s] -= D_b[q][r] * value;

	/* (qp|rs) */
	fock_a[q][p] += Dt[r][s] * value;
	fock_a[q][r] -= D_a[p][s] * value;
	fock_b[q][p] += Dt[r][s] * value;
	fock_b[q][r] -= D_b[p][s] * value;

	/* (qp|sr) */
	fock_a[q][p] += Dt[s][r] * value;
        fock_a[q][s] -= D_a[p][r] * value;
	fock_b[q][p] += Dt[s][r] * value;
        fock_b[q][s] -= D_b[p][r] * value;

	/* (rs|pq) */
	fock_a[r][s] += Dt[p][q] * value;
	fock_a[r][p] -= D_a[s][q] * value;
	fock_b[r][s] += Dt[p][q] * value;
	fock_b[r][p] -= D_b[s][q] * value;

	/* (rs|qp) */
	fock_a[r][s] += Dt[q][p] * value;
	fock_a[r][q] -= D_a[s][p] * value;
	fock_b[r][s] += Dt[q][p] * value;
	fock_b[r][q] -= D_b[s][p] * value;

	/* (sr|pq) */
	fock_a[s][r] += Dt[p][q] * value;
	fock_a[s][p] -= D_a[r][q] * value;
	fock_b[s][r] += Dt[p][q] * value;
	fock_b[s][p] -= D_b[r][q] * value;

	/* (sr|qp) */
	fock_a[s][r] += Dt[q][p] * value;
	fock_a[s][q] -= D_a[r][p] * value;
	fock_b[s][r] += Dt[q][p] * value;
	fock_b[s][q] -= D_b[r][p] * value;
      }
      else if(p!=q && r!=s && pq==rs) {

	/* (pq|sr) */
        fock_a[p][q] += Dt[s][r] * value;
	fock_a[p][s] -= D_a[q][r] * value;
        fock_b[p][q] += Dt[s][r] * value;
	fock_b[p][s] -= D_b[q][r] * value;

	/* (qp|rs) */
	fock_a[q][p] += Dt[r][s] * value;
	fock_a[q][r] -= D_a[p][s] * value;
	fock_b[q][p] += Dt[r][s] * value;
	fock_b[q][r] -= D_b[p][s] * value;

	/* (qp|sr) */
	fock_a[q][p] += Dt[s][r] * value;
        fock_a[q][s] -= D_a[p][r] * value;
	fock_b[q][p] += Dt[s][r] * value;
        fock_b[q][s] -= D_b[p][r] * value;

      }
      else if(p!=q && r==s) {

	/* (qp|rs) */
	fock_a[q][p] += Dt[r][s] * value;
	fock_a[q][r] -= D_a[p][s] * value;
	fock_b[q][p] += Dt[r][s] * value;
	fock_b[q][r] -= D_b[p][s] * value;

	/* (rs|pq) */
	fock_a[r][s] += Dt[p][q] * value;
	fock_a[r][p] -= D_a[s][q] * value;
	fock_b[r][s] += Dt[p][q] * value;
	fock_b[r][p] -= D_b[s][q] * value;

	/* (rs|qp) */
	fock_a[r][s] += Dt[q][p] * value;
	fock_a[r][q] -= D_a[s][p] * value;
	fock_b[r][s] += Dt[q][p] * value;
	fock_b[r][q] -= D_b[s][p] * value;

      }
      else if(p==q && r!=s) {

	/* (pq|sr) */
        fock_a[p][q] += Dt[s][r] * value;
	fock_a[p][s] -= D_a[q][r] * value;
        fock_b[p][q] += Dt[s][r] * value;
	fock_b[p][s] -= D_b[q][r] * value;

	/* (rs|pq) */
	fock_a[r][s] += Dt[p][q] * value;
	fock_a[r][p] -= D_a[s][q] * value;
	fock_b[r][s] += Dt[p][q] * value;
	fock_b[r][p] -= D_b[s][q] * value;

	/* (sr|pq) */
	fock_a[s][r] += Dt[p][q] * value;
	fock_a[s][p] -= D_a[r][q] * value;
	fock_b[s][r] += Dt[p][q] * value;
	fock_b[s][p] -= D_b[r][q] * value;

      }
      else if(p==q && r==s && pq!=rs) {

	/* (rs|pq) */
	fock_a[r][s] += Dt[p][q] * value;
	fock_a[r][p] -= D_a[s][q] * value;
	fock_b[r][s] += Dt[p][q] * value;
	fock_b[r][p] -= D_b[s][q] * value;

      }
    }

    if(!lastbuf) iwl_buf_fetch(&InBuf);

  } while (!lastbuf);
  iwl_buf_close(&InBuf, 1);

  free_block(Dt);
}
}} // namespace psi::ccenergy
