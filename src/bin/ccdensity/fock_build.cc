/*! \file
    \ingroup CCDENSITY
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

#define IOFF_MAX 32641

namespace psi { namespace ccdensity {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void rhf_fock_build(double **fock, double  **D)
{
  int i, j, ij;
  int *ioff;
  int nso, ntri, stat;
  double *scratch;
  int lastbuf, idx, p, q, r, s, pq, rs;
  double value;
  Value *valptr;
  Label *lblptr;
  struct iwlbuf InBuf;

  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;

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
}} // namespace psi::ccdensity
