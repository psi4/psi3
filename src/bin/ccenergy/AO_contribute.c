#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include <iwl.h>
#include <qt.h>
#include <dpd.h>

void AO_contribute(struct iwlbuf *InBuf, dpdbuf4 *tau1_AO, dpdbuf4 *tau2_AO)
{
  int idx, p, q, r, s;
  double value;
  Value *valptr;
  Label *lblptr;
  int Gp, Gq, Gr, Gs, Gpr, Gps, Gqr, Gqs, Grp, Gsp, Grq, Gsq;
  int pr, ps, qr, qs, rp, rq, sp, sq, pq, rs;

  lblptr = InBuf->labels;
  valptr = InBuf->values;

  for(idx=4*InBuf->idx; InBuf->idx < InBuf->inbuf; InBuf->idx++) {
    p = abs((int) lblptr[idx++]);
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    value = (double) valptr[InBuf->idx];

    Gp = tau1_AO->params->psym[p]; 
    Gq = tau1_AO->params->psym[q]; 
    Gr = tau1_AO->params->psym[r]; 
    Gs = tau1_AO->params->psym[s];

    Gpr = Grp = Gp^Gr;
    Gps = Gsp = Gp^Gs;
    Gqr = Grq = Gq^Gr;
    Gqs = Gsq = Gq^Gs;

    pq = tau1_AO->params->rowidx[p][q];  
    rs = tau1_AO->params->rowidx[r][s];

    pr = tau1_AO->params->rowidx[p][r];
    rp = tau1_AO->params->rowidx[r][p];
    ps = tau1_AO->params->rowidx[p][s];
    sp = tau1_AO->params->rowidx[s][p];
    qr = tau1_AO->params->rowidx[q][r];
    rq = tau1_AO->params->rowidx[r][q];
    qs = tau1_AO->params->rowidx[q][s];
    sq = tau1_AO->params->rowidx[s][q];

      /* (pq|rs) */
    C_DAXPY(tau1_AO->params->coltot[Gpr], value, tau1_AO->matrix[Gpr][qs], 1,
	    tau2_AO->matrix[Gpr][pr], 1);

    if(p!=q && r!=s && pq != rs) {

      /* (pq|sr) */
      C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1,
	      tau2_AO->matrix[Gps][ps], 1);

      /* (qp|rs) */
      C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1,
	      tau2_AO->matrix[Gqr][qr], 1);

      /* (qp|sr) */
      C_DAXPY(tau1_AO->params->coltot[Gqs], value, tau1_AO->matrix[Gqs][pr], 1,
	      tau2_AO->matrix[Gqs][qs], 1);

      /* (rs|pq) */
      C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
	      tau2_AO->matrix[Grp][rp], 1);

      /* (sr|pq) */
      C_DAXPY(tau1_AO->params->coltot[Gsp], value, tau1_AO->matrix[Gsp][rq], 1, 
	      tau2_AO->matrix[Gsp][sp], 1);

      /* (rs|qp) */
      C_DAXPY(tau1_AO->params->coltot[Grq], value, tau1_AO->matrix[Grq][sp], 1,
	      tau2_AO->matrix[Grq][rq], 1);

      /* (sr|qp) */
      C_DAXPY(tau1_AO->params->coltot[Gsq], value, tau1_AO->matrix[Gsq][rp], 1,
	      tau2_AO->matrix[Gsq][sq],1 );

    }
    else if(p!=q && r!=s && pq==rs) {

      /* (pq|sr) */
      C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1,
	      tau2_AO->matrix[Gps][ps], 1);

      /* (qp|rs) */
      C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1,
	      tau2_AO->matrix[Gqr][qr], 1);

      /* (qp|sr) */
      C_DAXPY(tau1_AO->params->coltot[Gqs], value, tau1_AO->matrix[Gqs][pr], 1,
	      tau2_AO->matrix[Gqs][qs], 1);

    }
    else if(p!=q && r==s) {

      /* (qp|rs) */
      C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1,
	      tau2_AO->matrix[Gqr][qr], 1);

      /* (rs|pq) */
      C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
	      tau2_AO->matrix[Grp][rp], 1);

      /* (rs|qp) */
      C_DAXPY(tau1_AO->params->coltot[Grq], value, tau1_AO->matrix[Grq][sp], 1,
	      tau2_AO->matrix[Grq][rq], 1);

    }

    else if(p==q && r!=s) {

      /* (pq|sr) */
      C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1,
	      tau2_AO->matrix[Gps][ps], 1);

      /* (rs|pq) */
      C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
	      tau2_AO->matrix[Grp][rp], 1);

      /* (sr|pq) */
      C_DAXPY(tau1_AO->params->coltot[Gsp], value, tau1_AO->matrix[Gsp][rq], 1, 
	      tau2_AO->matrix[Gsp][sp], 1);

    }

    else if(p==q && r==s && pq != rs) {

      /* (rs|pq) */
      C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
	      tau2_AO->matrix[Grp][rp], 1);

    }
  }

  return;
}
