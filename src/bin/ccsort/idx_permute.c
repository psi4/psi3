#include <stdio.h>
#include <iwl.h>
#include "dpd.h"

void idx_error(char *message, int p, int q, int r, int s, int pq, int rs,
	       int pq_sym, FILE *outfile);

void idx_permute(struct dpdfile *File, struct iwlbuf *OutBuf,
		 int **bucket_map, int p, int q, int r, int s,
		 int perm_pr, int perm_qs, int perm_prqs,
		 double value, FILE *outfile)
{
  int p_sym, q_sym, r_sym, s_sym;
  int pq_sym, rs_sym, rq_sym, ps_sym, qp_sym, sp_sym, sr_sym, qr_sym;
  int pq, rs, rq, ps, qp, sr, qr, sp;
  int perm_pq, perm_rs;
  struct dpdparams *Params;
  int this_bucket;

  Params = File->params;
  perm_pq = Params->perm_pq;
  perm_rs = Params->perm_rs;
  
  /* Get the orbital symmetries */
  p_sym = Params->psym[p]; q_sym = Params->qsym[q];
  r_sym = Params->rsym[r]; s_sym = Params->ssym[s];

  /* Go through the allowed permutations --- NB these are Dirac permutations */

  /* Get the main symmetry block */
  pq_sym = p_sym^q_sym;

  /* Get the row and column indices and assign the value */
  pq = Params->rowidx[p][q];
  rs = Params->colidx[r][s];
  if((pq >= Params->rowtot[pq_sym]) || (rs >= Params->coltot[pq_sym]))
      idx_error("Params_make: pq, rs", p,q,r,s,pq,rs,pq_sym,outfile);
  /*
  File->matrix[pq_sym][pq][rs] = value;
  */

  this_bucket = bucket_map[p][q];
  iwl_buf_wrt_val(&OutBuf[this_bucket], p, q, r, s, value, 0, outfile, 0);

  if(perm_pr) {
      rq_sym = r_sym^q_sym;
      rq = Params->rowidx[r][q];
      ps = Params->colidx[p][s];
      if((rq >= Params->rowtot[rq_sym]) || (ps >= Params->coltot[rq_sym]))
	  idx_error("Params_make: rq, ps", p,q,r,s,rq,ps,rq_sym,outfile);
      /*
      File->matrix[rq_sym][rq][ps] = value;
      */

      this_bucket = bucket_map[r][q];
      iwl_buf_wrt_val(&OutBuf[this_bucket], r, q, p, s, value, 0, outfile, 0);
    }

  if(perm_qs) {
      ps_sym = p_sym^s_sym;
      ps = Params->rowidx[p][s];
      rq = Params->colidx[r][q];
      if((ps >= Params->rowtot[ps_sym]) || (rq >= Params->coltot[ps_sym]))
	  idx_error("Params_make: ps, rq", p,q,r,s,ps,rq,ps_sym,outfile);
      /*
      File->matrix[ps_sym][ps][rq] = value;
      */

      this_bucket = bucket_map[p][s];
      iwl_buf_wrt_val(&OutBuf[this_bucket], p, s, r, q, value, 0, outfile, 0);
    }

  if(perm_pr && perm_qs) {
      rs_sym = r_sym^s_sym;
      rs = Params->rowidx[r][s];
      pq = Params->colidx[p][q];
      if((rs >= Params->rowtot[rs_sym]) || (pq >= Params->coltot[rs_sym]))
	  idx_error("Params_make: rs, pq", p,q,r,s,rs,pq,rs_sym,outfile);
      /*
      File->matrix[rs_sym][rs][pq] = value;
      */

      this_bucket = bucket_map[r][s];
      iwl_buf_wrt_val(&OutBuf[this_bucket], r, s, p, q, value, 0, outfile, 0);
    }

  if(perm_prqs) {
      qp_sym = q_sym^p_sym;
      qp = Params->rowidx[q][p];
      sr = Params->colidx[s][r];
      if((qp >= Params->rowtot[qp_sym]) || (sr >= Params->coltot[qp_sym]))
	  idx_error("Params_make: qp, sr", p,q,r,s,qp,sr,qp_sym,outfile);
      /*
      File->matrix[qp_sym][qp][sr] = value;
      */

      this_bucket = bucket_map[q][p];
      iwl_buf_wrt_val(&OutBuf[this_bucket], q, p, s, r, value, 0, outfile, 0);
      

      if(perm_pr) {
	  qr_sym = q_sym^r_sym;
	  qr = Params->rowidx[q][r];
	  sp = Params->colidx[s][p];
	  if((qr >= Params->rowtot[qr_sym])||(sp >= Params->coltot[qr_sym]))
	      idx_error("Params_make: qr, sp", p,q,r,s,qr,sp,qr_sym,
			       outfile);
	  /*
	  File->matrix[qr_sym][qr][sp] = value;
	  */

	  this_bucket = bucket_map[q][r];
	  iwl_buf_wrt_val(&OutBuf[this_bucket], q, r, s, p, value, 0, outfile, 0);
	}

      if(perm_qs) {
	  sp_sym = s_sym^p_sym;
	  sp = Params->rowidx[s][p];
	  qr = Params->colidx[q][r];
	  if((sp >= Params->rowtot[sp_sym])||(qr >= Params->coltot[sp_sym]))
	      idx_error("Params_make: sp, qr", p,q,r,s,sp,qr,sp_sym,
			       outfile);
	  /*
	  File->matrix[sp_sym][sp][qr] = value;
	  */

	  this_bucket = bucket_map[s][p];
	  iwl_buf_wrt_val(&OutBuf[this_bucket], s, p, q, r, value, 0, outfile, 0);
	}
      
      if(perm_pr && perm_qs) {
	  sr_sym = s_sym^r_sym;
	  sr = Params->rowidx[s][r];
	  qp = Params->colidx[q][p];
	  if((sr >= Params->rowtot[sr_sym])||(qp >= Params->coltot[sr_sym]))
	      idx_error("Params_make: sr, qp", p,q,r,s,sr,qp,qp_sym,
			       outfile);
	  /*
	  File->matrix[sr_sym][sr][qp] = value;
	  */

	  this_bucket = bucket_map[s][r];
	  iwl_buf_wrt_val(&OutBuf[this_bucket], s, r, q, p, value, 0, outfile, 0);
	}
    }
}
