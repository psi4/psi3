#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void precondition(dpdfile2 *RIA, dpdfile2 *Ria, 
  dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb,
  double eval)
{
  dpdfile2 DIA, Dia;
  dpdbuf4 DIJAB, Dijab, DIjAb;
  int h, nirreps, i, j, a, b, ij, ab;
  double tval;
  int irrep = 0; /* to be set by RIA irrep */

  nirreps = RIA->params->nirreps;

  dpd_file2_mat_init(RIA);
  dpd_file2_mat_rd(RIA);
  dpd_file2_init(&DIA, EOM_D, irrep, 0, 1, "DIA");
  dpd_file2_mat_init(&DIA);
  dpd_file2_mat_rd(&DIA);
  for(h=0; h < nirreps; h++)
     for(i=0; i < RIA->params->rowtot[h]; i++)
        for(a=0; a < RIA->params->coltot[h]; a++) {
           tval = eval - DIA.matrix[h][i][a];
           if (fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
        }
  dpd_file2_mat_wrt(RIA);
  dpd_file2_mat_close(RIA);
  dpd_file2_mat_close(&DIA);
  dpd_file2_close(&DIA);

  dpd_file2_mat_init(Ria);
  dpd_file2_mat_rd(Ria);
  dpd_file2_init(&Dia, EOM_D, irrep, 0, 1, "Dia");
  dpd_file2_mat_init(&Dia);
  dpd_file2_mat_rd(&Dia);
  for(h=0; h < nirreps; h++)
     for(i=0; i < Ria->params->rowtot[h]; i++)
        for(a=0; a < Ria->params->coltot[h]; a++) {
           tval = eval - Dia.matrix[h][i][a];
           if (fabs(tval) > 0.0001) Ria->matrix[h][i][a] /= tval;
        }
  dpd_file2_mat_wrt(Ria);
  dpd_file2_mat_close(Ria);
  dpd_file2_mat_close(&Dia);
  dpd_file2_close(&Dia);


  dpd_buf4_init(&DIJAB, EOM_D, irrep, 2, 7, 2, 7, 0, "DIJAB");
  for(h=0; h < RIJAB->params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(RIJAB, h);
    dpd_buf4_mat_irrep_init(&DIJAB, h);
    dpd_buf4_mat_irrep_rd(RIJAB, h);
    dpd_buf4_mat_irrep_rd(&DIJAB, h);
    for(ij=0; ij < RIJAB->params->rowtot[h]; ij++)
       for(ab=0; ab < RIJAB->params->coltot[h]; ab++) {
           tval = eval - DIJAB.matrix[h][ij][ab];
           if (fabs(tval) > 0.0001) RIJAB->matrix[h][ij][ab] /= tval;
      }
    dpd_buf4_mat_irrep_wrt(RIJAB, h);
    dpd_buf4_mat_irrep_close(RIJAB, h);
    dpd_buf4_mat_irrep_close(&DIJAB, h);
  }
  dpd_buf4_close(&DIJAB);


  dpd_buf4_init(&Dijab, EOM_D, irrep, 2, 7, 2, 7, 0, "Dijab");
  for(h=0; h < Rijab->params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(Rijab, h);
    dpd_buf4_mat_irrep_init(&Dijab, h);
    dpd_buf4_mat_irrep_rd(Rijab, h);
    dpd_buf4_mat_irrep_rd(&Dijab, h);
    for(ij=0; ij < Rijab->params->rowtot[h]; ij++)
       for(ab=0; ab < Rijab->params->coltot[h]; ab++) {
           tval = eval - Dijab.matrix[h][ij][ab];
           if (fabs(tval) > 0.0001) Rijab->matrix[h][ij][ab] /= tval;
      }
    dpd_buf4_mat_irrep_wrt(Rijab, h);
    dpd_buf4_mat_irrep_close(Rijab, h);
    dpd_buf4_mat_irrep_close(&Dijab, h);
  }
  dpd_buf4_close(&Dijab);


  dpd_buf4_init(&DIjAb, EOM_D, irrep, 0, 5, 0, 5, 0, "DIjAb");
  for(h=0; h < RIjAb->params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(RIjAb, h);
    dpd_buf4_mat_irrep_init(&DIjAb, h);
    dpd_buf4_mat_irrep_rd(RIjAb, h);
    dpd_buf4_mat_irrep_rd(&DIjAb, h);
    for(ij=0; ij < RIjAb->params->rowtot[h]; ij++)
       for(ab=0; ab < RIjAb->params->coltot[h]; ab++) {
           tval = eval - DIjAb.matrix[h][ij][ab];
           if (fabs(tval) > 0.0001) RIjAb->matrix[h][ij][ab] /= tval;
      }
    dpd_buf4_mat_irrep_wrt(RIjAb, h);
    dpd_buf4_mat_irrep_close(RIjAb, h);
    dpd_buf4_mat_irrep_close(&DIjAb, h);
  }
  dpd_buf4_close(&DIjAb);

  return;
}

void precondition_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, double eval)
{
  dpdfile2 DIA;
  dpdbuf4 DIjAb;
  int h, nirreps, i, j, a, b, ij, ab;
  double tval;
  int irrep = 0; /* to be set by RIA irrep */

  /* Local correlation decs */
  int nocc, nvir, natom, nao;
  int **pairdomain;
  double ***V, ***W, *eps_occ, **eps_vir;
  double *T1tilde, *T1bar, **T2tilde, **T2bar, **X1, **X2;

  nirreps = RIA->params->nirreps;

  if(params.local && !strcmp(local.method,"WERNER")) {

    nocc = local.nocc;
    nvir = local.nvir;
    natom = local.natom;
    pairdomain = local.pairdomain;
    pairdom_len = local.pairdom_len;
    pairdom_nrlen = local.pairdom_nrlen;
    eps_occ = local.eps_occ;
    eps_vir = local.eps_vir;

    dpd_file2_mat_init(RIA);
    dpd_file2_mat_rd(RIA);

    for(i=0; i < nocc; i++) {
      ii = i * nocc +i;

      T1tilde = init_array(pairdom_len[ii]);
      T1bar = init_array(pairdom_nrlen[ii]);

      /* Transform the virtuals to the redundant projected virtual basis */
      C_DGEMV('t', nvir, pairdom_len[ii], 1.0, &(V[ii][0][0]), pairdom_len[ii], 
	      &(RIA->matrix[0][i][0]), 1, 0.0, &(T1tilde[0]), 1);

      /* Transform the virtuals to the non-redundant virtual basis */
      C_DGEMV('t', pairdom_len[ii], pairdom_nrlen[ii], 1.0, &(W[ii][0][0]), pairdom_nrlen[ii], 
	      &(T1tilde[0]), 1, 0.0, &(T1bar[0]), 1);

      for(a=0; a < pairdom_nrlen[ii]; a++) 
	T1bar[a] /= (eval + eps_occ[i] - eps_vir[ij][a]);

      /* Transform the new T1's to the redundant projected virtual basis */
      C_DGEMV('n', pairdom_len[ii], pairdom_nrlen[ii], 1.0, &(W[ii][0][0]), pairdom_nrlen[ii],
	      &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);

      /* Transform the new T1's to the MO basis */
      C_DGEMV('n', nvir, pairdom_len[ii], 1.0, &(V[ii][0][0]), pairdom_len[ii], 
	      &(T1tilde[0]), 1, 0.0, &(RIA->matrix[0][i][0]), 1);

      free(T1bar);
      free(T1tilde);
    }

    dpd_file2_mat_wrt(RIA);
    dpd_file2_mat_close(RIA);

    dpd_buf4_mat_irrep_init(RIjAb, 0);
    dpd_buf4_mat_irrep_rd(RIjAb, 0);

    X1 = block_matrix(nao,nvir);
    X2 = block_matrix(nvir,nao);

    T2tilde = block_matrix(nao,nao);
    T2bar = block_matrix(nvir, nvir);
    for(i=0,ij=0; i < nocc; i++) {
      for(j=0; j < nocc; j++,ij++) {

	/* Transform the virtuals to the redundant projected virtual basis */
	C_DGEMM('t', 'n', pairdom_len[ij], nvir, nvir, 1.0, &(V[ij][0][0]), pairdom_len[ij],
		&(RIjAb->matrix[0][ij][0]), nvir, 0.0, &(X1[0][0]), nvir);
	C_DGEMM('n', 'n', pairdom_len[ij], pairdom_len[ij], nvir, 1.0, &(X1[0][0]), nvir,
		&(V[ij][0][0]), pairdom_len[ij], 0.0, &(T2tilde[0][0]), nao);

	/* Transform the virtuals to the non-redundant virtual basis */
	C_DGEMM('t', 'n', pairdom_nrlen[ij], pairdom_len[ij], pairdom_len[ij], 1.0, 
		&(W[ij][0][0]), pairdom_nrlen[ij], &(T2tilde[0][0]), nao, 0.0, &(X2[0][0]), nao);
	C_DGEMM('n', 'n', pairdom_nrlen[ij], pairdom_nrlen[ij], pairdom_len[ij], 1.0, 
		&(X2[0][0]), nao, &(W[ij][0][0]), pairdom_nrlen[ij], 0.0, &(T2bar[0][0]), nvir);

	/* Divide the new amplitudes by the denominators */
	for(a=0; a < pairdom_nrlen[ij]; a++) {
	  for(b=0; b < pairdom_nrlen[ij]; b++) {
	    T2bar[a][b] /= (eval + eps_occ[i] + eps_occ[j] - eps_vir[ij][a] - eps_vir[ij][b]);
	  }
	}

	/* Transform the new T2's to the redundant virtual basis */
	C_DGEMM('n', 'n', pairdom_len[ij], pairdom_nrlen[ij], pairdom_nrlen[ij], 1.0, 
		&(W[ij][0][0]), pairdom_nrlen[ij], &(T2bar[0][0]), nvir, 0.0, &(X1[0][0]), nvir);
	C_DGEMM('n','t', pairdom_len[ij], pairdom_len[ij], pairdom_nrlen[ij], 1.0, 
		&(X1[0][0]), nvir, &(W[ij][0][0]), pairdom_nrlen[ij], 0.0, &(T2tilde[0][0]), nao);

	/* Transform the new T2's to the MO basis */
	C_DGEMM('n', 'n', nvir, pairdom_len[ij], pairdom_len[ij], 1.0, 
		&(V[ij][0][0]), pairdom_len[ij], &(T2tilde[0][0]), nao, 0.0, &(X2[0][0]), nao);
	C_DGEMM('n', 't', nvir, nvir, pairdom_len[ij], 1.0, &(X2[0][0]), nao,
		&(V[ij][0][0]), pairdom_len[ij], 0.0, &(RIjAb->matrix[0][ij][0]), nvir);
      }
    }
    free_block(T2tilde);
    free_block(T2bar);

    free_block(X1);
    free_block(X2);

    dpd_buf4_mat_irrep_wrt(RIjAb, 0);
    dpd_buf4_mat_irrep_close(RIjAb, 0);
  }
  else {

    dpd_file2_mat_init(RIA);
    dpd_file2_mat_rd(RIA);
    dpd_file2_init(&DIA, EOM_D, irrep, 0, 1, "DIA");
    dpd_file2_mat_init(&DIA);
    dpd_file2_mat_rd(&DIA);
    for(h=0; h < nirreps; h++)
      for(i=0; i < RIA->params->rowtot[h]; i++)
        for(a=0; a < RIA->params->coltot[h]; a++) {
	  tval = eval - DIA.matrix[h][i][a];
	  if (fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
        }
    dpd_file2_mat_wrt(RIA);
    dpd_file2_mat_close(RIA);
    dpd_file2_mat_close(&DIA);
    dpd_file2_close(&DIA);

    dpd_buf4_init(&DIjAb, EOM_D, irrep, 0, 5, 0, 5, 0, "DIjAb");
    for(h=0; h < RIjAb->params->nirreps; h++) {
      dpd_buf4_mat_irrep_init(RIjAb, h);
      dpd_buf4_mat_irrep_init(&DIjAb, h);
      dpd_buf4_mat_irrep_rd(RIjAb, h);
      dpd_buf4_mat_irrep_rd(&DIjAb, h);
      for(ij=0; ij < RIjAb->params->rowtot[h]; ij++)
	for(ab=0; ab < RIjAb->params->coltot[h]; ab++) {
	  tval = eval - DIjAb.matrix[h][ij][ab];
	  if (fabs(tval) > 0.0001) RIjAb->matrix[h][ij][ab] /= tval;
	}
      dpd_buf4_mat_irrep_wrt(RIjAb, h);
      dpd_buf4_mat_irrep_close(RIjAb, h);
      dpd_buf4_mat_irrep_close(&DIjAb, h);
    }
    dpd_buf4_close(&DIjAb);
  }

  return;
}
