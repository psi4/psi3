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
  int h, nirreps, i, j, a, b, ij, ab, C_irr;
  double tval;

  nirreps = RIA->params->nirreps;
  C_irr = RIA->my_irrep;

  dpd_file2_mat_init(RIA);
  dpd_file2_mat_rd(RIA);
  dpd_file2_init(&DIA, EOM_D, C_irr, 0, 1, "DIA");
  dpd_file2_mat_init(&DIA);
  dpd_file2_mat_rd(&DIA);
  for(h=0; h < nirreps; h++)
     for(i=0; i < RIA->params->rowtot[h]; i++)
        for(a=0; a < RIA->params->coltot[h^C_irr]; a++) {
           tval = eval - DIA.matrix[h][i][a];
           if (fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
        }
  dpd_file2_mat_wrt(RIA);
  dpd_file2_mat_close(RIA);
  dpd_file2_mat_close(&DIA);
  dpd_file2_close(&DIA);

  dpd_file2_mat_init(Ria);
  dpd_file2_mat_rd(Ria);
  dpd_file2_init(&Dia, EOM_D, C_irr, 0, 1, "Dia");
  dpd_file2_mat_init(&Dia);
  dpd_file2_mat_rd(&Dia);
  for(h=0; h < nirreps; h++)
     for(i=0; i < Ria->params->rowtot[h]; i++)
        for(a=0; a < Ria->params->coltot[h^C_irr]; a++) {
           tval = eval - Dia.matrix[h][i][a];
           if (fabs(tval) > 0.0001) Ria->matrix[h][i][a] /= tval;
        }
  dpd_file2_mat_wrt(Ria);
  dpd_file2_mat_close(Ria);
  dpd_file2_mat_close(&Dia);
  dpd_file2_close(&Dia);


  dpd_buf4_init(&DIJAB, EOM_D, C_irr, 2, 7, 2, 7, 0, "DIJAB");
  for(h=0; h < RIJAB->params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(RIJAB, h);
    dpd_buf4_mat_irrep_init(&DIJAB, h);
    dpd_buf4_mat_irrep_rd(RIJAB, h);
    dpd_buf4_mat_irrep_rd(&DIJAB, h);
    for(ij=0; ij < RIJAB->params->rowtot[h]; ij++)
       for(ab=0; ab < RIJAB->params->coltot[h^C_irr]; ab++) {
           tval = eval - DIJAB.matrix[h][ij][ab];
           if (fabs(tval) > 0.0001) RIJAB->matrix[h][ij][ab] /= tval;
      }
    dpd_buf4_mat_irrep_wrt(RIJAB, h);
    dpd_buf4_mat_irrep_close(RIJAB, h);
    dpd_buf4_mat_irrep_close(&DIJAB, h);
  }
  dpd_buf4_close(&DIJAB);


  dpd_buf4_init(&Dijab, EOM_D, C_irr, 2, 7, 2, 7, 0, "Dijab");
  for(h=0; h < Rijab->params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(Rijab, h);
    dpd_buf4_mat_irrep_init(&Dijab, h);
    dpd_buf4_mat_irrep_rd(Rijab, h);
    dpd_buf4_mat_irrep_rd(&Dijab, h);
    for(ij=0; ij < Rijab->params->rowtot[h]; ij++)
       for(ab=0; ab < Rijab->params->coltot[h^C_irr]; ab++) {
           tval = eval - Dijab.matrix[h][ij][ab];
           if (fabs(tval) > 0.0001) Rijab->matrix[h][ij][ab] /= tval;
      }
    dpd_buf4_mat_irrep_wrt(Rijab, h);
    dpd_buf4_mat_irrep_close(Rijab, h);
    dpd_buf4_mat_irrep_close(&Dijab, h);
  }
  dpd_buf4_close(&Dijab);


  dpd_buf4_init(&DIjAb, EOM_D, C_irr, 0, 5, 0, 5, 0, "DIjAb");
  for(h=0; h < RIjAb->params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(RIjAb, h);
    dpd_buf4_mat_irrep_init(&DIjAb, h);
    dpd_buf4_mat_irrep_rd(RIjAb, h);
    dpd_buf4_mat_irrep_rd(&DIjAb, h);
    for(ij=0; ij < RIjAb->params->rowtot[h]; ij++)
       for(ab=0; ab < RIjAb->params->coltot[h^C_irr]; ab++) {
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
  int h, nirreps, i, j, a, b, ij, ab, C_irr;
  double tval;

  nirreps = RIA->params->nirreps;
  C_irr = RIA->my_irrep;

  dpd_file2_mat_init(RIA);
  dpd_file2_mat_rd(RIA);
  dpd_file2_init(&DIA, EOM_D, C_irr, 0, 1, "DIA");
  dpd_file2_mat_init(&DIA);
  dpd_file2_mat_rd(&DIA);
  for(h=0; h < nirreps; h++)
    for(i=0; i < RIA->params->rowtot[h]; i++)
      for(a=0; a < RIA->params->coltot[h^C_irr]; a++) {
	tval = eval - DIA.matrix[h][i][a];
	if (fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
      }
  dpd_file2_mat_wrt(RIA);
  dpd_file2_mat_close(RIA);
  dpd_file2_mat_close(&DIA);
  dpd_file2_close(&DIA);

  dpd_buf4_init(&DIjAb, EOM_D, C_irr, 0, 5, 0, 5, 0, "DIjAb");
  for(h=0; h < RIjAb->params->nirreps; h++) {
    dpd_buf4_mat_irrep_init(RIjAb, h);
    dpd_buf4_mat_irrep_init(&DIjAb, h);
    dpd_buf4_mat_irrep_rd(RIjAb, h);
    dpd_buf4_mat_irrep_rd(&DIjAb, h);
    for(ij=0; ij < RIjAb->params->rowtot[h]; ij++)
      for(ab=0; ab < RIjAb->params->coltot[h^C_irr]; ab++) {
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
