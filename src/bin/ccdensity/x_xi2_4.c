#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* compute terms 4 and 6 of xi2 amplitudes */
/* Xijab += P(ij) P(ab) (Rme Lia + Rmnef Linaf) <mj||eb> */

void x_xi2_4(void)
{
  dpdfile2 RIA, Ria, LIA, Lia;
  int L_irr, R_irr, G_irr, nirreps;
  int I, A, M, E, i, a, m, e, h, row, col, Isym, Esym, Asym, Msym;
  dpdbuf4 D, R2, L2, H2, I2, Z, Z2, XIJAB, Xijab, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  if(params.ref == 0 || params.ref == 1) {
    /* construct RL = Rme Lia + Rmnef Linaf */
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
    dpd_buf4_copy(&I2, EOM_TMP1, "RL_OVOV");
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
    dpd_buf4_copy(&I2, EOM_TMP1, "RL_ovov");
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
    dpd_buf4_copy(&I2, EOM_TMP1, "RL_OVov");
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
    dpd_buf4_copy(&I2, EOM_TMP1, "RL_ovOV");
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
    dpd_buf4_copy(&I2, EOM_TMP1, "RL_oVoV");
    dpd_buf4_close(&I2);
    dpd_buf4_init(&I2, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
    dpd_buf4_copy(&I2, EOM_TMP1, "RL_OvOv");
    dpd_buf4_close(&I2);

    /* RL_OVOV(me,ia) += Rme Lia */
    dpd_file2_init(&RIA, CC_GR, R_irr, 0, 1, "RIA");
    dpd_file2_init(&Ria, CC_GR, R_irr, 0, 1, "Ria");
    dpd_file2_init(&LIA, CC_GL, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_GL, L_irr, 0, 1, "Lia");
    dpd_file2_mat_init(&RIA);
    dpd_file2_mat_init(&Ria);
    dpd_file2_mat_init(&LIA);
    dpd_file2_mat_init(&Lia);
    dpd_file2_mat_rd(&RIA);
    dpd_file2_mat_rd(&Ria);
    dpd_file2_mat_rd(&LIA);
    dpd_file2_mat_rd(&Lia);

    dpd_buf4_init(&I2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVOV");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&I2, h);
      dpd_buf4_mat_irrep_rd(&I2, h);
      for(row=0; row < I2.params->rowtot[h]; row++) {
        m = I2.params->roworb[h][row][0];
        e = I2.params->roworb[h][row][1];
        M = RIA.params->rowidx[m]; Msym = RIA.params->psym[m];
        E = RIA.params->colidx[e]; Esym = RIA.params->qsym[e];
        for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
          i = I2.params->colorb[h^G_irr][col][0];
          a = I2.params->colorb[h^G_irr][col][1];
          I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
          A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
          if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
            I2.matrix[h][row][col] += RIA.matrix[Msym][M][E] * LIA.matrix[Isym][I][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&I2, h);
    dpd_buf4_mat_irrep_close(&I2, h);
    }
    dpd_buf4_close(&I2);

    dpd_buf4_init(&I2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovov");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&I2, h);
      dpd_buf4_mat_irrep_rd(&I2, h);
      for(row=0; row < I2.params->rowtot[h]; row++) {
        m = I2.params->roworb[h][row][0];
        e = I2.params->roworb[h][row][1];
        M = Ria.params->rowidx[m]; Msym = Ria.params->psym[m];
        E = Ria.params->colidx[e]; Esym = Ria.params->qsym[e];
        for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
          i = I2.params->colorb[h^G_irr][col][0];
          a = I2.params->colorb[h^G_irr][col][1];
          I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
          A = Lia.params->colidx[a]; Asym = Lia.params->qsym[a];
          if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
            I2.matrix[h][row][col] += Ria.matrix[Msym][M][E] * Lia.matrix[Isym][I][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&I2, h);
    dpd_buf4_mat_irrep_close(&I2, h);
    }
    dpd_buf4_close(&I2);

    dpd_buf4_init(&I2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVov");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&I2, h);
      dpd_buf4_mat_irrep_rd(&I2, h);
      for(row=0; row < I2.params->rowtot[h]; row++) {
        m = I2.params->roworb[h][row][0];
        e = I2.params->roworb[h][row][1];
        M = RIA.params->rowidx[m]; Msym = RIA.params->psym[m];
        E = RIA.params->colidx[e]; Esym = RIA.params->qsym[e];
        for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
          i = I2.params->colorb[h^G_irr][col][0];
          a = I2.params->colorb[h^G_irr][col][1];
          I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
          A = Lia.params->colidx[a]; Asym = Lia.params->qsym[a];
          if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
            I2.matrix[h][row][col] += RIA.matrix[Msym][M][E] * Lia.matrix[Isym][I][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&I2, h);
    dpd_buf4_mat_irrep_close(&I2, h);
    }
    dpd_buf4_close(&I2);

    dpd_buf4_init(&I2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovOV");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&I2, h);
      dpd_buf4_mat_irrep_rd(&I2, h);
      for(row=0; row < I2.params->rowtot[h]; row++) {
        m = I2.params->roworb[h][row][0];
        e = I2.params->roworb[h][row][1];
        M = Ria.params->rowidx[m]; Msym = Ria.params->psym[m];
        E = Ria.params->colidx[e]; Esym = Ria.params->qsym[e];
        for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
          i = I2.params->colorb[h^G_irr][col][0];
          a = I2.params->colorb[h^G_irr][col][1];
          I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
          A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
          if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
            I2.matrix[h][row][col] += Ria.matrix[Msym][M][E] * LIA.matrix[Isym][I][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&I2, h);
    dpd_buf4_mat_irrep_close(&I2, h);
    }
    dpd_buf4_close(&I2);

    dpd_file2_mat_close(&RIA);
    dpd_file2_mat_close(&Ria);
    dpd_file2_mat_close(&LIA);
    dpd_file2_mat_close(&Lia);
    dpd_file2_close(&RIA);
    dpd_file2_close(&Ria);
    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);

    /* compute quantity to be permuted: Z2(IJ,AB) */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (IA,JB)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVOV");
    dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovOV");
    dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 0, 5, "Z2 (IJ,AB)");
    dpd_buf4_close(&Z2);

    /* add in Z2(IJ,AB) permutations */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (IJ,AB)");
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 5, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z2, &XIJAB, 1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "XIJAB", -1.0);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "XIJAB", -1.0);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 2, 7, "XIJAB", 1.0);
    dpd_buf4_close(&Z2);

    /* compute quantity to be permuted: Z2(ij,ab) */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (ia,jb)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovov");
    dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVov");
    dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);
    dpd_buf4_sort(&Z2, EOM_TMP1, prqs, 0, 5, "Z2 (ij,ab)");
    dpd_buf4_close(&Z2);

    /* add in Z2(ij,ab) permutations */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 5, 2, 7, 0, "Xijab");
    dpd_buf4_axpy(&Z2, &XIJAB, 1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qprs, 2, 7, "Xijab", -1.0);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, pqsr, 2, 7, "Xijab", -1.0);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, qpsr, 2, 7, "Xijab", 1.0);
    dpd_buf4_close(&Z2);
    
    /* XIjAb += ZMIAE <Mj|Eb> + ZmIeA <mj||eb> */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (IA,jb)");
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVOV");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovOV");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);

    /* XIjAb += ZMjEb <MI|EA> + Zmjeb <mI|eA> */
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVov");
    dpd_contract444(&D, &Z, &Z2, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovov");
    dpd_contract444(&D, &Z, &Z2, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);
    dpd_buf4_sort_axpy(&Z2, EOM_XI, prqs, 0, 5, "XIjAb", 1.0);
    dpd_buf4_close(&Z2);

    /* XIjAb += ZmjEA <mI|bE> + ZMIeb <Mj|Ae> */
    dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (Ib,jA)");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_oVoV");
    dpd_contract444(&D, &Z, &Z2, 1, 1, 1.0, 0.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);

    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OvOv");
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&D);

    dpd_buf4_sort_axpy(&Z2, EOM_XI, prsq, 0, 5, "XIjAb", 1.0);
    dpd_buf4_close(&Z2);
  }
  return;
}

