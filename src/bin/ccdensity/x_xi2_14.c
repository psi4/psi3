#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* compute term 14 of Xi2 amplitudes */
/* Xijab += P(ij) P(ab) (Lmjeb) Rme Fia */

void x_xi2_14(void)
{
  dpdfile2 IIA, Iia, FME, Fme;
  int L_irr, R_irr, G_irr, nirreps;
  int I, A, J, B, i, a, j, b, h, row, col, Isym, Asym, Jsym, Bsym;
  int II, AA, JJ, BB, ii, aa, jj, bb, IIsym, AAsym, JJsym, BBsym;
  dpdbuf4 Z, XIJAB, Xijab, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  if(params.ref == 0 || params.ref == 1) {
    dpd_file2_init(&IIA, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_file2_init(&Iia, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    dpd_file2_mat_init(&IIA); dpd_file2_mat_init(&Iia);
    dpd_file2_mat_rd(&IIA); dpd_file2_mat_rd(&Iia);

    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
    dpd_file2_mat_init(&FME); dpd_file2_mat_init(&Fme);
    dpd_file2_mat_rd(&FME); dpd_file2_mat_rd(&Fme);

    /* build Z(IJ,AB) = FME(I,A) L2R1_OV(J,B) */
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (IJ,AB)");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&Z, h);
      for(row=0; row < Z.params->rowtot[h]; row++) {
        i = Z.params->roworb[h][row][0];
        j = Z.params->roworb[h][row][1];
        I = FME.params->rowidx[i]; Isym = FME.params->psym[i];
        J = IIA.params->rowidx[j]; Jsym = IIA.params->psym[j];
        for(col=0; col < Z.params->coltot[h^G_irr]; col++) {
          a = Z.params->colorb[h^G_irr][col][0];
          b = Z.params->colorb[h^G_irr][col][1];
          A = FME.params->colidx[a]; Asym = FME.params->qsym[a];
          B = IIA.params->colidx[b]; Bsym = IIA.params->qsym[b];
          if( (Isym==Asym) && ((Jsym^Bsym)==G_irr) )
            Z.matrix[h][row][col] += FME.matrix[Isym][I][A] * IIA.matrix[Jsym][J][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&Z, h);
    dpd_buf4_mat_irrep_close(&Z, h);
    }

    /* XIJAB += P(IJ) P(AB) Z(IJ,AB) */
    dpd_buf4_init(&XIJAB, EOM_XI, G_irr, 0, 5, 2, 7, 0, "XIJAB");
    dpd_buf4_axpy(&Z, &XIJAB, 1.0);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_sort_axpy(&Z, EOM_XI, qprs, 2, 7, "XIJAB", -1.0);
    dpd_buf4_sort_axpy(&Z, EOM_XI, pqsr, 2, 7, "XIJAB", -1.0);
    dpd_buf4_sort_axpy(&Z, EOM_XI, qpsr, 2, 7, "XIJAB", 1.0);
    dpd_buf4_close(&Z);

    /* build Z(ij,ab) = Fme(i,a) L2R1_ov(j,b) */
    dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (ij,ab)");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&Z, h);
      for(row=0; row < Z.params->rowtot[h]; row++) {
        i = Z.params->roworb[h][row][0];
        j = Z.params->roworb[h][row][1];
        I = Fme.params->rowidx[i]; Isym = Fme.params->psym[i];
        J = Iia.params->rowidx[j]; Jsym = Iia.params->psym[j];
        for(col=0; col < Z.params->coltot[h^G_irr]; col++) {
          a = Z.params->colorb[h^G_irr][col][0];
          b = Z.params->colorb[h^G_irr][col][1];
          A = Fme.params->colidx[a]; Asym = Fme.params->qsym[a];
          B = Iia.params->colidx[b]; Bsym = Iia.params->qsym[b];
          if( (Isym==Asym) && ((Jsym^Bsym)==G_irr) )
            Z.matrix[h][row][col] += Fme.matrix[Isym][I][A] * Iia.matrix[Jsym][J][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&Z, h);
    dpd_buf4_mat_irrep_close(&Z, h);
    }

    /* Xijab += P(ij) P(ab) Z(ij,ab) */
    dpd_buf4_init(&Xijab, EOM_XI, G_irr, 0, 5, 2, 7, 0, "Xijab");
    dpd_buf4_axpy(&Z, &Xijab, 1.0);
    dpd_buf4_close(&Xijab);
    dpd_buf4_sort_axpy(&Z, EOM_XI, qprs, 2, 7, "Xijab", -1.0);
    dpd_buf4_sort_axpy(&Z, EOM_XI, pqsr, 2, 7, "Xijab", -1.0);
    dpd_buf4_sort_axpy(&Z, EOM_XI, qpsr, 2, 7, "Xijab", 1.0);
    dpd_buf4_close(&Z);

    /* XIjAb += FME(I,A) L2R1_ov(j,b) + IIA(I,A) F(j,b) */
    dpd_buf4_init(&XIjAb, EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&XIjAb, h);
      dpd_buf4_mat_irrep_rd(&XIjAb, h);
      for(row=0; row < XIjAb.params->rowtot[h]; row++) {
        i = XIjAb.params->roworb[h][row][0];
        j = XIjAb.params->roworb[h][row][1];
        I = FME.params->rowidx[i]; Isym = FME.params->psym[i];
        J = Iia.params->rowidx[j]; Jsym = Iia.params->psym[j];
        II = IIA.params->rowidx[i]; IIsym = IIA.params->psym[i];
        JJ = Fme.params->rowidx[j]; JJsym = Fme.params->psym[j];
        for(col=0; col < XIjAb.params->coltot[h^G_irr]; col++) {
          a = XIjAb.params->colorb[h^G_irr][col][0];
          b = XIjAb.params->colorb[h^G_irr][col][1];
          A = FME.params->colidx[a]; Asym = FME.params->qsym[a];
          B = Iia.params->colidx[b]; Bsym = Iia.params->qsym[b];
          AA = IIA.params->colidx[a]; AAsym = IIA.params->qsym[a];
          BB = Fme.params->colidx[b]; BBsym = Fme.params->qsym[b];
          if( (Isym==Asym) && ((Jsym^Bsym)==G_irr) )
            XIjAb.matrix[h][row][col] += FME.matrix[Isym][I][A] * Iia.matrix[Jsym][J][B];
          if( ((IIsym^AAsym)==G_irr) && (JJsym==BBsym) )
            XIjAb.matrix[h][row][col] += IIA.matrix[IIsym][II][AA] * Fme.matrix[JJsym][JJ][BB];
      }
    }
    dpd_buf4_mat_irrep_wrt(&XIjAb, h);
    dpd_buf4_mat_irrep_close(&XIjAb, h);
    }
    dpd_buf4_close(&XIjAb);

    dpd_file2_mat_close(&IIA); dpd_file2_mat_close(&Iia);
    dpd_file2_close(&IIA); dpd_file2_close(&Iia);
    dpd_file2_mat_close(&FME); dpd_file2_mat_close(&Fme);
    dpd_file2_close(&FME); dpd_file2_close(&Fme);
  }
}
