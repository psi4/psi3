#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void x_te_intermediates(void)
{
  dpdfile2 R1, L1;
  dpdbuf4 V, L, R, T2;
  int G_irr, L_irr, R_irr;
  G_irr = params.G_irr; L_irr = params.L_irr; R_irr = params.R_irr;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    /* R2L2_OOOO = 0.5 Rijef Lklef */
    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_OOOO");
    dpd_buf4_init(&R, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&L, CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_oooo");
    dpd_buf4_init(&R, CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
    dpd_buf4_init(&L, CC_GL, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, 0, 0, 0, 0, 0, 0, "R2L2_OoOo");
    dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, 0, 0, 0, 0, 0, 0, "R2L2_oOoO");
    dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    /* Tau2L2_OOOO = 0.5 Tau_ijef Lklef */
    dpd_buf4_init(&V, EOM_TMP, L_irr, 2, 2, 2, 2, 0, "Tau2L2_OOOO");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_buf4_init(&L, CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, EOM_TMP, L_irr, 2, 2, 2, 2, 0, "Tau2L2_oooo");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_buf4_init(&L, CC_GL, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_OoOo");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_oOoO");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract444(&T2, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&V);

    /* R2L2_OVOV = Rimae Ljmbe */
    dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "Riajb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_contract444(&R, &L, &V, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "RIAJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_contract444(&R, &L, &V, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "RiaJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "Riajb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "RjAIb");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LIbjA");
    dpd_contract444(&R, &L, &V, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
    dpd_buf4_init(&R, CC_GR, R_irr, 10, 10, 10, 10, 0, "RIbjA");
    dpd_buf4_init(&L, CC_GL, L_irr, 10, 10, 10, 10, 0, "LjAIb");
    dpd_contract444(&R, &L, &V, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    /* L2R1_OOVO = Lijae Rke */
    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_OOVO");
    dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_oovo");
    dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjaB");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract424(&L, &R1, &V, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_OOVO");
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 2, 10, "L2R1_OOVO(pqsr)");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_oovo");
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 2, 10, "L2R1_oovo(pqsr)");
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 0, 10, "L2R1_OoVo(pqsr)");
    dpd_buf4_sort(&V, EOM_TMP, qprs, 0, 11, "L2R1_OoVo(qprs)");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OoVo(pqsr)");
    dpd_buf4_sort(&V, EOM_TMP, qprs, 0, 10, "L2R1_OoVo(qpsr)");
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO");
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 0, 10, "L2R1_OovO(pqsr)");
    dpd_buf4_sort(&V, EOM_TMP, qprs, 0, 11, "L2R1_OovO(qprs)");
    dpd_buf4_close(&V);

    /* L1R2_OOVO = Rijae Lke */
    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L1R2_OOVO");
    dpd_buf4_init(&R, CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
    dpd_contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L1R2_oovo");
    dpd_buf4_init(&R, CC_GR, R_irr, 2, 5, 2, 7, 0, "Rijab");
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "Lia");
    dpd_contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L1R2_OoVo");
    dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "Lia");
    dpd_contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L1R2_OovO");
    dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjaB");
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
    dpd_contract424(&R, &L1, &V, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L1R2_OOVO");
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 2, 10, "L1R2_OOVO(pqsr)");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L1R2_oovo");
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 2, 10, "L1R2_oovo(pqsr)");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L1R2_OoVo");
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 0, 10, "L1R2_OoVo(pqsr)");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L1R2_OoVo(pqsr)");
    dpd_buf4_sort(&V, EOM_TMP, qprs, 0, 10, "L1R2_OoVo(qpsr)");
    dpd_buf4_close(&V);
    dpd_buf4_init(&V, EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L1R2_OovO");
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 0, 10, "L1R2_OovO(pqsr)");
    dpd_buf4_close(&V);

    /* L2R1_VVOV = Limab Rmc */
    dpd_buf4_init(&V, EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_VVOV");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L);
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 7, 11, "L2R1_VVOV(pqsr)");
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_vvov");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 7, 2, 7, 0, "Lijab");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L);
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 7, 11, "L2R1_vvov(pqsr)");
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 5, 10, 5, 10, 0, "L2R1_VvOv");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
    dpd_contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L);
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 5, 11, "L2R1_VvOv(pqsr)");
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 5, 10, 5, 10, 0, "L2R1_VvoV");
    dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJAb");
    dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
    dpd_contract424(&L, &R1, &V, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&R1);
    dpd_buf4_close(&L);
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 5, 11, "L2R1_VvoV(pqsr)");
    dpd_buf4_close(&V);

    /* R2L1_VVOV = Rimab Lmc */
    dpd_buf4_init(&V, EOM_TMP, G_irr, 7, 10, 7, 10, 0, "R2L1_VVOV");
    dpd_buf4_init(&R, CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
    dpd_contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R);
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 7, 11, "R2L1_VVOV(pqsr)");
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 7, 10, 7, 10, 0, "R2L1_vvov");
    dpd_buf4_init(&R, CC_GR, R_irr, 0, 7, 2, 7, 0, "Rijab");
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "Lia");
    dpd_contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R);
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 7, 11, "R2L1_vvov(pqsr)");
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 5, 10, 5, 10, 0, "R2L1_VvOv");
    dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "Lia");
    dpd_contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R);
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 5, 11, "R2L1_VvOv(pqsr)");
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 5, 10, 5, 10, 0, "R2L1_VvoV");
    dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJAb");
    dpd_file2_init(&L1, CC_GL, L_irr, 0, 1, "LIA");
    dpd_contract424(&R, &L1, &V, 1, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&R);
    dpd_buf4_sort(&V, EOM_TMP, pqsr, 5, 11, "R2L1_VvoV(pqsr)");
    dpd_buf4_close(&V);

  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_OOOO");
    dpd_buf4_init(&R, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
    dpd_buf4_init(&L, CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 12, 12, 12, 12, 0, "R2L2_oooo");
    dpd_buf4_init(&R, CC_GR, R_irr, 12, 17, 12, 17, 0, "Rijab");
    dpd_buf4_init(&L, CC_GL, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 22, 22, 22, 22, 0, "R2L2_OoOo");
    dpd_buf4_init(&R, CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
    dpd_buf4_init(&L, CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
    dpd_buf4_init(&R, CC_GR, R_irr, 20, 20, 20, 20, 0, "RIAJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_init(&R, CC_GR, R_irr, 20, 30, 20, 30, 0, "RIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
    dpd_buf4_init(&R, CC_GR, R_irr, 30, 30, 30, 30, 0, "Riajb");
    dpd_buf4_init(&L, CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_init(&R, CC_GR, R_irr, 20, 30, 20, 30, 0, "RIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_contract444(&R, &L, &V, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
    dpd_buf4_init(&R, CC_GR, R_irr, 20, 30, 20, 30, 0, "RIAjb");
    dpd_buf4_init(&L, CC_GL, L_irr, 30, 30, 30, 30, 0, "Liajb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_init(&R, CC_GR, R_irr, 20, 20, 20, 20, 0, "RIAJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_contract444(&R, &L, &V, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
    dpd_buf4_init(&R, CC_GR, R_irr, 30, 20, 30, 20, 0, "RiaJB");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 20, 20, 20, 0, "LIAJB");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_init(&R, CC_GR, R_irr, 30, 30, 30, 30, 0, "Riajb");
    dpd_buf4_init(&L, CC_GL, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_contract444(&R, &L, &V, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&L);
    dpd_buf4_close(&R);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
    dpd_buf4_init(&L, CC_GL, L_irr, 24, 27, 24, 27, 0, "LIbjA");
    dpd_buf4_init(&R, CC_GR, R_irr, 27, 24, 27, 24, 0, "RjAIb");
    dpd_contract444(&R, &L, &V, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);

    dpd_buf4_init(&V, EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
    dpd_buf4_init(&L, CC_GL, L_irr, 27, 24, 27, 24, 0, "LjAIb");
    dpd_buf4_init(&R, CC_GR, R_irr, 24, 27, 24, 27, 0, "RIbjA");
    dpd_contract444(&R, &L, &V, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&L);
    dpd_buf4_close(&V);
  }
}
