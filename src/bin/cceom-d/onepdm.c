#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* onepdm(): Computes the non-R0 parts of the unrelaxed EOM 1pdm
* intermediates are defined in intermediates.c
* 
*   D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f]
*
*   D[a][b] = +LR_vv[a][b] + t1[n][b] * L2R1_ov[n][a]
*
*   D[a][i] = +L2R1_ov[i][a]
*
*   D[i][a] = + L1R2_ov[i][a] 
*            - t1[m][a] * LR_oo[m][i] 
*            - t1[i][e] * LR_vv[e][a]
*            - r1[m][a] * LT2_oo[m][i]
*            - r1[i][e] * LT2_vv[e][a]
*            + L2R1_ov[M][E] * (t2[i][m][a][e] - t1[i][e] * t1[m][a])
*
* RAK, July 2003
*/

void onepdm(void)
{
  dpdfile2 DIA, Dia, DIJ, DAB, Dij, Dab, TIA, Tia;
  dpdfile2 LIA, Lia, RIA, Ria, I, XIJ, Xij;
  dpdbuf4 T2, L2, R2, I2;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&TIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&Tia, CC_OEI, 0, 0, 1, "tia");

    dpd_file2_init(&RIA, CC_RAMPS, 0, 0, 1, "RIA");
    dpd_file2_init(&Ria, CC_RAMPS, 0, 0, 1, "Ria");

    dpd_file2_init(&LIA, CC_LAMPS, 0, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_LAMPS, 0, 0, 1, "Lia");

    /* D[i][j] = -LR_oo[j][i] - t1[i][f] * L2R1_ov[j][f] */

    dpd_file2_init(&DIJ, EOM_D, 0, 0, 0, "DIJ");
    dpd_file2_scm(&DIJ, 0.0);

    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR_OO");
    dpd_file2_axpy(&I, &DIJ, -1.0, 1);
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_OV");
    dpd_contract222(&TIA, &I, &DIJ, 0, 0, -1.0, 1.0);
    dpd_file2_close(&I);
    dpd_file2_close(&DIJ);

    dpd_file2_init(&Dij, EOM_D, 0, 0, 0, "Dij");
    dpd_file2_scm(&Dij, 0.0);

    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR_oo");
    dpd_file2_axpy(&I, &Dij, -1.0, 1);
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_ov");
    dpd_contract222(&Tia, &I, &Dij, 0, 0, -1.0, 1.0);
    dpd_file2_close(&I);
    dpd_file2_close(&Dij);

    /* D[a][b] = +LR_vv[a][b] + L2R1_ov[n][a] * t1[n][b] */
    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR_VV");
    dpd_file2_copy(&I, EOM_D, "DAB");
    dpd_file2_close(&I);

    dpd_file2_init(&DAB, EOM_D, 0, 1, 1, "DAB");
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_OV");
    dpd_contract222(&I, &TIA, &DAB, 1, 1, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_file2_close(&DAB);

    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR_vv");
    dpd_file2_copy(&I, EOM_D, "Dab");
    dpd_file2_close(&I);

    dpd_file2_init(&Dab, EOM_D, 0, 1, 1, "Dab");
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_ov");
    dpd_contract222(&I, &Tia, &Dab, 1, 1, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_file2_close(&Dab);

    /* D[a][i] = +L2R1_ov[i][a] */

    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_OV");
    dpd_file2_copy(&I, EOM_D, "DAI");
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_ov");
    dpd_file2_copy(&I, EOM_D, "Dai");
    dpd_file2_close(&I);

    /*
       D[I][A] = + L1R2_OV[I][A] 
                 - LR_OO[M][I] * t1[M][A]
                 - t1[I][E] * LR_vv[E][A]
                 - LT2_OO[M][I] * r1[M][A]
                 - r1[I][E] * LT2_vv[E][A]
                 + L2R1_ov[M][E] * (t2[i][m][a][e] - t1[i][e] * t1[m][a])
    */

    /* D[i][a] = L1R2_ov[i][a] */
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L1R2_OV");
    dpd_file2_copy(&I, EOM_D, "DIA");
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L1R2_ov");
    dpd_file2_copy(&I, EOM_D, "Dia");
    dpd_file2_close(&I);

    dpd_file2_init(&DIA, EOM_D, 0, 0, 1, "DIA");
    dpd_file2_init(&Dia, EOM_D, 0, 0, 1, "Dia");

    /* - LR_OO[M][I] * t1[M][A] */

    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR_OO");
    dpd_contract222(&I, &TIA, &DIA, 1, 1, -1.0, 1.0);
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LR_oo");
    dpd_contract222(&I, &Tia, &Dia, 1, 1, -1.0, 1.0);
    dpd_file2_close(&I);

    /* - t1[I][E] * LR_vv[E][A] */

    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR_VV");
    dpd_contract222(&TIA, &I, &DIA, 0, 1, -1.0, 1.0);
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LR_vv");
    dpd_contract222(&Tia, &I, &Dia, 0, 1, -1.0, 1.0);
    dpd_file2_close(&I);

    /* - LT2_OO[M][I] * r1[M][A] */

    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT2_OO");
    dpd_contract222(&I, &RIA, &DIA, 1, 1, -1.0, 1.0);
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP, 0, 0, 0, "LT2_oo");
    dpd_contract222(&I, &Ria, &Dia, 1, 1, -1.0, 1.0);
    dpd_file2_close(&I);

    /* - r1[I][E] * LT2_vv[E][A] */

    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT2_VV");
    dpd_contract222(&RIA, &I, &DIA, 0, 1, -1.0, 1.0);
    dpd_file2_close(&I);

    dpd_file2_init(&I, EOM_TMP, 0, 1, 1, "LT2_vv");
    dpd_contract222(&Ria, &I, &Dia, 0, 1, -1.0, 1.0);
    dpd_file2_close(&I);

    /* + L2R1_ov[M][E] * t2[i][m][a][e] */

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB"); 
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_OV");
    dpd_dot24(&I, &T2, &DIA, 0, 0, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb"); 
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_ov");
    dpd_dot24(&I, &T2, &DIA, 0, 0, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab"); 
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_ov");
    dpd_dot24(&I, &T2, &Dia, 0, 0, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_buf4_close(&T2);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB"); 
    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_OV");
    dpd_dot24(&I, &T2, &Dia, 0, 0, 1.0, 1.0);
    dpd_file2_close(&I);
    dpd_buf4_close(&T2);
    
    /* - (t1[i][e] * L2R1_ov[M][E]) * t1[m][a] */

    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_OV");
    dpd_file2_init(&XIJ, EOM_TMP, 0, 0, 0, "XIJ");
    dpd_contract222(&TIA, &I, &XIJ, 0, 0, 1.0, 0.0);
    dpd_file2_close(&I);

    dpd_file2_init(&XIJ, EOM_TMP, 0, 0, 0, "XIJ");
    dpd_contract222(&XIJ, &TIA, &DIA, 0, 1, -1.0, 1.0);
    dpd_file2_close(&XIJ);

    dpd_file2_init(&I, EOM_TMP, 0, 0, 1, "L2R1_ov");
    dpd_file2_init(&Xij, EOM_TMP, 0, 0, 0, "Xij");
    dpd_contract222(&Tia, &I, &Xij, 0, 0, 1.0, 0.0);
    dpd_file2_close(&I);

    dpd_file2_init(&Xij, EOM_TMP, 0, 0, 0, "Xij");
    dpd_contract222(&Xij, &Tia, &Dia, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Xij);

    dpd_file2_close(&DIA);
    dpd_file2_close(&Dia);

    dpd_file2_close(&TIA);
    dpd_file2_close(&Tia);
    dpd_file2_close(&RIA);
    dpd_file2_close(&Ria);
    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
  }
}
