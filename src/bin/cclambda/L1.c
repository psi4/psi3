#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void L1_build(void) {
  struct oe_dpdfile newLIA, newLia, LIA, Lia;
  struct oe_dpdfile dIA, dia, Fme, FME;
  struct oe_dpdfile LFaet2, LFAEt2, LFmit2, LFMIt2;
  struct oe_dpdfile GMI, Gmi, Gae;
  struct oe_dpdfile GAE;
  struct dpdbuf WMBEJ, Wmbej, WMbEj, WmBeJ;
  struct dpdbuf WMBIJ, Wmbij, WMbIj, WmBiJ;
  struct dpdbuf LIJAB, Lijab, LIjAb, LiJaB, L2;
  struct dpdbuf WMNIE, Wmnie, WMnIe, WmNiE;
  struct dpdbuf WAMEF, Wamef, WAmEf, WaMeF, W;


 /* L1 RHS = Fia */ 
  dpd_oe_file_init(&Fme,CC_OEI,0,1,"Fme",0,outfile);
  dpd_oe_file_init(&FME,CC_OEI,0,1,"FME",0,outfile);
  dpd_oe_copy(&Fme, CC_OEI, "New Lia", 0, outfile);
  dpd_oe_copy(&FME, CC_OEI, "New LIA", 0, outfile);
  dpd_oe_file_close(&Fme);
  dpd_oe_file_close(&FME);

  dpd_oe_file_init(&newLIA, CC_OEI, 0, 1, "New LIA", 0, outfile);
  dpd_oe_file_init(&newLia, CC_OEI, 0, 1, "New Lia", 0, outfile);


 /* L1 RHS += Lie*Fea */
  dpd_oe_file_init(&LIA, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_init(&Lia, CC_OEI, 0, 1, "Lia", 0, outfile);

  dpd_oe_file_init(&LFAEt2, CC_OEI, 1, 1, "FAEt", 0, outfile);
  dpd_oe_file_init(&LFaet2, CC_OEI, 1, 1, "Faet", 0, outfile);
  dpd_contract111(&Lia,&LFaet2,&newLia, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_contract111(&LIA,&LFAEt2,&newLIA, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&LFaet2);
  dpd_oe_file_close(&LFAEt2);


 /* L1 RHS += -Lma*Fim */
  dpd_oe_file_init(&LFMIt2,CC_OEI,0,0,"FMIt",0,outfile);
  dpd_oe_file_init(&LFmit2,CC_OEI,0,0,"Fmit",0,outfile);
  dpd_contract111(&LFmit2,&Lia,&newLia, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_contract111(&LFMIt2,&LIA,&newLIA, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&LFmit2);
  dpd_oe_file_close(&LFMIt2);


 /* L1 RHS += Lme*Wieam */
  dpd_buf_init(&WMBEJ, CC_HBAR, 10, 10, 10, 10, 0, "WMBEJ", 
               0, outfile);
  dpd_contract121(&WMBEJ, &LIA, &newLIA, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&WMBEJ);

  dpd_buf_init(&WMbEj, CC_HBAR, 10, 10, 10, 10, 0, "WMbEj", 
               0, outfile);
  dpd_contract121(&WMbEj, &Lia, &newLIA, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&WMbEj);

  dpd_buf_init(&Wmbej, CC_HBAR, 10, 10, 10, 10, 0, "Wmbej", 
               0, outfile);
  dpd_contract121(&Wmbej, &Lia, &newLia, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&Wmbej);

  dpd_buf_init(&WmBeJ, CC_HBAR, 10, 10, 10, 10, 0, "WmBeJ", 
               0, outfile);
  dpd_contract121(&WmBeJ, &LIA, &newLia, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&WmBeJ);

  dpd_oe_file_close(&LIA);
  dpd_oe_file_close(&Lia);


 /* L1 RHS += 1/2 Limef*Wefam */
  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "WEIAB", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_contract122(&L2, &W, &newLIA, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WEiAb", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_contract122(&L2, &W, &newLIA, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);

  dpd_buf_init(&W, CC_HBAR, 11, 7, 11, 7, 0, "Weiab", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_contract122(&L2, &W, &newLia, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);
  dpd_buf_init(&W, CC_HBAR, 11, 5, 11, 5, 0, "WeIaB", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_contract122(&L2, &W, &newLia, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&W);
  dpd_buf_close(&L2);

 /* L1 RHS += -1/2 Lmnae*Wiemn */
  dpd_buf_init(&WMBIJ, CC_HBAR, 10, 2, 10, 2, 0, "WMBIJ", 0, outfile);
  dpd_buf_init(&LIJAB, CC_LAMPS,  2, 5,  2, 7, 0, "LIJAB", 0, outfile);
  dpd_contract122(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&LIJAB);
  dpd_buf_close(&WMBIJ);

  dpd_buf_init(&WMbIj, CC_HBAR, 10, 0, 10, 0, 0, "WMbIj", 0, outfile);
  dpd_buf_init(&LIjAb, CC_LAMPS,  0, 5,  0, 5, 0, "LIjAb", 0, outfile);
  dpd_contract122(&WMbIj, &LIjAb, &newLIA, 0, 2, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&LIjAb);
  dpd_buf_close(&WMbIj);

  dpd_buf_init(&Wmbij, CC_HBAR, 10, 2, 10, 2, 0, "Wmbij", 0, outfile);
  dpd_buf_init(&Lijab, CC_LAMPS,  2, 5,  2, 7, 0, "Lijab", 0, outfile);
  dpd_contract122(&Wmbij, &Lijab, &newLia, 0, 2, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&Lijab);
  dpd_buf_close(&Wmbij);

  dpd_buf_init(&WmBiJ, CC_HBAR, 10, 0, 10, 0, 0, "WmBiJ", 0, outfile);
  dpd_buf_init(&LiJaB, CC_LAMPS,  0, 5,  0, 5, 0, "LiJaB", 0, outfile);
  dpd_contract122(&WmBiJ, &LiJaB, &newLia, 0, 2, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&LiJaB);
  dpd_buf_close(&WmBiJ);


 /* L1 RHS += -Gef*Weifa */

  dpd_oe_file_init(&GAE, CC_OEI, 1, 1, "GAE", 0, outfile);
  dpd_oe_file_init(&Gae, CC_OEI, 1, 1, "Gae", 0, outfile);

  dpd_buf_init(&WAMEF, CC_HBAR, 10, 5, 10, 7, 0, "WAMEF", 0, outfile);
  dpd_dot23(&GAE,&WAMEF,&newLIA, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&WAMEF);

  dpd_buf_init(&WaMeF, CC_HBAR, 10, 5, 10, 5, 0, "WaMeF", 0, outfile);
  dpd_dot23(&Gae,&WaMeF,&newLIA, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&WaMeF);

  dpd_buf_init(&Wamef, CC_HBAR, 10, 5, 10, 7, 0, "Wamef", 0, outfile);
  dpd_dot23(&Gae,&Wamef,&newLia, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&Wamef);

  dpd_buf_init(&WAmEf, CC_HBAR, 10, 5, 10, 5, 0, "WAmEf", 0, outfile);
  dpd_dot23(&GAE,&WAmEf,&newLia, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&WAmEf);

  dpd_oe_file_close(&Gae);
  dpd_oe_file_close(&GAE);


 /* L1 RHS += -Gmn*Wmina */

  dpd_oe_file_init(&GMI, CC_OEI, 0, 0, "GMI", 0, outfile);
  dpd_oe_file_init(&Gmi, CC_OEI, 0, 0, "Gmi", 0, outfile);

  dpd_buf_init(&WMNIE, CC_HBAR, 0, 11, 2, 11, 0, "WMNIE", 0, outfile);
  dpd_dot14(&GMI, &WMNIE, &newLIA, 0, 0, -1.0, 1.0, 0, outfile); 
  dpd_buf_close(&WMNIE);

  dpd_buf_init(&WmNiE, CC_HBAR, 0, 11, 0, 11, 0, "WmNiE", 0, outfile);
  dpd_dot14(&Gmi, &WmNiE, &newLIA, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&WmNiE);

  dpd_buf_init(&Wmnie, CC_HBAR, 0, 11, 2, 11, 0, "Wmnie", 0, outfile);
  dpd_dot14(&Gmi, &Wmnie, &newLia, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&Wmnie);

  dpd_buf_init(&WMnIe, CC_HBAR, 0, 11, 0, 11, 0, "WMnIe", 0, outfile);
  dpd_dot14(&GMI, &WMnIe, &newLia, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&WMnIe);

  dpd_oe_file_close(&Gmi);
  dpd_oe_file_close(&GMI);


  /* newLia * Dia */
  dpd_oe_file_init(&dIA, CC_OEI, 0, 1, "dIA", 0, outfile);
  dpd_oe_dirprd(&dIA, &newLIA, 0, outfile);
  dpd_oe_file_close(&dIA);

  dpd_oe_file_init(&dia, CC_OEI, 0, 1, "dia", 0, outfile);
  dpd_oe_dirprd(&dia, &newLia, 0, outfile);
  dpd_oe_file_close(&dia);

  dpd_oe_file_close(&newLIA);
  dpd_oe_file_close(&newLia);


  return;
}


