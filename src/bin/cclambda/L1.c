#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void L1_build(void) {
  dpdfile2 newLIA, newLia, LIA, Lia;
  dpdfile2 dIA, dia, Fme, FME;
  dpdfile2 LFaet2, LFAEt2, LFmit2, LFMIt2;
  dpdfile2 GMI, Gmi, Gae;
  dpdfile2 GAE;
  dpdbuf4 WMBEJ, Wmbej, WMbEj, WmBeJ;
  dpdbuf4 WMBIJ, Wmbij, WMbIj, WmBiJ;
  dpdbuf4 LIJAB, Lijab, LIjAb, LiJaB, L2;
  dpdbuf4 WMNIE, Wmnie, WMnIe, WmNiE;
  dpdbuf4 WAMEF, Wamef, WAmEf, WaMeF, W;


 /* L1 RHS = Fia */ 
  dpd_file2_init(&Fme,CC_OEI, 0, 0, 1, "Fme");
  dpd_file2_init(&FME,CC_OEI, 0, 0, 1, "FME");
  dpd_file2_copy(&Fme, CC_OEI, "New Lia");
  dpd_file2_copy(&FME, CC_OEI, "New LIA");
  dpd_file2_close(&Fme);
  dpd_file2_close(&FME);

  dpd_file2_init(&newLIA, CC_OEI, 0, 0, 1, "New LIA");
  dpd_file2_init(&newLia, CC_OEI, 0, 0, 1, "New Lia");


 /* L1 RHS += Lie*Fea */
  dpd_file2_init(&LIA, CC_OEI, 0, 0, 1, "LIA");
  dpd_file2_init(&Lia, CC_OEI, 0, 0, 1, "Lia");

  dpd_file2_init(&LFAEt2, CC_OEI, 0, 1, 1, "FAEt");
  dpd_file2_init(&LFaet2, CC_OEI, 0, 1, 1, "Faet");
  dpd_contract222(&Lia,&LFaet2,&newLia, 0, 1, 1.0, 1.0);
  dpd_contract222(&LIA,&LFAEt2,&newLIA, 0, 1, 1.0, 1.0);
  dpd_file2_close(&LFaet2);
  dpd_file2_close(&LFAEt2);


 /* L1 RHS += -Lma*Fim */
  dpd_file2_init(&LFMIt2,CC_OEI, 0, 0, 0, "FMIt");
  dpd_file2_init(&LFmit2,CC_OEI, 0, 0, 0, "Fmit");
  dpd_contract222(&LFmit2,&Lia,&newLia, 0, 1, -1.0, 1.0);
  dpd_contract222(&LFMIt2,&LIA,&newLIA, 0, 1, -1.0, 1.0);
  dpd_file2_close(&LFmit2);
  dpd_file2_close(&LFMIt2);


 /* L1 RHS += Lme*Wieam */
  dpd_buf4_init(&WMBEJ, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMBEJ");
  dpd_contract422(&WMBEJ, &LIA, &newLIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&WMBEJ);

  dpd_buf4_init(&WMbEj, CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  dpd_contract422(&WMbEj, &Lia, &newLIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&WMbEj);

  dpd_buf4_init(&Wmbej, CC_HBAR, 0, 10, 10, 10, 10, 0, "Wmbej");
  dpd_contract422(&Wmbej, &Lia, &newLia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&Wmbej);

  dpd_buf4_init(&WmBeJ, CC_HBAR, 0, 10, 10, 10, 10, 0, "WmBeJ");
  dpd_contract422(&WmBeJ, &LIA, &newLia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&WmBeJ);

  dpd_file2_close(&LIA);
  dpd_file2_close(&Lia);

 /* L1 RHS += 1/2 Limef*Wefam */
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 7, 11, 7, 0, "WEIAB");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB");
  dpd_contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);

  dpd_buf4_init(&W, CC_HBAR, 0, 11, 7, 11, 7, 0, "Weiab");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 7, 2, 7, 0, "Lijab");
  dpd_contract442(&L2, &W, &newLia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LiJaB");
  dpd_contract442(&L2, &W, &newLia, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&W);
  dpd_buf4_close(&L2);

 /* L1 RHS += -1/2 Lmnae*Wiemn */
  dpd_buf4_init(&WMBIJ, CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
  dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "LIJAB");
  dpd_contract442(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1.0, 1.0);
  dpd_buf4_close(&LIJAB);
  dpd_buf4_close(&WMBIJ);

  dpd_buf4_init(&WMbIj, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
  dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1.0, 1.0);
  dpd_buf4_close(&LIjAb);
  dpd_buf4_close(&WMbIj);

  dpd_buf4_init(&Wmbij, CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
  dpd_buf4_init(&Lijab, CC_LAMPS, 0, 2, 5, 2, 7, 0, "Lijab");
  dpd_contract442(&Wmbij, &Lijab, &newLia, 0, 2, -1.0, 1.0);
  dpd_buf4_close(&Lijab);
  dpd_buf4_close(&Wmbij);

  dpd_buf4_init(&WmBiJ, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
  dpd_buf4_init(&LiJaB, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LiJaB");
  dpd_contract442(&WmBiJ, &LiJaB, &newLia, 0, 2, -1.0, 1.0);
  dpd_buf4_close(&LiJaB);
  dpd_buf4_close(&WmBiJ);


 /* L1 RHS += -Gef*Weifa */

  dpd_file2_init(&GAE, CC_OEI, 0, 1, 1, "GAE");
  dpd_file2_init(&Gae, CC_OEI, 0, 1, 1, "Gae");

  dpd_buf4_init(&WAMEF, CC_HBAR, 0, 10, 5, 10, 7, 0, "WAMEF");
  dpd_dot23(&GAE,&WAMEF,&newLIA, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&WAMEF);

  dpd_buf4_init(&WaMeF, CC_HBAR, 0, 10, 5, 10, 5, 0, "WaMeF");
  dpd_dot23(&Gae,&WaMeF,&newLIA, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&WaMeF);

  dpd_buf4_init(&Wamef, CC_HBAR, 0, 10, 5, 10, 7, 0, "Wamef");
  dpd_dot23(&Gae,&Wamef,&newLia, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&Wamef);

  dpd_buf4_init(&WAmEf, CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf");
  dpd_dot23(&GAE,&WAmEf,&newLia, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&WAmEf);

  dpd_file2_close(&Gae);
  dpd_file2_close(&GAE);

 /* L1 RHS += -Gmn*Wmina */

  dpd_file2_init(&GMI, CC_OEI, 0, 0, 0, "GMI");
  dpd_file2_init(&Gmi, CC_OEI, 0, 0, 0, "Gmi");

  dpd_buf4_init(&WMNIE, CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE");
  dpd_dot14(&GMI, &WMNIE, &newLIA, 0, 0, -1.0, 1.0); 
  dpd_buf4_close(&WMNIE);

  dpd_buf4_init(&WmNiE, CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE");
  dpd_dot14(&Gmi, &WmNiE, &newLIA, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&WmNiE);

  dpd_buf4_init(&Wmnie, CC_HBAR, 0, 0, 11, 2, 11, 0, "Wmnie");
  dpd_dot14(&Gmi, &Wmnie, &newLia, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&Wmnie);

  dpd_buf4_init(&WMnIe, CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe");
  dpd_dot14(&GMI, &WMnIe, &newLia, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&WMnIe);

  dpd_file2_close(&Gmi);
  dpd_file2_close(&GMI);


  /* newLia * Dia */
  dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&dIA, &newLIA);
  dpd_file2_close(&dIA);

  dpd_file2_init(&dia, CC_OEI, 0, 0, 1, "dia");
  dpd_file2_dirprd(&dia, &newLia);
  dpd_file2_close(&dia);


  dpd_file2_close(&newLIA);
  dpd_file2_close(&newLia);

  return;
}


