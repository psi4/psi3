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

  if (params.ground) {
    /* L1 RHS = Fia */ 
    if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
      dpd_file2_init(&Fme,CC_OEI, 0, 0, 1, "Fme");
      dpd_file2_init(&FME,CC_OEI, 0, 0, 1, "FME");
      dpd_file2_copy(&Fme, CC_OEI, "New Lia");
      dpd_file2_copy(&FME, CC_OEI, "New LIA");
      dpd_file2_close(&Fme);
      dpd_file2_close(&FME);
    }
    else if(params.ref == 2) { /** UHF **/
      dpd_file2_init(&Fme,CC_OEI, 0, 2, 3, "Fme");
      dpd_file2_init(&FME,CC_OEI, 0, 0, 1, "FME");
      dpd_file2_copy(&Fme, CC_OEI, "New Lia");
      dpd_file2_copy(&FME, CC_OEI, "New LIA");
      dpd_file2_close(&Fme);
      dpd_file2_close(&FME);
    }
  }

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_file2_init(&newLIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&newLia, CC_OEI, L_irr, 0, 1, "New Lia");
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&newLIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&newLia, CC_OEI, L_irr, 2, 3, "New Lia");
  }

  /* include no connected terms for excited states */
  if (!params.ground) {
    dpd_file2_scm(&newLIA, 0.0);
    dpd_file2_scm(&newLia, 0.0);
  }

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    /* L1 RHS += Lie*Fea */
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "Lia");

    dpd_file2_init(&LFAEt2, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&LFaet2, CC_OEI, 0, 1, 1, "Fae");
    dpd_contract222(&Lia,&LFaet2,&newLia, 0, 1, 1.0, 1.0);
    dpd_contract222(&LIA,&LFAEt2,&newLIA, 0, 1, 1.0, 1.0);
    dpd_file2_close(&LFaet2);
    dpd_file2_close(&LFAEt2);


    /* L1 RHS += -Lma*Fim */
    dpd_file2_init(&LFMIt2,CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_init(&LFmit2,CC_OEI, 0, 0, 0, "Fmi");
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
  }
  else if(params.ref == 2) { /** UHF **/

    /* L1 RHS += Lie*Fea */
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 2, 3, "Lia");

    dpd_file2_init(&LFAEt2, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_init(&LFaet2, CC_OEI, 0, 3, 3, "Faet");
    dpd_contract222(&Lia,&LFaet2,&newLia, 0, 1, 1, 1);
    dpd_contract222(&LIA,&LFAEt2,&newLIA, 0, 1, 1, 1);
    dpd_file2_close(&LFaet2);
    dpd_file2_close(&LFAEt2);

    /* L1 RHS += -Lma*Fim */
    dpd_file2_init(&LFMIt2,CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_init(&LFmit2,CC_OEI, 0, 2, 2, "Fmit");
    dpd_contract222(&LFmit2,&Lia,&newLia, 0, 1, -1, 1);
    dpd_contract222(&LFMIt2,&LIA,&newLIA, 0, 1, -1, 1);
    dpd_file2_close(&LFmit2);
    dpd_file2_close(&LFMIt2);

    /* L1 RHS += Lme*Wieam */
    dpd_buf4_init(&WMBEJ, CC_HBAR, 0, 20, 20, 20, 20, 0, "WMBEJ");
    dpd_contract422(&WMBEJ, &LIA, &newLIA, 0, 0, 1, 1);
    dpd_buf4_close(&WMBEJ);

    dpd_buf4_init(&WMbEj, CC_HBAR, 0, 20, 30, 20, 30, 0, "WMbEj");
    dpd_contract422(&WMbEj, &Lia, &newLIA, 0, 0, 1, 1);
    dpd_buf4_close(&WMbEj);

    dpd_buf4_init(&Wmbej, CC_HBAR, 0, 30, 30, 30, 30, 0, "Wmbej");
    dpd_contract422(&Wmbej, &Lia, &newLia, 0, 0, 1, 1);
    dpd_buf4_close(&Wmbej);

    dpd_buf4_init(&WmBeJ, CC_HBAR, 0, 30, 20, 30, 20, 0, "WmBeJ");
    dpd_contract422(&WmBeJ, &LIA, &newLia, 0, 0, 1, 1);
    dpd_buf4_close(&WmBeJ);

    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);

  }

  /* L1 RHS += 1/2 Limef*Wefam */
  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_buf4_init(&W, CC_HBAR, 0, 11, 7, 11, 7, 0, "WEIAB");
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WEiAb");
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&L2, &W, &newLIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&W, CC_HBAR, 0, 11, 7, 11, 7, 0, "Weiab");
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 7, 2, 7, 0, "Lijab");
    dpd_contract442(&L2, &W, &newLia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&W, CC_HBAR, 0, 11, 5, 11, 5, 0, "WeIaB");
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&L2, &W, &newLia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&W, CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WEiAb");
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&L2, &W, &newLIA, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&W, CC_HBAR, 0, 31, 17, 31, 17, 0, "Weiab");
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 10, 17, 12, 17, 0, "Lijab");
    dpd_contract442(&L2, &W, &newLia, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&W, CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&L2, &W, &newLia, 0, 0, 1, 1);
    dpd_buf4_close(&W);
    dpd_buf4_close(&L2);

  }

  /* L1 RHS += -1/2 Lmnae*Wiemn */
  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_buf4_init(&WMBIJ, CC_HBAR, 0, 10, 2, 10, 2, 0, "WMBIJ");
    dpd_buf4_init(&LIJAB, CC_LAMPS, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1.0, 1.0);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&WMBIJ);

    dpd_buf4_init(&WMbIj, CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1.0, 1.0);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&WMbIj);

    dpd_buf4_init(&Wmbij, CC_HBAR, 0, 10, 2, 10, 2, 0, "Wmbij");
    dpd_buf4_init(&Lijab, CC_LAMPS, L_irr, 2, 5, 2, 7, 0, "Lijab");
    dpd_contract442(&Wmbij, &Lijab, &newLia, 0, 2, -1.0, 1.0);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&Wmbij);

    dpd_buf4_init(&WmBiJ, CC_HBAR, 0, 10, 0, 10, 0, 0, "WmBiJ");
    dpd_buf4_init(&LiJaB, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&WmBiJ, &LiJaB, &newLia, 0, 2, -1.0, 1.0);
    dpd_buf4_close(&LiJaB);
    dpd_buf4_close(&WmBiJ);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&WMBIJ, CC_HBAR, 0, 20, 2, 20, 2, 0, "WMBIJ");
    dpd_buf4_init(&LIJAB, CC_LAMPS, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&WMBIJ, &LIJAB, &newLIA, 0, 2, -1, 1);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&WMBIJ);

    dpd_buf4_init(&WMbIj, CC_HBAR, 0, 24, 22, 24, 22, 0, "WMbIj");
    dpd_buf4_init(&LIjAb, CC_LAMPS, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&WMbIj, &LIjAb, &newLIA, 0, 2, -1, 1);
    dpd_buf4_close(&LIjAb);
    dpd_buf4_close(&WMbIj);

    dpd_buf4_init(&Wmbij, CC_HBAR, 0, 30, 12, 30, 12, 0, "Wmbij");
    dpd_buf4_init(&Lijab, CC_LAMPS, L_irr, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract442(&Wmbij, &Lijab, &newLia, 0, 2, -1, 1);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&Wmbij);

    dpd_buf4_init(&WmBiJ, CC_HBAR, 0, 27, 23, 27, 23, 0, "WmBiJ");
    dpd_buf4_init(&LiJaB, CC_LAMPS, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&WmBiJ, &LiJaB, &newLia, 0, 2, -1, 1);
    dpd_buf4_close(&LiJaB);
    dpd_buf4_close(&WmBiJ);

  }


  /* L1 RHS += -Gef*Weifa */
  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&GAE, CC_OEI, L_irr, 1, 1, "GAE");
    dpd_file2_init(&Gae, CC_OEI, L_irr, 1, 1, "Gae");

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
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&GAE, CC_OEI, L_irr, 1, 1, "GAE");
    dpd_file2_init(&Gae, CC_OEI, L_irr, 3, 3, "Gae");

    dpd_buf4_init(&W, CC_HBAR, 0, 21, 5, 21, 7, 0, "WAMEF");
    dpd_dot13(&GAE,&W,&newLIA, 0, 0, -1, 1);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, 0, 25, 29, 25, 29, 0, "WaMeF");
    dpd_dot13(&Gae,&W,&newLIA, 0, 0, -1, 1);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, 0, 31, 15, 31, 17, 0, "Wamef");
    dpd_dot13(&Gae,&W,&newLia, 0, 0, -1, 1);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, 0, 26, 28, 26, 28, 0, "WAmEf");
    dpd_dot13(&GAE,&W,&newLia, 0, 0, -1, 1);
    dpd_buf4_close(&W);

    dpd_file2_close(&Gae);
    dpd_file2_close(&GAE);

  }

  /* L1 RHS += -Gmn*Wmina */
  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&GMI, CC_OEI, L_irr, 0, 0, "GMI");
    dpd_file2_init(&Gmi, CC_OEI, L_irr, 0, 0, "Gmi");

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

  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&GMI, CC_OEI, L_irr, 0, 0, "GMI");
    dpd_file2_init(&Gmi, CC_OEI, L_irr, 2, 2, "Gmi");

    dpd_buf4_init(&W, CC_HBAR, 0, 0, 21, 2, 21, 0, "WMNIE");
    dpd_dot14(&GMI, &W, &newLIA, 0, 0, -1, 1); 
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE");
    dpd_dot14(&Gmi, &W, &newLIA, 0, 0, -1, 1);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, 0, 10, 31, 12, 31, 0, "Wmnie");
    dpd_dot14(&Gmi, &W, &newLia, 0, 0, -1, 1);
    dpd_buf4_close(&W);

    dpd_buf4_init(&W, CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe");
    dpd_dot14(&GMI, &W, &newLia, 0, 0, -1, 1);
    dpd_buf4_close(&W);

    dpd_file2_close(&Gmi);
    dpd_file2_close(&GMI);
  }

  dpd_file2_close(&newLIA);
  dpd_file2_close(&newLia);


  /* newLia * Dia */
  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&newLIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_copy(&newLIA, CC_OEI, "New LIA Increment");
    dpd_file2_close(&newLIA);

    dpd_file2_init(&newLIA, CC_OEI, L_irr, 0, 1, "New LIA Increment");
    if(params.local && local.filter_singles) local_filter_T1(&newLIA);
    else {
      dpd_file2_init(&dIA, CC_OEI, L_irr, 0, 1, "dIA");
      dpd_file2_dirprd(&dIA, &newLIA);
      dpd_file2_close(&dIA);
    }
    dpd_file2_close(&newLIA);

    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&LIA, CC_OEI, "New LIA");
    dpd_file2_close(&LIA);
    dpd_file2_init(&newLIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "New LIA Increment");
    dpd_file2_axpy(&LIA, &newLIA, 1, 0);
    dpd_file2_close(&LIA);

    dpd_file2_copy(&newLIA, CC_OEI, "New Lia");  /* spin-adaptation for RHF */
    dpd_file2_close(&newLIA);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&newLIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_copy(&newLIA, CC_OEI, "New LIA Increment");
    dpd_file2_close(&newLIA);

    dpd_file2_init(&newLIA, CC_OEI, L_irr, 0, 1, "New LIA Increment");
    dpd_file2_init(&dIA, CC_OEI, L_irr, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newLIA);
    dpd_file2_close(&dIA);
    dpd_file2_close(&newLIA);

    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "LIA");
    dpd_file2_copy(&LIA, CC_OEI, "New LIA");
    dpd_file2_close(&LIA);
    dpd_file2_init(&newLIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&LIA, CC_OEI, L_irr, 0, 1, "New LIA Increment");
    dpd_file2_axpy(&LIA, &newLIA, 1, 0);
    dpd_file2_close(&LIA);
    dpd_file2_close(&newLIA);

    dpd_file2_init(&newLia, CC_OEI, L_irr, 0, 1, "New Lia");
    dpd_file2_copy(&newLia, CC_OEI, "New Lia Increment");
    dpd_file2_close(&newLia);

    dpd_file2_init(&newLia, CC_OEI, L_irr, 0, 1, "New Lia Increment");
    dpd_file2_init(&dia, CC_OEI, L_irr, 0, 1, "dia");
    dpd_file2_dirprd(&dia, &newLia);
    dpd_file2_close(&dia);
    dpd_file2_close(&newLia);

    dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "Lia");
    dpd_file2_copy(&Lia, CC_OEI, "New Lia");
    dpd_file2_close(&Lia);
    dpd_file2_init(&newLia, CC_OEI, L_irr, 0, 1, "New Lia");
    dpd_file2_init(&Lia, CC_OEI, L_irr, 0, 1, "New Lia Increment");
    dpd_file2_axpy(&Lia, &newLia, 1, 0);
    dpd_file2_close(&Lia);
    dpd_file2_close(&newLia);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&newLIA, CC_OEI, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&dIA, CC_OEI, L_irr, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newLIA);
    dpd_file2_close(&dIA);
    dpd_file2_close(&newLIA);

    dpd_file2_init(&newLia, CC_OEI, L_irr, 2, 3, "New Lia");
    dpd_file2_init(&dia, CC_OEI, L_irr, 2, 3, "dia");
    dpd_file2_dirprd(&dia, &newLia);
    dpd_file2_close(&dia);
    dpd_file2_close(&newLia);
  }

  return;
}


