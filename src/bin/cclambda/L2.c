#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void DL2(struct L_Params L_params);
void FaeL2(int L_irr);
void FmiL2(int L_irr);
void WijmnL2(int L_irr);
void WefabL2(int L_irr);
void WejabL2(int L_irr);
void WijmbL2(int L_irr);
void L1FL2(int L_irr);
void WmbejL2(int L_irr);
void GaeL2(int L_irr);
void GmiL2(int L_irr);
void dijabL2(int L_irr);

void BL2_AO(int L_irr);

void L2_build(struct L_Params L_params) {
  int L_irr;
  L_irr = L_params.irrep;

  DL2(L_params);
#ifdef EOM_DEBUG
check_sum("DL2", L_irr);
#endif

  FaeL2(L_irr);
#ifdef EOM_DEBUG
check_sum("FaeL2", L_irr);
#endif

  FmiL2(L_irr);
#ifdef EOM_DEBUG
check_sum("FmiL2", L_irr);
#endif

  WijmnL2(L_irr);
#ifdef EOM_DEBUG
check_sum("WijmnL2", L_irr);
#endif

  WefabL2(L_irr);
#ifdef EOM_DEBUG
check_sum("WefabL2", L_irr);
#endif

  WejabL2(L_irr);
#ifdef EOM_DEBUG
check_sum("WejabL2", L_irr);
#endif

  WijmbL2(L_irr);
#ifdef EOM_DEBUG
check_sum("WijmbL2", L_irr);
#endif

  WmbejL2(L_irr); 
#ifdef EOM_DEBUG
check_sum("WmbejL2", L_irr);
#endif

  L1FL2(L_irr);
#ifdef EOM_DEBUG
check_sum("L1FL2", L_irr);
#endif

  GaeL2(L_irr);
#ifdef EOM_DEBUG
check_sum("GaeL2", L_irr);
#endif

  GmiL2(L_irr);
#ifdef EOM_DEBUG
check_sum("GmiL2", L_irr);
#endif

  dijabL2(L_irr);
#ifdef EOM_DEBUG
  check_sum("after D2s", L_irr);
#endif
}

