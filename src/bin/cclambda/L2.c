#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void DL2(void);
void FaeL2(void);
void FmiL2(void);
void WijmnL2(void);
void WefabL2(void);
void WejabL2(void);
void WijmbL2(void);
void L1FL2(void);
void WmbejL2(void);
void GaeL2(void);
void GmiL2(void);
void dijabL2(void);

void BL2_AO(void);

void L2_build(void) {

  DL2();
  FaeL2();
  FmiL2();
  WijmnL2();
  if(params.aobasis) BL2_AO();
  WefabL2();
  WejabL2();
  WijmbL2();
  WmbejL2(); 
  L1FL2();
  GaeL2();
  GmiL2();
  dijabL2();
}

