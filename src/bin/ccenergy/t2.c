#include <stdio.h>
#define EXTERN
#include "globals.h"

void DT2(void), FaetT2(void), FmitT2(void), WmnijT2(void), WmbejT2(void);
void BT2(void), ZT2(void), FT2(void), ET2(void), CT2(void), dijabT2(void);

void t2_build(void)
{
  DT2();
  FaetT2(); 
  FmitT2();
  timer_on("WmnijT2", outfile);
  WmnijT2();
  timer_off("WmnijT2", outfile);
  timer_on("BT2", outfile);
  BT2();
  timer_off("BT2", outfile);
  timer_on("ZT2", outfile);
  ZT2();
  timer_off("ZT2", outfile);
  timer_on("FT2", outfile);
  FT2();
  timer_off("FT2", outfile);
  timer_on("ET2", outfile);
  ET2();
  timer_off("ET2", outfile);
  timer_on("WmbejT2", outfile);
  WmbejT2();
  timer_off("WmbejT2", outfile);
  timer_on("CT2", outfile);
  CT2();
  timer_off("CT2", outfile);
  dijabT2();
}
