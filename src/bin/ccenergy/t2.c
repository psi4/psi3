#include <stdio.h>
#include <stdlib.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void DT2(void), FaetT2(void), FmitT2(void), WmnijT2(void), WmbejT2(void);
void BT2(void), ZT2(void), FT2(void), ET2(void), CT2(void), dijabT2(void);
void BT2_AO(void);
void status(char *, FILE *);

void t2_build(void)
{

  DT2();
  if(params.print & 2) status("<ij||ab> -> T2", outfile);

  FaetT2();
  FmitT2();
  if(params.print & 2) status("F -> T2", outfile);

  timer_on("WmnijT2", outfile);
  WmnijT2();
  if(params.print & 2) status("Wmnij -> T2", outfile);
  timer_off("WmnijT2", outfile);

  timer_on("BT2", outfile);
  if(params.aobasis) BT2_AO();
  else BT2();
  if(params.print & 2) status("<ab||cd> -> T2", outfile);
  timer_off("BT2", outfile);

  timer_on("ZT2", outfile);
  ZT2();
  if(params.print & 2) status("Z -> T2", outfile);
  timer_off("ZT2", outfile);

  timer_on("FT2", outfile);
  FT2();
  if(params.print & 2) status("<ia||bc> -> T2", outfile);
  timer_off("FT2", outfile);

  timer_on("ET2", outfile);
  ET2();
  if(params.print & 2) status("<ij||ka> -> T2", outfile);
  timer_off("ET2", outfile);

  timer_on("WmbejT2", outfile);
  WmbejT2();
  if(params.print & 2) status("Wmbej -> T2", outfile);
  timer_off("WmbejT2", outfile);

  timer_on("CT2", outfile);
  CT2();
  if(params.print & 2) status("<ia||jb> -> T2", outfile);
  timer_off("CT2", outfile);

  dijabT2();
}
