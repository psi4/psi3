#include<stdio.h>
#include<stdlib.h>
#include<libint.h>
#include<psifiles.h>
#include"defines.h"
#define EXTERN
#include"global.h"


void init_globals()
{
  int_stack = NULL;
  Shell_Data = NULL;
  int_fjttable.n1 = 0;
  int_fjttable.d = NULL;

  IOUnits.itap30 = 30;
  IOUnits.itap33 = PSIF_AO_TEI;
  IOUnits.itapS = PSIF_AO_S;
  IOUnits.itapT = PSIF_AO_T;
  IOUnits.itapV = PSIF_AO_V;
  IOUnits.itapDSCF = PSIF_DSCF;
  IOUnits.itapD = PSIF_AO_OPDM;
  IOUnits.itapG = PSIF_AO_TPDM;
  IOUnits.itapR12 = PSIF_AO_R12;
  IOUnits.itapT1  = PSIF_AO_R12T1;
  IOUnits.itapERI_MO = PSIF_MO_TEI;
  IOUnits.itapR12_MO = PSIF_MO_R12;
  IOUnits.itapR12T2_MO = PSIF_MO_R12T1;

  return;
}
