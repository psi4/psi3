#include<stdio.h>
#include<stdlib.h>
#include<libint.h>
#include<psifiles.h>
#include"defines.h"
#define EXTERN
#include"global.h"


void init_globals()
{
  IOUnits.itap30 = 30;
  IOUnits.itap33 = PSIF_SO_TEI;
  IOUnits.itapS = PSIF_SO_S;
  IOUnits.itapT = PSIF_SO_T;
  IOUnits.itapV = PSIF_SO_V;
  IOUnits.itapS_AO = PSIF_AO_S;
  IOUnits.itapMX_AO = PSIF_AO_MX;
  IOUnits.itapMY_AO = PSIF_AO_MY;
  IOUnits.itapMZ_AO = PSIF_AO_MZ;
  IOUnits.itapDSCF = PSIF_DSCF;
  IOUnits.itapD = PSIF_AO_OPDM;
  IOUnits.itapG = PSIF_AO_TPDM;
  IOUnits.itapR12 = PSIF_SO_R12;
  IOUnits.itapT1  = PSIF_SO_R12T1;
  IOUnits.itapERI_MO = PSIF_MO_TEI;
  IOUnits.itapR12_MO = PSIF_MO_R12;
  IOUnits.itapR12T2_MO = PSIF_MO_R12T1;

  return;
}
