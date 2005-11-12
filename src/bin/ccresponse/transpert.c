#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#define EXTERN
#include <psifiles.h>
#include "globals.h"

void transpert(char *pert, double sign)
{
  int nao, nso, nmo, noei_ao;
  int alpha;
  int i, j, ij;
  double *scratch, **TMP, **X, **target;
  char *name;
  double prefactor, anti;

  nao = moinfo.nao;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  noei_ao = moinfo.noei_ao;

  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);
  scratch = init_array(noei_ao);

  if(!strcmp(pert,"Mu")) { prefactor = 1.0; anti = 1.0; sign = 1.0; }
  else if(!strcmp(pert, "L")) { prefactor = -0.5; anti = -1.0; }
  else if(!strcmp(pert, "P")) { prefactor = -1.0; anti = -1.0; }

  for(alpha=0; alpha < 3; alpha++) {

    target = block_matrix(nmo,nmo);

    if(!strcmp(pert,"Mu")) {
      if(alpha == 0) { name = PSIF_AO_MX; moinfo.MUX = target; }
      else if(alpha == 1) { name = PSIF_AO_MY; moinfo.MUY = target; }
      else if(alpha == 2) { name = PSIF_AO_MZ; moinfo.MUZ = target; }
    }
    else if(!strcmp(pert,"L")) {
      if(alpha == 0) { name = PSIF_AO_LX; moinfo.LX = target; }
      else if(alpha == 1) { name = PSIF_AO_LY; moinfo.LY = target; }
      else if(alpha == 2) { name = PSIF_AO_LZ; moinfo.LZ = target; }
    }
    else if(!strcmp(pert,"P")) {
      if(alpha == 0) { name = PSIF_AO_NablaX; moinfo.PX = target; }
      else if(alpha == 1) { name = PSIF_AO_NablaY; moinfo.PY = target; }
      else if(alpha == 2) { name = PSIF_AO_NablaZ; moinfo.PZ = target; }
    }

    iwl_rdone(PSIF_OEI, name, scratch, noei_ao, 0, 0, outfile);
    for(i=0,ij=0; i < nao; i++)
      for(j=0; j <= i; j++,ij++) {
	TMP[i][j] = prefactor * sign * scratch[ij];
	TMP[j][i] = anti * prefactor * sign * scratch[ij];
      }

    C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(moinfo.usotao[0][0]),nao,
	    0,&(X[0][0]),nao);
    C_DGEMM('n','n',nso,nso,nao,1,&(moinfo.usotao[0][0]),nao,&(X[0][0]),nao,
	    0,&(TMP[0][0]),nao);

    C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
	    0,&(X[0][0]),nao);
    C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
	    0,&(target[0][0]),nmo);

    zero_arr(scratch,noei_ao);

  }

  free(scratch);
  free_block(TMP);
  free_block(X);
}
