#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include "MOInfo.h"
#define EXTERN
#include <psifiles.h>
#include "globals.h"

void transL(void)
{
  int nao, noei_ao, nso, nmo;
  int i, j;
  double **TMP, **X;
  double **LX, **LY, **LZ;  /* MO-basis dipole integrals */

  nao = moinfo.nao;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  noei_ao = moinfo.nao*(moinfo.nao+1)/2;

  /**** Transform the magnetic dipole integrals to the MO basis ****/

  LX = block_matrix(nmo,nmo);
  LY = block_matrix(nmo,nmo);
  LZ = block_matrix(nmo,nmo);

  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);

  /* NB: (a|m|b) = -1/2 (a|L|b) */
  /* NB: The angular momentum integrals are antisymmetric! */
  psio_open(PSIF_OEI, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_OEI, PSIF_AO_LX, (char *) &(TMP[0][0]), nao*nao*sizeof(double));
  psio_close(PSIF_OEI, 1);
  for(i=0; i < nao; i++)
    for(j=0; j < nao; j++)
      TMP[i][j] *= -0.5;

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(moinfo.usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(moinfo.usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(LX[0][0]),nmo);

  psio_open(PSIF_OEI, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_OEI, PSIF_AO_LY, (char *) &(TMP[0][0]), nao*nao*sizeof(double));
  psio_close(PSIF_OEI, 1);
  for(i=0; i < nao; i++)
    for(j=0; j < nao; j++)
      TMP[i][j] *= -0.5;

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(moinfo.usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(moinfo.usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(LY[0][0]),nmo);

  psio_open(PSIF_OEI, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_OEI, PSIF_AO_LZ, (char *) &(TMP[0][0]), nao*nao*sizeof(double));
  psio_close(PSIF_OEI, 1);
  for(i=0; i < nao; i++)
    for(j=0; j < nao; j++)
      TMP[i][j] *= -0.5;

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(moinfo.usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(moinfo.usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(moinfo.scf[0][0]),nmo,
	  0,&(X[0][0]),nao);
  C_DGEMM('t','n',nmo,nmo,nso,1,&(moinfo.scf[0][0]),nmo,&(X[0][0]),nao,
	  0,&(LZ[0][0]),nmo);

  free_block(TMP);
  free_block(X);

/*   printf("MO-Basis LX Integrals:\n"); */
/*   mat_print(LX, nmo, nmo, stdout); */
/*   printf("MO-Basis LY Integrals:\n"); */
/*   mat_print(LY, nmo, nmo, stdout); */
/*   printf("MO-Basis LZ Integrals:\n"); */
/*   mat_print(LZ, nmo, nmo, stdout); */

  moinfo.LX = LX;
  moinfo.LY = LY;
  moinfo.LZ = LZ;
}
