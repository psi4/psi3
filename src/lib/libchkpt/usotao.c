/*!
  \file usotao.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include "chkpt.h"
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_usotao(): Read in the SO to AO transformation matrix 
**
** takes no arguments.
**
** returns: usotao = A num_so by num_ao matrix of doubles
** 
** \ingroup (CHKPT)
*/

double **chkpt_rd_usotao(void)
{
  double **usotao;
  int num_ao, num_so, i;
  psio_address ptr;
  char *keyword;
  keyword = chkpt_build_keyword("SO->AO transmat");

  num_ao = chkpt_rd_nao();
  num_so = chkpt_rd_nso();

  usotao = block_matrix(num_so,num_ao);
  ptr = PSIO_ZERO;

  for(i=0;i<num_so;i++)
    psio_read(PSIF_CHKPT, keyword, (char *) usotao[i], (int) num_ao*sizeof(double), ptr, &ptr);

  free(keyword);
  return usotao;
}


/*!
** chkpt_wt_usotao(): Writes out the SO to AO transformation matrix 
**
** \param usotao = A num_so by num_ao matrix of doubles
**
** returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_usotao(double **usotao)
{
  int num_ao, num_so, i;
  psio_address ptr;
  char *keyword;
  keyword = chkpt_build_keyword("SO->AO transmat");

  num_ao = chkpt_rd_nao();
  num_so = chkpt_rd_nso();

  ptr = PSIO_ZERO;
  for(i=0;i<num_so;i++)
    psio_write(PSIF_CHKPT, keyword, (char *) usotao[i], (int) num_ao*sizeof(double), ptr, &ptr);

  free(keyword);
}
