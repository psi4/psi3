#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_usotao_new(): Read in an SO to AO transformation matrix 
**
**  takes no arguments.
**
**  returns: double **usotao Read in a num_so by num_ao matrix of doubles
*/


double **file30_rd_usotao_new(void)
{
  double **usotao;
  int num_ao, num_so, i;
  PSI_FPTR usotao_ptr;

  num_ao = file30_rd_nao();
  num_so = file30_rd_nso();
  usotao_ptr = (PSI_FPTR) (info30_.mpoint[40] - 1)*sizeof(int);

  usotao = block_matrix(num_so,num_ao);

  for(i=0;i<num_so;i++)
    wreadw(info30_.filenum, (char *) usotao[i], (int) num_ao*sizeof(double),
	   usotao_ptr, &usotao_ptr);

  return usotao;
}
