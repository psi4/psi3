/*!
  \file rd_usotbf.c
*/

#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_usotbf(): Read in an SO to basis functions transformation matrix 
**
**  takes no arguments.
**
**  returns: double **usotbf Read in a num_so by num_so matrix of doubles
*/


double **file30_rd_usotbf(void)
{
  double **usotbf;
  int num_so, i;
  PSI_FPTR usotbf_ptr;

  num_so = file30_rd_nso();
  usotbf_ptr = (PSI_FPTR) (info30_.mpoint[41] - 1)*sizeof(int);

  usotbf = block_matrix(num_so,num_so);

  for(i=0;i<num_so;i++)
    wreadw(info30_.filenum, (char *) usotbf[i], (int) num_so*sizeof(double),
	   usotbf_ptr, &usotbf_ptr);

  return usotbf;
}
