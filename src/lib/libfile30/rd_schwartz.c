/*!
  \file rd_schwartz.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_rd_schwartz(): Read in the table of maxima of Schwartz integrals (ij|ij) for
**                       each shell doublet;
**
**  takes no arguments.
**
**  returns: NULL if no table is present in the file, double** (num_shells by num_shells) otherwise
**
**               ** MAY THE SCHWARTZ BE WITH YOU!!! **
*/


double **file30_rd_schwartz(void)
{
  double **schwartz;
  int num_shells,i;
  PSI_FPTR schwartz_ptr;

  if (info30_.mpoint[45] == 0)
    return NULL;
  else {
    schwartz_ptr = (PSI_FPTR) (info30_.mpoint[45] - 1)*sizeof(int);
    num_shells = file30_rd_nshell();
    schwartz = block_matrix(num_shells, num_shells);

    for(i=0;i<num_shells;i++) {
      wreadw(info30_.filenum, (char *) schwartz[i], num_shells*sizeof(double),
	     schwartz_ptr, &schwartz_ptr);
    }

    return schwartz;
  }
  
}
