/*!
  \file wt_schwartz.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_wt_schwartz(): Write out the table of maxima of Schwartz integrals (ij|ij) for
**                       each shell doublet;
**
**  takess:
** \param double** (num_shells by num_shells)
**
**  returns: void
*/


void file30_wt_schwartz(double **schwartz)
{
  int num_shells,i;
  PSI_FPTR schwartz_ptr;

  schwartz_ptr = (PSI_FPTR) (info30_.mpoint[45] - 1)*sizeof(int);

  if (info30_.mpoint[45] == 0)
    /*--- Append at the end of the file ---*/
    schwartz_ptr == (PSI_FPTR) (info30_.mconst[0] - 1)*sizeof(int);
  else {
    schwartz_ptr = (PSI_FPTR) (info30_.mpoint[45] - 1)*sizeof(int);
  }
  num_shells = file30_rd_nshell();

  for(i=0;i<num_shells;i++) {
      wwritw(info30_.filenum, (char *) schwartz[i], num_shells*sizeof(double),
	     schwartz_ptr, &schwartz_ptr);
  }

  return;
}
