/*!
  \file rd_felement.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** rd_felement():  Reads in element labels including dummy atoms
**
**   takes no arguments.
**
**   returns: char **label element label matrix
*/


char **file30_rd_felement(void)
{
  PSI_FPTR label_ptr;
  char **label;
  int nentry, i;

  nentry = file30_rd_nentry();

  label_ptr = (PSI_FPTR) (info30_.mpoint[28] - 1)*sizeof(int);

  label = (char **)malloc(nentry*sizeof(char*));
  for(i=0;i<nentry;++i)
	label[i] = (char *)malloc(8*sizeof(char));

  for(i=0;i<nentry;++i) 
    wreadw(info30_.filenum, (char *) label[i], (int) 2*sizeof(int), label_ptr, &label_ptr);

  return label;  
}
