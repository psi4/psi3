#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*
** rd_sym_label():  Reads in the symmetry label.
**
**   takes no arguments.
**
**   returns: char *symmetry   symmetry label.
*/


char *file30_rd_sym_label(void)
{
  PSI_FPTR sym_label_ptr;
  char *sym_label;

  sym_label_ptr = (PSI_FPTR) (info30_.mpoint[37] - 1)*sizeof(int);

  sym_label = (char *)malloc(sizeof(int));
  wreadw(info30_.filenum, (char *) sym_label, (int) sizeof(int), sym_label_ptr, &sym_label_ptr);
  sym_label[3] = '\0';

  return sym_label;  
}
