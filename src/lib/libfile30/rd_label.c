/*!
  \file rd_label.c
*/

#include <string.h>
#include "file30.h"
#include "file30.gbl"

/*!
** file30_rd_label():  Reads the main file30 label.
**
**   takes no arguments.
**
**   returns: char *
*/

char *file30_rd_label(void)
{
  int i;
  char *label;
  label = strdup(info30_.label);
  return label;
}
