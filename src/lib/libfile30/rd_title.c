#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_title(): Reads in the title from file30. 
**
**   takes no arguments.
**
**   returns a char *title 
*/

char *file30_rd_title(void)
{
  return info30_.label;
}
