#include "file30.h"
#include "file30.gbl"

/*
** file30_rd_num_unique_shell():  Reads in the number of symmetry unique shells. 
**
**   takes no arguments.
**
**   returns: int num_unique_shells   number of symmetry unique shells.
*/


int file30_rd_num_unique_shell(void)
{
  return info30_.mconst[7];
}
